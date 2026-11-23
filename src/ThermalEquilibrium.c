#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "ChemicalEquilibriumCooling.h"
#include "CoolingHeatingRates.h"
#include "EquilibriumAbundance.h"
#include "GrainSurfaceChemistry.h"
#include "PhotoelectronDeposition.h"
#ifdef EQC_LEGACY_XRAY_CUTOFF
#include "Xray.h"
#endif

typedef struct {
    const StructEQCData *initial;
    const CEQUVBRates *uvb;
    CEQUVBEquilibriumStats *stats;
    CEQSecondaryIonizationSeed *secondarySeed;
    double timestep;
    double accuracy;
    double chemistryAccuracy;
    bool fixedMolecularHydrogen;
#ifdef EQC_LEGACY_XRAY_CUTOFF
    bool hasXrayIonization;
#endif
} ThermalProblem;

typedef struct {
    double temperature;
    double residual;
    bool converged;
#ifdef EQC_LEGACY_XRAY_CUTOFF
    bool aboveXrayCutoff;
#endif
} ThermalPoint;

static double NetHeatingRate(const StructEQCData *state)
{
    return (state->H_tot - state->L_tot) / state->Rho;
}

static EQCLibStatus EvaluateEquilibriumRates(const ThermalProblem *problem, StructEQCData *trial,
                                             CEQSecondaryIonizationSeed *seed)
{
    if(problem->fixedMolecularHydrogen)
        return CEQCalculateEquilibriumRatesHoldingH2(trial,problem->chemistryAccuracy,problem->uvb,seed);
    return CEQCalculateChemicalAbundanceCoolingHeatingRatesWithSeed(
        trial,problem->chemistryAccuracy,problem->uvb,problem->stats,seed);
}

static EQCLibStatus EvaluateThermalPoint(const ThermalProblem *problem, double temperature,
                                         StructEQCData *trial, ThermalPoint *point)
{
    if(!isfinite(temperature) || temperature <= 0.0)
        return EQC_NO_CONVERGENCE;

    /* Keep the accepted composition as the branch reference for every trial. */
    *trial = *problem->initial;
    trial->T = temperature;
    trial->Td = EffectiveDustTemperature(trial);
    if(problem->stats) problem->stats->temperatureTrials++;
    CEQSecondaryIonizationSeed seed = {0};
    if(problem->secondarySeed) seed = *problem->secondarySeed;
    EQCLibStatus status = EvaluateEquilibriumRates(problem,trial,problem->secondarySeed ? &seed : NULL);
    if(status != EQC_SUCCESS && problem->secondarySeed && (seed.hydrogen > 0.0 || seed.helium > 0.0)){
        seed = (CEQSecondaryIonizationSeed){0};
        status = EvaluateEquilibriumRates(problem,trial,&seed);
    }
    if(status != EQC_SUCCESS)
        return status;

    trial->U = EQCLibConvertTtoU(trial);
    double rate = NetHeatingRate(trial);
    double initialU = problem->initial->U;
    double energyChange = trial->U - initialU;
    double sourceChange = problem->timestep * rate;
    double residual = energyChange - sourceChange;
    if(!isfinite(trial->U) || trial->U <= 0.0 || !isfinite(rate) || !isfinite(residual))
        return EQC_NO_CONVERGENCE;

    double scale = problem->timestep == 0.0 ? initialU : fmax(fabs(energyChange), fabs(sourceChange));
    double roundoff = 8.0 * DBL_EPSILON * fmax(initialU, trial->U);
    if(problem->timestep > 0.0){
        /* H_tot and L_tot can be much larger than their difference. */
        double rateRoundoff = (8.0 * DBL_EPSILON * fmax(fabs(trial->H_tot), fabs(trial->L_tot))) / trial->Rho;
        roundoff += problem->timestep * rateRoundoff;
    }
    if(!isfinite(roundoff))
        return EQC_NO_CONVERGENCE;
    point->temperature = temperature;
#ifdef EQC_LEGACY_XRAY_CUTOFF
    point->aboveXrayCutoff = problem->hasXrayIonization &&
                            trial->nel/trial->n[ipHYDROGEN] > EQC_XRAY_IONIZATION_MAX_XE;
#endif
    point->residual = residual;
    point->converged = fabs(residual) <= problem->accuracy * scale + roundoff;
    if(problem->secondarySeed) *problem->secondarySeed = seed;
    return EQC_SUCCESS;
}

static EQCLibStatus RefineThermalBracket(const ThermalProblem *problem,
                                         ThermalPoint near, ThermalPoint far, EQCLibStatus boundaryStatus,
                                         StructEQCData *trial)
{
    bool continueInterpolation = false;
    ThermalPoint displaced = {0};
    for(int iteration = 0; iteration < 128; iteration++){
        double fraction = 0.5;
        const bool interpolate = boundaryStatus == EQC_SUCCESS &&
            (iteration % 2 == 0 || continueInterpolation);
        double bestResidual = 0.0;
        if(interpolate){
            /* Use an interpolated point, with regular bisection to shrink the bracket. */
            bestResidual = fmin(fabs(near.residual),fabs(far.residual));
            double scale = fmax(fabs(near.residual), fabs(far.residual));
            double nearWeight = fabs(near.residual) / scale;
            double farWeight = fabs(far.residual) / scale;
            fraction = nearWeight / (nearWeight + farWeight);
        }
        double temperature = near.temperature + (far.temperature - near.temperature) * fraction;
        double lower = fmin(near.temperature, far.temperature);
        double upper = fmax(near.temperature, far.temperature);
        if(continueInterpolation && interpolate){
            /* Inverse quadratic interpolation at zero residual, using the
               last displaced endpoint. Keep the secant if it leaves the bracket. */
            const double scale = fmax(fabs(displaced.residual),fmax(fabs(near.residual),fabs(far.residual)));
            const double a = near.residual/scale, b = far.residual/scale, c = displaced.residual/scale;
            if(a != b && a != c && b != c){
                const double farWeight = (a/(b-a))*(c/(b-c));
                const double thirdWeight = (a/(c-a))*(b/(c-b));
                const double estimate = near.temperature + farWeight*(far.temperature-near.temperature)
                                                     + thirdWeight*(displaced.temperature-near.temperature);
                if(isfinite(estimate) && estimate > lower && estimate < upper) temperature = estimate;
            }
        }
        if(temperature <= lower || temperature >= upper)
            temperature = lower + 0.5 * (upper - lower);
        if(temperature <= lower || temperature >= upper)
            break;

        ThermalPoint next;
        EQCLibStatus status = EvaluateThermalPoint(problem, temperature, trial, &next);
        /* Keep the enclosure and retry interpolation while it makes rapid
           residual progress. Otherwise return to the original bisections. */
        continueInterpolation = problem->secondarySeed && interpolate && status == EQC_SUCCESS &&
            fabs(next.residual) <= 0.25*bestResidual;
        if(status != EQC_SUCCESS){
            /* Search up to the first invalid interval from the same side as the initial state. */
            far.temperature = temperature;
            boundaryStatus = status;
        }else if(next.converged){
            return EQC_SUCCESS;
        }else if((next.residual < 0.0) == (near.residual < 0.0)){
            displaced = near;
            near = next;
        }else{
            displaced = far;
            far = next;
            boundaryStatus = EQC_SUCCESS;
        }
    }
    return boundaryStatus == EQC_SUCCESS ? EQC_NO_CONVERGENCE : boundaryStatus;
}

static double PredictThermalTemperature(ThermalPoint previous, ThermalPoint near, double factor,
                                        int remaining, int *stride)
{
    *stride = 1;
    const double defaultTemperature = near.temperature*factor;
    const double residual = fabs(near.residual);
    const double decrease = fabs(previous.residual)-residual;
    if(decrease <= 0.0)
        return defaultTemperature;

    const double logStep = 0.8*fabs(log(near.temperature/previous.temperature))*(residual/decrease);
    if(!isfinite(logStep) || logStep <= fabs(log(factor)))
        return defaultTemperature;

    const double prediction = near.temperature*exp(copysign(fmin(logStep,log(1.5)),factor-1.0));
    if(!isfinite(prediction) || prediction <= 0.0)
        return defaultTemperature;
    /* Keep the usual spacing when approaching the CO chemistry switch at 10^4 K. */
    const double carbonBoundary = 1.e4;
    const double minimumFactor = sqrt(1.1);
    if(EQCLibRunParameters.UseCOFormation &&
       fmin(near.temperature,prediction) <= carbonBoundary*minimumFactor &&
       fmax(near.temperature,prediction) >= carbonBoundary/minimumFactor)
        return defaultTemperature;
    /* Keep the original sample lattice so narrow branches near a guard
       do not depend on the phase of a preceding wide step. */
    double temperature = defaultTemperature;
    for(; *stride < remaining; (*stride)++){
        double nextTemperature = temperature*factor;
        if(nextTemperature == temperature ||
           (factor > 1.0 ? nextTemperature > prediction : nextTemperature < prediction))
            break;
        temperature = nextTemperature;
    }
    return temperature;
}

static EQCLibStatus VerifyThermalBracket(const ThermalProblem *problem, ThermalPoint near,
                                         ThermalPoint far, double factor, StructEQCData *trial)
{
    /* Probe the original spacing until the wide trial so the first detected crossing wins. */
    StructEQCData probe;
    for(double temperature = near.temperature*factor;
        factor > 1.0 ? temperature < far.temperature : temperature > far.temperature;
        temperature = near.temperature*factor){
        if(temperature == near.temperature) break;
        ThermalPoint point = {.temperature = temperature};
        EQCLibStatus status = EvaluateThermalPoint(problem,temperature,&probe,&point);
        if(status != EQC_SUCCESS)
            return RefineThermalBracket(problem,near,point,status,trial);
        if(point.converged){
            *trial = probe;
            return EQC_SUCCESS;
        }
        if((point.residual < 0.0) != (near.residual < 0.0))
            return RefineThermalBracket(problem,near,point,EQC_SUCCESS,trial);
        near = point;
    }
    if(far.converged) return EQC_SUCCESS;
    return RefineThermalBracket(problem,near,far,EQC_SUCCESS,trial);
}

static EQCLibStatus SearchThermalState(const ThermalProblem *problem, double temperature, StructEQCData *trial)
{
    ThermalPoint near;
    EQCLibStatus status = EvaluateThermalPoint(problem, temperature, trial, &near);
    if(status != EQC_SUCCESS || near.converged)
        return status;

    const double factor = near.residual < 0.0 ? sqrt(1.1) : 1.0 / sqrt(1.1);
    ThermalPoint previous = near;
    /* Bound the temperature range by the original number of lattice steps. */
    for(int remaining = 1024; remaining > 0;){
        const double defaultTemperature = near.temperature*factor;
        int stride;
        temperature = PredictThermalTemperature(previous,near,factor,remaining,&stride);
        if(!isfinite(temperature) || temperature <= 0.0 || temperature == near.temperature)
            break;

        ThermalPoint far = {.temperature = temperature};
        status = EvaluateThermalPoint(problem,temperature,trial,&far);
        if(stride > 1 &&
           (status != EQC_SUCCESS ||
#ifdef EQC_LEGACY_XRAY_CUTOFF
            far.aboveXrayCutoff != near.aboveXrayCutoff ||
#endif
            ((near.residual < 0.0) == (far.residual < 0.0) && fabs(far.residual) > fabs(near.residual)))){
            /* Recover the fine grid at invalid intervals, legacy jumps, and on a reversal. */
            far.temperature = defaultTemperature;
            status = EvaluateThermalPoint(problem,defaultTemperature,trial,&far);
            stride = 1;
        }
        if(status != EQC_SUCCESS)
            return RefineThermalBracket(problem,near,far,status,trial);
        const bool crossing = (near.residual < 0.0) != (far.residual < 0.0);
        if(stride > 1 && (far.converged || crossing))
            return VerifyThermalBracket(problem,near,far,factor,trial);
        if(far.converged)
            return EQC_SUCCESS;
        if(crossing)
            return RefineThermalBracket(problem,near,far,EQC_SUCCESS,trial);
        remaining -= stride;
        previous = near;
        near = far;
    }
    return EQC_NO_CONVERGENCE;
}

static EQCLibStatus SolveThermalState(const StructEQCData *initial, double timestep, double accuracy,
                                      StructEQCData *trial, CEQUVBEquilibriumStats *stats,
                                      bool fixedMolecularHydrogen)
{
    if(initial == NULL || !isfinite(initial->U) || initial->U <= 0.0 ||
       !isfinite(initial->Rho) || initial->Rho <= 0.0)
        return EQC_INVALID_INPUT;

    double temperature = initial->T;
    if(!isfinite(temperature) || temperature <= 0.0){
        *trial = *initial;
        temperature = EQCLibConvertUtoT(trial);
    }
    if(!isfinite(temperature) || temperature <= 0.0)
        return EQC_INVALID_INPUT;

    CEQUVBRates uvb;
    EQCLibStatus status = CEQEvaluateUVBackgroundWithMetallicity(initial->NHIuvb,initial->Z,&uvb);
    if(status != EQC_SUCCESS) return status;

    const double minimumChemistryAccuracy = 16.0 * DBL_EPSILON;
    /* Only the secondary coefficients are carried between temperatures. The
       accepted cell remains the reference for charge and molecular branches. */
    const bool hasSecondaryPower = CEQHasUVBHighPower(&uvb,initial->n[ipHELIUM]);
    CEQSecondaryIonizationSeed initialSeed = {0};
    if(hasSecondaryPower){
        const CEQPhotoelectronDeposition deposit = CEQDepositUVBPhotoelectrons(initial,&uvb);
        initialSeed = (CEQSecondaryIonizationSeed){deposit.ionizationHI,deposit.ionizationHeI};
    }
    CEQSecondaryIonizationSeed secondarySeed = initialSeed;
    ThermalProblem problem = {
        .initial = initial,
        .uvb = &uvb,
        .stats = stats,
        .secondarySeed = hasSecondaryPower ? &secondarySeed : NULL,
        .timestep = timestep,
        .accuracy = accuracy,
        .fixedMolecularHydrogen = fixedMolecularHydrogen,
#ifdef EQC_LEGACY_XRAY_CUTOFF
        .hasXrayIonization = EQCLibRunParameters.UseSoftXrayIonization && initial->X0 > 0.0,
#endif
        .chemistryAccuracy = fmax(minimumChemistryAccuracy,
                                  fmin(EQCLibRunParameters.ChemistryAccuracy, 0.1 * accuracy))
    };
    status = SearchThermalState(&problem, temperature, trial);
    if(status == EQC_NO_CONVERGENCE && problem.chemistryAccuracy > minimumChemistryAccuracy){
        /* Near thermal balance, small changes in the chemistry can dominate the net rate. */
        double intermediateAccuracy = 0.1 * problem.chemistryAccuracy;
        problem.chemistryAccuracy = minimumChemistryAccuracy;
        secondarySeed = initialSeed;
        status = SearchThermalState(&problem, temperature, trial);
        /* An almost machine-precision chemistry solve can itself stall. Try
           intermediate accuracies while retaining the thermal residual test. */
        for(; status == EQC_NO_CONVERGENCE && intermediateAccuracy > minimumChemistryAccuracy;
            intermediateAccuracy *= 0.1){
            problem.chemistryAccuracy = intermediateAccuracy;
            secondarySeed = initialSeed;
            status = SearchThermalState(&problem, temperature, trial);
        }
    }
    return status;
}

static EQCLibStatus NetCoolingRateAtEnergy(StructEQCData *state, double *rate, bool fixedMolecularHydrogen)
{
    if(rate == NULL)
        return EQC_INVALID_INPUT;

    StructEQCData trial;
    EQCLibStatus status = SolveThermalState(state,0.0,1.e-8,&trial,NULL,fixedMolecularHydrogen);
    if(status == EQC_SUCCESS){
        trial.U = state->U;
        *state = trial;
        *rate = NetHeatingRate(state);
    }
    return status;
}

EQCLibStatus EQCLibTryNetCoolingRateChemcalEquilibrium(StructEQCData *state, double *rate)
{
    return NetCoolingRateAtEnergy(state,rate,false);
}

EQCLibStatus EQCLibTryNetCoolingRateFixedH2(StructEQCData *state, double *rate)
{
    return NetCoolingRateAtEnergy(state,rate,true);
}

EQCLibStatus EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructEQCData *state,
                                                                          double accuracy, double *Du)
{
    return CEQSolveEquilibriumEnergyWithStats(state,accuracy,Du,NULL);
}

EQCLibStatus CEQSolveEquilibriumEnergyWithStats(StructEQCData *state, double accuracy, double *Du,
                                              CEQUVBEquilibriumStats *stats)
{
    if(state == NULL || Du == NULL || !isfinite(state->dt) || state->dt <= 0.0 ||
       !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0)
        return EQC_INVALID_INPUT;

    StructEQCData trial;
    EQCLibStatus status = SolveThermalState(state,state->dt,accuracy,&trial,stats,false);
    if(status == EQC_SUCCESS){
        double rate = (trial.U - state->U) / state->dt;
        if(!isfinite(rate))
            return EQC_NO_CONVERGENCE;
        /* The caller still applies U += dt * Du. */
        trial.U = state->U;
        *state = trial;
        *Du = rate;
    }
    return status;
}

/* At fixed T and charge, an exponential step integrates both formation and
   loss with frozen rate coefficients. Sampling again at the new composition
   resolves H2 self-shielding and composition-dependent collisional loss. */
static double MolecularFractionStep(StructEQCData *trial, const CEQMolecularHydrogenCoefficients *coefficients,
                                    double neutralFraction, double fraction, double timestep)
{
    const double nH = trial->n[ipHYDROGEN];
    trial->nmol[ipH2] = (0.5*nH)*fraction;
    trial->natom[ipHYDROGEN] = (nH-2.0*trial->nmol[ipH2])*neutralFraction;
    const double formation = 2.0*coefficients->formation*neutralFraction;
    const double destruction = CEQH2DestructionRate(trial,coefficients);
    const double total = formation + destruction;
    if(!isfinite(total) || !(formation >= 0.0 && destruction >= 0.0)) return NAN;
    if(total == 0.0 || timestep == 0.0) return fraction;
    const double exponent = -total*timestep;
    return fmin(1.0,fraction*exp(exponent) + (formation/total)*(-expm1(exponent)));
}

static EQCLibStatus AdvanceMolecularHydrogen(StructEQCData *state, double timestep, double accuracy)
{
    const CEQMolecularHydrogenCoefficients coefficients = CEQPrepareMolecularHydrogenCoefficients(state);
    const double nH = state->n[ipHYDROGEN];
    const double atomicDensity = nH - 2.0*state->nmol[ipH2];
    double neutralFraction = 1.0;
    if(atomicDensity > 0.0){
        neutralFraction = fmin(1.0,state->natom[ipHYDROGEN]/atomicDensity);
    }else if(state->nel > 0.0 && coefficients.formation > 0.0){
        const CEQHydrogenRates rates = CEQHydrogenIonizationRates(state);
        double ionized, unused;
        const int status = CEQSolveIonLadder(1.0,rates.ionization,rates.recombination,0.0,0.0,
                                            &neutralFraction,&ionized,&unused);
        if(status != EQC_SUCCESS) return (EQCLibStatus)status;
    }

    const double initialFraction = 2.0*(state->nmol[ipH2]/nH);
    double fraction = initialFraction;
    double remaining = timestep, step = timestep;
    StructEQCData trial = *state;
    for(int attempts = 0; remaining > 0.0 && attempts < 4096; attempts++){
        const double coarse = MolecularFractionStep(&trial,&coefficients,neutralFraction,fraction,step);
        const double firstHalf = 0.5*step;
        const double half = MolecularFractionStep(&trial,&coefficients,neutralFraction,fraction,firstHalf);
        const double fine = MolecularFractionStep(&trial,&coefficients,neutralFraction,half,step-firstHalf);
        if(!isfinite(coarse) || !isfinite(fine)) return EQC_NO_CONVERGENCE;
        const double tolerance = accuracy*(1.e-4 + fmax(fraction,fine));
        const double error = fabs(fine-coarse);
        if(error <= tolerance){
            if(remaining-step == remaining) return EQC_NO_CONVERGENCE;
            fraction = fine;
            remaining -= step;
            if(error < 0.25*tolerance) step *= 2.0;
        }else{
            step *= 0.5;
        }
        step = fmin(step,remaining);
    }
    if(remaining > 0.0) return EQC_NO_CONVERGENCE;
    if(fraction != initialFraction) state->nmol[ipH2] = (0.5*nH)*fraction;
    return EQC_SUCCESS;
}

EQCLibStatus CEQAdvanceH2Thermochemistry(StructEQCData *state, double accuracy)
{
    if(!EQCLibRunParameters.UseH2Formation || EQCLibRunParameters.UseCOFormation ||
       state->nmol[ipCO] != 0.0) return EQC_INVALID_INPUT;

    /* Equilibrate the charge and map U to T without changing H2. Evolve H2
       for half a step on each side of the implicit energy update. */
    StructEQCData composition, thermal;
    EQCLibStatus status = SolveThermalState(state,0.0,accuracy,&composition,NULL,true);
    if(status != EQC_SUCCESS) return status;
    composition.U = state->U;
    const double half = 0.5*state->dt;
    status = AdvanceMolecularHydrogen(&composition,half,accuracy);
    if(status != EQC_SUCCESS) return status;
    status = SolveThermalState(&composition,state->dt,accuracy,&thermal,NULL,true);
    if(status != EQC_SUCCESS) return status;

    const double before = thermal.nmol[ipH2];
    status = AdvanceMolecularHydrogen(&thermal,state->dt-half,accuracy);
    if(status != EQC_SUCCESS) return status;
    if(thermal.nmol[ipH2] == before){
        *state = thermal;
    }else{
        /* Refresh the equation of state and rates at the accepted U and H2. */
        status = SolveThermalState(&thermal,0.0,accuracy,&composition,NULL,true);
        if(status != EQC_SUCCESS) return status;
        composition.U = thermal.U;
        *state = composition;
    }
    return EQC_SUCCESS;
}

double EQCLibNetCoolingRateChemcalEquilibrium(StructEQCData *state)
{
    double rate;
    return EQCLibTryNetCoolingRateChemcalEquilibrium(state, &rate) == EQC_SUCCESS ? rate : NAN;
}

double EQCLibEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructEQCData *state, double accuracy)
{
    double rate;
    return EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(state, accuracy, &rate) == EQC_SUCCESS ? rate : NAN;
}
