#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "ChemicalEquilibriumCooling.h"
#include "GrainSurfaceChemistry.h"

typedef struct {
    const StructCEQData *initial;
    double timestep;
    double accuracy;
    double chemistryAccuracy;
    bool hasXrayIonization;
} ThermalProblem;

typedef struct {
    double temperature;
    double residual;
    bool converged;
    bool aboveXrayCutoff;
} ThermalPoint;

static double NetHeatingRate(const StructCEQData *state)
{
    return (state->H_tot - state->L_tot) / state->Rho;
}

static CEQLibStatus EvaluateThermalPoint(const ThermalProblem *problem, double temperature,
                                         StructCEQData *trial, ThermalPoint *point)
{
    if(!isfinite(temperature) || temperature <= 0.0)
        return CEQ_NO_CONVERGENCE;

    /* Keep the accepted composition as the branch reference for every trial. */
    *trial = *problem->initial;
    trial->T = temperature;
    trial->Td = EffectiveDustTemperature(trial);
    CEQLibStatus status = CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(trial, problem->chemistryAccuracy);
    if(status != CEQ_SUCCESS)
        return status;

    trial->U = CEQLibConvertTtoU(trial);
    double rate = NetHeatingRate(trial);
    double initialU = problem->initial->U;
    double energyChange = trial->U - initialU;
    double sourceChange = problem->timestep * rate;
    double residual = energyChange - sourceChange;
    if(!isfinite(trial->U) || trial->U <= 0.0 || !isfinite(rate) || !isfinite(residual))
        return CEQ_NO_CONVERGENCE;

    double scale = problem->timestep == 0.0 ? initialU : fmax(fabs(energyChange), fabs(sourceChange));
    double roundoff = 8.0 * DBL_EPSILON * fmax(initialU, trial->U);
    if(problem->timestep > 0.0){
        /* H_tot and L_tot can be much larger than their difference. */
        double rateRoundoff = (8.0 * DBL_EPSILON * fmax(fabs(trial->H_tot), fabs(trial->L_tot))) / trial->Rho;
        roundoff += problem->timestep * rateRoundoff;
    }
    if(!isfinite(roundoff))
        return CEQ_NO_CONVERGENCE;
    point->temperature = temperature;
    point->aboveXrayCutoff = problem->hasXrayIonization &&
                             trial->nel/trial->n[ipHYDROGEN] > 0.5;
    point->residual = residual;
    point->converged = fabs(residual) <= problem->accuracy * scale + roundoff;
    return CEQ_SUCCESS;
}

static CEQLibStatus RefineThermalBracket(const ThermalProblem *problem,
                                         ThermalPoint near, ThermalPoint far, CEQLibStatus boundaryStatus,
                                         StructCEQData *trial)
{
    for(int iteration = 0; iteration < 128; iteration++){
        double fraction = 0.5;
        if(boundaryStatus == CEQ_SUCCESS && iteration % 2 == 0){
            /* Use an interpolated point, with regular bisection to shrink the bracket. */
            double scale = fmax(fabs(near.residual), fabs(far.residual));
            double nearWeight = fabs(near.residual) / scale;
            double farWeight = fabs(far.residual) / scale;
            fraction = nearWeight / (nearWeight + farWeight);
        }
        double temperature = near.temperature + (far.temperature - near.temperature) * fraction;
        double lower = fmin(near.temperature, far.temperature);
        double upper = fmax(near.temperature, far.temperature);
        if(temperature <= lower || temperature >= upper)
            temperature = lower + 0.5 * (upper - lower);
        if(temperature <= lower || temperature >= upper)
            break;

        ThermalPoint next;
        CEQLibStatus status = EvaluateThermalPoint(problem, temperature, trial, &next);
        if(status != CEQ_SUCCESS){
            /* Search up to the first invalid interval from the same side as the initial state. */
            far.temperature = temperature;
            boundaryStatus = status;
        }else if(next.converged){
            return CEQ_SUCCESS;
        }else if((next.residual < 0.0) == (near.residual < 0.0)){
            near = next;
        }else{
            far = next;
            boundaryStatus = CEQ_SUCCESS;
        }
    }
    return boundaryStatus == CEQ_SUCCESS ? CEQ_NO_CONVERGENCE : boundaryStatus;
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
    if(CEQLibRunParameters.UseCOFormation &&
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

static CEQLibStatus VerifyThermalBracket(const ThermalProblem *problem, ThermalPoint near,
                                         ThermalPoint far, double factor, StructCEQData *trial)
{
    /* Probe the original spacing until the wide trial so the first detected crossing wins. */
    StructCEQData probe;
    for(double temperature = near.temperature*factor;
        factor > 1.0 ? temperature < far.temperature : temperature > far.temperature;
        temperature = near.temperature*factor){
        if(temperature == near.temperature) break;
        ThermalPoint point = {.temperature = temperature};
        CEQLibStatus status = EvaluateThermalPoint(problem,temperature,&probe,&point);
        if(status != CEQ_SUCCESS)
            return RefineThermalBracket(problem,near,point,status,trial);
        if(point.converged){
            *trial = probe;
            return CEQ_SUCCESS;
        }
        if((point.residual < 0.0) != (near.residual < 0.0))
            return RefineThermalBracket(problem,near,point,CEQ_SUCCESS,trial);
        near = point;
    }
    if(far.converged) return CEQ_SUCCESS;
    return RefineThermalBracket(problem,near,far,CEQ_SUCCESS,trial);
}

static CEQLibStatus SearchThermalState(const ThermalProblem *problem, double temperature, StructCEQData *trial)
{
    ThermalPoint near;
    CEQLibStatus status = EvaluateThermalPoint(problem, temperature, trial, &near);
    if(status != CEQ_SUCCESS || near.converged)
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
           (status != CEQ_SUCCESS || far.aboveXrayCutoff != near.aboveXrayCutoff ||
            ((near.residual < 0.0) == (far.residual < 0.0) && fabs(far.residual) > fabs(near.residual)))){
            /* Recover the fine grid at invalid or discontinuous intervals and on a reversal. */
            far.temperature = defaultTemperature;
            status = EvaluateThermalPoint(problem,defaultTemperature,trial,&far);
            stride = 1;
        }
        if(status != CEQ_SUCCESS)
            return RefineThermalBracket(problem,near,far,status,trial);
        const bool crossing = (near.residual < 0.0) != (far.residual < 0.0);
        if(stride > 1 && (far.converged || crossing))
            return VerifyThermalBracket(problem,near,far,factor,trial);
        if(far.converged)
            return CEQ_SUCCESS;
        if(crossing)
            return RefineThermalBracket(problem,near,far,CEQ_SUCCESS,trial);
        remaining -= stride;
        previous = near;
        near = far;
    }
    return CEQ_NO_CONVERGENCE;
}

static CEQLibStatus SolveThermalState(const StructCEQData *initial, double timestep, double accuracy,
                                      StructCEQData *trial)
{
    if(initial == NULL || !isfinite(initial->U) || initial->U <= 0.0 ||
       !isfinite(initial->Rho) || initial->Rho <= 0.0)
        return CEQ_INVALID_INPUT;

    double temperature = initial->T;
    if(!isfinite(temperature) || temperature <= 0.0){
        *trial = *initial;
        temperature = CEQLibConvertUtoT(trial);
    }
    if(!isfinite(temperature) || temperature <= 0.0)
        return CEQ_INVALID_INPUT;

    const double minimumChemistryAccuracy = 16.0 * DBL_EPSILON;
    ThermalProblem problem = {
        .initial = initial,
        .timestep = timestep,
        .accuracy = accuracy,
        .hasXrayIonization = CEQLibRunParameters.UseSoftXrayIonization && initial->X0 > 0.0,
        .chemistryAccuracy = fmax(minimumChemistryAccuracy,
                                  fmin(CEQLibRunParameters.ChemistryAccuracy, 0.1 * accuracy))
    };
    CEQLibStatus status = SearchThermalState(&problem, temperature, trial);
    if(status == CEQ_NO_CONVERGENCE && problem.chemistryAccuracy > minimumChemistryAccuracy){
        /* Near thermal balance, small changes in the chemistry can dominate the net rate. */
        problem.chemistryAccuracy = minimumChemistryAccuracy;
        status = SearchThermalState(&problem, temperature, trial);
    }
    return status;
}

CEQLibStatus CEQLibTryNetCoolingRateChemcalEquilibrium(StructCEQData *state, double *rate)
{
    if(rate == NULL)
        return CEQ_INVALID_INPUT;

    StructCEQData trial;
    CEQLibStatus status = SolveThermalState(state, 0.0, 1.e-8, &trial);
    if(status == CEQ_SUCCESS){
        trial.U = state->U;
        *state = trial;
        *rate = NetHeatingRate(state);
    }
    return status;
}

CEQLibStatus CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructCEQData *state,
                                                                          double accuracy, double *Du)
{
    if(state == NULL || Du == NULL || !isfinite(state->dt) || state->dt <= 0.0 ||
       !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0)
        return CEQ_INVALID_INPUT;

    StructCEQData trial;
    CEQLibStatus status = SolveThermalState(state, state->dt, accuracy, &trial);
    if(status == CEQ_SUCCESS){
        double rate = (trial.U - state->U) / state->dt;
        if(!isfinite(rate))
            return CEQ_NO_CONVERGENCE;
        /* The caller still applies U += dt * Du. */
        trial.U = state->U;
        *state = trial;
        *Du = rate;
    }
    return status;
}

double CEQLibNetCoolingRateChemcalEquilibrium(StructCEQData *state)
{
    double rate;
    return CEQLibTryNetCoolingRateChemcalEquilibrium(state, &rate) == CEQ_SUCCESS ? rate : NAN;
}

double CEQLibEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructCEQData *state, double accuracy)
{
    double rate;
    return CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(state, accuracy, &rate) == CEQ_SUCCESS ? rate : NAN;
}
