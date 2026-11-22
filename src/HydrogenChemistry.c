#include <float.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "ChemicalEquilibriumCooling.h"
#include "CoolingHeatingRates.h"
#include "EquilibriumAbundance.h"
#include "PhotoelectronDeposition.h"

enum {
    Ionized, LogEnergy, FirstLadder,
    MaxVariables = FirstLadder+7, /* Three ion ladders and optional H2. */
};

typedef struct {
    StructEQCData state;
    const CEQUVBRates *uvb;
    CEQRateCoefficients coefficients;
    double initialEnergy, siliconElectrons;
    int variables, ladderStart;
    bool evolveEnergy;
    bool chargeFree;
    bool hasFastPhotoelectrons;
    int molecularOffset; /* H2, when active, follows the ion ladders. */
    int oxygenOffset, heliumOffset, carbonOffset; /* Zero when the element is absent. */
    bool hasLastDerivative;
    double lastY[MaxVariables], lastDerivative[MaxVariables];
} ThermochemistryProblem;

static double NeutralFraction(const double fractions[])
{
    return (1.0-fmax(fractions[0],fractions[1]))-fmin(fractions[0],fractions[1]);
}

/* Implicit solver stages can leave the physical domain. Keep the rate
   coefficients bounded, preserving the smaller fraction at the neutral boundary. */
static void BoundFractions(double fractions[])
{
    fractions[0] = fmax(0.0,fmin(1.0,fractions[0]));
    fractions[1] = fmax(0.0,fmin(1.0,fractions[1]));
    const int larger = fractions[0] >= fractions[1] ? 0 : 1;
    const int smaller = 1-larger;
    fractions[larger] = fmin(fractions[larger],1.0-fractions[smaller]);
}

static bool FractionsWithinBounds(const double fractions[], double tolerance)
{
    return isfinite(fractions[0]) && isfinite(fractions[1]) &&
           fractions[0] >= -tolerance && fractions[1] >= -tolerance &&
           NeutralFraction(fractions) >= -tolerance;
}

static void SetLadderAbundances(double total, const double fractions[], int offset,
                                double *neutral, double *singly, double *doubly)
{
    *neutral = total;
    *singly = *doubly = 0.0;
    if(offset){
        const double *ions = fractions+offset;
        *neutral = total*fmax(0.0,NeutralFraction(ions));
        *singly = total*ions[0];
        *doubly = total*ions[1];
    }
}

static void SetAbundances(StructEQCData *state, const double hydrogen[], const double fractions[],
                          const ThermochemistryProblem *problem)
{
    const double nH = state->n[ipHYDROGEN];
    state->natom[ipHYDROGEN] = nH*fmax(0.0,NeutralFraction(hydrogen));
    state->nion[ipHYDROGEN] = nH*hydrogen[0];
    state->nmol[ipH2] = nH*(0.5*hydrogen[1]);
    SetLadderAbundances(state->n[ipOXYGEN],fractions,problem->oxygenOffset,
                        &state->natom[ipOXYGEN],&state->nion[ipOXYGEN],&state->nOIII);
    SetLadderAbundances(state->n[ipHELIUM],fractions,problem->heliumOffset,
                        &state->natom[ipHELIUM],&state->nion[ipHELIUM],&state->nHeIII);
    SetLadderAbundances(state->n[ipCARBON],fractions,problem->carbonOffset,
                        &state->natom[ipCARBON],&state->nion[ipCARBON],&state->nCIII);
    state->nel = EQCLibGetElectronNumberDensity(state);
}

static bool SetThermalState(StructEQCData *state, double energy)
{
    if(!isfinite(energy) || energy <= 0.0) return false;
    state->U = energy;
    state->MeanMolecularWeight = EQCLibGetMeanMolecularWeight(state);
    state->DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(state);
    if(!isfinite(state->MeanMolecularWeight) || state->MeanMolecularWeight <= 0.0 ||
       !isfinite(state->DegreeOfFreedom) || state->DegreeOfFreedom <= 0.0) return false;

    state->Gamma = (1.0+state->DegreeOfFreedom)/state->DegreeOfFreedom;
    state->Gm1 = state->Gamma-1.0;
    state->GGm1 = state->Gamma*state->Gm1;
    state->T = EQCLibConvertUtoT(state);
    if(!isfinite(state->T) || state->T <= 0.0) return false;
    state->invT = 1.0/state->T;
    state->logT = log10(state->T);
    return isfinite(state->invT);
}

/* Reflect negative trial ion fractions only in the electron-dependent rate
   coefficients, so the signed recombination flux restores the boundary. */
static double TrialIonFraction(double stage, double bounded)
{
    return stage < 0.0 ? fmin(-stage,1.0) : bounded;
}

static double TrialLadderCharge(const double stage[], const double bounded[])
{
    return TrialIonFraction(stage[0],bounded[0]) + 2.0*TrialIonFraction(stage[1],bounded[1]);
}

static bool ValidRate(double rate)
{
    return isfinite(rate) && rate >= 0.0;
}

/* The two net upward fluxes per nucleus on a three-state ion ladder. */
static bool LadderFluxes(CEQIonLadderRates rates, const double ions[], double *lower, double *upper)
{
    if(!ValidRate(rates.ionization0) || !ValidRate(rates.recombination1) ||
       !ValidRate(rates.ionization1) || !ValidRate(rates.recombination2)) return false;
    *lower = rates.ionization0*NeutralFraction(ions) - rates.recombination1*ions[0];
    *upper = rates.ionization1*ions[0] - rates.recombination2*ions[1];
    return true;
}

static bool LadderDerivative(CEQIonLadderRates rates, const double ions[], double derivative[])
{
    double lower, upper;
    if(!LadderFluxes(rates,ions,&lower,&upper)) return false;
    derivative[0] = lower-upper;
    derivative[1] = upper;
    return true;
}

static bool OxygenDerivative(const ThermochemistryProblem *problem, const double y[], double neutral,
                             double derivative[])
{
    const StructEQCData *state = &problem->state;
    const CEQRateCoefficients *coefficients = &problem->coefficients;
    const int oxygen = problem->oxygenOffset;
    const double *ions = y+oxygen;
    double toOII, toOIII;
    if(!LadderFluxes(CEQOxygenIonizationRates(state,coefficients),ions,&toOII,&toOIII)) return false;
    if(EQCLibRunParameters.UseChargeTransfer){
        const double nH = state->n[ipHYDROGEN];
        const double nHI = nH*neutral, nHII = nH*y[Ionized];
        if(!ValidRate(coefficients->toOxygenII) || !ValidRate(coefficients->toHydrogenI) ||
           !ValidRate(coefficients->oxygenIIIToHydrogenI)) return false;

        /* Per oxygen nucleus: O I -> O II and O III -> O II. Each changes
           the H II population by the opposite change in oxygen charge. */
        const double exchangeOI = coefficients->toOxygenII*nHII*NeutralFraction(ions)
                                - coefficients->toHydrogenI*nHI*ions[0];
        const double exchangeOIII = coefficients->oxygenIIIToHydrogenI*nHI*ions[1];
        toOII += exchangeOI;
        toOIII -= exchangeOIII;
        derivative[Ionized] += (state->n[ipOXYGEN]/nH)*(exchangeOIII-exchangeOI);
    }
    derivative[oxygen] = toOII-toOIII;
    derivative[oxygen+1] = toOIII;
    return true;
}

static int ChemistryDerivative(double time, const double y[], double derivative[], void *parameters)
{
    (void)time;
    ThermochemistryProblem *problem = parameters;
    /* The system is autonomous at fixed radiation. The driver and Jacobian
       can request the same stage, so reuse an exactly matching evaluation. */
    const size_t bytes = (size_t)problem->variables*sizeof(*y);
    if(problem->hasLastDerivative && memcmp(y,problem->lastY,bytes) == 0){
        memcpy(derivative,problem->lastDerivative,bytes);
        return GSL_SUCCESS;
    }
    StructEQCData *state = &problem->state;
    const CEQRateCoefficients *coefficients = &problem->coefficients;
    const double nH = state->n[ipHYDROGEN];
    const int oxygen = problem->oxygenOffset, helium = problem->heliumOffset, carbon = problem->carbonOffset;
    const int molecularOffset = problem->molecularOffset;
    const int ladderEnd = molecularOffset ? molecularOffset : problem->variables;

    double fractions[MaxVariables] = {0};
    for(int i=0;i<problem->variables;i++){
        if(!isfinite(y[i])) return GSL_FAILURE;
        fractions[i] = y[i];
    }
    const double trialHydrogen[] = {problem->chargeFree ? 0.0 : y[Ionized],
                                   molecularOffset ? y[molecularOffset] : 0.0};
    double hydrogen[] = {trialHydrogen[0],trialHydrogen[1]};
    BoundFractions(hydrogen);
    for(int i=problem->ladderStart;i<ladderEnd;i+=2) BoundFractions(fractions+i);
    if(problem->evolveEnergy || problem->hasFastPhotoelectrons)
        SetAbundances(state,hydrogen,fractions,problem);
    CEQPhotoelectronDeposition deposition = {0};
    if(problem->evolveEnergy || problem->hasFastPhotoelectrons)
        deposition = CEQDepositUVBPhotoelectrons(state,problem->uvb);
    if(problem->evolveEnergy){
        /* The energy variable is ln(U/U_initial), so all thermal stages have
           positive energy and the tolerance does not depend on its CGS scale. */
        if(!SetThermalState(state,problem->initialEnergy*exp(y[LogEnergy]))) return GSL_FAILURE;
        CEQTemperaturePoint temperaturePoint;
        const CEQTemperaturePoint *point = CEQLocateLogTemperature(state->logT,&temperaturePoint)
                                         ? &temperaturePoint : NULL;
        if(CEQCalculateNonEquilibriumEnergyRates(state,point,deposition.heating) != EQC_SUCCESS) return GSL_FAILURE;
        derivative[LogEnergy] = (state->H_tot-state->L_tot)/state->Rho/state->U;
        if(!isfinite(derivative[LogEnergy])) return GSL_FAILURE;
        CEQRefreshNonEquilibriumRateCoefficients(state,&problem->coefficients,point);
    }

    const double ionized = hydrogen[0], molecular = hydrogen[1];
    const double neutral = fmax(0.0,NeutralFraction(hydrogen));
    state->nel = nH*TrialIonFraction(trialHydrogen[0],ionized)+problem->siliconElectrons;
    if(oxygen) state->nel += state->n[ipOXYGEN]*TrialLadderCharge(y+oxygen,fractions+oxygen);
    if(helium) state->nel += state->n[ipHELIUM]*TrialLadderCharge(y+helium,fractions+helium);
    if(carbon) state->nel += state->n[ipCARBON]*TrialLadderCharge(y+carbon,fractions+carbon);
    if(!isfinite(state->nel)) return GSL_FAILURE;
    const int invalidRates = problem->evolveEnergy ? GSL_FAILURE : GSL_EBADFUNC;
    const double neutralStage = NeutralFraction(trialHydrogen);
    derivative[Ionized] = 0.0;
    if(!problem->chargeFree){
        const CEQHydrogenRates rates = CEQHydrogenIonizationRatesWithCoefficients(state,coefficients);
        if(!ValidRate(rates.ionization) || !ValidRate(rates.recombination)) return invalidRates;
        /* Signed stage populations restore the boundary after an overshoot. */
        derivative[Ionized] = rates.ionization*neutralStage - rates.recombination*trialHydrogen[0];
        /* The secondary-electron budget uses the physical stage populations. */
        if(deposition.secondaryHI != 0.0) derivative[Ionized] += deposition.secondaryHI/nH;
    }

    if(molecularOffset){
        const CEQMolecularHydrogenCoefficients *h2 = &coefficients->molecularHydrogen;
        const double photo = H2PhotoDissociationRateWithDustFactor(state->G0,(0.5*molecular)*state->Nt,
                                                                   h2->doppler,h2->dustShielding);
        const double collision = H2DissociationFromCoefficients(h2->collision,state->nel,nH*neutral,
                                                               (0.5*molecular)*nH);
        if(!ValidRate(h2->formation) || !ValidRate(photo) || !ValidRate(collision)) return invalidRates;
        derivative[molecularOffset] = 2.0*(h2->formation*neutralStage) - (photo+collision)*trialHydrogen[1];
    }
    if(oxygen && !OxygenDerivative(problem,y,neutralStage,derivative)) return GSL_FAILURE;
    if(helium && !LadderDerivative(CEQHeliumIonizationRates(state,coefficients),y+helium,derivative+helium))
        return GSL_FAILURE;
    if(helium && deposition.secondaryHeI != 0.0)
        derivative[helium] += deposition.secondaryHeI/state->n[ipHELIUM];
    if(carbon && !LadderDerivative(CEQCarbonIonizationRates(state,coefficients),y+carbon,derivative+carbon))
        return GSL_FAILURE;
    for(int i=0;i<problem->variables;i++)
        if(!isfinite(derivative[i])) return GSL_FAILURE;
    memcpy(problem->lastY,y,bytes);
    memcpy(problem->lastDerivative,derivative,bytes);
    problem->hasLastDerivative = true;
    return GSL_SUCCESS;
}

static int ChemistryJacobian(double time, const double y[], double matrix[], double dfdt[], void *parameters)
{
    const ThermochemistryProblem *problem = parameters;
    const int variables = problem->variables;
    const int molecular = problem->molecularOffset;
    const double hydrogen[] = {y[Ionized],molecular ? y[molecular] : 0.0};
    double base[MaxVariables];
    int status = ChemistryDerivative(time,y,base,parameters);
    if(status != GSL_SUCCESS) return status;

    for(int j=0;j<variables;j++){
        double nearby[MaxVariables];
        for(int i=0;i<variables;i++) nearby[i] = y[i];
        const bool energy = problem->evolveEnergy && j == LogEnergy;
        const double step = sqrt(DBL_EPSILON)*fmax(fabs(y[j]),energy ? 1.0 : 1.e-4);
        /* Approach the full-population boundary from the interior. */
        bool interior = false;
        if(!energy && y[j] > step){
            const double *pair = j == Ionized || j == molecular ? hydrogen : y+j-(j-problem->ladderStart)%2;
            interior = NeutralFraction(pair) < step;
        }
        nearby[j] += interior ? -step : step;
        double shifted[MaxVariables];
        status = ChemistryDerivative(time,nearby,shifted,parameters);
        if(status != GSL_SUCCESS) return status;
        for(int i=0;i<variables;i++){
            const double element = (shifted[i]-base[i])/(nearby[j]-y[j]);
            if(!isfinite(element)) return GSL_FAILURE;
            matrix[i*variables+j] = element;
        }
        dfdt[j] = 0.0;
    }
    return GSL_SUCCESS;
}

static bool InitializeLadder(double total, double singly, double doubly, double fractions[], double tolerance)
{
    if(!ValidRate(total) || !ValidRate(singly) || singly > total || !ValidRate(doubly) || doubly > total)
        return false;
    if(total == 0.0) return true;
    fractions[0] = singly/total;
    fractions[1] = doubly/total;
    return FractionsWithinBounds(fractions,tolerance);
}

typedef enum { FixedHydrogen, FixedHydrogenHelium, FixedHydrogenHeliumOxygen } FixedTemperatureNetwork;

static EQCLibStatus AdvanceChemistry(StructEQCData *state, double accuracy, bool evolveEnergy, FixedTemperatureNetwork network)
{
    if(!state || !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0 ||
       !isfinite(state->dt) || state->dt < 0.0)
        return EQC_INVALID_INPUT;
    if(evolveEnergy ? (!isfinite(state->U) || state->U <= 0.0 ||
                       !isfinite(state->Rho) || state->Rho <= 0.0 || !ValidRate(state->Z)) :
                     (!isfinite(state->T) || state->T <= 0.0))
        return EQC_INVALID_INPUT;
    if(evolveEnergy && state->dt > 0.0 && EQCLibRunParameters.UseCOFormation)
        return EQC_INVALID_INPUT;

    const double conditions[] = {state->G0,state->X0,state->NLyn,state->Nt,state->Zd};
    for(unsigned int i=0;i<sizeof conditions/sizeof conditions[0];i++)
        if(!isfinite(conditions[i]) || conditions[i] < 0.0) return EQC_INVALID_INPUT;
    for(int i=ipHYDROGEN+1;i<NumberChemicalElement;i++){
        if(evolveEnergy && (i == ipOXYGEN || i == ipHELIUM || i == ipCARBON || i == ipSILICON)) continue;
        if((i == ipHELIUM && network != FixedHydrogen) ||
           (i == ipOXYGEN && network == FixedHydrogenHeliumOxygen)) continue;
        if(state->n[i] != 0.0 || state->natom[i] != 0.0 || state->nion[i] != 0.0)
            return EQC_INVALID_INPUT;
    }
    if((!evolveEnergy && ((network == FixedHydrogen && state->nHeIII != 0.0) ||
                         (network != FixedHydrogenHeliumOxygen && state->nOIII != 0.0) ||
                         state->nCIII != 0.0 || state->nSiIII != 0.0)) ||
       state->nmol[ipCO] != 0.0) return EQC_INVALID_INPUT;

    const double nH = state->n[ipHYDROGEN];
    const double nHII = state->nion[ipHYDROGEN], nH2 = state->nmol[ipH2];
    if(!isfinite(nH) || nH <= 0.0 || !isfinite(nHII) || nHII < 0.0 || nHII > nH ||
       !isfinite(nH2) || nH2 < 0.0 || nH2 > 0.5*nH ||
       (!EQCLibRunParameters.UseH2Formation && nH2 != 0.0)) return EQC_INVALID_INPUT;

    const double nO = state->n[ipOXYGEN], nHe = state->n[ipHELIUM], nC = state->n[ipCARBON];
    const double nSi = state->n[ipSILICON], nSiII = state->nion[ipSILICON], nSiIII = state->nSiIII;
    const double roundoff = 8.0*DBL_EPSILON;
    double silicon[2] = {0};
    if(!InitializeLadder(nSi,nSiII,nSiIII,silicon,roundoff)) return EQC_INVALID_INPUT;
    const double siliconElectrons = nSiII+2.0*nSiIII;
    CEQUVBRates uvb;
    const EQCLibStatus uvbStatus = CEQEvaluateUVBackgroundWithMetallicity(state->NHIuvb,state->Z,&uvb);
    if(uvbStatus != EQC_SUCCESS) return uvbStatus;
    const bool hasUVBIonization = uvb.ionization[CEQ_UVB_HI] > 0.0 ||
                                 uvb.ionization[CEQ_UVB_HEI] > 0.0 ||
                                 uvb.ionization[CEQ_UVB_HEII] > 0.0 ||
                                 (nO > 0.0 && uvb.metalIonization[CEQ_UVB_OI] > 0.0) ||
                                 (nC > 0.0 && uvb.metalIonization[CEQ_UVB_CI] > 0.0);
    /* Without an ion source, zero charge is an invariant of this network.
       Keep finite-difference perturbations from seeding electron reactions. */
    const bool chargeFree = !hasUVBIonization && uvb.metalPowerHigh == 0.0 && state->NLyn == 0.0 && state->X0 == 0.0 &&
                            (nC == 0.0 || state->G0 == 0.0) && nHII == 0.0 &&
                            state->nion[ipOXYGEN] == 0.0 && state->nOIII == 0.0 &&
                            state->nion[ipHELIUM] == 0.0 && state->nHeIII == 0.0 &&
                            state->nion[ipCARBON] == 0.0 && state->nCIII == 0.0 &&
                            siliconElectrons == 0.0;
    double hydrogen[] = {nHII/nH,2.0*(nH2/nH)};
    double y[MaxVariables] = {0};
    const int ladderStart = evolveEnergy ? FirstLadder : 1;
    int variables = ladderStart;
    const int oxygen = !chargeFree && nO > 0.0 ? variables : 0;
    if(oxygen) variables += 2;
    const int helium = !chargeFree && nHe > 0.0 ? variables : 0;
    if(helium) variables += 2;
    const int carbon = !chargeFree && nC > 0.0 ? variables : 0;
    if(carbon) variables += 2;
    if(!FractionsWithinBounds(hydrogen,roundoff) || !isfinite(nH+2.0*(nO+nHe+nC+nSi)) || !isfinite(nO/nH) ||
       !InitializeLadder(nO,state->nion[ipOXYGEN],state->nOIII,y+oxygen,roundoff) ||
       !InitializeLadder(nHe,state->nion[ipHELIUM],state->nHeIII,y+helium,roundoff) ||
       !InitializeLadder(nC,state->nion[ipCARBON],state->nCIII,y+carbon,roundoff))
        return EQC_INVALID_INPUT;
    if(state->dt == 0.0) return EQC_SUCCESS;
    BoundFractions(hydrogen);
    y[Ionized] = hydrogen[0];
    for(int i=ladderStart;i<variables;i+=2) BoundFractions(y+i);
    const int ladderEnd = variables;
    const int molecular = EQCLibRunParameters.UseH2Formation ? variables++ : 0;
    if(molecular) y[molecular] = hydrogen[1];

    StructEQCData updated = *state;
    /* Silicon is a fixed trace species. Only its neutral density is derived. */
    if(evolveEnergy) updated.natom[ipSILICON] = nSi*fmax(0.0,NeutralFraction(silicon));
    ThermochemistryProblem problem = {.state = updated, .uvb = &uvb,
                                     .initialEnergy = state->U, .siliconElectrons = siliconElectrons,
                                     .variables = variables, .ladderStart = ladderStart,
                                     .evolveEnergy = evolveEnergy, .chargeFree = chargeFree, .molecularOffset = molecular,
                                     .hasFastPhotoelectrons = CEQHasUVBHighPower(&uvb,nHe),
                                     .oxygenOffset = oxygen, .heliumOffset = helium, .carbonOffset = carbon};
    CEQInitializeNonEquilibriumRateCoefficients(state,&problem.coefficients,&uvb);
    if(!evolveEnergy){
        CEQTemperaturePoint temperaturePoint;
        const CEQTemperaturePoint *point = CEQLocateTemperature(state->T,&temperaturePoint) ? &temperaturePoint : NULL;
        CEQRefreshNonEquilibriumRateCoefficients(state,&problem.coefficients,point);
    }
    double derivative[MaxVariables];
    if(ChemistryDerivative(0.0,y,derivative,&problem) != GSL_SUCCESS) return EQC_INVALID_INPUT;

    const double absoluteTolerance = fmax(DBL_MIN,1.e-3*accuracy);
    double speed = 0.0;
    for(int i=0;i<variables;i++) speed = fmax(speed,fabs(derivative[i]));
    if(speed > 0.0){
        const gsl_odeiv2_system system = {ChemistryDerivative,ChemistryJacobian,variables,&problem};
        const double firstStep = fmin(state->dt,1.e-3/speed);
        gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&system,gsl_odeiv2_step_msbdf,firstStep,
                                                                 absoluteTolerance,accuracy);
        if(!driver) return EQC_NO_CONVERGENCE;
        gsl_odeiv2_driver_set_nmax(driver,100000);
        double elapsed = 0.0;
        const int status = gsl_odeiv2_driver_apply(driver,&elapsed,state->dt,y);
        gsl_odeiv2_driver_free(driver);
        if(status != GSL_SUCCESS) return EQC_NO_CONVERGENCE;
    }

    hydrogen[0] = chargeFree ? 0.0 : y[Ionized];
    hydrogen[1] = molecular ? y[molecular] : 0.0;
    const double boundaryTolerance = roundoff+accuracy;
    if(!FractionsWithinBounds(hydrogen,boundaryTolerance)) return EQC_NO_CONVERGENCE;
    BoundFractions(hydrogen);
    for(int i=ladderStart;i<ladderEnd;i+=2){
        if(!FractionsWithinBounds(y+i,boundaryTolerance)) return EQC_NO_CONVERGENCE;
        BoundFractions(y+i);
    }
    SetAbundances(&updated,hydrogen,y,&problem);
    if(evolveEnergy && (!SetThermalState(&updated,state->U*exp(y[LogEnergy])) ||
                       CEQCalculateNonEquilibriumCoolingHeatingRatesWithUVB(&updated,&uvb) != EQC_SUCCESS))
        return EQC_NO_CONVERGENCE;
    *state = updated;
    return EQC_SUCCESS;
}

EQCLibStatus EQCLibTryAdvanceHydrogenChemistry(StructEQCData *state, double accuracy)
{
    return AdvanceChemistry(state,accuracy,false,FixedHydrogen);
}

EQCLibStatus EQCLibTryAdvanceHydrogenThermochemistry(StructEQCData *state, double accuracy)
{
    return AdvanceChemistry(state,accuracy,true,FixedHydrogen);
}

EQCLibStatus CEQAdvanceFixedTemperatureHydrogenHeliumChemistry(StructEQCData *state, double accuracy)
{
    return AdvanceChemistry(state,accuracy,false,FixedHydrogenHelium);
}

EQCLibStatus CEQAdvanceFixedTemperatureHHeOChemistry(StructEQCData *state, double accuracy)
{
    return AdvanceChemistry(state,accuracy,false,FixedHydrogenHeliumOxygen);
}
