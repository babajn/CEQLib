#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "CosmicRay.h"
#include "EquilibriumAbundance.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#include "PhotoelectronDeposition.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "Xray.h"

/* Wolfire et al. (1995), Appendix A, gives the accuracy of the soft X-ray
   ionization fit down to xe = 1e-4. Keep local probes at or above that limit. */
static const double XrayFitMinimumElectronFraction = 1.0e-4;

typedef struct {
    double density;
    double residual;
    double charge;
    int status;
} ChargePoint;

typedef struct {
    const StructEQCData *initial;
    double previousElectronDensity;
    double previousCO;
    double accuracy;
    const CEQCarbonShieldingInterval *carbonIntervals;
    const CEQRateCoefficients *coefficients;
    const CEQUVBRates *uvb;
    CEQUVBEquilibriumStats *uvbStats;
    int carbonIntervalCount;
    StructEQCData solution;
    CEQSecondaryIonizationSeed solutionSecondary;
    bool found;
    double seedH2;   /* H2 density of the latest successful trial, to seed the next. */
    double seedSecondaryHI, seedSecondaryHeI;
    bool hasUVBSecondaries;
    bool hasXrayIonization;
    bool hasMolecularChargeCoupling;
    bool oneWayChargeTransfer;
    unsigned int solutions;
    int failure;
} ChargeSearch;

static bool UsesTwoWayChargeTransfer(void)
{
    return EQCLibRunParameters.UseHIIRegion && EQCLibRunParameters.UseChargeTransfer;
}

/* In UVB-only mode, the H solution provides the collision partners for O. */
static int SolveOxygen(StructEQCData *state, const CEQHydrogenRates *hydrogen,
                       const CEQIonLadderRates *oxygen, const ChargeSearch *search)
{
    const CEQRateCoefficients *c = search->coefficients;
    if(UsesTwoWayChargeTransfer())
        return CEQSolveHydrogenOxygenChargeTransfer(state,hydrogen,oxygen,c,search->accuracy);
    CEQIonLadderRates rates = *oxygen;
    if(search->oneWayChargeTransfer){
        const double nHI = state->natom[ipHYDROGEN];
        rates.ionization0 += state->nion[ipHYDROGEN]*c->toOxygenII;
        rates.recombination1 += nHI*c->toHydrogenI;
        rates.recombination2 += nHI*c->oxygenIIIToHydrogenI;
    }
    return CEQSolveIonLadder(state->n[ipOXYGEN]-state->nmol[ipCO],
                             rates.ionization0,rates.recombination1,rates.ionization1,rates.recombination2,
                             &state->natom[ipOXYGEN],&state->nion[ipOXYGEN],&state->nOIII);
}

static int SolveHelium(StructEQCData *state, const CEQIonLadderRates *helium)
{
    return CEQSolveIonLadder(state->n[ipHELIUM],
        helium->ionization0, helium->recombination1, helium->ionization1, helium->recombination2,
        &state->natom[ipHELIUM], &state->nion[ipHELIUM], &state->nHeIII);
}

static int SolveOxygenAndSilicon(StructEQCData *state, const CEQHydrogenRates *hydrogen,
                                CEQIonLadderRates *oxygen, const ChargeSearch *search)
{
    const CEQRateCoefficients *c = search->coefficients;
    const double ne = state->nel;
    const double nH = state->n[ipHYDROGEN];
    const double xe = ne/nH;
    const double T = state->T;
    const double X0 = state->X0;
    const double G0 = state->G0;
    const double Nt = state->Nt;
    const double Zd = state->Zd;

    *oxygen = CEQOxygenIonizationRates(state,c);
    const int status = SolveOxygen(state,hydrogen,oxygen,search);
    if(status != EQC_SUCCESS) return status;

    return CEQSolveIonLadder(state->n[ipSILICON],
        c->silicon.photo0 + ne*c->silicon.collision0
            + X0*(CosmicRayIonizationRateSiI(xe) + SoftXrayIonizationRateSiI(xe,Nt)),
        ne*c->silicon.recombination1 + nH*SiIIRecombinationRateOnDust(T,G0,ne,Zd),
        c->silicon.photo1 + ne*c->silicon.collision1
            + X0*(CosmicRayIonizationRateSiII(xe) + SoftXrayIonizationRateSiII(xe,Nt)),
        ne*c->silicon.recombination2,
        &state->natom[ipSILICON], &state->nion[ipSILICON], &state->nSiIII);
}

static int SolveIonsAtFixedCharge(StructEQCData *state, const ChargeSearch *search,
                                  const CEQHydrogenRates *adjustedHydrogen, const CEQIonLadderRates *adjustedHelium)
{
    const CEQRateCoefficients *coefficients = search->coefficients;
    const CEQHydrogenRates hydrogen = adjustedHydrogen ?
        *adjustedHydrogen : CEQHydrogenIonizationRatesWithCoefficients(state,coefficients);
    int status = CEQSolveHydrogenEquilibriumWithRates(state,&hydrogen,coefficients,search->accuracy);
    if(status == EQC_SUCCESS)
        status = CEQSolveCarbonWithIntervals(state,search->previousCO,search->accuracy,
                                             search->carbonIntervals,search->carbonIntervalCount,
                                             coefficients,false);
    double carbonH2 = state->nmol[ipH2];
    CEQIonLadderRates oxygen;
    if(status == EQC_SUCCESS){
        const CEQIonLadderRates helium = adjustedHelium ?
            *adjustedHelium : CEQHeliumIonizationRates(state,coefficients);
        status = SolveHelium(state,&helium);
        if(status == EQC_SUCCESS)
            status = SolveOxygenAndSilicon(state,&hydrogen,&oxygen,search);
    }
    /* Check the CO fluxes whenever charge transfer changes H2. Even a small
       change can cross a shielding table boundary. Update H and O again only
       when the CO solve changes the oxygen available for charge transfer. */
    if(status == EQC_SUCCESS && UsesTwoWayChargeTransfer() && EQCLibRunParameters.UseH2Formation &&
       EQCLibRunParameters.UseCOFormation && state->T <= 1.e4 && state->n[ipCARBON] > 0.0){
        for(int pass = 0; state->nmol[ipH2] != carbonH2; pass++){
            if(pass == 16){
                status = EQC_NO_CONVERGENCE;
                break;
            }
            carbonH2 = state->nmol[ipH2];
            const double previousCO = state->nmol[ipCO];
            status = CEQSolveCarbonWithIntervals(state,previousCO,search->accuracy,
                                                 search->carbonIntervals,search->carbonIntervalCount,
                                                 coefficients,true);
            if(status != EQC_SUCCESS || state->nmol[ipCO] == previousCO) break;
            status = SolveOxygen(state,&hydrogen,&oxygen,search);
            if(status != EQC_SUCCESS) break;
        }
    }
    return status;
}

static double SecondaryRateResidual(double previous, double next, double primary)
{
    const double scale = primary + fmax(previous,next);
    return scale > 0.0 ? fabs(next-previous)/scale : 0.0;
}

static int SolveIonsWithUVBSecondaries(StructEQCData *state, ChargeSearch *search)
{
    const CEQRateCoefficients *primary = search->coefficients;
    const CEQFixedChargeLadder hydrogenLadder = CEQHydrogenLadderAtFixedCharge(state,primary);
    const CEQFixedChargeLadder heliumLadder = CEQHeliumLadderAtFixedCharge(state,primary);
    const double hydrogen = CEQEvaluateFixedChargeLadder(hydrogenLadder,primary->hydrogen.photo0).ionization0;
    const double helium = CEQEvaluateFixedChargeLadder(heliumLadder,primary->helium.photo0).ionization0;
    /* One-way exchange lets H and He converge before solving O and CO. Keep
       the existing CO search order in modes without one-way exchange. */
    const bool deferMetals = !UsesTwoWayChargeTransfer() &&
        (search->oneWayChargeTransfer || !EQCLibRunParameters.UseCOFormation ||
         state->T > 1.e4 || state->n[ipCARBON] == 0.0);
    double secondaryHI = search->seedSecondaryHI, secondaryHeI = search->seedSecondaryHeI;
    const double tolerance = fmax(4.0*DBL_EPSILON,0.25*search->accuracy);
    double damping = 1.0, previousResidual = INFINITY;

    /* At fixed ne, a trial pair of secondary rates permits the existing H/He
       and molecular balances. Recompute the photon absorption, the H ion
       fraction and the available neutral targets from the resulting state. */
    for(int iteration = 0; iteration < 128; iteration++){
        if(search->uvbStats) search->uvbStats->iterations++;
        const CEQIonLadderRates h = CEQEvaluateFixedChargeLadder(hydrogenLadder,primary->hydrogen.photo0+secondaryHI);
        const CEQHydrogenRates hydrogenRates = {.ionization = h.ionization0, .recombination = h.recombination1};
        const CEQIonLadderRates heliumRates = CEQEvaluateFixedChargeLadder(heliumLadder,primary->helium.photo0+secondaryHeI);
        int status;
        if(deferMetals){
            status = CEQSolveHydrogenEquilibriumWithRates(state,&hydrogenRates,primary,search->accuracy);
            if(status == EQC_SUCCESS) status = SolveHelium(state,&heliumRates);
        }else{
            status = SolveIonsAtFixedCharge(state,search,&hydrogenRates,&heliumRates);
        }
        if(status != EQC_SUCCESS) return status;

        const double power = CEQHighPhotoelectronPower(state,search->uvb);
        const CEQPhotoelectronIonization rates = CEQPhotoelectronIonizationRates(state,power);
        const double nextHeI = state->n[ipHELIUM] > 0.0 ? rates.helium : 0.0;
        if(!isfinite(rates.hydrogen) || !isfinite(nextHeI))
            return EQC_NO_CONVERGENCE;
        const double residual = fmax(SecondaryRateResidual(secondaryHI,rates.hydrogen,hydrogen),
                                     SecondaryRateResidual(secondaryHeI,nextHeI,helium));
        if(residual <= tolerance){
            if(power == 0.0 && (secondaryHI > 0.0 || secondaryHeI > 0.0)){
                /* An exactly vanishing source must recover the primary-only
                   state, even when the incoming rates fit within tolerance. */
                secondaryHI = secondaryHeI = 0.0;
                state->nmol[ipH2] = search->seedH2;
                continue;
            }
            if(deferMetals){
                status = CEQSolveCarbonWithIntervals(state,search->previousCO,search->accuracy,
                                                      search->carbonIntervals,search->carbonIntervalCount,
                                                      primary,false);
                if(status != EQC_SUCCESS) return status;
                CEQIonLadderRates oxygen;
                status = SolveOxygenAndSilicon(state,&hydrogenRates,&oxygen,search);
                if(status != EQC_SUCCESS) return status;
            }
            search->seedSecondaryHI = rates.hydrogen;
            search->seedSecondaryHeI = nextHeI;
            return EQC_SUCCESS;
        }

        /* Damp stalled or alternating iterates while keeping both rates
           nonnegative. The accepted molecular branch seeds the next pass. */
        if(residual >= 0.8*previousResidual) damping = fmax(0.25,0.5*damping);
        secondaryHI = (1.0-damping)*secondaryHI + damping*rates.hydrogen;
        secondaryHeI = (1.0-damping)*secondaryHeI + damping*nextHeI;
        previousResidual = residual;
    }
    return EQC_NO_CONVERGENCE;
}

static ChargePoint EvaluateCharge(ChargeSearch *search, double electronDensity)
{
    if(search->uvbStats) search->uvbStats->chargeTrials++;
    StructEQCData state = *search->initial;
    state.nel = electronDensity;
    /* Continue the H2 branch from the latest successful charge trial. */
    state.nmol[ipH2] = search->seedH2;
    const int status = search->hasUVBSecondaries ? SolveIonsWithUVBSecondaries(&state,search) :
                                                  SolveIonsAtFixedCharge(&state,search,NULL,NULL);
    if(status == EQC_SUCCESS) search->seedH2 = state.nmol[ipH2];
    ChargePoint point = {.density = electronDensity, .status = status};
    if(status != EQC_SUCCESS){
        if(status == EQC_NO_CONVERGENCE || search->failure == EQC_NO_EQUILIBRIUM)
            search->failure = status;
        return point;
    }

    const double charge = EQCLibGetElectronNumberDensity(&state);
    if(!isfinite(charge) || charge < 0.0){
        point.status = search->failure = EQC_NO_CONVERGENCE;
        return point;
    }
    double scale = fmax(charge,electronDensity);
    point.charge = charge;
    point.residual = scale == 0.0 ? 0.0 : (charge-electronDensity)/scale;

    /* The soft X-ray fit tends to zero as xe approaches zero in the neutral
       gas regime. Its formal zero should not select a neutral branch while
       H I or He I are exposed to an active X-ray field. */
    bool neutralXrayGas = electronDensity == 0.0 && search->hasXrayIonization &&
        (state.natom[ipHYDROGEN] > 0.0 || state.natom[ipHELIUM] > 0.0);
    if(!neutralXrayGas && fabs(point.residual) <= search->accuracy){
        double distance = fabs(electronDensity-search->previousElectronDensity);
        double bestDistance = fabs(search->solution.nel-search->previousElectronDensity);
        bool preferCharged = search->previousElectronDensity > 0.0 &&
            (electronDensity > 0.0) != (search->solution.nel > 0.0);
        if(!search->found || (preferCharged ? electronDensity > 0.0 :
                distance < bestDistance ||
                (distance == bestDistance && electronDensity < search->solution.nel))){
            search->solution = state;
            search->solutionSecondary = (CEQSecondaryIonizationSeed){search->seedSecondaryHI,search->seedSecondaryHeI};
            search->found = true;
        }
        search->solutions++;
    }
    return point;
}

static void SearchPhotoionizedCharge(ChargeSearch *search, ChargePoint hint, double maximumCharge)
{
    /* The ion charge at a nearly neutral seed gives a much closer first
       estimate than expanding the search by small factors around that seed. */
    const int failure = search->failure;
    double density = hint.charge;
    ChargePoint previous = hint;
    for(int i = 0; i < 6 && density > 0.0 && density <= maximumCharge && density != previous.density; i++){
        const ChargePoint next = EvaluateCharge(search,density);
        if(search->found) return;
        if(next.status != EQC_SUCCESS) break;

        /* Use a constant scale for the secant so it sees the slope of q-ne.
           The convergence check still uses the relative charge residual. */
        const double imbalance = (next.charge-next.density)/maximumCharge;
        const double previousImbalance = (previous.charge-previous.density)/maximumCharge;
        const double difference = imbalance-previousImbalance;
        const double estimate = difference == 0.0 ? next.charge :
            next.density - (imbalance/difference)*(next.density-previous.density);
        density = estimate > 0.0 && estimate <= maximumCharge ? estimate : next.charge;
        previous = next;
    }
    search->failure = failure;
}

typedef enum {
    LocalChargeMiss,
    LocalChargeRoot,
    LocalChargeBracket
} LocalChargeResult;

/* Place new probes within [minimumDensity, maximumDensity). */
static LocalChargeResult ExpandChargeBracket(ChargeSearch *search, ChargePoint hint,
                                             double minimumDensity, double maximumDensity,
                                             ChargePoint *lower, ChargePoint *upper)
{
    const bool increase = hint.residual > 0.0;
    const int failure = search->failure;
    ChargePoint previous = hint;
    double factor = 1.05;
    for(int step=0; step<40; step++){
        const double density = increase ? previous.density * factor : previous.density / factor;
        if(density <= 0.0 || density < minimumDensity || density >= maximumDensity ||
                density == previous.density)
            break;

        const ChargePoint next = EvaluateCharge(search,density);
        if(next.status != EQC_SUCCESS)
            break;
        if(fabs(next.residual) <= search->accuracy)
            return LocalChargeRoot;
        if((next.residual < 0.0) != (previous.residual < 0.0)){
            *lower = increase ? previous : next;
            *upper = increase ? next : previous;
            return LocalChargeBracket;
        }
        previous = next;
        factor = fmin(10.0,factor*factor);
    }
    search->failure = failure;
    return LocalChargeMiss;
}

static bool NarrowChargeBracket(ChargeSearch *search, ChargePoint *lower, ChargePoint *upper)
{
    /* For independent atomic ladders in this regime, charge divided by ne
       decreases with ne. The new O rates also depend on the H ion fraction,
       so use the original bracket when one-way transfer is active. */
    if(search->oneWayChargeTransfer || search->hasMolecularChargeCoupling || search->hasUVBSecondaries ||
            (search->hasXrayIonization &&
             lower->density/search->initial->n[ipHYDROGEN] <= EQC_XRAY_IONIZATION_MAX_XE))
        return false;

    const double seed = search->previousElectronDensity;
    const bool increase = lower->density == seed && lower->residual > 0.0;
    if(seed <= 0.0 || (!increase && (upper->density != seed || upper->residual >= 0.0)))
        return false;

    const ChargePoint hint = increase ? *lower : *upper;
    return ExpandChargeBracket(search,hint,nextafter(lower->density,upper->density),
                               upper->density,lower,upper) == LocalChargeRoot;
}

static void BisectCharge(ChargeSearch *search, ChargePoint lower, ChargePoint upper,
                        bool interpolate)
{
    if(lower.status != EQC_SUCCESS || upper.status != EQC_SUCCESS ||
            fabs(lower.residual) <= search->accuracy || fabs(upper.residual) <= search->accuracy ||
            (lower.residual < 0.0) == (upper.residual < 0.0))
        return;

    if(search->hasXrayIonization &&
            lower.density < XrayFitMinimumElectronFraction*search->initial->n[ipHYDROGEN])
        interpolate = false;
    if(interpolate && NarrowChargeBracket(search,&lower,&upper)) return;

    for(int iteration=0; iteration<160; iteration++){
        double middle;
        if(lower.density > 0.0 && lower.density < 0.25*upper.density){
            middle = sqrt(lower.density)*sqrt(upper.density);
        }else{
            middle = lower.density + 0.5*(upper.density-lower.density);
        }
        if(middle <= lower.density || middle >= upper.density) break;

        double density = middle;
        if(interpolate && iteration % 2 == 0){
            const double fraction = lower.residual/(lower.residual-upper.residual);
            double estimate = lower.density + fraction*(upper.density-lower.density);
            if(lower.density > 0.0){
                const double logLower = log(lower.density);
                estimate = exp(logLower + fraction*(log(upper.density)-logLower));
            }
            if(estimate > lower.density && estimate < upper.density) density = estimate;
        }

        const int failure = search->failure;
        ChargePoint trial = EvaluateCharge(search,density);
        if(trial.status != EQC_SUCCESS && density != middle){
            search->failure = failure;
            trial = EvaluateCharge(search,middle);
        }
        if(trial.status != EQC_SUCCESS) return;
        if(fabs(trial.residual) <= search->accuracy) return;
        if((trial.residual < 0.0) == (lower.residual < 0.0)) lower = trial;
        else upper = trial;
    }
    search->failure = EQC_NO_CONVERGENCE;
}

static void SearchChargeAroundHint(ChargeSearch *search, ChargePoint hint,
                                   double minimumDensity, double maximumDensity)
{
    if(hint.status != EQC_SUCCESS || hint.density <= 0.0 || hint.density < minimumDensity)
        return;

    ChargePoint lower, upper;
    const int failure = search->failure;
    if(ExpandChargeBracket(search,hint,minimumDensity,maximumDensity,&lower,&upper) == LocalChargeBracket)
        BisectCharge(search,lower,upper,true);
    if(!search->found) search->failure = failure;
}

static void SearchChargeInterval(ChargeSearch *search, double lowerDensity, double upperDensity,
                                 ChargePoint hint)
{
    unsigned int initialSolutions = search->solutions;
    ChargePoint lower = EvaluateCharge(search,lowerDensity);
    ChargePoint upper = EvaluateCharge(search,upperDensity);
    double xrayRangeStart = XrayFitMinimumElectronFraction*search->initial->n[ipHYDROGEN];
    bool probeBelowXrayRange = lowerDensity == 0.0 && search->hasXrayIonization &&
        hint.density < xrayRangeStart && lower.status == EQC_SUCCESS && lower.residual == 0.0;

    if(lowerDensity == 0.0){
        /* Test the positive side separately; pure collisional gas can also have
           an exactly neutral solution at zero electron density. */
        double positiveDensity = sqrt(DBL_MIN)*sqrt(upperDensity);
        // Probe the lower end of the fitted X-ray electron fraction range.
        if(search->hasXrayIonization)
            positiveDensity = fmin(xrayRangeStart,0.5*upperDensity);
        if(positiveDensity > lowerDensity && positiveDensity < upperDensity){
            ChargePoint positive = EvaluateCharge(search,positiveDensity);
            if(positive.status == EQC_SUCCESS){
                /* Keep the neutral-limit probes on the existing bisection grid. */
                if(hint.status == EQC_SUCCESS && hint.density > lower.density && hint.density < positive.density){
                    BisectCharge(search,lower,hint,false);
                    BisectCharge(search,hint,positive,false);
                }else{
                    BisectCharge(search,lower,positive,false);
                }
                lower = positive;
            }
        }
    }

    unsigned int beforeBisection = search->solutions;
    if(hint.status == EQC_SUCCESS && hint.density > lower.density && hint.density < upper.density){
        BisectCharge(search,lower,hint,true);
        BisectCharge(search,hint,upper,true);
    }else{
        BisectCharge(search,lower,upper,true);
    }

    if(probeBelowXrayRange && search->solutions == initialSolutions){
        ChargePoint previous = lower;
        for(int i=0;i<64 && search->solutions == initialSolutions;i++){
            double density = 0.5*previous.density;
            if(density == 0.0) break;
            ChargePoint current = EvaluateCharge(search,density);
            BisectCharge(search,current,previous,false);
            previous = current;
        }
    }

    /* Molecular shielding can remove a CO root over part of the interval.
       Probe neighboring valid subintervals without bracketing across failures. */
    if(EQCLibRunParameters.UseCOFormation && search->initial->T <= 1.0e4 &&
            search->solutions == beforeBisection &&
            (search->solutions == initialSolutions || lower.status != EQC_SUCCESS || upper.status != EQC_SUCCESS)){
        ChargePoint previous = lower;
        const int subdivisions = 32;
        for(int i=1;i<=subdivisions;i++){
            ChargePoint current = i == subdivisions ? upper :
                EvaluateCharge(search,lower.density + (upper.density-lower.density)*i/subdivisions);
            BisectCharge(search,previous,current,true);
            previous = current;
        }
        if(search->solutions == beforeBisection){
            /* Linear samples cover high electron fractions; logarithmic samples
               also locate admissible intervals in almost neutral gas. */
            previous = upper;
            for(int i=0;i<64;i++){
                double density = 0.5*previous.density;
                if(density <= lower.density) break;
                ChargePoint current = EvaluateCharge(search,density);
                BisectCharge(search,current,previous,true);
                previous = current;
            }
        }
    }
}

static bool ConservesElements(const StructEQCData *state)
{
    const double higher[NumberChemicalElement] = {
        0.0, state->nHeIII, state->nCIII, state->nOIII, state->nSiIII
    };
    const double molecules[NumberChemicalElement] = {
        2.0*state->nmol[ipH2], 0.0, state->nmol[ipCO], state->nmol[ipCO], 0.0
    };
    for(int i=0;i<NumberChemicalElement;i++){
        double sum = state->natom[i]+state->nion[i]+higher[i]+molecules[i];
        if(!isfinite(sum) || fabs(sum-state->n[i]) > 32.0*DBL_EPSILON*state->n[i])
            return false;
    }
    return true;
}

int EQCLibChemicalReactionEquilibrium(StructEQCData *state, double accuracy)
{
    if(state == NULL) return EQC_INVALID_INPUT;
    CEQUVBRates uvb;
    const EQCLibStatus status = CEQEvaluateUVBackgroundWithMetallicity(state->NHIuvb,state->Z,&uvb);
    if(status != EQC_SUCCESS) return status;
    return CEQChemicalReactionEquilibriumWithUVB(state,accuracy,&uvb);
}

int CEQChemicalReactionEquilibriumWithUVB(StructEQCData *state, double accuracy, const CEQUVBRates *uvb)
{
    return CEQChemicalReactionEquilibriumWithUVBStats(state,accuracy,uvb,NULL);
}

int CEQChemicalReactionEquilibriumWithUVBStats(StructEQCData *state, double accuracy, const CEQUVBRates *uvb,
                                                CEQUVBEquilibriumStats *stats)
{
    return CEQChemicalReactionEquilibriumWithUVBSeed(state,accuracy,uvb,stats,NULL);
}

int CEQChemicalReactionEquilibriumWithUVBSeed(StructEQCData *state, double accuracy, const CEQUVBRates *uvb,
                                               CEQUVBEquilibriumStats *stats, CEQSecondaryIonizationSeed *seed)
{
    if(!state || !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0 ||
            !isfinite(state->T) || state->T <= 0.0) return EQC_INVALID_INPUT;

    const double conditions[] = {state->G0, state->X0, state->NLyn, state->Nt, state->Z, state->Zd};
    for(unsigned int i=0;i<sizeof(conditions)/sizeof(conditions[0]);i++){
        if(!isfinite(conditions[i]) || conditions[i] < 0.0) return EQC_INVALID_INPUT;
    }
    const double nH = state->n[ipHYDROGEN];
    if(nH <= 0.0) return EQC_INVALID_INPUT;
    double maximumCharge = nH;
    for(int i=0;i<NumberChemicalElement;i++){
        if(!isfinite(state->n[i]) || state->n[i] < 0.0) return EQC_INVALID_INPUT;
        if(i != ipHYDROGEN) maximumCharge += 2.0*state->n[i];
    }
    if(!isfinite(maximumCharge)) return EQC_INVALID_INPUT;

    const int shieldingPoints = EQCLibRunParameters.UseCOFormation && state->T <= 1.0e4 ?
                                CEQCOShieldingPointCount() : 0;
    CEQCarbonShieldingInterval carbonIntervals[shieldingPoints+1];
    const int carbonIntervalCount = CEQPrepareCarbonShieldingIntervals(state,carbonIntervals,shieldingPoints);
    /* T, stellar flux, dust and the attenuated UVB rates are fixed during the search. */
    CEQRateCoefficients coefficients;
    CEQPrepareRateCoefficients(state,&coefficients,uvb);
    ChargeSearch search = {
        .initial = state,
        .coefficients = &coefficients,
        .uvb = uvb,
        .uvbStats = stats,
        .seedH2 = state->nmol[ipH2],
        .previousElectronDensity = isfinite(state->nel) ? fmin(fmax(state->nel,0.0),maximumCharge) : 0.0,
        .previousCO = EQCLibRunParameters.UseCOFormation ? state->nmol[ipCO] : 0.0,
        .accuracy = accuracy,
        .carbonIntervals = carbonIntervals,
        .carbonIntervalCount = carbonIntervalCount,
        .hasUVBSecondaries = CEQHasUVBHighPower(uvb,state->n[ipHELIUM]),
        .hasXrayIonization = EQCLibRunParameters.UseSoftXrayIonization && state->X0 > 0.0,
        .hasMolecularChargeCoupling = EQCLibRunParameters.UseH2Formation ||
                                     (EQCLibRunParameters.UseCOFormation && state->T <= 1.0e4),
        .oneWayChargeTransfer = !EQCLibRunParameters.UseHIIRegion && EQCLibRunParameters.UseChargeTransfer &&
                                state->n[ipOXYGEN] > 0.0 && CEQUVBackgroundIsEnabled(),
        .failure = EQC_NO_EQUILIBRIUM
    };
    if(search.hasUVBSecondaries && seed && isfinite(seed->hydrogen) && seed->hydrogen >= 0.0 &&
       isfinite(seed->helium) && seed->helium >= 0.0){
        search.seedSecondaryHI = seed->hydrogen;
        search.seedSecondaryHeI = state->n[ipHELIUM] > 0.0 ? seed->helium : 0.0;
    }
    ChargePoint hint = EvaluateCharge(&search,search.previousElectronDensity);
    /* With only photoionization, electron collisions and electron
       recombination, every ionization/recombination ratio is of the form
       A + B/ne at fixed T. The ion charge is nonincreasing, giving a unique
       positive root. Leave branch selection to the original search when
       molecular coupling, secondary UVB ionization, grain recombination,
       cosmic rays or X-rays apply.
       H-O charge transfer conserves charge and couples the H and O ladders
       only through O/H. The argument above does not cover that coupling, so
       HIIRegionTest compares it with solves from an ionized seed. */
    if(!search.found && hint.status == EQC_SUCCESS && hint.residual > 0.5 &&
            EQCLibRunParameters.UseHIIRegion && state->NLyn > 0.0 &&
            !search.hasMolecularChargeCoupling && !search.hasUVBSecondaries && state->X0 == 0.0 &&
            (!EQCLibRunParameters.UseGrainSurfaceProcesses || state->Zd == 0.0))
        SearchPhotoionizedCharge(&search,hint,maximumCharge);
    if(!search.found){
        double rangeStart = 0.0;
        if(search.hasXrayIonization){
            /* Leave the extrapolation below the fit range to the neutral-limit probes. */
            rangeStart = XrayFitMinimumElectronFraction*nH;
            if(rangeStart/nH < XrayFitMinimumElectronFraction)
                rangeStart = nextafter(rangeStart,INFINITY);
        }
#ifdef EQC_LEGACY_XRAY_CUTOFF
        const bool preferXrayIntervals = search.hasXrayIonization;
#else
        /* Carrying H2 between probes can change the molecular branch if the
           probe order changes. Try the usual intervals for molecular gas
           before searching the full, continuous charge range. */
        const bool preferXrayIntervals = search.hasXrayIonization && search.hasMolecularChargeCoupling;
#endif
        if(preferXrayIntervals){
            double boundary = EQC_XRAY_IONIZATION_MAX_XE*nH;
            /* The rate uses ne/nH, so begin the second segment at the first
               representable density whose computed fraction exceeds 0.5. */
            while(boundary/nH > EQC_XRAY_IONIZATION_MAX_XE) boundary = nextafter(boundary,0.0);
            double afterBoundary = nextafter(boundary,INFINITY);
            while(afterBoundary/nH <= EQC_XRAY_IONIZATION_MAX_XE) afterBoundary = nextafter(afterBoundary,INFINITY);
            if(hint.density <= boundary){
                SearchChargeAroundHint(&search,hint,rangeStart,afterBoundary);
            }else{
                SearchChargeAroundHint(&search,hint,afterBoundary,maximumCharge);
            }
            if(!search.found){
                SearchChargeInterval(&search,0.0,boundary,hint);
                if(afterBoundary <= maximumCharge)
                    SearchChargeInterval(&search,afterBoundary,maximumCharge,hint);
            }
        }
#ifdef EQC_LEGACY_XRAY_CUTOFF
        if(!preferXrayIntervals)
#else
        if(!search.found)
#endif
        {
            SearchChargeAroundHint(&search,hint,rangeStart,maximumCharge);
            if(!search.found) SearchChargeInterval(&search,0.0,maximumCharge,hint);
        }
    }

    if(!search.found) return search.failure;
    StructEQCData *solution = &search.solution;
    if(!ConservesElements(solution)) return EQC_NO_CONVERGENCE;
    solution->MeanMolecularWeight = EQCLibGetMeanMolecularWeight(solution);
    solution->DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(solution);
    if(!isfinite(solution->MeanMolecularWeight) || solution->MeanMolecularWeight <= 0.0 ||
            !isfinite(solution->DegreeOfFreedom) || solution->DegreeOfFreedom <= 0.0)
        return EQC_NO_CONVERGENCE;
    solution->Gamma = (1.0+solution->DegreeOfFreedom)/solution->DegreeOfFreedom;
    solution->Gm1 = solution->Gamma-1.0;
    solution->GGm1 = solution->Gamma*solution->Gm1;
    if(seed) *seed = search.solutionSecondary;
    *state = *solution;
    return EQC_SUCCESS;
}
