#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "CosmicRay.h"
#include "EquilibriumAbundance.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
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
    const StructCEQData *initial;
    double previousElectronDensity;
    double previousCO;
    double accuracy;
    const CEQCarbonShieldingInterval *carbonIntervals;
    int carbonIntervalCount;
    StructCEQData solution;
    bool found;
    bool hasXrayIonization;
    bool hasMolecularChargeCoupling;
    unsigned int solutions;
    int failure;
} ChargeSearch;

static int SolveOtherIons(StructCEQData *state, const CEQHydrogenRates *hydrogen, double accuracy)
{
    const double ne = state->nel;
    const double nH = state->n[ipHYDROGEN];
    const double xe = ne/nH;
    const double T = state->T;
    const double NLyn = state->NLyn;
    const double X0 = state->X0;
    const double G0 = state->G0;
    const double Nt = state->Nt;
    const double Zd = state->Zd;

    int status = CEQSolveIonLadder(state->n[ipHELIUM],
        HeIPhotoIonizationRate(NLyn) + ne*HeICollisionalIonizationRate(T)
            + X0*(CosmicRayIonizationRateHeI(xe) + SoftXrayIonizationRateHeI(xe,Nt)),
        ne*(HeIIRecombinationRate(T) + HeIIDielectronicRecombinationRate(T))
            + nH*HeIIRecombinationRateOnDust(T,G0,ne,Zd),
        HeIIPhotoIonizationRate(NLyn) + ne*HeIICollisionalIonizationRate(T),
        ne*HeIIIRecombinationRate(T),
        &state->natom[ipHELIUM], &state->nion[ipHELIUM], &state->nHeIII);
    if(status != CEQ_SUCCESS) return status;

    const CEQIonLadderRates oxygen = {
        .ionization0 = OIPhotoionizationRate(NLyn) + ne*OICollisionalIonizationRate(T)
                     + X0*CosmicRayIonizationRateOI(xe),
        .recombination1 = ne*OIIRecombinationRate(T) + nH*OIIRecombinationRateOnDust(T,G0,ne,Zd),
        .ionization1 = OIIPhotoionizationRate(NLyn) + ne*OIICollisionalIonizationRate(T),
        .recombination2 = ne*OIIIRecombinationRate(T),
    };
    /* In H II region mode, couple O and H through charge transfer. */
    status = CEQLibRunParameters.UseHIIRegion ?
        CEQSolveHydrogenOxygenChargeTransfer(state,hydrogen,&oxygen,accuracy) :
        CEQSolveIonLadder(state->n[ipOXYGEN]-state->nmol[ipCO],
                          oxygen.ionization0, oxygen.recombination1,
                          oxygen.ionization1, oxygen.recombination2,
                          &state->natom[ipOXYGEN], &state->nion[ipOXYGEN], &state->nOIII);
    if(status != CEQ_SUCCESS) return status;

    return CEQSolveIonLadder(state->n[ipSILICON],
        SiIPhotoionizationRate(NLyn) + ne*SiICollisionalIonizationRate(T)
            + X0*(CosmicRayIonizationRateSiI(xe) + SoftXrayIonizationRateSiI(xe,Nt)),
        ne*SiIIRecombinationRate(T) + nH*SiIIRecombinationRateOnDust(T,G0,ne,Zd),
        SiIIPhotoionizationRate(NLyn) + ne*SiIICollisionalIonizationRate(T)
            + X0*(CosmicRayIonizationRateSiII(xe) + SoftXrayIonizationRateSiII(xe,Nt)),
        ne*SiIIIRecombinationRate(T),
        &state->natom[ipSILICON], &state->nion[ipSILICON], &state->nSiIII);
}

static ChargePoint EvaluateCharge(ChargeSearch *search, double electronDensity)
{
    StructCEQData state = *search->initial;
    state.nel = electronDensity;
    const CEQHydrogenRates hydrogen = CEQHydrogenIonizationRates(&state);
    int status = CEQSolveHydrogenEquilibriumWithRates(&state,&hydrogen,search->accuracy);
    if(status == CEQ_SUCCESS)
        status = CEQSolveCarbonWithIntervals(&state,search->previousCO,search->accuracy,
                                             search->carbonIntervals,search->carbonIntervalCount);
    if(status == CEQ_SUCCESS)
        status = SolveOtherIons(&state,&hydrogen,search->accuracy);

    ChargePoint point = {.density = electronDensity, .status = status};
    if(status != CEQ_SUCCESS){
        if(status == CEQ_NO_CONVERGENCE || search->failure == CEQ_NO_EQUILIBRIUM)
            search->failure = status;
        return point;
    }

    const double charge = CEQLibGetElectronNumberDensity(&state);
    if(!isfinite(charge) || charge < 0.0){
        point.status = search->failure = CEQ_NO_CONVERGENCE;
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
        if(next.status != CEQ_SUCCESS) break;

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
        if(next.status != CEQ_SUCCESS)
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
    /* Outside the X-ray and molecular regimes, the charge returned by each
       atomic ladder divided by ne decreases with ne. Its ionization rates
       divided by ne are nonincreasing and its recombination rates increase.
       A local bracket within this interval therefore retains the same root. */
    if(search->hasMolecularChargeCoupling ||
            (search->hasXrayIonization && lower->density/search->initial->n[ipHYDROGEN] <= 0.5))
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
    if(lower.status != CEQ_SUCCESS || upper.status != CEQ_SUCCESS ||
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
        if(trial.status != CEQ_SUCCESS && density != middle){
            search->failure = failure;
            trial = EvaluateCharge(search,middle);
        }
        if(trial.status != CEQ_SUCCESS) return;
        if(fabs(trial.residual) <= search->accuracy) return;
        if((trial.residual < 0.0) == (lower.residual < 0.0)) lower = trial;
        else upper = trial;
    }
    search->failure = CEQ_NO_CONVERGENCE;
}

static void SearchChargeAroundHint(ChargeSearch *search, ChargePoint hint,
                                   double minimumDensity, double maximumDensity)
{
    if(hint.status != CEQ_SUCCESS || hint.density <= 0.0 || hint.density < minimumDensity)
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
        hint.density < xrayRangeStart && lower.status == CEQ_SUCCESS && lower.residual == 0.0;

    if(lowerDensity == 0.0){
        /* Test the positive side separately; pure collisional gas can also have
           an exactly neutral solution at zero electron density. */
        double positiveDensity = sqrt(DBL_MIN)*sqrt(upperDensity);
        // Probe the lower end of the fitted X-ray electron fraction range.
        if(search->hasXrayIonization)
            positiveDensity = fmin(xrayRangeStart,0.5*upperDensity);
        if(positiveDensity > lowerDensity && positiveDensity < upperDensity){
            ChargePoint positive = EvaluateCharge(search,positiveDensity);
            if(positive.status == CEQ_SUCCESS){
                /* Keep the neutral-limit probes on the existing bisection grid. */
                if(hint.status == CEQ_SUCCESS && hint.density > lower.density && hint.density < positive.density){
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
    if(hint.status == CEQ_SUCCESS && hint.density > lower.density && hint.density < upper.density){
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
    if(CEQLibRunParameters.UseCOFormation && search->initial->T <= 1.0e4 &&
            search->solutions == beforeBisection &&
            (search->solutions == initialSolutions || lower.status != CEQ_SUCCESS || upper.status != CEQ_SUCCESS)){
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

static bool ConservesElements(const StructCEQData *state)
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

int CEQLibChemicalReactionEquilibrium(StructCEQData *state, const double accuracy)
{
    if(!state || !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0 ||
            !isfinite(state->T) || state->T <= 0.0) return CEQ_INVALID_INPUT;

    const double conditions[] = {state->G0, state->X0, state->NLyn, state->Nt, state->Z, state->Zd};
    for(unsigned int i=0;i<sizeof(conditions)/sizeof(conditions[0]);i++){
        if(!isfinite(conditions[i]) || conditions[i] < 0.0) return CEQ_INVALID_INPUT;
    }
    const double nH = state->n[ipHYDROGEN];
    if(nH <= 0.0) return CEQ_INVALID_INPUT;
    double maximumCharge = nH;
    for(int i=0;i<NumberChemicalElement;i++){
        if(!isfinite(state->n[i]) || state->n[i] < 0.0) return CEQ_INVALID_INPUT;
        if(i != ipHYDROGEN) maximumCharge += 2.0*state->n[i];
    }
    if(!isfinite(maximumCharge)) return CEQ_INVALID_INPUT;

    const int shieldingPoints = CEQLibRunParameters.UseCOFormation && state->T <= 1.0e4 ?
                                CEQCOShieldingPointCount() : 0;
    CEQCarbonShieldingInterval carbonIntervals[shieldingPoints+1];
    const int carbonIntervalCount = CEQPrepareCarbonShieldingIntervals(state,carbonIntervals,shieldingPoints);
    ChargeSearch search = {
        .initial = state,
        .previousElectronDensity = isfinite(state->nel) ? fmin(fmax(state->nel,0.0),maximumCharge) : 0.0,
        .previousCO = CEQLibRunParameters.UseCOFormation ? state->nmol[ipCO] : 0.0,
        .accuracy = accuracy,
        .carbonIntervals = carbonIntervals,
        .carbonIntervalCount = carbonIntervalCount,
        .hasXrayIonization = CEQLibRunParameters.UseSoftXrayIonization && state->X0 > 0.0,
        .hasMolecularChargeCoupling = CEQLibRunParameters.UseH2Formation ||
                                     (CEQLibRunParameters.UseCOFormation && state->T <= 1.0e4),
        .failure = CEQ_NO_EQUILIBRIUM
    };
    ChargePoint hint = EvaluateCharge(&search,search.previousElectronDensity);
    /* With only photoionization, electron collisions and electron
       recombination, every ionization/recombination ratio is of the form
       A + B/ne at fixed T. The ion charge is nonincreasing, giving a unique
       positive root. Leave branch selection to the original search when
       molecular coupling, grain recombination, cosmic rays or X-rays apply.
       H-O charge transfer conserves charge and couples the H and O ladders
       only through O/H. The argument above does not cover that coupling, so
       HIIRegionTest compares it with solves from an ionized seed. */
    if(!search.found && hint.status == CEQ_SUCCESS && hint.residual > 0.5 &&
            CEQLibRunParameters.UseHIIRegion && state->NLyn > 0.0 &&
            !search.hasMolecularChargeCoupling && state->X0 == 0.0 &&
            (!CEQLibRunParameters.UseGrainSurfaceProcesses || state->Zd == 0.0))
        SearchPhotoionizedCharge(&search,hint,maximumCharge);
    if(!search.found){
        if(search.hasXrayIonization){
            double boundary = 0.5*nH;
            /* The rate uses ne/nH, so begin the second segment at the first
               representable density whose computed fraction exceeds 0.5. */
            while(boundary/nH > 0.5) boundary = nextafter(boundary,0.0);
            double afterBoundary = nextafter(boundary,INFINITY);
            while(afterBoundary/nH <= 0.5) afterBoundary = nextafter(afterBoundary,INFINITY);
            if(hint.density <= boundary){
                double rangeStart = XrayFitMinimumElectronFraction*nH;
                if(rangeStart/nH < XrayFitMinimumElectronFraction)
                    rangeStart = nextafter(rangeStart,INFINITY);
                SearchChargeAroundHint(&search,hint,rangeStart,afterBoundary);
            }else{
                SearchChargeAroundHint(&search,hint,afterBoundary,maximumCharge);
            }
            if(!search.found){
                SearchChargeInterval(&search,0.0,boundary,hint);
                if(afterBoundary <= maximumCharge)
                    SearchChargeInterval(&search,afterBoundary,maximumCharge,hint);
            }
        }else{
            SearchChargeAroundHint(&search,hint,0.0,maximumCharge);
            if(!search.found) SearchChargeInterval(&search,0.0,maximumCharge,hint);
        }
    }

    if(!search.found) return search.failure;
    StructCEQData *solution = &search.solution;
    if(!ConservesElements(solution)) return CEQ_NO_CONVERGENCE;
    solution->MeanMolecularWeight = CEQLibGetMeanMolecularWeight(solution);
    solution->DegreeOfFreedom = CEQLibGetMeanDegreeOfFreedom(solution);
    if(!isfinite(solution->MeanMolecularWeight) || solution->MeanMolecularWeight <= 0.0 ||
            !isfinite(solution->DegreeOfFreedom) || solution->DegreeOfFreedom <= 0.0)
        return CEQ_NO_CONVERGENCE;
    solution->Gamma = (1.0+solution->DegreeOfFreedom)/solution->DegreeOfFreedom;
    solution->Gm1 = solution->Gamma-1.0;
    solution->GGm1 = solution->Gamma*solution->Gm1;
    *state = *solution;
    return CEQ_SUCCESS;
}
