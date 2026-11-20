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

typedef struct {
    double density;
    double residual;
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

static int SolveOtherIons(StructCEQData *state)
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

    status = CEQSolveIonLadder(state->n[ipOXYGEN]-state->nmol[ipCO],
        OIPhotoionizationRate(NLyn) + ne*OICollisionalIonizationRate(T)
            + X0*CosmicRayIonizationRateOI(xe),
        ne*OIIRecombinationRate(T) + nH*OIIRecombinationRateOnDust(T,G0,ne,Zd),
        OIIPhotoionizationRate(NLyn) + ne*OIICollisionalIonizationRate(T),
        ne*OIIIRecombinationRate(T),
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
    int status = CEQSolveHydrogenEquilibrium(&state,search->accuracy);
    if(status == CEQ_SUCCESS)
        status = CEQSolveCarbonWithIntervals(&state,search->previousCO,search->accuracy,
                                             search->carbonIntervals,search->carbonIntervalCount);
    if(status == CEQ_SUCCESS)
        status = SolveOtherIons(&state);

    ChargePoint point = {.density = electronDensity, .status = status};
    if(status != CEQ_SUCCESS){
        if(status == CEQ_NO_CONVERGENCE || search->failure == CEQ_NO_EQUILIBRIUM)
            search->failure = status;
        return point;
    }

    double charge = CEQLibGetElectronNumberDensity(&state);
    if(!isfinite(charge) || charge < 0.0){
        point.status = search->failure = CEQ_NO_CONVERGENCE;
        return point;
    }
    double scale = fmax(charge,electronDensity);
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

    ChargePoint previous = increase ? *lower : *upper;
    const int failure = search->failure;
    double factor = 1.05;
    for(int step=0; step<40; step++){
        const double density = increase ? previous.density * factor : previous.density / factor;
        if(density <= lower->density || density >= upper->density || density == previous.density)
            break;

        const ChargePoint next = EvaluateCharge(search,density);
        if(next.status != CEQ_SUCCESS)
            break;
        if(fabs(next.residual) <= search->accuracy)
            return true;
        if((next.residual < 0.0) != (previous.residual < 0.0)){
            *lower = increase ? previous : next;
            *upper = increase ? next : previous;
            return false;
        }
        previous = next;
        factor = fmin(10.0,factor*factor);
    }
    search->failure = failure;
    return false;
}

static void BisectCharge(ChargeSearch *search, ChargePoint lower, ChargePoint upper,
                        bool interpolate)
{
    if(lower.status != CEQ_SUCCESS || upper.status != CEQ_SUCCESS ||
            fabs(lower.residual) <= search->accuracy || fabs(upper.residual) <= search->accuracy ||
            (lower.residual < 0.0) == (upper.residual < 0.0))
        return;

    if(search->hasXrayIonization && lower.density < 1.0e-4*search->initial->n[ipHYDROGEN])
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

static void SearchChargeInterval(ChargeSearch *search, double lowerDensity, double upperDensity,
                                 ChargePoint hint)
{
    unsigned int initialSolutions = search->solutions;
    ChargePoint lower = EvaluateCharge(search,lowerDensity);
    ChargePoint upper = EvaluateCharge(search,upperDensity);
    double xrayRangeStart = 1.0e-4*search->initial->n[ipHYDROGEN];
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
    if(!search.found){
        if(search.hasXrayIonization){
            double boundary = 0.5*nH;
            /* The rate uses ne/nH, so begin the second segment at the first
               representable density whose computed fraction exceeds 0.5. */
            while(boundary/nH > 0.5) boundary = nextafter(boundary,0.0);
            double afterBoundary = nextafter(boundary,INFINITY);
            while(afterBoundary/nH <= 0.5) afterBoundary = nextafter(afterBoundary,INFINITY);
            SearchChargeInterval(&search,0.0,boundary,hint);
            if(afterBoundary <= maximumCharge)
                SearchChargeInterval(&search,afterBoundary,maximumCharge,hint);
        }else{
            SearchChargeInterval(&search,0.0,maximumCharge,hint);
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
