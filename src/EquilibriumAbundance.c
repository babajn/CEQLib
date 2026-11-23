#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include "Astro.h"
#include "Constants.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "EquilibriumAbundance.h"
#include "CollisionalGasPhaseChemistry.h"
#include "PhotoelectronDeposition.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#include "CosmicRay.h"
#include "Xray.h"
    
#define TURBULENCE_VELOCITY (3.0)   // km/s
double b_turb = TURBULENCE_VELOCITY;    // km/s

typedef enum {
    TemperatureIonizeHI, TemperatureRecombineHIICaseA, TemperatureRecombineHIICaseB,
    TemperatureIonizeHeI, TemperatureIonizeHeII,
    TemperatureRecombineHeIICaseA, TemperatureRecombineHeIICaseB,
    TemperatureRecombineHeIIICaseA, TemperatureRecombineHeIIICaseB,
    TemperatureIonizeCI, TemperatureIonizeCII, TemperatureRecombineCII, TemperatureRecombineCIII,
    TemperatureIonizeOI, TemperatureIonizeOII, TemperatureRecombineOII, TemperatureRecombineOIII,
    TemperatureTransferToOII, TemperatureTransferToHI, TemperatureTransferOIII,
    TemperatureDissociateH2Electron, TemperatureDissociateH2Atom, TemperatureDissociateH2Molecule,
    ChemistryTemperatureTableCount,
} ChemistryTemperatureRate;

static CEQTemperatureTable chemistryTemperatureTables[ChemistryTemperatureTableCount];

static double RecombineHeIICaseA(double T){
    return HeIIRecombinationRateCaseA(T) + HeIIDielectronicRecombinationRate(T);
}

static double RecombineHeIICaseB(double T){
    return HeIIRecombinationRateCaseB(T) + HeIIDielectronicRecombinationRate(T);
}

static double DissociateH2Electron(double T){
    return H2DissociationCoefficients(T).electron;
}

static double DissociateH2Atom(double T){
    return H2DissociationCoefficients(T).atomicHydrogen;
}

static double DissociateH2Molecule(double T){
    return H2DissociationCoefficients(T).molecularHydrogen;
}

static double (*const chemistryTemperatureRates[ChemistryTemperatureTableCount])(double) = {
    [TemperatureIonizeHI] = HICollisionalIonizationRate,
    [TemperatureRecombineHIICaseA] = HIIRecombinationRateCaseA,
    [TemperatureRecombineHIICaseB] = HIIRecombinationRateCaseB,
    [TemperatureIonizeHeI] = HeICollisionalIonizationRate,
    [TemperatureIonizeHeII] = HeIICollisionalIonizationRate,
    [TemperatureRecombineHeIICaseA] = RecombineHeIICaseA,
    [TemperatureRecombineHeIICaseB] = RecombineHeIICaseB,
    [TemperatureRecombineHeIIICaseA] = HeIIIRecombinationRateCaseA,
    [TemperatureRecombineHeIIICaseB] = HeIIIRecombinationRateCaseB,
    [TemperatureIonizeCI] = CICollisionalIonizationRate,
    [TemperatureIonizeCII] = CIICollisionalIonizationRate,
    [TemperatureRecombineCII] = CIIRecombinationRate,
    [TemperatureRecombineCIII] = CIIIRecombinationRate,
    [TemperatureIonizeOI] = OICollisionalIonizationRate,
    [TemperatureIonizeOII] = OIICollisionalIonizationRate,
    [TemperatureRecombineOII] = OIIRecombinationRate,
    [TemperatureRecombineOIII] = OIIIRecombinationRate,
    [TemperatureTransferToOII] = OIChargeTransferRateToHII,
    [TemperatureTransferToHI] = OIIChargeTransferRateToHI,
    [TemperatureTransferOIII] = OIIIChargeTransferRateToHI,
    [TemperatureDissociateH2Electron] = DissociateH2Electron,
    [TemperatureDissociateH2Atom] = DissociateH2Atom,
    [TemperatureDissociateH2Molecule] = DissociateH2Molecule,
};

void CEQInitializeChemistryTemperatureTables(void){
    /* C II branches apply to the optional Glover & Jappsen fit. Protect the
       boundaries for either build, and the clamps in O III charge transfer. */
    static const double carbonBreaks[] = {7950.0,21140.0};
    static const double transferBreaks[] = {100.0,1.e5};
    for(int rate=0;rate<ChemistryTemperatureTableCount;rate++){
        const double *breaks = NULL;
        size_t count = 0;
        if(rate == TemperatureRecombineCII){
            breaks = carbonBreaks;
            count = sizeof carbonBreaks/sizeof carbonBreaks[0];
        }else if(rate == TemperatureTransferOIII){
            breaks = transferBreaks;
            count = sizeof transferBreaks/sizeof transferBreaks[0];
        }
        CEQFillTemperatureTable(&chemistryTemperatureTables[rate],chemistryTemperatureRates[rate],breaks,count);
    }
}

static inline double ChemistryCoefficient(ChemistryTemperatureRate rate, double T, const CEQTemperaturePoint *point){
    if(!point) return chemistryTemperatureRates[rate](T);
    return CEQReadTemperatureTable(&chemistryTemperatureTables[rate],point,T,chemistryTemperatureRates[rate]);
}

H2CollisionalDissociationCoefficients CEQH2DissociationCoefficientsWithTemperatureTable(
        double T, const CEQTemperaturePoint *point){
    if(!point) return H2DissociationCoefficients(T);
    H2CollisionalDissociationCoefficients coefficients = {
        .electron = CEQReadTemperatureTable(&chemistryTemperatureTables[TemperatureDissociateH2Electron],point,T,NULL),
        .atomicHydrogen = CEQReadTemperatureTable(&chemistryTemperatureTables[TemperatureDissociateH2Atom],point,T,NULL),
        .molecularHydrogen = CEQReadTemperatureTable(&chemistryTemperatureTables[TemperatureDissociateH2Molecule],point,T,NULL),
    };
    if(!isfinite(coefficients.electron) || !isfinite(coefficients.atomicHydrogen) ||
       !isfinite(coefficients.molecularHydrogen)){
        const H2CollisionalDissociationCoefficients direct = H2DissociationCoefficients(T);
        if(!isfinite(coefficients.electron)) coefficients.electron = direct.electron;
        if(!isfinite(coefficients.atomicHydrogen)) coefficients.atomicHydrogen = direct.atomicHydrogen;
        if(!isfinite(coefficients.molecularHydrogen)) coefficients.molecularHydrogen = direct.molecularHydrogen;
    }
    return coefficients;
}


/* Select the stationary class accessible from the neutral species whenever a
 * vanishing transition disconnects the ion ladder. */
int CEQSolveIonLadder(double density, double ionization0, double recombination1,
                      double ionization1, double recombination2,
                      double *neutral, double *ion, double *doubleIon){
    if(!isfinite(density) || density < 0.0 ||
       !isfinite(ionization0) || ionization0 < 0.0 ||
       !isfinite(recombination1) || recombination1 < 0.0 ||
       !isfinite(ionization1) || ionization1 < 0.0 ||
       !isfinite(recombination2) || recombination2 < 0.0) return EQC_INVALID_INPUT;

    if(density == 0.0 || ionization0 == 0.0){
        *neutral = density;
        *ion = *doubleIon = 0.0;
        return EQC_SUCCESS;
    }
    if(ionization1 > 0.0 && recombination2 == 0.0){
        *neutral = *ion = 0.0;
        *doubleIon = density;
        return EQC_SUCCESS;
    }

    const double lowerScale = fmax(ionization0,recombination1);
    const double upperScale = fmax(ionization1,recombination2);
    const double upperDown = ionization1 == 0.0 ? 1.0 : recombination2/upperScale;
    const double upperUp = ionization1 == 0.0 ? 0.0 : ionization1/upperScale;
    double weights[] = {
        (recombination1/lowerScale)*upperDown,
        (ionization0/lowerScale)*upperDown,
        (ionization0/lowerScale)*upperUp,
    };
    double sum = weights[0] + weights[1] + weights[2];
    if(sum < DBL_MIN){
        // Extreme rate ratios can underflow both scaled products.
        double logs[] = {log(recombination1) + log(recombination2),
                         log(ionization0) + log(recombination2),
                         log(ionization0) + log(ionization1)};
        const double highest = fmax(logs[0],fmax(logs[1],logs[2]));
        for(int i=0;i<3;i++) weights[i] = exp(logs[i] - highest);
        sum = weights[0] + weights[1] + weights[2];
    }
    *neutral = density*(weights[0]/sum);
    *ion = density*(weights[1]/sum);
    *doubleIon = density*(weights[2]/sum);
    return EQC_SUCCESS;
}

/* A bounded residual keeps the tolerance relative even for trace abundances. */
static double FluxResidual(double destroyed, double formed){
    if(!isfinite(destroyed) || destroyed < 0.0 || !isfinite(formed) || formed < 0.0) return NAN;
    const double scale = fmax(destroyed,formed);
    return scale > 0.0 ? (destroyed-formed)/scale : 0.0;
}

typedef struct {
    double fraction;
    double residual;
} FluxPoint;

typedef struct {
    const StructEQCData *state;
    const CEQMolecularHydrogenCoefficients *h2;
    double neutralFraction;
    bool useAtomicFraction;
} HydrogenBalance;

typedef struct {
    const HydrogenBalance *balance;
    double nH2, fixedDestruction;
} TraceAtomicFlux;

/* Under binary64 round-to-nearest, 1 - atomic rounds to 1 in this range. */
static const double TraceAtomicFraction = 0x1p-54;

static double HydrogenFluxResidual(double fraction, const void *parameters);
static bool CacheTraceAtomicFlux(const HydrogenBalance *balance, TraceAtomicFlux *cache);
static double CachedTraceAtomicResidual(double fraction, const TraceAtomicFlux *cache);
static int SkipSaturatedHydrogenBisections(FluxPoint *upper, const void *parameters, int budget, bool wide);

static int SolveFluxBracket(double (*residual)(double,const void *), const void *parameters,
                            FluxPoint lower, FluxPoint upper, double accuracy,
                            double initialFraction, double *root){
    if(!isfinite(lower.residual) || !isfinite(upper.residual)) return EQC_INVALID_INPUT;
    if(fabs(lower.residual) <= accuracy){ *root = lower.fraction; return EQC_SUCCESS; }
    if(fabs(upper.residual) <= accuracy){ *root = upper.fraction; return EQC_SUCCESS; }
    if((lower.residual < 0.0) == (upper.residual < 0.0)) return EQC_NO_EQUILIBRIUM;

    // Interpolate between periodic bisections and retain the sign bracket.
    const int limit = 2*(DBL_MAX_EXP+DBL_MANT_DIG);
    bool tryTraceSkip = residual == HydrogenFluxResidual;
    bool tryWideSkip = tryTraceSkip;
    bool tryTraceCache = tryTraceSkip;
    bool useTraceCache = false;
    TraceAtomicFlux traceCache;
    for(int i=0;i<limit;i++){
        const bool useInitialFraction = i == 0 && initialFraction > lower.fraction && initialFraction < upper.fraction;
        const bool saturated = !useInitialFraction && lower.fraction == 0.0
                            && lower.residual == 1.0 && upper.residual == -1.0;
        if(tryWideSkip && saturated && upper.fraction > TraceAtomicFraction){
            tryWideSkip = false;
            i += SkipSaturatedHydrogenBisections(&upper,parameters,limit-i,true);
            if(i == limit) break;
        }
        if(tryTraceSkip && saturated && upper.fraction <= TraceAtomicFraction){
            /* The opposite saturated residuals make the secant a midpoint, too.
               Count the omitted halvings so later interpolation keeps its parity. */
            tryTraceSkip = false;
            i += SkipSaturatedHydrogenBisections(&upper,parameters,limit-i,false);
            if(i == limit) break;
        }
        if(tryTraceCache && upper.fraction <= TraceAtomicFraction){
            tryTraceCache = false;
            useTraceCache = CacheTraceAtomicFlux(parameters,&traceCache);
        }
        double middle = lower.fraction + 0.5*(upper.fraction-lower.fraction);
        if(useInitialFraction){
            middle = initialFraction;
        }else if(i % 2 == 0){
            const double fraction = lower.residual/(lower.residual-upper.residual);
            const double estimate = lower.fraction + fraction*(upper.fraction-lower.fraction);
            if(estimate > lower.fraction && estimate < upper.fraction) middle = estimate;
        }else if(lower.fraction > 0.0 && lower.fraction < 0.25*upper.fraction){
            // Bisect in log space for a bracket spanning trace abundances.
            const double estimate = sqrt(lower.fraction)*sqrt(upper.fraction);
            if(estimate > lower.fraction && estimate < upper.fraction) middle = estimate;
        }
        if(middle == lower.fraction || middle == upper.fraction) break;
        const double value = useTraceCache ? CachedTraceAtomicResidual(middle,&traceCache)
                                           : residual(middle,parameters);
        if(!isfinite(value)) return EQC_INVALID_INPUT;
        if(fabs(value) <= accuracy){ *root = middle; return EQC_SUCCESS; }
        if((value < 0.0) == (lower.residual < 0.0))
            lower = (FluxPoint){middle,value};
        else
            upper = (FluxPoint){middle,value};
    }
    return EQC_NO_CONVERGENCE;
}

typedef struct {
    double destroyed;
    double formed;
} HydrogenFluxes;

static HydrogenFluxes MolecularHydrogenFluxes(const HydrogenBalance *balance, double atomic, double molecular){
    const StructEQCData *state = balance->state;
    const CEQMolecularHydrogenCoefficients *h2 = balance->h2;
    const double nH = state->n[ipHYDROGEN];
    const double nHI = (nH*atomic)*balance->neutralFraction;
    const double nH2 = (0.5*nH)*molecular;
    const double photo = H2PhotoDissociationRateWithDustFactor(state->G0,(0.5*molecular)*state->Nt,
                                                               h2->doppler,h2->dustShielding);
    const double collision = H2DissociationFromCoefficients(h2->collision,state->nel,nHI,nH2);
    return (HydrogenFluxes){(photo+h2->cosmicRayDissociation+collision)*nH2,h2->formation*nHI};
}

double CEQH2DestructionRate(const StructEQCData *state, const CEQMolecularHydrogenCoefficients *h2){
    const double nH2 = state->nmol[ipH2];
    const double column = (nH2/state->n[ipHYDROGEN])*state->Nt;
    const double photo = H2PhotoDissociationRateWithDustFactor(state->G0,column,h2->doppler,h2->dustShielding);
    const double collision = H2DissociationFromCoefficients(h2->collision,state->nel,state->natom[ipHYDROGEN],nH2);
    return photo+h2->cosmicRayDissociation+collision;
}

static double HydrogenFluxResidual(double fraction, const void *parameters){
    const HydrogenBalance *balance = parameters;
    const double molecular = balance->useAtomicFraction ? 1.0-fraction : fraction;
    const double atomic = balance->useAtomicFraction ? fraction : 1.0-fraction;
    const HydrogenFluxes flux = MolecularHydrogenFluxes(balance,atomic,molecular);
    return FluxResidual(flux.destroyed,flux.formed);
}

static bool CacheTraceAtomicFlux(const HydrogenBalance *balance, TraceAtomicFlux *cache){
    if(FLT_RADIX != 2 || DBL_MANT_DIG != 53 || FLT_EVAL_METHOD != 0 ||
       !balance->useAtomicFraction || fegetround() != FE_TONEAREST) return false;

    const StructEQCData *state = balance->state;
    cache->balance = balance;
    cache->nH2 = 0.5*state->n[ipHYDROGEN];
    const CEQMolecularHydrogenCoefficients *h2 = balance->h2;
    const double photo = H2PhotoDissociationRateWithDustFactor(state->G0,0.5*state->Nt,
                                                               h2->doppler,h2->dustShielding);
    cache->fixedDestruction = photo+h2->cosmicRayDissociation;
    return true;
}

static double CachedTraceAtomicResidual(double fraction, const TraceAtomicFlux *cache){
    const HydrogenBalance *balance = cache->balance;
    const StructEQCData *state = balance->state;
    const CEQMolecularHydrogenCoefficients *h2 = balance->h2;
    const double nHI = (state->n[ipHYDROGEN]*fraction)*balance->neutralFraction;
    const double collision = H2DissociationFromCoefficients(h2->collision,state->nel,nHI,cache->nH2);
    return FluxResidual((cache->fixedDestruction+collision)*cache->nH2,h2->formation*nHI);
}

static int SkipSaturatedHydrogenBisections(FluxPoint *upper, const void *parameters, int budget, bool wide){
    const HydrogenBalance *balance = parameters;
    /* Bound every skipped residual with an upper bound on destruction and
       the formation flux at the new, smaller atomic fraction. */
    if(FLT_RADIX != 2 || DBL_MANT_DIG != 53 || FLT_EVAL_METHOD != 0 ||
       !balance->useAtomicFraction || !(upper->fraction >= DBL_MIN) ||
       !(balance->state->n[ipHYDROGEN] > 0.0) || !(balance->neutralFraction > 0.0) ||
       !(balance->h2->formation > 0.0) || !(balance->h2->collision.atomicHydrogen >= 0.0) ||
       fegetround() != FE_TONEAREST) return 0;

    HydrogenFluxes flux;
    if(!wide){
        /* With constant H2 density and shielding, destruction increases
           with H I, so the old upper endpoint bounds the whole range. */
        flux = MolecularHydrogenFluxes(balance,upper->fraction,1.0);
    }else{
        const StructEQCData *s = balance->state;
        const CEQMolecularHydrogenCoefficients *h2 = balance->h2;
        if(!(s->G0 >= 0.0) || !(h2->dustShielding >= 0.0) ||
           !(h2->doppler > 0.0) || !isfinite(h2->doppler) ||
           !(s->Nt >= 0.0) || !isfinite(s->Nt) || !(s->nel >= 0.0) ||
           !(h2->collision.electron >= 0.0) || !(h2->collision.molecularHydrogen >= 0.0)) return 0;
        /* Bound the photo rate by the optically thin value with a floating-
           point safety margin. In both shielding fits the factor is at most
           1/0.965 even ignoring the decrease in the exponential. Require a
           normal optically thin result unless the photo source is absent. */
        const double thin = H2PhotoDissociationRateWithDustFactor(s->G0,0.0,h2->doppler,h2->dustShielding);
        if(s->G0 > 0.0 && h2->dustShielding > 0.0 && !(thin >= DBL_MIN)) return 0;
        const double nH = s->n[ipHYDROGEN];
        const double nHImax = (nH*upper->fraction)*balance->neutralFraction;
        const double nH2max = 0.5*nH;
        const double maximumCollision = H2DissociationFromCoefficients(h2->collision,s->nel,nHImax,nH2max);
        const double destructionBound = 4.0*thin+h2->cosmicRayDissociation+maximumCollision;
        flux = (HydrogenFluxes){destructionBound*nH2max,h2->formation*nHImax};
    }
    if(!(flux.destroyed > 0.0) || !isfinite(flux.destroyed) ||
       !(flux.formed > 0.0) || !isfinite(flux.formed)) return 0;

    /* Estimate how far we can halve while staying below half an ulp of the
       formation term. The bound below verifies the estimate with rounded
       fluxes. Keep the endpoints normal so all skipped halvings are exact. */
    int steps = ilogb(flux.formed)-ilogb(flux.destroyed)-58;
    const int normalSteps = ilogb(upper->fraction)-(DBL_MIN_EXP-1);
    if(steps > normalSteps) steps = normalSteps;
    if(steps > budget) steps = budget;
    if(steps < 8) return 0;

    const double candidate = scalbn(upper->fraction,-steps);
    const double nHI = (balance->state->n[ipHYDROGEN]*candidate)*balance->neutralFraction;
    const double bound = (balance->h2->formation*nHI)*0x1p-56;
    if(!(bound >= DBL_MIN) || flux.destroyed > bound) return 0;

    /* At each skipped point the subtraction rounds to -formation, so the
       original residual is exactly -1 and the lower endpoint does not move. */
    upper->fraction = candidate;
    return steps;
}

#ifndef EQC_H2_LEGACY_SHIELDING
/* Stability of a balanced state against changes in the H2 fraction.
   The derivative is with respect to ln(molecular/atomic). */
static bool StableMolecularHydrogen(const HydrogenBalance *balance, double atomic, double molecular){
    if(atomic == 0.0 || molecular == 0.0) return true;
    const StructEQCData *state = balance->state;
    const CEQMolecularHydrogenCoefficients *h2 = balance->h2;
    const double nH = state->n[ipHYDROGEN];
    const double byAtom = nH*balance->neutralFraction*h2->collision.atomicHydrogen;
    const double byMolecule = 0.5*nH*h2->collision.molecularHydrogen;
    double shieldSlope;
    const double photo = CEQH2PhotoRateAndShieldingSlope(state->G0,(0.5*molecular)*state->Nt,
                                                         h2->doppler,h2->dustShielding,&shieldSlope);
    const double destruction = photo+h2->cosmicRayDissociation+state->nel*h2->collision.electron
                               + byAtom*atomic+byMolecule*molecular;
    return destruction+atomic*(photo*shieldSlope+molecular*(byMolecule-byAtom)) >= 0.0;
}
#endif

static void RefreshMolecularHydrogenCoefficients(const StructEQCData *state,
                                                 CEQMolecularHydrogenCoefficients *h2,
                                                 const CEQTemperaturePoint *point){
    const double dustTemperature = EffectiveDustTemperature(state);
    h2->formation = state->n[ipHYDROGEN]*H2FormationRate(state->T,state->Zd,dustTemperature);
    h2->doppler = H2DopplerWidth(state->T);
    h2->collision = CEQH2DissociationCoefficientsWithTemperatureTable(state->T,point);
}

static CEQMolecularHydrogenCoefficients MolecularHydrogenRadiationCoefficients(const StructEQCData *state){
    return (CEQMolecularHydrogenCoefficients){
        .dustShielding = DustShieldingFactorFUV(state->Nt,state->Zd),
        .cosmicRayDissociation = CEQH2CosmicRayDissociationRate(state->X0),
    };
}

CEQMolecularHydrogenCoefficients CEQPrepareMolecularHydrogenCoefficients(const StructEQCData *state){
    CEQMolecularHydrogenCoefficients h2 = MolecularHydrogenRadiationCoefficients(state);
    RefreshMolecularHydrogenCoefficients(state,&h2,NULL);
    return h2;
}

#ifndef EQC_H2_LEGACY_SHIELDING
/* Below this H2 column the line shielding changes the photodissociation rate
   by less than 0.2%, so the balance has a single root there. */
static const double TraceShieldingColumn = 1.e12;

/* Searched fraction at which the H2 column reaches TraceShieldingColumn,
   or a negative value when the whole range stays below it. */
static double TraceShieldingFraction(const HydrogenBalance *balance){
    const double Nt = balance->state->Nt;
    if(!(0.5*Nt > TraceShieldingColumn)) return -1.0;
    const double molecular = 2.0*TraceShieldingColumn/Nt;
    return balance->useAtomicFraction ? 1.0-molecular : molecular;
}

static bool BelowTraceShielding(const HydrogenBalance *balance, double boundary, double fraction){
    if(boundary < 0.0) return true;
    return balance->useAtomicFraction ? fraction > boundary : fraction < boundary;
}
#endif

/* Follow the sign of net H2 production away from the incoming fraction.
   Returns EQC_NO_EQUILIBRIUM when there is no bracket in this half of the range. */
static int SolveMolecularHydrogenNear(const HydrogenBalance *balance, double start, double half,
                                      double accuracy, double *root){
    FluxPoint near = {start,HydrogenFluxResidual(start,balance)};
    if(!isfinite(near.residual)) return EQC_INVALID_INPUT;
    if(fabs(near.residual) <= accuracy){ *root = start; return EQC_SUCCESS; }
    const bool upward = balance->useAtomicFraction ? near.residual > 0.0 : near.residual < 0.0;
    double step = fmax(2.0*fabs(near.residual),1.e-12);
#ifndef EQC_H2_LEGACY_SHIELDING
    /* Trace H2 may sit many decades from the incoming fraction. Grow the step
       freely below the shielding column and keep it local above. */
    const double boundary = TraceShieldingFraction(balance);
    if(!BelowTraceShielding(balance,boundary,near.fraction)) step = fmin(step,0.25);
#endif
    for(int i=0;i<64;i++){
        FluxPoint far;
        if(upward){
            double fraction = near.fraction*(1.0+step);
#ifndef EQC_H2_LEGACY_SHIELDING
            /* Stop at the shielding column so that no pair of shielded roots is skipped. */
            if(boundary > near.fraction && fraction > boundary) fraction = boundary;
#endif
            far = fraction >= 0.5 ? (FluxPoint){0.5,half} : (FluxPoint){fraction,HydrogenFluxResidual(fraction,balance)};
        }else{
            double fraction = near.fraction/(1.0+step);
#ifndef EQC_H2_LEGACY_SHIELDING
            if(boundary > 0.0 && boundary < near.fraction && fraction < boundary) fraction = boundary;
#endif
            if(!(fraction > 0.0)) return EQC_NO_EQUILIBRIUM;
            far = (FluxPoint){fraction,HydrogenFluxResidual(fraction,balance)};
        }
        if(!isfinite(far.residual)) return EQC_INVALID_INPUT;
        if(fabs(far.residual) <= accuracy){ *root = far.fraction; return EQC_SUCCESS; }
        if((far.residual < 0.0) != (near.residual < 0.0)){
            /* Illinois regula falsi: the residual is smooth inside this narrow bracket.
               The weights wa and wb are the end residuals scaled down for
               interpolation only; a and b keep the evaluated residuals. */
            FluxPoint a = near, b = far;
            double wa = a.residual, wb = b.residual;
            int side = 0;
            for(int j=0;j<100;j++){
                double x = (a.fraction*wb - b.fraction*wa)/(wb - wa);
                const double low = fmin(a.fraction,b.fraction), high = fmax(a.fraction,b.fraction);
                if(!(x > low && x < high)) x = low + 0.5*(high-low);
                if(!(x > low && x < high)) break;
                const FluxPoint c = {x,HydrogenFluxResidual(x,balance)};
                if(!isfinite(c.residual)) return EQC_INVALID_INPUT;
                if(fabs(c.residual) <= accuracy){ *root = x; return EQC_SUCCESS; }
                if((c.residual < 0.0) == (b.residual < 0.0)){
                    b = c; wb = c.residual;
                    if(side == 1) wa *= 0.5;
                    side = 1;
                }else{
                    a = c; wa = c.residual;
                    if(side == -1) wb *= 0.5;
                    side = -1;
                }
            }
            const FluxPoint lower = a.fraction < b.fraction ? a : b;
            const FluxPoint upper = a.fraction < b.fraction ? b : a;
            return SolveFluxBracket(HydrogenFluxResidual,balance,lower,upper,accuracy,-1.0,root);
        }
        if(far.fraction >= 0.5) return EQC_NO_EQUILIBRIUM;
        near = far;
#ifdef EQC_H2_LEGACY_SHIELDING
        step *= 4.0;
#else
        /* The line-shielding fit can support multiple equilibria. Keep the
           bracket search local as the step grows. */
        step = BelowTraceShielding(balance,boundary,near.fraction) ? 4.0*step : fmin(4.0*step,0.25);
#endif
    }
    return EQC_NO_EQUILIBRIUM;
}

/* Unshielded balance 0.5 y (A + B a + C y) = F a, quadratic in y. */
static double UnshieldedMolecularFraction(const StructEQCData *state, const CEQMolecularHydrogenCoefficients *h2,
                                          double neutral, double byElectron, double byAtom, double byMolecule){
    const double photo = H2PhotoDissociationRateWithDustFactor(state->G0,0.0,h2->doppler,h2->dustShielding);
    const double formation = h2->formation*neutral;
    const double q2 = 0.5*(byMolecule-byAtom);
    const double q1 = 0.5*(photo+h2->cosmicRayDissociation+byElectron+byAtom) + formation;
    const double discriminant = q1*q1 + 4.0*q2*formation;
    const double guess = discriminant >= 0.0 ? (2.0*formation)/(q1+sqrt(discriminant)) : 0.5;
    return guess > 0.0 && guess < 1.0 ? guess : 0.5;
}

/* Newton iteration on r = ln(destroyed/formed) in the logit z = ln(y/a) of
   the H2 fraction y and atomic fraction a = 1-y. The logit step updates y and
   a without forming 1-y, so trace H2 and trace H I both stay accurate. The
   convergence test is the bounded flux residual of the bracket solver.
   Returns EQC_NO_CONVERGENCE to leave the bracket solver to the caller. */
static int SolveMolecularHydrogenNewton(const StructEQCData *state, const CEQMolecularHydrogenCoefficients *h2,
                                        double neutral, double previousH2, double accuracy,
                                        double *atomic, double *molecular){
    const double nH = state->n[ipHYDROGEN];
    const double byElectron = state->nel*h2->collision.electron;
    const double byAtom = (nH*neutral)*h2->collision.atomicHydrogen;     /* times a */
    const double byMolecule = (0.5*nH)*h2->collision.molecularHydrogen;  /* times y */
    const double formation = h2->formation*(nH*neutral);                 /* times a */
    if(!(formation > 0.0) || !isfinite(formation)) return EQC_NO_CONVERGENCE;

    double y = 2.0*(previousH2/nH);
    if(!(y > 0.0 && y < 1.0))
        y = UnshieldedMolecularFraction(state,h2,neutral,byElectron,byAtom,byMolecule);
    double a = 1.0-y;
    for(int iteration=0; iteration<16; iteration++){
        const double column = (0.5*y)*state->Nt;
        /* Match the bracket residual's evaluation order. */
        const double nHI = (nH*a)*neutral;
        const double nH2 = (0.5*nH)*y;
        double shieldSlope;
        const double photo = CEQH2PhotoRateAndShieldingSlope(state->G0,column,h2->doppler,h2->dustShielding,
                                                             &shieldSlope);
        const double collision = H2DissociationFromCoefficients(h2->collision,state->nel,nHI,nH2);
        const double destruction = photo+h2->cosmicRayDissociation+collision;
        const double destroyed = destruction*nH2;
        const double formed = h2->formation*nHI;
        if(!(destroyed >= DBL_MIN && formed >= DBL_MIN)) return EQC_NO_CONVERGENCE;
        const double residual = FluxResidual(destroyed,formed);
        if(!isfinite(residual)) return EQC_NO_CONVERGENCE;
        /* dr/dz = 1 + y a (d destruction/dy)/destruction.
           The positive derivative selects a stable chemical balance. */
        const double slope = 1.0 + (a*photo*shieldSlope + (y*a)*(byMolecule-byAtom))/destruction;
        if(!(slope > 0.0)) return EQC_NO_CONVERGENCE;
        if(fabs(residual) <= accuracy){
            *atomic = a;
            *molecular = y;
            return EQC_SUCCESS;
        }
        const double step = fmax(-8.0,fmin(8.0,-log(destroyed/formed)/slope));
        const double grown = y*exp(step);
        const double total = a+grown;
        y = grown/total;
        a = a/total;
        if(!(y > 0.0 && a > 0.0)) return EQC_NO_CONVERGENCE;
    }
    return EQC_NO_CONVERGENCE;
}

/* Split the H nuclei between H2 and atomic H for a given neutral fraction of
   the atomic H. previousH2 seeds the search. */
static int SolveMolecularHydrogen(const StructEQCData *state, const CEQMolecularHydrogenCoefficients *h2,
                                  double neutral, double previousH2, double accuracy,
                                  double *atomic, double *molecular){
    /* Keep the original iteration at tolerances sensitive to roundoff. */
    if(accuracy > 1.e-12 &&
       SolveMolecularHydrogenNewton(state,h2,neutral,previousH2,accuracy,atomic,molecular) == EQC_SUCCESS)
        return EQC_SUCCESS;
    const double nH = state->n[ipHYDROGEN];
    HydrogenBalance balance = {.state = state, .h2 = h2, .neutralFraction = neutral};
    const double half = HydrogenFluxResidual(0.5,&balance);
    if(!isfinite(half)) return EQC_INVALID_INPUT;
    double fraction;
    int status;
    if(half == 0.0 && HydrogenFluxResidual(0.0,&balance) == 0.0 &&
       HydrogenFluxResidual(1.0,&balance) == 0.0){
        // With both reaction fluxes zero, continue the previous abundance.
        fraction = 2.0*(previousH2/nH);
        if(!isfinite(fraction)) return EQC_INVALID_INPUT;
        fraction = fmin(1.0,fmax(0.0,fraction));
        status = EQC_SUCCESS;
    }else{
        // Use the smaller fraction so that nearly molecular gas retains its trace H I.
        const bool bracketAtomic = half < 0.0;
        double molecularGuess = 2.0*(previousH2/nH);
        balance.useAtomicFraction = bracketAtomic;
#ifndef EQC_H2_LEGACY_SHIELDING
        /* A branch can lie in either half, irrespective of the sign at 0.5.
           Search the incoming branch before bracketing from the endpoints. */
        if(molecularGuess == 0.0){
            molecularGuess = UnshieldedMolecularFraction(state,h2,neutral,state->nel*h2->collision.electron,
                (nH*neutral)*h2->collision.atomicHydrogen,(0.5*nH)*h2->collision.molecularHydrogen);
        }
        if(molecularGuess > 0.0 && molecularGuess < 1.0)
            balance.useAtomicFraction = molecularGuess > 0.5;
#endif
        double initialFraction = 0.0;
#ifdef EQC_H2_LEGACY_SHIELDING
        if(previousH2 > 0.0 && previousH2 < nH)
#else
        if(molecularGuess > 0.0 && molecularGuess < 1.0)
#endif
            initialFraction = balance.useAtomicFraction ? 1.0-molecularGuess : molecularGuess;
        status = EQC_NO_EQUILIBRIUM;
        if(initialFraction > 0.0 && initialFraction < 0.5)
            status = SolveMolecularHydrogenNear(&balance,initialFraction,half,accuracy,&fraction);
#ifndef EQC_H2_LEGACY_SHIELDING
        if(status == EQC_SUCCESS &&
           !StableMolecularHydrogen(&balance,balance.useAtomicFraction ? fraction : 1.0-fraction,
                                    balance.useAtomicFraction ? 1.0-fraction : fraction)){
            status = EQC_NO_EQUILIBRIUM;
            initialFraction = 0.0;
        }
        if(balance.useAtomicFraction != bracketAtomic && status != EQC_SUCCESS){
            balance.useAtomicFraction = bracketAtomic;
            initialFraction = 0.0;
        }
#endif
        if(status != EQC_SUCCESS){
            const FluxPoint lower = {0.0,HydrogenFluxResidual(0.0,&balance)};
            const FluxPoint upper = {0.5,half};
            status = SolveFluxBracket(HydrogenFluxResidual,&balance,lower,upper,accuracy,initialFraction,&fraction);
        }
    }
    if(status != EQC_SUCCESS) return status;
    *atomic = balance.useAtomicFraction ? fraction : 1.0-fraction;
    *molecular = balance.useAtomicFraction ? 1.0-fraction : fraction;
#ifndef EQC_H2_LEGACY_SHIELDING
    if(!StableMolecularHydrogen(&balance,*atomic,*molecular)) return EQC_NO_CONVERGENCE;
#endif
    return EQC_SUCCESS;
}

static CEQLadderCoefficients CarbonPhotoCoefficients(const StructEQCData *s, const CEQUVBRates *uvb){
    /* Without H II region mode, NLyn is ignored and C keeps its FUV photoionization. */
    const double NLyn = EQCLibRunParameters.UseHIIRegion ? s->NLyn : 0.0;
    return (CEQLadderCoefficients){
        .photo0 = (NLyn == 0.0 ? CIFUVPhotoionizationRate(s->G0,VisualExtinction(s->Nt,s->Zd))
                               : CIPhotoionizationRate(NLyn)) + uvb->metalIonization[CEQ_UVB_CI],
        .photo1 = CIIPhotoionizationRate(NLyn) + uvb->metalIonization[CEQ_UVB_CII],
    };
}

static void RefreshCarbonCoefficients(double T, CEQLadderCoefficients *carbon, const CEQTemperaturePoint *point){
    carbon->collision0 = ChemistryCoefficient(TemperatureIonizeCI,T,point);
    carbon->recombination1 = ChemistryCoefficient(TemperatureRecombineCII,T,point);
    carbon->collision1 = ChemistryCoefficient(TemperatureIonizeCII,T,point);
    carbon->recombination2 = ChemistryCoefficient(TemperatureRecombineCIII,T,point);
}

static CEQLadderCoefficients CarbonCoefficients(const StructEQCData *s, const CEQUVBRates *uvb){
    CEQLadderCoefficients carbon = CarbonPhotoCoefficients(s,uvb);
    RefreshCarbonCoefficients(s->T,&carbon,NULL);
    return carbon;
}

static double COFormationCoefficient(double T){
    /* Radiative association forming CH2+, k148 in Glover et al. (2010). */
    return (4.e-16*pow(T/300.0,-0.2))*EQCLibRunParameters.CpForCO;
}

static CEQMolecularCarbonCoefficients MolecularCarbonCoefficients(const StructEQCData *s, const CEQUVBRates *uvb){
    return (CEQMolecularCarbonCoefficients){
        .formationPerH2 = COFormationCoefficient(s->T),
        .photoCI = 3.1e-10*s->G0*DustShieldingFactor(s->Nt,s->Zd,3.0) + uvb->metalIonization[CEQ_UVB_CI],
        .photoCII = uvb->metalIonization[CEQ_UVB_CII],
        .photoCO = COPhotoDissociationRate(s->G0,s->Nt,s->Zd),
    };
}

static void InitializeRateCoefficients(const StructEQCData *s, CEQRateCoefficients *c, bool equilibrium,
                                        const CEQUVBRates *uvb){
    const double NLyn = s->NLyn;
    *c = (CEQRateCoefficients){
        .hydrogen.photo0 = HIPhotoIonizationRate(NLyn) + uvb->ionization[CEQ_UVB_HI],
        .hasMolecularHydrogen = EQCLibRunParameters.UseH2Formation,
    };
    if(equilibrium || s->n[ipHELIUM] > 0.0){
        c->helium.photo0 = HeIPhotoIonizationRate(NLyn) + uvb->ionization[CEQ_UVB_HEI];
        c->helium.photo1 = HeIIPhotoIonizationRate(NLyn) + uvb->ionization[CEQ_UVB_HEII];
    }
    if(equilibrium || s->n[ipCARBON] > 0.0) c->carbon = CarbonPhotoCoefficients(s,uvb);
    if(equilibrium || s->n[ipOXYGEN] > 0.0){
        c->oxygen.photo0 = OIPhotoionizationRate(NLyn) + uvb->metalIonization[CEQ_UVB_OI];
        c->oxygen.photo1 = OIIPhotoionizationRate(NLyn) + uvb->metalIonization[CEQ_UVB_OII];
    }
    if(equilibrium){
        c->silicon.photo0 = SiIPhotoionizationRate(NLyn) + uvb->metalIonization[CEQ_UVB_SI_I];
        c->silicon.photo1 = SiIIPhotoionizationRate(NLyn) + uvb->metalIonization[CEQ_UVB_SI_II];
    }
    if(c->hasMolecularHydrogen) c->molecularHydrogen = MolecularHydrogenRadiationCoefficients(s);
}

static void RefreshHydrogenCoefficients(double T, CEQLadderCoefficients *hydrogen, const CEQTemperaturePoint *point){
    const bool caseB = EQCLibRunParameters.UseRecombinationCaseB;
    hydrogen->collision0 = ChemistryCoefficient(TemperatureIonizeHI,T,point);
    hydrogen->recombination1 = ChemistryCoefficient(caseB ? TemperatureRecombineHIICaseB : TemperatureRecombineHIICaseA,
                                                    T,point);
}

static void RefreshRateCoefficients(const StructEQCData *s, CEQRateCoefficients *c, bool equilibrium,
                                    const CEQTemperaturePoint *point){
    const double T = s->T;
    RefreshHydrogenCoefficients(T,&c->hydrogen,point);
    if(equilibrium || s->n[ipHELIUM] > 0.0){
        const bool caseB = EQCLibRunParameters.UseRecombinationCaseB;
        c->helium.collision0 = ChemistryCoefficient(TemperatureIonizeHeI,T,point);
        c->helium.recombination1 = ChemistryCoefficient(caseB ? TemperatureRecombineHeIICaseB : TemperatureRecombineHeIICaseA,
                                                        T,point);
        c->helium.collision1 = ChemistryCoefficient(TemperatureIonizeHeII,T,point);
        c->helium.recombination2 = ChemistryCoefficient(caseB ? TemperatureRecombineHeIIICaseB : TemperatureRecombineHeIIICaseA,
                                                        T,point);
    }
    if(equilibrium || s->n[ipCARBON] > 0.0) RefreshCarbonCoefficients(T,&c->carbon,point);
    if(equilibrium || s->n[ipOXYGEN] > 0.0){
        c->oxygen.collision0 = ChemistryCoefficient(TemperatureIonizeOI,T,point);
        c->oxygen.recombination1 = ChemistryCoefficient(TemperatureRecombineOII,T,point);
        c->oxygen.collision1 = ChemistryCoefficient(TemperatureIonizeOII,T,point);
        c->oxygen.recombination2 = ChemistryCoefficient(TemperatureRecombineOIII,T,point);
        if(equilibrium || EQCLibRunParameters.UseChargeTransfer){
            c->toOxygenII = ChemistryCoefficient(TemperatureTransferToOII,T,point);
            c->toHydrogenI = ChemistryCoefficient(TemperatureTransferToHI,T,point);
            c->oxygenIIIToHydrogenI = ChemistryCoefficient(TemperatureTransferOIII,T,point);
        }
    }
    if(equilibrium){
        c->silicon.collision0 = SiICollisionalIonizationRate(T);
        c->silicon.recombination1 = SiIIRecombinationRate(T);
        c->silicon.collision1 = SiIICollisionalIonizationRate(T);
        c->silicon.recombination2 = SiIIIRecombinationRate(T);
    }
    if(c->hasMolecularHydrogen) RefreshMolecularHydrogenCoefficients(s,&c->molecularHydrogen,point);
}

void CEQPrepareRateCoefficients(const StructEQCData *s, CEQRateCoefficients *c, const CEQUVBRates *uvb){
    InitializeRateCoefficients(s,c,true,uvb);
    RefreshRateCoefficients(s,c,true,NULL);
    if(EQCLibRunParameters.UseCOFormation && s->T <= 1.e4 && s->n[ipCARBON] > 0.0)
        c->molecularCarbon = MolecularCarbonCoefficients(s,uvb);
}

void CEQInitializeNonEquilibriumRateCoefficients(const StructEQCData *s, CEQRateCoefficients *c,
                                                 const CEQUVBRates *uvb){
    InitializeRateCoefficients(s,c,false,uvb);
}

void CEQRefreshNonEquilibriumRateCoefficients(const StructEQCData *s, CEQRateCoefficients *c,
                                              const CEQTemperaturePoint *point){
    RefreshRateCoefficients(s,c,false,point);
}

CEQFixedChargeLadder CEQHydrogenLadderAtFixedCharge(const StructEQCData *s, const CEQRateCoefficients *c,
                                                   double softXrayRate){
    const double xe = s->nel/s->n[ipHYDROGEN];
    return (CEQFixedChargeLadder){
        .collision0 = s->nel*c->hydrogen.collision0,
        .external0 = s->X0*(CosmicRayIonizationRateHI(xe) + softXrayRate),
        .recombination1 = s->nel*c->hydrogen.recombination1
                        + s->n[ipHYDROGEN]*HIIRecombinationRateOnDust(s->T,s->G0,s->nel,s->Zd),
    };
}

CEQHydrogenRates CEQHydrogenIonizationRatesWithCoefficients(const StructEQCData *s,
                                                            const CEQRateCoefficients *c){
    const double softXrayRate = SoftXrayIonizationRateHI(s->nel/s->n[ipHYDROGEN],s->Nt);
    return CEQEvaluateFixedChargeHydrogenLadder(CEQHydrogenLadderAtFixedCharge(s,c,softXrayRate),
                                                c->hydrogen.photo0);
}

CEQHydrogenRates CEQHydrogenIonizationRates(const StructEQCData *s){
    CEQUVBRates uvb;
    if(CEQEvaluateUVBackgroundWithMetallicity(s->NHIuvb,s->Z,&uvb) != EQC_SUCCESS)
        return (CEQHydrogenRates){.ionization = NAN, .recombination = NAN};
    CEQRateCoefficients c = {.hydrogen.photo0 = HIPhotoIonizationRate(s->NLyn) + uvb.ionization[CEQ_UVB_HI]};
    const double power = CEQHighPhotoelectronPower(s,&uvb);
    if(power > 0.0)
        c.hydrogen.photo0 += CEQDepositPhotoelectrons(s,0.0,power).ionizationHI;
    RefreshHydrogenCoefficients(s->T,&c.hydrogen,NULL);
    return CEQHydrogenIonizationRatesWithCoefficients(s,&c);
}

CEQFixedChargeLadder CEQHeliumLadderAtFixedCharge(const StructEQCData *state, const CEQRateCoefficients *c,
                                                 double softXrayRate){
    const double ne = state->nel;
    const double nH = state->n[ipHYDROGEN];
    const double xe = ne/nH;
    return (CEQFixedChargeLadder){
        .collision0 = ne*c->helium.collision0,
        .external0 = state->X0*(CosmicRayIonizationRateHeI(xe) + softXrayRate),
        .recombination1 = ne*c->helium.recombination1
                        + nH*HeIIRecombinationRateOnDust(state->T,state->G0,ne,state->Zd),
        .ionization1 = c->helium.photo1 + ne*c->helium.collision1,
        .recombination2 = ne*c->helium.recombination2,
    };
}

CEQIonLadderRates CEQHeliumIonizationRates(const StructEQCData *state, const CEQRateCoefficients *c){
    const double softXrayRate = SoftXrayIonizationRateHeI(state->nel/state->n[ipHYDROGEN],state->Nt);
    return CEQEvaluateFixedChargeLadder(CEQHeliumLadderAtFixedCharge(state,c,softXrayRate),c->helium.photo0);
}

CEQIonLadderRates CEQCarbonIonizationRates(const StructEQCData *state, const CEQRateCoefficients *coefficients){
    CEQLadderCoefficients c;
    if(coefficients != NULL) c = coefficients->carbon;
    else{
        CEQUVBRates uvb;
        if(CEQEvaluateUVBackgroundWithMetallicity(state->NHIuvb,state->Z,&uvb) != EQC_SUCCESS)
            return (CEQIonLadderRates){NAN,NAN,NAN,NAN};
        c = CarbonCoefficients(state,&uvb);
    }
    const double ne = state->nel;
    const double nH = state->n[ipHYDROGEN];
    const double xe = ne/nH;
    return (CEQIonLadderRates){
        .ionization0 = c.photo0 + ne*c.collision0
                     + state->X0*(CosmicRayIonizationRateCI(xe) + SoftXrayIonizationRateCI(xe,state->Nt)),
        .recombination1 = ne*c.recombination1
                        + nH*CIIRecombinationRateOnDust(state->T,state->G0,ne,state->Zd),
        .ionization1 = c.photo1 + ne*c.collision1,
        .recombination2 = ne*c.recombination2,
    };
}

/* Oxygen ladder rates without charge transfer, in s^-1 per ion or atom. */
CEQIonLadderRates CEQOxygenIonizationRates(const StructEQCData *state, const CEQRateCoefficients *c){
    const double ne = state->nel;
    const double nH = state->n[ipHYDROGEN];
    return (CEQIonLadderRates){
        .ionization0 = c->oxygen.photo0 + ne*c->oxygen.collision0
                     + state->X0*CosmicRayIonizationRateOI(ne/nH),
        .recombination1 = ne*c->oxygen.recombination1
                        + nH*OIIRecombinationRateOnDust(state->T,state->G0,ne,state->Zd),
        .ionization1 = c->oxygen.photo1 + ne*c->oxygen.collision1,
        .recombination2 = ne*c->oxygen.recombination2,
    };
}

typedef struct {
    StructEQCData *state;
    const CEQIonLadderRates *oxygen;
    const CEQMolecularHydrogenCoefficients *h2;   /* NULL without H2 formation. */
    double atomicHydrogen, freeOxygen, accuracy;
    double ionization, recombination;   /* H rates per atom without charge transfer. */
    double toOxygenII, toHydrogenI, oxygenIIIToHydrogenI;
    double startS, startAtomic, startMolecular;   /* H2 split already solved at startS. */
    int status;
} ChargeTransferBalance;

/* Set H I/H II for s = ln(n_HII/n_HI), with the H2 split re-solved for that
   neutral fraction, then solve the O ladder with charge transfer and return
   s - ln(I/R) for the resulting per-atom H rates. */
static double ChargeTransferResidual(double s, ChargeTransferBalance *b){
    StructEQCData *state = b->state;
    double atomicHydrogen = b->atomicHydrogen;
    if(b->h2 != NULL){
        double atomic = b->startAtomic, molecular = b->startMolecular;
        const double neutral = 1.0/(1.0+exp(s));
        if(s != b->startS){
            const int status = SolveMolecularHydrogen(state,b->h2,neutral,state->nmol[ipH2],
                                                      b->accuracy,&atomic,&molecular);
            if(status != EQC_SUCCESS){
                b->status = status;
                return NAN;
            }
        }
        const double nH = state->n[ipHYDROGEN];
        atomicHydrogen = nH*atomic;
        state->nmol[ipH2] = (0.5*nH)*molecular;
    }
    state->nion[ipHYDROGEN] = atomicHydrogen/(1.0+exp(-s));
    state->natom[ipHYDROGEN] = atomicHydrogen/(1.0+exp(s));
    const double nHI = state->natom[ipHYDROGEN], nHII = state->nion[ipHYDROGEN];
    const CEQIonLadderRates *o = b->oxygen;
    if(CEQSolveIonLadder(b->freeOxygen, o->ionization0 + nHII*b->toOxygenII,
                         o->recombination1 + nHI*b->toHydrogenI, o->ionization1,
                         o->recombination2 + nHI*b->oxygenIIIToHydrogenI,
                         &state->natom[ipOXYGEN], &state->nion[ipOXYGEN], &state->nOIII) != EQC_SUCCESS)
        return NAN;
    const double ionization = b->ionization + state->nion[ipOXYGEN]*b->toHydrogenI
                            + state->nOIII*b->oxygenIIIToHydrogenI;
    const double recombination = b->recombination + state->natom[ipOXYGEN]*b->toOxygenII;
    return s - (log(ionization) - log(recombination));
}

static double AvailableAtomicHydrogen(const StructEQCData *state){
    const double nH = state->n[ipHYDROGEN], nH2 = state->nmol[ipH2];
    if(!isfinite(nH) || nH <= 0.0 || !isfinite(nH2) || nH2 < 0.0 || nH2 > 0.5*nH) return NAN;
    const double atomic = nH-2.0*nH2;
    return atomic >= 0.0 ? atomic : NAN;
}

static int SolveHydrogenOxygenChargeTransfer(StructEQCData *state, const CEQHydrogenRates *hydrogen,
                                             const CEQIonLadderRates *oxygen, const CEQRateCoefficients *coefficients,
                                             double accuracy, bool fixedH2){
    const double freeOxygen = state->n[ipOXYGEN]-state->nmol[ipCO];
    const double atomicHydrogen = fixedH2 ? AvailableAtomicHydrogen(state)
                                          : state->natom[ipHYDROGEN]+state->nion[ipHYDROGEN];
    ChargeTransferBalance b = {
        .state = state, .oxygen = oxygen, .atomicHydrogen = atomicHydrogen, .freeOxygen = freeOxygen,
        .ionization = hydrogen->ionization, .recombination = hydrogen->recombination,
        .toOxygenII = coefficients->toOxygenII,
        .toHydrogenI = coefficients->toHydrogenI,
        .oxygenIIIToHydrogenI = coefficients->oxygenIIIToHydrogenI,
        .accuracy = accuracy,
        .status = EQC_INVALID_INPUT,
    };
    if(!isfinite(freeOxygen) || freeOxygen < 0.0 || !isfinite(atomicHydrogen) ||
       !isfinite(b.ionization) || !isfinite(b.recombination)) return EQC_INVALID_INPUT;

    if(freeOxygen == 0.0 || atomicHydrogen == 0.0 || b.ionization == 0.0 || b.recombination == 0.0){
        /* Without both H I and H II, charge transfer only acts on oxygen; keep the H ladder. */
        const double nHI = state->natom[ipHYDROGEN], nHII = state->nion[ipHYDROGEN];
        return CEQSolveIonLadder(freeOxygen, oxygen->ionization0 + nHII*b.toOxygenII,
                                 oxygen->recombination1 + nHI*b.toHydrogenI, oxygen->ionization1,
                                 oxygen->recombination2 + nHI*b.oxygenIIIToHydrogenI,
                                 &state->natom[ipOXYGEN], &state->nion[ipOXYGEN], &state->nOIII);
    }

    /* Charge transfer adds at most these rates to the H ladder, so the root
       lies in [lower, upper], where the residual is <= 0 and >= 0. */
    double lower = log(b.ionization) - log(b.recombination + freeOxygen*b.toOxygenII);
    double upper = log(b.ionization + freeOxygen*fmax(b.toHydrogenI,b.oxygenIIIToHydrogenI))
                 - log(b.recombination);
    double lowerResidual = NAN, upperResidual = NAN;

    /* Start from the root without charge transfer. The map s -> ln(I/R) changes
       slowly when O/H is small, so a fixed-point step followed by secant steps
       converges in a few evaluations; the bracket keeps every step bounded. */
    double s = log(b.ionization) - log(b.recombination);
    if(EQCLibRunParameters.UseH2Formation && !fixedH2){
        /* The incoming state holds an H2 split already solved for its own H
           ionization ratio: the uncoupled root after the H solve, or the
           previous coupled root when the caller repeats this after updating CO.
           Start there so the first evaluation reuses that split. */
        const double nH = state->n[ipHYDROGEN];
        b.h2 = coefficients->hasMolecularHydrogen ? &coefficients->molecularHydrogen : NULL;
        if(b.h2 == NULL) return EQC_INVALID_INPUT;
        b.startS = log(state->nion[ipHYDROGEN]) - log(state->natom[ipHYDROGEN]);
        b.startAtomic = atomicHydrogen/nH;
        b.startMolecular = 2.0*(state->nmol[ipH2]/nH);
        if(b.startS > lower && b.startS < upper) s = b.startS;
    }
    double previous = NAN, previousResidual = NAN;
    for(int iteration = 0; iteration < 200; iteration++){
        if(!(s > lower && s < upper)){
            s = isfinite(lowerResidual) && isfinite(upperResidual) ?
                (lower*upperResidual - upper*lowerResidual)/(upperResidual - lowerResidual) :
                lower + 0.5*(upper-lower);
            if(!(s > lower && s < upper)) s = lower + 0.5*(upper-lower);
        }
        if(iteration > 0 && (s <= lower || s >= upper)) break;
        const double residual = ChargeTransferResidual(s,&b);
        if(!isfinite(residual)) return b.status;
        if(fabs(residual) <= accuracy) return EQC_SUCCESS;
        if(residual < 0.0){ lower = s; lowerResidual = residual; }
        else              { upper = s; upperResidual = residual; }

        double next = s - residual;   /* Fixed-point step s <- ln(I/R). */
        if(isfinite(previous) && residual != previousResidual)
            next = s - residual*(s-previous)/(residual-previousResidual);
        if(iteration % 4 == 3) next = lower + 0.5*(upper-lower);
        previous = s; previousResidual = residual;
        s = next;
    }
    return EQC_NO_CONVERGENCE;
}

int CEQSolveHydrogenOxygenChargeTransfer(StructEQCData *state, const CEQHydrogenRates *hydrogen,
                                         const CEQIonLadderRates *oxygen,
                                         const CEQRateCoefficients *coefficients, double accuracy){
    return SolveHydrogenOxygenChargeTransfer(state,hydrogen,oxygen,coefficients,accuracy,false);
}

int CEQSolveHydrogenOxygenChargeTransferWithFixedH2(StructEQCData *state, const CEQHydrogenRates *hydrogen,
                                                    const CEQIonLadderRates *oxygen,
                                                    const CEQRateCoefficients *coefficients, double accuracy){
    return SolveHydrogenOxygenChargeTransfer(state,hydrogen,oxygen,coefficients,accuracy,true);
}

int CEQSolveHydrogenWithFixedH2(StructEQCData *state, const CEQHydrogenRates *rates){
    const double atomic = AvailableAtomicHydrogen(state);
    double unused;
    return CEQSolveIonLadder(atomic,rates->ionization,rates->recombination,0.0,0.0,
                            &state->natom[ipHYDROGEN],&state->nion[ipHYDROGEN],&unused);
}

int CEQSolveHydrogenEquilibrium(StructEQCData *state, double accuracy){
    const double nH = state->n[ipHYDROGEN];
    if(!isfinite(nH) || nH <= 0.0) return EQC_INVALID_INPUT;
    const CEQHydrogenRates rates = CEQHydrogenIonizationRates(state);
    return CEQSolveHydrogenEquilibriumWithRates(state,&rates,NULL,accuracy);
}

int CEQSolveHydrogenEquilibriumWithRates(StructEQCData *state, const CEQHydrogenRates *rates,
                                         const CEQRateCoefficients *coefficients, double accuracy){
    const double nH = state->n[ipHYDROGEN];
    if(!isfinite(nH) || nH <= 0.0 || !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0){
        return EQC_INVALID_INPUT;
    }
    double neutral, ion, unused;
    int status = CEQSolveIonLadder(1.0,rates->ionization,rates->recombination,0.0,0.0,&neutral,&ion,&unused);
    if(status != EQC_SUCCESS) return status;

    double atomic = 1.0;
    double molecular = 0.0;
    if(EQCLibRunParameters.UseH2Formation){
        const CEQMolecularHydrogenCoefficients h2 = coefficients != NULL && coefficients->hasMolecularHydrogen ?
            coefficients->molecularHydrogen : CEQPrepareMolecularHydrogenCoefficients(state);
        status = SolveMolecularHydrogen(state,&h2,neutral,state->nmol[ipH2],accuracy,&atomic,&molecular);
        if(status != EQC_SUCCESS) return status;
    }
    state->natom[ipHYDROGEN] = (nH*atomic)*neutral;
    state->nion[ipHYDROGEN] = (nH*atomic)*ion;
    state->nmol[ipH2] = (0.5*nH)*molecular;
    return EQC_SUCCESS;
}

void EquilibriumHydrogenFractionalAbundance(StructEQCData *PhydroPP){
    CEQUVBRates uvb;
    if(CEQEvaluateUVBackgroundWithMetallicity(PhydroPP->NHIuvb,PhydroPP->Z,&uvb) != EQC_SUCCESS) return;

    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nHI    = PhydroPP->natom[ipHYDROGEN];
    double nH2    = PhydroPP->nmol [ipH2];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn;
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Nh2    = Ntot*(nH2/nH);
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;
    double Tdust  = EffectiveDustTemperature(PhydroPP);
    //double b_turb = TURBULENCE_VELOCITY;    // km/s
    double b5     = H2DopplerWidth(T);

    double dnHI_phot  = HIPhotoIonizationRate(NLyn) + uvb.ionization[CEQ_UVB_HI];
    const double power = CEQHighPhotoelectronPower(PhydroPP,&uvb);
    if(power > 0.0)
        dnHI_phot += CEQDepositPhotoelectrons(PhydroPP,0.0,power).ionizationHI;
    double dnHI_coli  = ne*HICollisionalIonizationRate(T);
    double dnHI_cr    = X0*CosmicRayIonizationRateHI(xe);
    double dnHI_xr    = X0*SoftXrayIonizationRateHI(xe,Ntot);
    double dnHII_rec  = ne*HIIRecombinationRate(T);
    double dnHII_dust = nH*HIIRecombinationRateOnDust(T,G0,ne,Zd); 
    double HIItoHI    = (dnHI_phot + dnHI_coli + dnHI_cr + dnHI_xr)/(dnHII_rec + dnHII_dust);
    double fHII, fHI, fH2;
    if(EQCLibRunParameters.UseH2Formation==true){
        double dnH2_form  = nH*H2FormationRate(T,Zd,Tdust);
        double dnH2_diss  = H2PhotoDissociationRate(G0,Ntot,Nh2,b5,Zd);
        double dnH2_dest  = H2CollisionalDissociationRate(T,nH,nHI,nH2,ne);
        double H2toHI     = dnH2_form/(dnH2_diss + CEQH2CosmicRayDissociationRate(X0) + dnH2_dest);
        double HIItoHt    = MIN(1.0,HIItoHI/(HIItoHI+2*H2toHI+1.0));
        double HItoHt     = MIN(1.0,1.0    /(HIItoHI+2*H2toHI+1.0));
        fHII = HIItoHt;
        fHI  = HItoHt;
        fH2  = 0.5*MAX(0.0,1.e0 - fHII - fHI);
    }else{
        double HIItoHt    = MIN(1.0,HIItoHI/(HIItoHI+1.0));
        fHII = HIItoHt;
        fHI  = MAX(0.0,1.e0 - fHII);
        fH2  = 0.0;
    }
    if(fHII<0||fHI<0||fH2<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"  fHII=%e fHI=%e fH2=%e\n",fHII,fHI,fH2);
        exit(1);
    }
	PhydroPP->nion[ipHYDROGEN]  = nH*fHII;
    PhydroPP->natom[ipHYDROGEN] = nH*fHI ;
    PhydroPP->nmol[ipH2]        = nH*fH2 ;
    return ;
}

//void EquilibriumHeliumFractionalAbundance(double T, double nH, double ne, double G0, double Ntot, double Z, double *fHeII, double *fHeI){
//    double xe          = ne/nH;
//    double X0          = G0;
//    double dnHeI_coli  = ne*HeICollisionalIonizationRate(T);
//    double dnHeI_cr    = X0*CosmicRayIonizationRateHeI(xe);
//    double dnHeI_xr    = X0*SoftXrayIonizationRateHeI(xe,Ntot);
//    double dnHeII_rec  = ne*(HeIIRecombinationRateCaseA(T) + HeIIDielectronicRecombinationRate(T));
//    double dnHeII_dust = nH*HeIIRecombinationRateOnDust(T,G0,ne,Z); 
//
//    double HeIItoHeI   = (dnHeI_coli + dnHeI_cr + dnHeI_xr)/(dnHeII_rec + dnHeII_dust);
//    *fHeII = HeIItoHeI/(1.0+HeIItoHeI);
//    *fHeI  = MAX(0.0,1.0 - (*fHeII));
//    if(*fHeII<0||*fHeI<0){
//        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
//        fprintf(stderr,"  fHeII=%e fHeI=%e\n",*fHeII,*fHeI);
//        exit(1);
//    }
//    
//    return ;
//}

void EquilibriumHeliumFractionalAbundance(StructEQCData *PhydroPP){
    CEQUVBRates uvb;
    if(CEQEvaluateUVBackgroundWithMetallicity(PhydroPP->NHIuvb,PhydroPP->Z,&uvb) != EQC_SUCCESS) return;
    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nHe    = PhydroPP->n[ipHELIUM];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn; 
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;
    
    double dnHeI_phti   = HeIPhotoIonizationRate(NLyn) + uvb.ionization[CEQ_UVB_HEI];
    const double power = CEQHighPhotoelectronPower(PhydroPP,&uvb);
    if(power > 0.0)
        dnHeI_phti += CEQDepositPhotoelectrons(PhydroPP,0.0,power).ionizationHeI;
    double dnHeI_coli   = ne*HeICollisionalIonizationRate(T);
    double dnHeI_cr     = X0*CosmicRayIonizationRateHeI(xe);
    double dnHeI_xr     = X0*SoftXrayIonizationRateHeI(xe,Ntot);
    double dnHeII_rec   = ne*(HeIIRecombinationRate(T) + HeIIDielectronicRecombinationRate(T));
    double dnHeII_dust  = nH*HeIIRecombinationRateOnDust(T,G0,ne,Zd); 

    double dnHeII_phti  = HeIIPhotoIonizationRate(NLyn) + uvb.ionization[CEQ_UVB_HEII];
    double dnHeII_coli  = ne*HeIICollisionalIonizationRate(T);
    double dnHeII_cr    = 0.e0; // OK
    double dnHeII_xr    = 0.e0; // 0K
    double dnHeIII_rec  = ne*HeIIIRecombinationRate(T);
    double dnHeIII_dust = 0.e0; // OK
    
    double HeIItoHeI   = (dnHeI_phti  + dnHeI_coli  + dnHeI_cr  + dnHeI_xr )/(dnHeII_rec  + dnHeII_dust );
    double HeIIItoHeII = (dnHeII_phti + dnHeII_coli + dnHeII_cr + dnHeII_xr)/(dnHeIII_rec + dnHeIII_dust);
    double fHeIII = HeIIItoHeII/(HeIIItoHeII + 1.e0 + 1.e0/HeIItoHeI);
    double fHeII  =        1.e0/(HeIIItoHeII + 1.e0 + 1.e0/HeIItoHeI);
    double fHeI   = MAX(0.0,1.e0 - fHeII - fHeIII);
    if(fHeIII<0||fHeII<0||fHeI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        //fprintf(stderr,"  ne=%e nH=%e T=%e\n",ne,nH,T);
        //fprintf(stderr,"  G0=%e X0=%e NLyn=%e\n",G0,X0,PhydroPP->NLyn);
        //fprintf(stderr,"  Nt=%e Z=%e\n",Ntot,Z);
        //fprintf(stderr,"  dnHeI_phti=%e dnHeI_coli=%e dnHeI_cr=%e dnHeI_xr=%e\n",dnHeI_phti,dnHeI_coli,dnHeI_cr,dnHeI_xr);
        fprintf(stderr,"  fHeIII=%e fHeII=%e fHeI=%e\n",fHeIII,fHeII,fHeI);
        exit(1);
    }
    PhydroPP->nHeIII          = nHe*fHeIII;
    PhydroPP->nion[ipHELIUM]  = nHe*fHeII ;
    PhydroPP->natom[ipHELIUM] = nHe*fHeI  ;

    return ;
}


//// ionization state of oxygen (charge transfer)
//// ref. eq(2.38) in Osterbrock & Ferland (2006)
////      eq(14.32) in Draine (2011)
#define dE_OI_HI (0.19*EV_TO_ERG)
void EquilibriumOxygenFractionalAbundanceChargeTransfer(StructEQCData *PhydroPP){
    double kT   = BOLTZMANN_CONSTANT_CGS*PhydroPP->T;
    double nHI  = PhydroPP->nion [ipHYDROGEN];
    double nHII = PhydroPP->natom[ipHYDROGEN];
    double nO   = PhydroPP->n[ipOXYGEN];
    double nCO  = PhydroPP->nmol [ipCO];
    double RateRatio = 9.0*exp(dE_OI_HI/kT)/8.0;
    double fOIII = 0.e0;
    double fOII  = nHII/(nHII + RateRatio*nHI);
    double fOI   = MAX(0.0,1.0 - fOII);
    if(fOII<0||fOI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"fOII = %e fOI = %e T = %e nHI = %e nHII = %e\n",
                fOII,fOI,kT/BOLTZMANN_CONSTANT_CGS,nHI,nHII);
        exit(1);
    }
    PhydroPP->nOIII           = MAX(0.0,nO-nCO)*fOIII;
    PhydroPP->nion[ipOXYGEN]  = MAX(0.0,nO-nCO)*fOII ;
    PhydroPP->natom[ipOXYGEN] = MAX(0.0,nO-nCO)*fOI  ;
    return ;
}

void EquilibriumOxygenFractionalAbundance(StructEQCData *PhydroPP){
    CEQUVBRates uvb;
    if(CEQEvaluateUVBackgroundWithMetallicity(PhydroPP->NHIuvb,PhydroPP->Z,&uvb) != EQC_SUCCESS) return;

    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nO     = PhydroPP->n[ipOXYGEN];
    double nCO    = PhydroPP->nmol [ipCO];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn;
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;
    
    double dnOI_phti   = OIPhotoionizationRate(NLyn) + uvb.metalIonization[CEQ_UVB_OI];
    double dnOI_coli   = ne*OICollisionalIonizationRate(T);
    double dnOI_cr     = X0*CosmicRayIonizationRateOI(xe);
    double dnOI_xr     = 0.e0;
    double dnOII_rec   = ne*OIIRecombinationRate(T);
    double dnOII_dust  = nH*OIIRecombinationRateOnDust(T,G0,ne,Zd);
    double dnOI_chtr   = 0.e0; //nHII*OIChargeTransferRateToHII(T);
    double dnOII_chtr  = 0.e0; //nHI *OIIChargeTransferRateToHI(T);

    double dnOII_phti  = OIIPhotoionizationRate(NLyn) + uvb.metalIonization[CEQ_UVB_OII];
    double dnOII_coli  = ne*OIICollisionalIonizationRate(T);
    double dnOII_cr    = 0.e0;
    double dnOII_xr    = 0.e0;
    double dnOIII_rec  = ne*OIIIRecombinationRate(T);
    double dnOIII_dust = 0.e0;
    
    double OIItoOI   = (dnOI_phti  + dnOI_coli  + dnOI_cr  + dnOI_xr  + dnOI_chtr)/(dnOII_rec  + dnOII_dust + dnOII_chtr);
    double OIIItoOII = (dnOII_phti + dnOII_coli + dnOII_cr + dnOII_xr)/(dnOIII_rec + dnOIII_dust);
    double fOIII = OIIItoOII/(OIIItoOII + 1.e0 + 1.e0/OIItoOI);
    double fOII  =      1.e0/(OIIItoOII + 1.e0 + 1.e0/OIItoOI);
    double fOI   = MAX(0.0,1.e0 - fOII - fOIII);
    if(fOIII<0||fOII<0||fOI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"fOIII = %e fOII = %e fOI = %e \n",fOIII,fOII,fOI);
        exit(1);
    }
    PhydroPP->nOIII           = MAX(0.0,nO-nCO)*fOIII;
    PhydroPP->nion[ipOXYGEN]  = MAX(0.0,nO-nCO)*fOII ;
    PhydroPP->natom[ipOXYGEN] = MAX(0.0,nO-nCO)*fOI  ;

    return ;
}


//#ifdef CHEMISTRY_CO_KETO_CASELLI_2008
// equilibrium carbon chemistry
// ref. Keto & Caselli, ApJ, 683, 238, (2008)
//      see footnote 2 of Glover & Clark, MNRAS, 421, 116, (2012)
static double COFormationRateCH2(double nH2, double T){
    return COFormationCoefficient(T)*nH2;
}

static double COtoCIIratio(double nH2, double T, double G0, double Ntot, double N2, double Nco, double Z){
    // eq.(14) in Glover & Clark, MNRAS, 421, 116, (2012)
    //double b_turb = TURBULENCE_VELOCITY;    // km/s
    double b5   = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
    double r_CO = COPhotoDissociationRate(G0,Ntot,Z)*COShieldingFactor(N2,Nco,b5);
    return COFormationRateCH2(nH2,T)/r_CO;
}

static double CItoCIIratio(double nH2, double T, double G0, double Ntot, double Z){
    // eq.(15) in Glover & Clark, MNRAS, 421, 116, (2012)
    double r_C = 3.1e-10*G0*DustShieldingFactor(Ntot,Z,3.0); // R171 in Table B2 of Glover et al. 2010
    return COFormationRateCH2(nH2,T)/r_C;
}

static int SolveAtomicCarbon(StructEQCData *state, double nC, const CEQRateCoefficients *coefficients){
    const CEQIonLadderRates rates = CEQCarbonIonizationRates(state,coefficients);
    return CEQSolveIonLadder(nC,
             rates.ionization0,rates.recombination1,rates.ionization1,rates.recombination2,
             &state->natom[ipCARBON],&state->nion[ipCARBON],&state->nCIII);
}

typedef struct {
    double formation;
    double unshieldedPhoto;
    double photoWithH2;
    double h2Column;
    double carbonColumn;
    double logCarbonColumn;
    double doppler;
    double referenceCO;
} CarbonBalance;

typedef struct {
    const CarbonBalance *balance;
    const CEQCarbonShieldingInterval *interval;
    bool useFreeCarbon;
} CarbonSegment;

typedef struct {
    double molecularFraction;
    double freeFraction;
    bool found;
} CarbonRoot;

static double CarbonSelfShielding(const CEQCarbonShieldingInterval *interval, double logCarbonColumn,
                                  double molecular){
    double logShielding = interval->logShieldingAtAnchor;
    if(interval->slope != 0.0){
        logShielding += interval->slope*(log10(molecular)+logCarbonColumn-interval->logColumnAtAnchor);
    }
    return pow(10.0,logShielding);
}

static double CarbonFluxResidual(const CarbonBalance *balance, double molecular, double freeCarbon,
                                 double selfShielding){
    return FluxResidual(balance->photoWithH2*selfShielding*molecular,
                        balance->formation*freeCarbon);
}

static double CarbonSegmentResidual(double fraction, const void *parameters){
    const CarbonSegment *segment = parameters;
    const CarbonBalance *balance = segment->balance;
    const double molecular = segment->useFreeCarbon ? 1.0-fraction : fraction;
    const double freeCarbon = segment->useFreeCarbon ? fraction : 1.0-fraction;
    const double shielding = CarbonSelfShielding(segment->interval,balance->logCarbonColumn,molecular);
    return CarbonFluxResidual(balance,molecular,freeCarbon,shielding);
}

static int FindCarbonRoot(CarbonSegment *segment, FluxPoint lower, FluxPoint upper,
                          double accuracy, CarbonRoot *best){
    const CarbonBalance *balance = segment->balance;
    const CEQCarbonShieldingInterval *interval = segment->interval;
    // In the CO-dominated gas the remaining carbon is the better-conditioned unknown.
    segment->useFreeCarbon = lower.fraction >= 0.5;
    if(segment->useFreeCarbon){
        // Both complements are exact on [0.5, 1], and the endpoint residuals carry over.
        const FluxPoint freeUpper = {1.0-lower.fraction,lower.residual};
        lower = (FluxPoint){1.0-upper.fraction,upper.residual};
        upper = freeUpper;
    }
    double root;
    // Resolve approximate roots inside a segment before comparing their locations.
    double initialFraction = 0.0;
    if(interval->slope == 0.0){
        const double destruction = balance->photoWithH2*interval->lowerSelfShielding;
        const double rate = destruction+balance->formation;
        if(rate > 0.0) initialFraction = (segment->useFreeCarbon ? destruction : balance->formation)/rate;
    }
    const int status = SolveFluxBracket(CarbonSegmentResidual,segment,lower,upper,
                                       fmin(accuracy,1.e-10),initialFraction,&root);
    if(status == EQC_NO_EQUILIBRIUM) return EQC_SUCCESS;
    if(status != EQC_SUCCESS) return status;

    const double molecular = segment->useFreeCarbon ? 1.0-root : root;
    const double freeCarbon = segment->useFreeCarbon ? root : 1.0-root;
    const double shield = COShieldingFactor(balance->h2Column,balance->carbonColumn*molecular,balance->doppler);
    const double residual = FluxResidual(balance->unshieldedPhoto*shield*molecular,
                                         balance->formation*freeCarbon);
    if(!isfinite(residual)) return EQC_INVALID_INPUT;
    // A root of a one-sided formula at an open endpoint is admissible only
    // if it also satisfies the actual shielding function.
    if(fabs(residual) > accuracy) return EQC_SUCCESS;

    const double distance = fabs(molecular-balance->referenceCO);
    const double bestDistance = fabs(best->molecularFraction-balance->referenceCO);
    if(!best->found || distance < bestDistance ||
       (distance == bestDistance && molecular < best->molecularFraction)){
        best->molecularFraction = molecular;
        best->freeFraction = freeCarbon;
        best->found = true;
    }
    return EQC_SUCCESS;
}

int CEQPrepareCarbonShieldingIntervals(const StructEQCData *state,
                                       CEQCarbonShieldingInterval *intervals, int points){
    const double nC = state->n[ipCARBON];
    if(points < 2 || nC <= 0.0 || state->n[ipOXYGEN] <= 0.0) return 0;
    const double carbonColumn = state->Nt*(nC/state->n[ipHYDROGEN]);
    if(!isfinite(carbonColumn) || carbonColumn < 0.0) return 0;
    const double logCarbonColumn = carbonColumn > 0.0 ? log10(carbonColumn) : -INFINITY;
    const double upperBound = fmin(nC,state->n[ipOXYGEN])/nC;
    double lower = 0.0;
    double anchor = 0.0;
    double logShield = 0.0;
    int count = 0;
    for(int i=0;i<=points && lower < upperBound;i++){
        CEQCarbonShieldingInterval interval = {
            .upperFraction = upperBound,
            .logColumnAtAnchor = anchor,
            .logShieldingAtAnchor = i == points ? -7.0 : logShield,
        };
        if(i < points){
            double nextAnchor, nextShield;
            CEQCOShieldingPoint(i,&nextAnchor,&nextShield);
            interval.upperFraction = fmin(upperBound,pow(10.0,nextAnchor-logCarbonColumn));
            if(i > 0) interval.slope = (nextShield-logShield)/(nextAnchor-anchor);
            anchor = nextAnchor;
            logShield = nextShield;
        }
        if(interval.upperFraction <= lower) continue;
        // Evaluate each endpoint with this interval's shielding law so table-edge jumps retain both limits.
        interval.lowerSelfShielding = CarbonSelfShielding(&interval,logCarbonColumn,lower);
        interval.upperSelfShielding = CarbonSelfShielding(&interval,logCarbonColumn,interval.upperFraction);
        intervals[count++] = interval;
        lower = interval.upperFraction;
    }
    return count;
}

static int SolveCarbonSegments(const CarbonBalance *balance, const CEQCarbonShieldingInterval *intervals,
                               int count, double accuracy, CarbonRoot *best){
    const bool activeFluxes = balance->formation > 0.0 && balance->photoWithH2 > 0.0;
    const double logRateRatio = activeFluxes ? log10(balance->formation) - log10(balance->photoWithH2) : 0.0;
    double lower = 0.0;
    for(int i=0;i<count;i++){
        const CEQCarbonShieldingInterval *interval = &intervals[i];
        const double upper = interval->upperFraction;
        CarbonSegment segment = {.balance = balance, .interval = interval};
        // Split at the possible minimum of D*S(f)*f - F*(1-f).
        double breaks[4] = {lower,upper};
        int size = 2;
        if(interval->slope < -1.0 && activeFluxes){
            const double ratio = logRateRatio - log10(-1.0-interval->slope) - interval->logShieldingAtAnchor;
            const double critical = pow(10.0,interval->logColumnAtAnchor-balance->logCarbonColumn
                                             + ratio/interval->slope);
            if(critical > lower && critical < upper){
                breaks[1] = critical;
                breaks[size++] = upper;
            }
        }
        // The change of root coordinate at 1/2 keeps both trace limits accurate.
        if(lower < 0.5 && upper > 0.5){
            int slot = size;
            while(slot > 0 && breaks[slot-1] > 0.5){
                breaks[slot] = breaks[slot-1];
                slot--;
            }
            breaks[slot] = 0.5;
            size++;
        }
        FluxPoint lowerEndpoint = {lower,CarbonFluxResidual(balance,lower,1.0-lower,interval->lowerSelfShielding)};
        for(int j=1;j<size;j++){
            const double fraction = breaks[j];
            const double shielding = j == size-1 ? interval->upperSelfShielding :
                CarbonSelfShielding(interval,balance->logCarbonColumn,fraction);
            const FluxPoint upperEndpoint = {fraction,CarbonFluxResidual(balance,fraction,1.0-fraction,shielding)};
            const int status = FindCarbonRoot(&segment,lowerEndpoint,upperEndpoint,accuracy,best);
            if(status != EQC_SUCCESS) return status;
            lowerEndpoint = upperEndpoint;
        }
        lower = upper;
    }
    return best->found ? EQC_SUCCESS : EQC_NO_EQUILIBRIUM;
}

int CEQSolveCarbonWithIntervals(StructEQCData *state, double previousCO, double accuracy,
                               const CEQCarbonShieldingInterval *intervals, int count,
                               const CEQRateCoefficients *coefficients, bool checkExistingCO){
    const double nC = state->n[ipCARBON];
    const double nO = state->n[ipOXYGEN];
    if(!isfinite(nC) || nC < 0.0 || !isfinite(nO) || nO < 0.0) return EQC_INVALID_INPUT;
    if(nC == 0.0){
        state->natom[ipCARBON] = state->nion[ipCARBON] = state->nCIII = state->nmol[ipCO] = 0.0;
        return EQC_SUCCESS;
    }
    const bool formCO = EQCLibRunParameters.UseCOFormation && state->T <= 1.e4;
    CEQRateCoefficients local;
    if(coefficients == NULL){
        CEQUVBRates uvb;
        const EQCLibStatus status = CEQEvaluateUVBackgroundWithMetallicity(state->NHIuvb,state->Z,&uvb);
        if(status != EQC_SUCCESS) return status;
        local.carbon = CarbonCoefficients(state,&uvb);
        if(formCO) local.molecularCarbon = MolecularCarbonCoefficients(state,&uvb);
        coefficients = &local;
    }
    if(!formCO){
        const int status = SolveAtomicCarbon(state,nC,coefficients);
        if(status == EQC_SUCCESS) state->nmol[ipCO] = 0.0;
        return status;
    }

    const double nH = state->n[ipHYDROGEN];
    const double hydroRatio = state->nmol[ipH2]/nH;
    const CEQMolecularCarbonCoefficients rates = coefficients->molecularCarbon;
    const double formation = rates.formationPerH2*state->nmol[ipH2];
    const double photoCI = rates.photoCI, photoCO = rates.photoCO;
    if(!isfinite(formation) || formation < 0.0 || !isfinite(photoCI) || photoCI < 0.0 ||
       !isfinite(photoCO) || photoCO < 0.0 || !isfinite(previousCO) || previousCO < 0.0 ||
       !isfinite(hydroRatio) || hydroRatio < 0.0 || !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0){
        return EQC_INVALID_INPUT;
    }
    const double rateScale = fmax(formation,photoCI);
    const bool stellarIonization = EQCLibRunParameters.UseHIIRegion && state->NLyn > 0.0;
    double carbonIon = 0.0;
    double neutralCarbon = 1.0;
    double doubleCarbonIon = 0.0;
    double effectiveFormation = 0.0;
    if(!stellarIonization && rates.photoCII > 0.0){
        /* Keep the molecular C I/C II balance and include the drain from
           C II to C III. Only C II can enter the CO formation channel. */
        const double recombinationCIII = state->nel*coefficients->carbon.recombination2;
        const int status = CEQSolveIonLadder(1.0,photoCI,formation,rates.photoCII,recombinationCIII,
                                             &neutralCarbon,&carbonIon,&doubleCarbonIon);
        if(status != EQC_SUCCESS) return status;
        effectiveFormation = formation*carbonIon;
    }else if(rateScale > 0.0){
        const double scaledFormation = formation/rateScale;
        const double scaledPhoto = photoCI/rateScale;
        carbonIon = scaledPhoto/(scaledFormation+scaledPhoto);
        neutralCarbon = scaledFormation/(scaledFormation+scaledPhoto);
        effectiveFormation = fmin(formation,photoCI)/(scaledFormation+scaledPhoto);
    }
    const double reference = fmin(previousCO,fmin(nC,nO))/nC;
    CarbonRoot root = {.molecularFraction = reference, .freeFraction = 1.0-reference};
    if(photoCO == 0.0 && effectiveFormation > 0.0){
        if(nO < nC) return EQC_NO_EQUILIBRIUM;
        root.molecularFraction = 1.0;
        root.freeFraction = 0.0;
    }else if(photoCO > 0.0 && effectiveFormation == 0.0){
        root.molecularFraction = 0.0;
        root.freeFraction = 1.0;
    }else if(photoCO > 0.0){
        if(CEQCOShieldingPointCount() < 2) return EQC_INVALID_INPUT;
        const double h2Column = state->Nt*hydroRatio;
        const double carbonColumn = state->Nt*(nC/nH);
        if(!isfinite(h2Column) || !isfinite(carbonColumn) || h2Column < 0.0 || carbonColumn < 0.0){
            return EQC_INVALID_INPUT;
        }
        const double doppler = sqrt(SQ(12.90*sqrt(1.e-4*state->T))+SQ(b_turb));
        const double scale = fmax(effectiveFormation,photoCO);
        if(checkExistingCO){
            /* Use the sum of the carbon ions to retain trace free carbon. */
            const double freeCarbon = state->natom[ipCARBON] + state->nion[ipCARBON] + state->nCIII;
            const double shield = COShieldingFactor(h2Column,state->Nt*(state->nmol[ipCO]/nH),doppler);
            const double residual = FluxResidual((photoCO/scale)*shield*reference,
                                                 (effectiveFormation/scale)*(freeCarbon/nC));
            if(fabs(residual) <= accuracy){
                /* The atomic ladder at fixed n_e does not depend on H2. The
                   molecular model's C I/C II split still needs updating. */
                if(!stellarIonization){
                    state->natom[ipCARBON] = freeCarbon*neutralCarbon;
                    state->nion[ipCARBON] = freeCarbon*carbonIon;
                    if(rates.photoCII > 0.0) state->nCIII = freeCarbon*doubleCarbonIon;
                }
                return EQC_SUCCESS;
            }
        }
        const double photoWithH2 = photoCO*COShieldingFactor(h2Column,0.0,doppler);
        CarbonBalance balance = {
            .formation = effectiveFormation/scale,
            .unshieldedPhoto = photoCO/scale,
            .photoWithH2 = photoWithH2/scale,
            .h2Column = h2Column,
            .carbonColumn = carbonColumn,
            .logCarbonColumn = carbonColumn > 0.0 ? log10(carbonColumn) : -INFINITY,
            .doppler = doppler,
            .referenceCO = reference,
        };
        const int status = SolveCarbonSegments(&balance,intervals,count,accuracy,&root);
        if(status != EQC_SUCCESS) return status;
    }

    state->nmol[ipCO] = fmin(nO,nC*root.molecularFraction);
    /* Ionizing radiation also acts on the carbon outside CO. The molecular
       model's C I/C II ratio only accounts for FUV and molecular reactions. */
    if(stellarIonization)
        return SolveAtomicCarbon(state,nC*root.freeFraction,coefficients);
    state->nion[ipCARBON] = (nC*root.freeFraction)*carbonIon;
    state->natom[ipCARBON] = (nC*root.freeFraction)*neutralCarbon;
    state->nCIII = (nC*root.freeFraction)*doubleCarbonIon;
    return EQC_SUCCESS;
}

static int SolveCarbonEquilibrium(StructEQCData *state, double previousCO, double accuracy,
                                  const CEQRateCoefficients *coefficients){
    const int points = CEQCOShieldingPointCount();
    CEQCarbonShieldingInterval intervals[points+1];
    const int count = CEQPrepareCarbonShieldingIntervals(state,intervals,points);
    return CEQSolveCarbonWithIntervals(state,previousCO,accuracy,intervals,count,coefficients,false);
}

int CEQSolveCarbonEquilibrium(StructEQCData *state, double previousCO, double accuracy){
    return SolveCarbonEquilibrium(state,previousCO,accuracy,NULL);
}

// eq(11)-(13) in Glover & Clark, MNRAS, 421, 116, (2012)
void EquilibriumCarbonFractionalAbundance(StructEQCData *PhydroPP){
    CEQUVBRates uvb;
    if(CEQEvaluateUVBackgroundWithMetallicity(PhydroPP->NHIuvb,PhydroPP->Z,&uvb) != EQC_SUCCESS) return;
    if(EQCLibRunParameters.UseCOFormation && PhydroPP->T <= 1.e4 &&
       (uvb.metalIonization[CEQ_UVB_CI] > 0.0 || uvb.metalIonization[CEQ_UVB_CII] > 0.0)){
        StructEQCData result = *PhydroPP;
        const CEQRateCoefficients coefficients = {
            .carbon = CarbonCoefficients(&result,&uvb), .molecularCarbon = MolecularCarbonCoefficients(&result,&uvb)
        };
        if(SolveCarbonEquilibrium(&result,result.nmol[ipCO],EQCLibRunParameters.ChemistryAccuracy,&coefficients) == EQC_SUCCESS)
            *PhydroPP = result;
        return;
    }
    double ne   = PhydroPP->nel;
    double nH   = PhydroPP->n[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    double nC   = PhydroPP->n[ipCARBON];
    double T    = PhydroPP->T;
    double G0   = PhydroPP->G0;
    double X0   = PhydroPP->X0;
    double Ntot = PhydroPP->Nt;
    double Zd   = PhydroPP->Zd;
    double fCIII = 0.e0, fCII = 0.e0, fCI = 0.e0, fCO = 0.e0;
    if(EQCLibRunParameters.UseCOFormation==true){
        double nCO = PhydroPP->nmol [ipCO];
        double Nh2 = Ntot*(nH2/nH);
        double Nco = Ntot*nCO/nH;
        double COtoCII = COtoCIIratio(nH2,T,G0,Ntot,Nh2,Nco,Zd); // eq.(14)
        double CItoCII = CItoCIIratio(nH2,T,G0,Ntot,Zd);         // eq.(15)
        double CIItoCtotal = 1.e0/(1.0 + COtoCII + CItoCII);    // eq.(11)
        fCIII = 0.e0;
        fCII  = CIItoCtotal;
        fCI   = CItoCII*CIItoCtotal;
        fCO   = MAX(0.0,1.e0 - fCII - fCI);
        if(fCIII<0||fCII<0||fCI<0||fCO<0){
            //double b_turb = 0.e0;    // km/s
            double b5     = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
            fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
            fprintf(stderr,"  COtoCII=%e CItoCII=%e \n",COtoCII,CItoCII);
            fprintf(stderr,"  nH2=%e T=%e G0=%e Ntot=%e Nh2=%e Nco=%e Zd=%e \n",nH2,T,G0,Ntot,Nh2,Nco,Zd);
            fprintf(stderr,"  COSelfShielding=%e b5=%e\n",COShieldingFactor(Nh2,Nco,b5),b5);
            fprintf(stderr,"  fCIII=%e fCII=%e fCI=%e fCO=%e\n",fCIII,fCII,fCI,fCO);
            exit(1);
        }
    //} else{
    //} if(EQCLibRunParameters.UseCOFormation==false||fCO<1.e-9){
    } if(EQCLibRunParameters.UseCOFormation==false||T>1.e+4){
        double NLyn       = PhydroPP->NLyn;
        double xe         = ne/nH;
        //double Av         = 0.0; //VisualExtinction(Ntot,Z);    for Wolfire et al. 1995
        double Av         = VisualExtinction(Ntot,Zd);
        
        double dnCI_phti  = 0.e0;
        if(NLyn==0.e0 || !EQCLibRunParameters.UseHIIRegion)  dnCI_phti = CIFUVPhotoionizationRate(G0,Av);
        else            dnCI_phti = CIPhotoionizationRate(NLyn);
        dnCI_phti += uvb.metalIonization[CEQ_UVB_CI];
        double dnCI_coli  = ne*CICollisionalIonizationRate(T);
        double dnCI_cr    = X0*CosmicRayIonizationRateCI(xe);
        double dnCI_xr    = X0*SoftXrayIonizationRateCI(xe,Ntot);
        double dnCII_rec  = ne*CIIRecombinationRate(T);
        double dnCII_dust = nH*CIIRecombinationRateOnDust(T,G0,ne,Zd); 

        double dnCII_phti  = CIIPhotoionizationRate(NLyn) + uvb.metalIonization[CEQ_UVB_CII];
        double dnCII_coli  = ne*CIICollisionalIonizationRate(T);
        double dnCII_cr    = 0.e0;
        double dnCII_xr    = 0.e0;
        double dnCIII_rec  = ne*CIIIRecombinationRate(T);
        double dnCIII_dust = 0.e0;
       
        double CIItoCI   = (dnCI_phti  + dnCI_coli  + dnCI_cr  + dnCI_xr )/(dnCII_rec  + dnCII_dust );
        double CIIItoCII = (dnCII_phti + dnCII_coli + dnCII_cr + dnCII_xr)/(dnCIII_rec + dnCIII_dust);

        //fCII = CIItoCI/(1.0+CIItoCI);
        //fCI  = MAX(0.0,1.0 - fCII);
        //fCO  = 0.e0;
        fCIII = CIIItoCII/(CIIItoCII + 1.e0 + 1.e0/CIItoCI);
        fCII  =      1.e0/(CIIItoCII + 1.e0 + 1.e0/CIItoCI);
        fCI   = MAX(0.0,1.e0 - fCII - fCIII);
        fCO   = 0.e0;
        if(fCIII<0||fCII<0||fCI<0||fCO<0){
            fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
            fprintf(stderr,"  fCIII=%e fCII=%e fCI=%e fCO=%e\n",fCIII,fCII,fCI,fCO);
            exit(1);
        }
    }
	PhydroPP->nCIII           = nC*fCIII;
	PhydroPP->nion[ipCARBON]  = nC*fCII ;
    PhydroPP->natom[ipCARBON] = nC*fCI  ;
    PhydroPP->nmol[ipCO]      = nC*fCO  ;
    return ;
}

void EquilibriumSiliconFractionalAbundance(StructEQCData *PhydroPP){
    CEQUVBRates uvb;
    if(CEQEvaluateUVBackgroundWithMetallicity(PhydroPP->NHIuvb,PhydroPP->Z,&uvb) != EQC_SUCCESS) return;
    double ne     = PhydroPP->nel;
    double nH     = PhydroPP->n[ipHYDROGEN];
    double nSi    = PhydroPP->n[ipSILICON];
    double T      = PhydroPP->T;
    double NLyn   = PhydroPP->NLyn;
    double G0     = PhydroPP->G0;
    double X0     = PhydroPP->X0;
    double Ntot   = PhydroPP->Nt;
    double Zd     = PhydroPP->Zd;
    double xe     = ne/nH;

    double dnSiI_phti   = SiIPhotoionizationRate(NLyn) + uvb.metalIonization[CEQ_UVB_SI_I];
    double dnSiI_coli   = ne*SiICollisionalIonizationRate(T);
    double dnSiI_cr     = X0*CosmicRayIonizationRateSiI(xe);
    double dnSiI_xr     = X0*SoftXrayIonizationRateSiI(xe,Ntot);
    double dnSiII_rec   = ne*SiIIRecombinationRate(T);
    double dnSiII_dust  = nH*SiIIRecombinationRateOnDust(T,G0,ne,Zd); 

    double dnSiII_phti  = SiIIPhotoionizationRate(NLyn) + uvb.metalIonization[CEQ_UVB_SI_II];
    double dnSiII_coli  = ne*SiIICollisionalIonizationRate(T);
    double dnSiII_cr    = X0*CosmicRayIonizationRateSiII(xe);
    double dnSiII_xr    = X0*SoftXrayIonizationRateSiII(xe,Ntot);
    double dnSiIII_rec  = ne*SiIIIRecombinationRate(T);
    double dnSiIII_dust = 0.e0; // OK
    
    double SiIItoSiI   = (dnSiI_phti  + dnSiI_coli  + dnSiI_cr  + dnSiI_xr )/(dnSiII_rec  + dnSiII_dust );
    double SiIIItoSiII = (dnSiII_phti + dnSiII_coli + dnSiII_cr + dnSiII_xr)/(dnSiIII_rec + dnSiIII_dust);
   
    double fSiIII = SiIIItoSiII/(SiIIItoSiII + 1.e0 + 1.e0/SiIItoSiI);
    double fSiII  =        1.e0/(SiIIItoSiII + 1.e0 + 1.e0/SiIItoSiI);
    double fSiI   = MAX(0.0,1.e0 - fSiII - fSiIII);
    if(fSiIII<0||fSiII<0||fSiI<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"  fSiIII=%e fSiII=%e fSiI=%e\n",fSiIII,fSiII,fSiI);
        exit(1);
    }
    PhydroPP->nSiIII            = nSi*fSiIII;
    PhydroPP->nion[ipSILICON]   = nSi*fSiII ;
    PhydroPP->natom[ipSILICON]  = nSi*fSiI  ;
    return ;
}

