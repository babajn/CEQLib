#include <math.h>

#include "Constants.h"
#include "PhotoelectronDeposition.h"

double CEQPhotoelectronHeatingFraction(double ionizedFraction)
{
    const double remaining = 1.0-pow(ionizedFraction,0.2663);
    return 0.9971*(1.0-pow(remaining,1.3163));
}

static double IonizedHydrogenFraction(const StructEQCData *state)
{
    return fmin(fmax(state->nion[ipHYDROGEN]/state->n[ipHYDROGEN],0.0),1.0);
}

static double FitIonizedFraction(double ionizedFraction)
{
    /* Shull & van Steenberg (1985), Table 2: high-energy fits for x >= 10^-4. */
    return fmax(ionizedFraction,1.0e-4);
}

static CEQPhotoelectronIonization PhotoelectronIonizationRates(const StructEQCData *state, double highPower,
                                                              double ionizedFraction, double fitFraction)
{
    CEQPhotoelectronIonization rates = {0};
    if(ionizedFraction == 1.0) return rates;
    const double hydrogenNeutral = state->natom[ipHYDROGEN] + 2.0*state->nmol[ipH2];
    /* The fit uses He/H = 0.1. Keep a reference He abundance for trace-He
       cells so the total number of He ionizations vanishes with n(He I). */
    const double referenceHelium = fmax(state->n[ipHELIUM],0.1*state->n[ipHYDROGEN]);
    const double heliumNeutral = fmax(state->natom[ipHELIUM],(1.0-ionizedFraction)*referenceHelium);
    const double logFitFraction = log(fitFraction);
    const double hydrogenShare = 0.3908*pow(-expm1(0.4092*logFitFraction),1.7592);
    const double heliumShare = 0.0554*pow(-expm1(0.4614*logFitFraction),1.666);
    if(hydrogenNeutral > 0.0)
        rates.hydrogen = (highPower/hydrogenNeutral)*(hydrogenShare/(13.6*EV_TO_ERG));
    if(heliumNeutral > 0.0)
        rates.helium = (highPower/heliumNeutral)*(heliumShare/(24.59*EV_TO_ERG));
    return rates;
}

CEQPhotoelectronIonization CEQPhotoelectronIonizationRates(const StructEQCData *state, double highPower)
{
    if(highPower == 0.0) return (CEQPhotoelectronIonization){0};
    if(!(state->n[ipHYDROGEN] > 0.0) || !(highPower > 0.0) || !isfinite(highPower))
        return (CEQPhotoelectronIonization){NAN,NAN};
    const double ionizedFraction = IonizedHydrogenFraction(state);
    return PhotoelectronIonizationRates(state,highPower,ionizedFraction,FitIonizedFraction(ionizedFraction));
}

CEQPhotoelectronDeposition CEQDepositPhotoelectrons(const StructEQCData *state, double lowPower, double highPower)
{
    if(highPower == 0.0) return (CEQPhotoelectronDeposition){.heating = lowPower};
    if(!(state->n[ipHYDROGEN] > 0.0) || !(lowPower >= 0.0) || !(highPower > 0.0))
        return (CEQPhotoelectronDeposition){.heating = NAN};

    const double ionizedFraction = IonizedHydrogenFraction(state);
    const double fitFraction = FitIonizedFraction(ionizedFraction);
    const double heatShare = CEQPhotoelectronHeatingFraction(fitFraction);
    CEQPhotoelectronDeposition result = {
        .heating = lowPower + heatShare*highPower,
        .otherPower = (1.0-heatShare)*highPower
    };
    if(ionizedFraction == 1.0) return result;

    const CEQPhotoelectronIonization rates = PhotoelectronIonizationRates(state,highPower,ionizedFraction,fitFraction);
    const double hydrogenThreshold = 13.6*EV_TO_ERG;
    const double heliumThreshold = 24.59*EV_TO_ERG;
    result.ionizationHI = rates.hydrogen;
    result.ionizationHeI = rates.helium;
    if(result.ionizationHI != 0.0)
        result.secondaryHI = state->natom[ipHYDROGEN]*result.ionizationHI;
    if(result.ionizationHeI != 0.0)
        result.secondaryHeI = state->natom[ipHELIUM]*result.ionizationHeI;
    result.otherPower -= hydrogenThreshold*result.secondaryHI + heliumThreshold*result.secondaryHeI;
    return result;
}

static double AbsorbedPhotoelectronPower(const StructEQCData *state, const double power[CEQ_UVB_SPECIES])
{
    return state->natom[ipHYDROGEN]*power[CEQ_UVB_HI]
         + state->natom[ipHELIUM]*power[CEQ_UVB_HEI]
         + state->nion[ipHELIUM]*power[CEQ_UVB_HEII];
}

double CEQHighPhotoelectronPower(const StructEQCData *state, const CEQUVBRates *uvb)
{
    const double atomicPower = AbsorbedPhotoelectronPower(state,uvb->heatingHigh);
    return uvb->metalPowerHigh == 0.0 ? atomicPower : atomicPower + state->n[ipHYDROGEN]*uvb->metalPowerHigh;
}

CEQPhotoelectronDeposition CEQDepositUVBPhotoelectrons(const StructEQCData *state, const CEQUVBRates *uvb)
{
    if(uvb->metalPowerLow == 0.0 && uvb->metalPowerHigh == 0.0 &&
       uvb->heating[CEQ_UVB_HI] == 0.0 && uvb->heating[CEQ_UVB_HEI] == 0.0 && uvb->heating[CEQ_UVB_HEII] == 0.0)
        return (CEQPhotoelectronDeposition){0};

    const double power = AbsorbedPhotoelectronPower(state,uvb->heating);
    const double highPower = AbsorbedPhotoelectronPower(state,uvb->heatingHigh);
    if(uvb->metalPowerLow == 0.0 && uvb->metalPowerHigh == 0.0)
        return CEQDepositPhotoelectrons(state,power-highPower,highPower);
    const double nH = state->n[ipHYDROGEN];
    const double lowPower = power-highPower + nH*uvb->metalPowerLow;
    return CEQDepositPhotoelectrons(state,lowPower,highPower + nH*uvb->metalPowerHigh);
}
