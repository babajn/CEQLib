#ifndef EQC_PHOTOELECTRON_DEPOSITION_H
#define EQC_PHOTOELECTRON_DEPOSITION_H

#include "UVBackground.h"

typedef struct {
    double heating, otherPower;           /* erg cm^-3 s^-1 */
    double ionizationHI, ionizationHeI;   /* s^-1 per target */
    double secondaryHI, secondaryHeI;     /* cm^-3 s^-1 */
} CEQPhotoelectronDeposition;

typedef struct {
    double hydrogen, helium;             /* s^-1 per target */
} CEQPhotoelectronIonization;

/* High-energy heating fit. The caller supplies the ionized fraction. */
double CEQPhotoelectronHeatingFraction(double ionizedFraction);

/* Only the per-target rates are needed during fixed-temperature equilibrium. */
CEQPhotoelectronIonization CEQPhotoelectronIonizationRates(const StructEQCData *state, double highPower);

/* Use bounded trial-state abundances and volume powers in erg cm^-3 s^-1. */
CEQPhotoelectronDeposition CEQDepositPhotoelectrons(const StructEQCData *state, double lowPower, double highPower);
double CEQHighPhotoelectronPower(const StructEQCData *state, const CEQUVBRates *uvb);
CEQPhotoelectronDeposition CEQDepositUVBPhotoelectrons(const StructEQCData *state, const CEQUVBRates *uvb);

#endif
