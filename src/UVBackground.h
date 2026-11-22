#ifndef EQC_UV_BACKGROUND_H
#define EQC_UV_BACKGROUND_H

#include <stddef.h>
#include "ChemicalEquilibriumCooling.h"
#include "MetalPhotoionization.h"

enum {
    CEQ_UVB_HI,
    CEQ_UVB_HEI,
    CEQ_UVB_HEII,
    CEQ_UVB_SPECIES
};

enum {
    CEQ_UVB_CI = CEQ_METAL_PI_CI,
    CEQ_UVB_CII = CEQ_METAL_PI_CII,
    CEQ_UVB_OI = CEQ_METAL_PI_OI,
    CEQ_UVB_OII = CEQ_METAL_PI_OII,
    CEQ_UVB_SI_I = CEQ_METAL_PI_SI_I,
    CEQ_UVB_SI_II = CEQ_METAL_PI_SI_II,
    CEQ_UVB_METAL_SPECIES = CEQ_METAL_PI_SPECIES
};

typedef struct {
    double ionization[CEQ_UVB_SPECIES];  /* s^-1 per absorber */
    double heating[CEQ_UVB_SPECIES];     /* Total photoelectron power, erg s^-1 per absorber */
    double heatingHigh[CEQ_UVB_SPECIES]; /* Power from electrons at or above 100 eV */
    double metalPowerLow, metalPowerHigh; /* Metal photoelectron power, erg s^-1 per H nucleus, at cell Z */
    double metalIonization[CEQ_UVB_METAL_SPECIES]; /* s^-1 per C, O or Si absorber */
} CEQUVBRates;

typedef struct {
    size_t redshiftCount, energyCount;
    double *redshifts;  /* Strictly increasing. */
    double *energyEV;   /* Strictly increasing. */
    double *intensity;  /* [redshift][energy], erg s^-1 cm^-2 Hz^-1 sr^-1. */
} CEQUVBSpectrum;

/* Transfer ownership of the arrays on success. Leave the output intact on failure. */
EQCLibStatus CEQReadUVBackgroundSpectrum(const char *directory, CEQUVBSpectrum *spectrum);
EQCLibStatus CEQInitializeUVBackground(int check);
bool CEQUVBackgroundIsEnabled(void);
EQCLibStatus CEQEvaluateUVBackgroundHemisphere(double column, CEQUVBRates *rates);
EQCLibStatus CEQEvaluateUVBackground(const double columns[2], CEQUVBRates *rates);
EQCLibStatus CEQEvaluateUVBackgroundHemisphereWithMetallicity(double column, double Z, CEQUVBRates *rates);
EQCLibStatus CEQEvaluateUVBackgroundWithMetallicity(const double columns[2], double Z, CEQUVBRates *rates);
size_t CEQUVBPersistentBytes(void);

static inline bool CEQHasUVBHighPower(const CEQUVBRates *rates, double nHe)
{
    return rates->metalPowerHigh > 0.0 || rates->heatingHigh[CEQ_UVB_HI] > 0.0 ||
        (nHe > 0.0 && (rates->heatingHigh[CEQ_UVB_HEI] > 0.0 || rates->heatingHigh[CEQ_UVB_HEII] > 0.0));
}

#endif
