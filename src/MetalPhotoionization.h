#ifndef EQC_METAL_PHOTOIONIZATION_H
#define EQC_METAL_PHOTOIONIZATION_H

#include <stddef.h>

enum {
    CEQ_METAL_PI_CI,
    CEQ_METAL_PI_CII,
    CEQ_METAL_PI_OI,
    CEQ_METAL_PI_OII,
    CEQ_METAL_PI_SI_I,
    CEQ_METAL_PI_SI_II,
    CEQ_METAL_PI_SPECIES
};

/* Each entry is cm^2 per ion. The optional diagnostic cross section counts
 * absorptions that eject more than one electron and is a subset of total. */
void CEQMetalPhotoionizationCrossSections(double energyEV,
                                         double total[CEQ_METAL_PI_SPECIES],
                                         double multiplyIonizing[CEQ_METAL_PI_SPECIES]);

/* Pass NULL to obtain the number of integration boundaries, in eV. */
size_t CEQMetalPhotoionizationBoundaries(double *energiesEV);

#endif
