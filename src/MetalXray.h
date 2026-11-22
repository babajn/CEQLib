#ifndef EQC_METAL_XRAY_H
#define EQC_METAL_XRAY_H

#include <stddef.h>

enum {
    CEQ_METAL_XRAY_C,
    CEQ_METAL_XRAY_N,
    CEQ_METAL_XRAY_O,
    CEQ_METAL_XRAY_NE,
    CEQ_METAL_XRAY_NA,
    CEQ_METAL_XRAY_MG,
    CEQ_METAL_XRAY_AL,
    CEQ_METAL_XRAY_SI,
    CEQ_METAL_XRAY_S,
    CEQ_METAL_XRAY_CL,
    CEQ_METAL_XRAY_AR,
    CEQ_METAL_XRAY_CA,
    CEQ_METAL_XRAY_CR,
    CEQ_METAL_XRAY_FE,
    CEQ_METAL_XRAY_NI,
    CEQ_METAL_XRAY_ELEMENTS
};

double CEQMetalXrayCrossSection(size_t element, double energyEV); /* cm^2 per atom */
double CEQMetalXraySolarAbundance(size_t element); /* n_i/n_H */
double CEQMetalXrayEnergyLoss(size_t element, double energyEV); /* eV per absorption */

/* Solar mixture cross section and cross-section-weighted electron energies.
 * Units are cm^2 per H nucleus and cm^2 eV per H nucleus, respectively. */
void CEQMetalXrayMixture(double energyEV, double *crossSection,
                         double *lowWeight, double *highWeight);

/* Pass NULL to obtain the array size. Otherwise write that many energies
 * in eV, without ordering or removing duplicates. */
size_t CEQMetalXrayBoundaries(double *energiesEV);

#endif
