#ifndef EQC_CHIMES_THERMAL_DIAGNOSTICS_H
#define EQC_CHIMES_THERMAL_DIAGNOSTICS_H

#include "CompareBackend.h"

typedef struct {
    double fH2, xe, xHII;            /* Ratios to the number of H nuclei. */
    double xCII, xOI, xSiII;         /* Fractions of the parent element. */
} ThermalComposition;

typedef struct {
    ThermalComposition composition;
    double xHI;
    double pe_g0, pe_psi;            /* Attenuated G0 and grain charging parameter. */
    double heat_pe, heat_cr, heat_h2, heat_photo, heat_other, heat;
    double cool_h2, cool_h2_diss, cool_cii, cool_oi, cool_ci, cool_si2, cool_co;
    double cool_dust_rec, cool_gas_dust, cool_other, cool;
    double cool_cii_e_only;         /* Isolated C II evaluation outside the budget sum. */
} ThermalRates;                     /* Heating and cooling in erg cm^-3 s^-1. */

/* Sample the native cell when fixed is NULL. Otherwise use the given
   colliders, neutral He and complementary C I, O II and Si I. */
int BackendThermalRates(const Cell *cell, double nH, double T,
                        const ThermalComposition *fixed, ThermalRates *rates);

#endif
