#ifndef EQC_CHIMES_H2_DIAGNOSTICS_H
#define EQC_CHIMES_H2_DIAGNOSTICS_H

#include "CompareBackend.h"

typedef struct {
    double fH2, xe, xHII;
    double form_grain, form_other;                /* H2 molecules cm^-3 s^-1 */
    double destroy_photo, destroy_coll, destroy_other;
    double k_form, k_photo0;                     /* cm^3 s^-1 and s^-1 */
    double dust_shield, self_shield, shield_other;
    double cool_h2, cool_h2_diss, heat_h2;         /* erg cm^-3 s^-1 */
    double heat, cool, energy;                   /* erg cm^-3 s^-1, energy erg cm^-3 */
    double net_oh, net_chplus;                   /* net H2 loss in reversible channels, cm^-3 s^-1 */
} H2Rates;

/* Negative fH2 samples the cell. A nonnegative fH2 samples a prescribed
   mixture of H I, H2, He I, C II, O I, Si I and electrons. */
int BackendH2Rates(Cell *cell, double nH, double T, double fH2, H2Rates *rates);

#endif
