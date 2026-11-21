#ifndef EQC_CHIMES_COMPARE_BACKEND_H
#define EQC_CHIMES_COMPARE_BACKEND_H

/* Common physical setup shared by the EQCLib and CHIMES back ends.
 * Solar gas with the EQCLib gas-phase abundances, an ISRF of G0 = 1.7 Habing,
 * a primary cosmic-ray HI ionization rate of 1.8e-17 s^-1, MW dust,
 * a shielding column N_H = 1e20 (nH/cm^-3)^(1/2) cm^-2, no X-rays and
 * no ionizing (> 13.6 eV) radiation. */
#define CMP_G0            1.7
#define CMP_CR_RATE       1.8e-17
#define CMP_DUST_TEMP     10.0
#define CMP_HE            0.1
#define CMP_C             3.019952e-4
#define CMP_O             4.570882e-4
#define CMP_SI            3.548134e-6
#define CMP_YEAR          3.15576e7

static inline double CmpColumn(double nH) { return 1.0e20*__builtin_sqrt(nH); }

typedef struct {
    double T, xe, xHII, fH2;
    double xHeII, xHeIII, xCII, xCIII, xCO, xOII, xOIII, xSiII;
    double heat, cool;     /* erg cm^-3 s^-1 */
} Snapshot;

typedef struct Cell Cell;

/* Each back end chooses its variant ("default"/"co" for EQCLib,
 * "matched"/"full" for CHIMES) at start-up. */
int BackendInit(const char *variant, const char *dataDirectory);
const char *BackendName(void);
int BackendMethodCount(void);
const char *BackendMethodName(int method);

Cell *BackendNewCell(void);
void BackendCopyCell(Cell *dst, const Cell *src);
/* Fixed-temperature chemical equilibrium; returns 0 on success and fills
 * *iterations with a back-end specific effort counter. */
int BackendEquilibrium(Cell *cell, double nH, double T, int *iterations);
/* Isochoric thermochemical update over dt seconds with the given method. */
int BackendAdvance(Cell *cell, double dt, int method);
/* Composition and rates at the current state (not timed). */
void BackendSnapshot(Cell *cell, Snapshot *snapshot, int computeRates);

#endif
