#ifndef CEQ_HII_REGION_OXYGEN_DATA_H
#define CEQ_HII_REGION_OXYGEN_DATA_H

#include <stddef.h>
#include "HIIregionLineCooling.h"

enum { OxygenTransitionCount = 10 };

typedef struct {
    double wavenumber;        /* Level energy divided by hc, in cm^-1. */
    double statisticalWeight;
} OxygenLevel;

typedef struct {
    unsigned int upper;
    unsigned int lower;
    double spontaneousRate;   /* Einstein A, in s^-1. */
    const double *collisionStrength;  /* Maxwellian effective collision strength. */
} OxygenTransition;

typedef struct {
    const OxygenLevel *levels;
    const OxygenTransition *transitions;
    const double *logTemperature;     /* log10(T/K). */
    size_t temperatureCount;
} OxygenIonData;

/* Five lowest fine-structure levels in increasing energy order from the
   NIST Atomic Spectra Database. The transition and collision data are
   tabulated in PyNeb (Morisset et al. 2020, doi:10.3390/atoms8040066). */

/* O II radiative rates: Zeippen (1982, MNRAS 198, 111), o_ii_atom_Z82-WFD96.
   Effective collision strengths: Kisielius et al. (2009, MNRAS 397, 903),
   o_ii_coll_Kal09, doi:10.1111/j.1365-2966.2009.14989.x. */
static const OxygenLevel OxygenIILevels[OxygenLevelCount] = {
    {0, 4},
    {26810.55, 6},
    {26830.57, 4},
    {40468.01, 4},
    {40470, 2},
};

static const double OxygenIILogTemperature[] = {
    2.00, 2.18, 2.30, 2.48, 2.70, 2.88, 3.00, 3.18, 3.30, 3.48,
    3.70, 3.88, 4.00, 4.18, 4.30, 4.48, 4.70, 4.88, 5.00
};

static const OxygenTransition OxygenIITransitions[OxygenTransitionCount] = {
    {1, 0, 3.82e-5, (const double[]){
        0.796, 0.797, 0.798, 0.801, 0.808, 0.817, 0.823, 0.830,
        0.832, 0.832, 0.831, 0.833, 0.834, 0.839, 0.844, 0.856,
        0.881, 0.905, 0.919
    }},
    {2, 0, 1.65e-4, (const double[]){
        0.531, 0.533, 0.533, 0.535, 0.540, 0.546, 0.550, 0.554,
        0.555, 0.554, 0.553, 0.553, 0.554, 0.557, 0.561, 0.569,
        0.585, 0.601, 0.611
    }},
    {2, 1, 1.2e-7, (const double[]){
        1.095, 1.086, 1.078, 1.072, 1.097, 1.151, 1.194, 1.239,
        1.254, 1.256, 1.241, 1.221, 1.203, 1.183, 1.179, 1.193,
        1.229, 1.257, 1.270
    }},
    {3, 0, 5.6400001e-2, (const double[]){
        0.244, 0.245, 0.245, 0.245, 0.245, 0.246, 0.246, 0.247,
        0.247, 0.249, 0.251, 0.253, 0.256, 0.260, 0.265, 0.274,
        0.290, 0.304, 0.312
    }},
    {3, 1, 1.17e-1, (const double[]){
        0.791, 0.793, 0.793, 0.794, 0.796, 0.797, 0.799, 0.801,
        0.804, 0.809, 0.820, 0.834, 0.851, 0.891, 0.930, 0.997,
        1.084, 1.144, 1.178
    }},
    {3, 2, 6.14e-2, (const double[]){
        0.439, 0.440, 0.440, 0.440, 0.441, 0.442, 0.443, 0.444,
        0.445, 0.448, 0.454, 0.462, 0.472, 0.494, 0.516, 0.551,
        0.595, 0.624, 0.639
    }},
    {4, 0, 2.32e-2, (const double[]){
        0.126, 0.126, 0.126, 0.126, 0.127, 0.127, 0.127, 0.127,
        0.128, 0.128, 0.129, 0.131, 0.132, 0.134, 0.136, 0.141,
        0.149, 0.155, 0.159
    }},
    {4, 1, 6.1500002e-2, (const double[]){
        0.315, 0.316, 0.316, 0.316, 0.317, 0.318, 0.318, 0.319,
        0.320, 0.322, 0.326, 0.332, 0.339, 0.356, 0.371, 0.396,
        0.427, 0.447, 0.458
    }},
    {4, 2, 1.02e-1, (const double[]){
        0.308, 0.308, 0.309, 0.309, 0.310, 0.310, 0.311, 0.312,
        0.313, 0.315, 0.319, 0.324, 0.331, 0.345, 0.360, 0.386,
        0.421, 0.445, 0.459
    }},
    {4, 3, 2.0800001e-11, (const double[]){
        0.273, 0.274, 0.274, 0.274, 0.274, 0.275, 0.275, 0.276,
        0.276, 0.277, 0.279, 0.282, 0.285, 0.294, 0.305, 0.327,
        0.361, 0.388, 0.405
    }},
};

static const OxygenIonData OxygenII = {
    .levels = OxygenIILevels,
    .transitions = OxygenIITransitions,
    .logTemperature = OxygenIILogTemperature,
    .temperatureCount = sizeof(OxygenIILogTemperature)/sizeof(OxygenIILogTemperature[0])
};

/* O III radiative rates: Froese Fischer & Tachiev (2004, ADNDT 87, 1)
   and Storey & Zeippen (2000, MNRAS 312, 813), o_iii_atom_FFT04-SZ00.
   Effective collision strengths: Storey, Sochi & Badnell (2014,
   MNRAS 441, 3028), o_iii_coll_SSB14, doi:10.1093/mnras/stu777. */
static const OxygenLevel OxygenIIILevels[OxygenLevelCount] = {
    {0, 1},
    {113.178, 3},
    {306.174, 5},
    {20273.27, 5},
    {43185.74, 1},
};

static const double OxygenIIILogTemperature[] = {
    2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,
    3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,
    4.00, 4.10, 4.20, 4.30, 4.40
};

static const OxygenTransition OxygenIIITransitions[OxygenTransitionCount] = {
    {1, 0, 2.596e-5, (const double[]){
        0.63500, 0.62649, 0.61540, 0.60216, 0.58741, 0.57193, 0.55668, 0.54292,
        0.53194, 0.52440, 0.51992, 0.51730, 0.51544, 0.51396, 0.51320, 0.51369,
        0.51585, 0.51993, 0.52604, 0.53378, 0.54210, 0.54956, 0.55511, 0.55866,
        0.56091
    }},
    {2, 0, 3.03e-11, (const double[]){
        0.22586, 0.22561, 0.22529, 0.22490, 0.22447, 0.22411, 0.22408, 0.22475,
        0.22633, 0.22869, 0.23126, 0.23343, 0.23489, 0.23583, 0.23671, 0.23794,
        0.23978, 0.24246, 0.24618, 0.25101, 0.25678, 0.26315, 0.26978, 0.27654,
        0.28349
    }},
    {2, 1, 9.632e-5, (const double[]){
        1.11212, 1.10993, 1.10711, 1.10362, 1.09952, 1.09530, 1.09235, 1.09280,
        1.09839, 1.10890, 1.12183, 1.13386, 1.14303, 1.14971, 1.15583, 1.16333,
        1.17355, 1.18745, 1.20566, 1.22793, 1.25256, 1.27703, 1.29937, 1.31924,
        1.33776
    }},
    {3, 0, 2.322e-6, (const double[]){
        0.23177, 0.23211, 0.23259, 0.23335, 0.23449, 0.23600, 0.23761, 0.23888,
        0.23934, 0.23872, 0.23701, 0.23441, 0.23123, 0.22786, 0.22485, 0.22306,
        0.22342, 0.22653, 0.23236, 0.24025, 0.24917, 0.25817, 0.26652, 0.27384,
        0.27983
    }},
    {3, 1, 6.7909998e-3, (const double[]){
        0.69748, 0.69849, 0.69994, 0.70220, 0.70562, 0.71013, 0.71497, 0.71876,
        0.72014, 0.71829, 0.71317, 0.70539, 0.69591, 0.68588, 0.67700, 0.67179,
        0.67301, 0.68249, 0.70010, 0.72382, 0.75062, 0.77758, 0.80262, 0.82453,
        0.84246
    }},
    {3, 2, 2.046e-2, (const double[]){
        1.17021, 1.17188, 1.17427, 1.17802, 1.18370, 1.19123, 1.19932, 1.20571,
        1.20807, 1.20503, 1.19654, 1.18365, 1.16800, 1.15154, 1.13712, 1.12887,
        1.13139, 1.14765, 1.17740, 1.21720, 1.26198, 1.30687, 1.34841, 1.38465,
        1.41423
    }},
    {4, 0, 0.0, (const double[]){
        0.02989, 0.02989, 0.02989, 0.02988, 0.02987, 0.02985, 0.02983, 0.02980,
        0.02976, 0.02970, 0.02963, 0.02952, 0.02937, 0.02918, 0.02900, 0.02891,
        0.02906, 0.02954, 0.03036, 0.03144, 0.03268, 0.03395, 0.03517, 0.03625,
        0.03710
    }},
    {4, 1, 2.255e-1, (const double[]){
        0.09003, 0.09003, 0.09001, 0.08998, 0.08994, 0.08989, 0.08982, 0.08972,
        0.08960, 0.08942, 0.08919, 0.08887, 0.08842, 0.08787, 0.08732, 0.08709,
        0.08756, 0.08901, 0.09148, 0.09475, 0.09849, 0.10233, 0.10601, 0.10927,
        0.11183
    }},
    {4, 2, 6.998e-4, (const double[]){
        0.15117, 0.15116, 0.15112, 0.15106, 0.15099, 0.15089, 0.15075, 0.15058,
        0.15036, 0.15007, 0.14969, 0.14917, 0.14846, 0.14756, 0.14670, 0.14636,
        0.14719, 0.14967, 0.15386, 0.15940, 0.16570, 0.17219, 0.17837, 0.18385,
        0.18816
    }},
    {4, 3, 1.685e+0, (const double[]){
        0.38267, 0.38291, 0.38316, 0.38346, 0.38382, 0.38428, 0.38486, 0.38561,
        0.38660, 0.38792, 0.38974, 0.39242, 0.39684, 0.40497, 0.41997, 0.44484,
        0.47990, 0.52127, 0.56212, 0.59558, 0.61740, 0.62672, 0.62540, 0.61639,
        0.60225
    }},
};

static const OxygenIonData OxygenIII = {
    .levels = OxygenIIILevels,
    .transitions = OxygenIIITransitions,
    .logTemperature = OxygenIIILogTemperature,
    .temperatureCount = sizeof(OxygenIIILogTemperature)/sizeof(OxygenIIILogTemperature[0])
};

#endif
