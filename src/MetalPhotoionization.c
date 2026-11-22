#include <math.h>

#include "MetalPhotoionization.h"
#include "PhotochemicalGasPhaseChemistry.h"

typedef struct {
    unsigned angularMomentum;
    double thresholdEV, scaleEV, crossSectionMb, ya, exponent, width;
    double multipleProbability;
} SubshellCrossSection;

typedef struct {
    double (*outerCrossSection)(double);
    size_t innerCount, shellCount;
    SubshellCrossSection shells[5];
} IonCrossSections;

/* Verner & Yakovlev (1995), A&AS 109, 125, Table 1 (CDS J/A+AS/109/125).
 * The Si II K threshold follows the authors' updated phfit2.f (1848 eV).
 * Shells are ordered 1s, 2s, 2p, 3s, 3p. The optional multiple-ejection
 * probabilities are from Kaastra & Mewe (1993), A&AS 97, 443, as grouped by
 * subshell in Cloudy c25.00 data/mewe_nelectron.dat. */
static const IonCrossSections Ions[CEQ_METAL_PI_SPECIES] = {
    [CEQ_METAL_PI_CI] = { PhotoionizationCrossSectionCI, 1, 3, {
        { 0, 291.0, 86.55, 74.21, 54.98, 1.503, 0, .9974 },
        { 0, 19.39, 10.26, 4564, 1.568, 10.85, 0, 0 },
        { 1, 11.26, 9.435, 1152, 5.687, 6.336, .4474, 0 },
    }},
    [CEQ_METAL_PI_CII] = { PhotoionizationCrossSectionCII, 1, 3, {
        { 0, 307.6, 91.13, 66.49, 96.09, 1.338, 0, .9981 },
        { 0, 30.47, 2.991, 1184, 3.085, 14.80, 0, 0 },
        { 1, 24.38, 10.94, 179.2, 33.08, 4.150, .5276, 0 },
    }},
    [CEQ_METAL_PI_OI] = { PhotoionizationCrossSectionOI, 1, 3, {
        { 0, 538.0, 177.4, 32.37, 381.2, 1.083, 0, .9906 },
        { 0, 28.48, 19.94, 241.5, 3.241, 8.037, 0, 0 },
        { 1, 13.62, 13.91, 122000, 1.364, 11.40, .4103, 0 },
    }},
    [CEQ_METAL_PI_OII] = { PhotoionizationCrossSectionOII, 1, 3, {
        { 0, 558.1, 169.0, 35.84, 189.4, 1.185, 0, .9891 },
        { 0, 45.99, 17.59, 196.2, 4.020, 7.999, 0, 0 },
        { 1, 35.12, 17.45, 518.6, 17.28, 4.995, .02182, 0 },
    }},
    [CEQ_METAL_PI_SI_I] = { PhotoionizationCrossSectionSiI, 3, 5, {
        { 0, 1846, 532.2, 11.84, 258.0, 1.102, 0, .9988 },
        { 0, 156.0, 70.17, 11.66, 47.42, 3.933, 0, 1 },
        { 1, 106.0, 78.08, 153.2, 5.765e6, 2.639, .0002774, .9973334 },
        { 0, 15.17, 14.13, 11.66, 22.88, 5.334, 0, 0 },
        { 1, 8.152, 22.12, 184.5, 3.849, 9.721, .2921, 0 },
    }},
    [CEQ_METAL_PI_SI_II] = { PhotoionizationCrossSectionSiII, 3, 5, {
        { 0, 1848, 458.0, 16.28, 77.06, 1.385, 0, .9993 },
        { 0, 161.9, 67.38, 12.63, 36.23, 4.172, 0, 1 },
        { 1, 118.6, 71.54, 183.2, 353.7, 3.133, .0002870, .9958334 },
        { 0, 22.40, 13.64, 8.454, 74.89, 4.676, 0, 0 },
        { 1, 16.35, 21.23, 69.75, 4.907, 9.525, .3169, 0 },
    }},
};

static double PartialCrossSection(const SubshellCrossSection *shell, double energyEV)
{
    if (energyEV < shell->thresholdEV) return 0.0;

    /* Verner & Yakovlev Eq. (1), arranged in terms of E_0/E to keep the
     * high-energy limit finite without multiplying very large and small powers. */
    const double inverseY = shell->scaleEV / energyEV;
    const double offset = 1.0 - inverseY;
    const double width = shell->width * inverseY;
    const double scale = sqrt(inverseY) + 1.0 / sqrt(shell->ya);
    return 1.e-18 * shell->crossSectionMb * (offset * offset + width * width) *
           pow(inverseY, 3.5 + shell->angularMomentum) * pow(scale, -shell->exponent);
}

void CEQMetalPhotoionizationCrossSections(double energyEV,
                                         double total[CEQ_METAL_PI_SPECIES],
                                         double multiplyIonizing[CEQ_METAL_PI_SPECIES])
{
    for (size_t i = 0; i < CEQ_METAL_PI_SPECIES; ++i) {
        const IonCrossSections *ion = &Ions[i];
        const double firstInnerEdge = ion->shells[ion->innerCount - 1].thresholdEV;
        double crossSection = 0.0, multiple = 0.0;
        if (energyEV > 0.0 && isfinite(energyEV)) {
            /* Match the authors' phfit2.f switch: use the aggregate outer
             * fit up to the first inner edge, then switch to all partial fits. */
            if (energyEV < firstInnerEdge) {
                crossSection = ion->outerCrossSection(energyEV);
            } else {
                for (size_t j = 0; j < ion->shellCount; ++j) {
                    const double partial = PartialCrossSection(&ion->shells[j], energyEV);
                    crossSection += partial;
                    multiple += partial * ion->shells[j].multipleProbability;
                }
            }
        }
        total[i] = crossSection;
        if (multiplyIonizing) multiplyIonizing[i] = multiple;
    }
}

size_t CEQMetalPhotoionizationBoundaries(double *energiesEV)
{
    size_t count = 0;
    for (size_t i = 0; i < CEQ_METAL_PI_SPECIES; ++i) {
        const IonCrossSections *ion = &Ions[i];
        for (size_t j = 0; j < ion->innerCount; ++j) {
            if (energiesEV) energiesEV[count] = ion->shells[j].thresholdEV;
            ++count;
        }
        if (energiesEV) energiesEV[count] = ion->shells[ion->shellCount - 1].thresholdEV;
        ++count;
    }
    return count;
}
