#include <math.h>
#include "MetalXray.h"

static const double MinimumEnergyEV = 30.0;
static const double MaximumEnergyEV = 10000.0;
static const double HighElectronEnergyEV = 100.0;
static const double BCMAvogadro = 6.022045e23;

typedef struct {
    double upperEV;
    size_t coefficientCount;
    double coefficient[9];
} MetalXraySegment;

typedef struct {
    double abundance, atomicMass;
    double firstIonizationEV, secondIonizationEV;
    double kYield, kAlphaEV, lEdgeEV, kEdgeEV;
    size_t segmentCount;
    MetalXraySegment segments[4];
} MetalXrayElement;

#include "MetalXrayData.h"

static double AtomicCrossSection(const MetalXrayElement *element, double energyEV,
                                 double logFitEnergyEV, double inverseEnergyCubed)
{
    const MetalXraySegment *segment = element->segments;
    while (segment < element->segments + element->segmentCount - 1 &&
           energyEV >= segment->upperEV) {
        ++segment;
    }
    size_t i = segment->coefficientCount - 1;
    double polynomial = segment->coefficient[i];
    while (i > 0) {
        --i;
        polynomial = polynomial * logFitEnergyEV + segment->coefficient[i];
    }
    return exp(polynomial) * inverseEnergyCubed * (element->atomicMass / BCMAvogadro);
}

static double EnergyLoss(const MetalXrayElement *element, double energyEV)
{
    if (energyEV >= element->kEdgeEV) {
        /* Radiative decay leaves one ionization and lets the photon escape.
         * Auger decay leaves two ionizations. Average over the decay branches. */
        return element->firstIonizationEV +
               (1.0 - element->kYield) * element->secondIonizationEV +
               element->kYield * element->kAlphaEV;
    }
    if (element->lEdgeEV > 0.0 && energyEV >= element->lEdgeEV) {
        /* L fluorescence yields are small for these elements. */
        return element->firstIonizationEV + element->secondIonizationEV;
    }
    return element->firstIonizationEV;
}

double CEQMetalXrayCrossSection(size_t element, double energyEV)
{
    if (element >= CEQ_METAL_XRAY_ELEMENTS || !(energyEV >= MinimumEnergyEV)) {
        return 0.0;
    }
    const double inverseEnergy = 1.0 / energyEV;
    /* Keep the 10 keV endpoint continuous and extrapolate as E^-3. */
    return AtomicCrossSection(&Elements[element], energyEV,
                              log(fmin(energyEV, MaximumEnergyEV)),
                              inverseEnergy * inverseEnergy * inverseEnergy);
}

double CEQMetalXraySolarAbundance(size_t element)
{
    return element < CEQ_METAL_XRAY_ELEMENTS ? Elements[element].abundance : 0.0;
}

double CEQMetalXrayEnergyLoss(size_t element, double energyEV)
{
    return element < CEQ_METAL_XRAY_ELEMENTS && energyEV >= MinimumEnergyEV ?
           EnergyLoss(&Elements[element], energyEV) : 0.0;
}

void CEQMetalXrayMixture(double energyEV, double *crossSection,
                         double *lowWeight, double *highWeight)
{
    *crossSection = 0.0;
    *lowWeight = 0.0;
    *highWeight = 0.0;
    if (!(energyEV >= MinimumEnergyEV) || !isfinite(energyEV)) {
        return;
    }
    const double inverseEnergy = 1.0 / energyEV;
    const double inverseEnergyCubed = inverseEnergy * inverseEnergy * inverseEnergy;
    const double logFitEnergyEV = log(fmin(energyEV, MaximumEnergyEV));
    for (size_t i = 0; i < CEQ_METAL_XRAY_ELEMENTS; ++i) {
        const MetalXrayElement *element = &Elements[i];
        const double opacity = element->abundance *
            AtomicCrossSection(element, energyEV, logFitEnergyEV, inverseEnergyCubed);
        const double kinetic = energyEV - EnergyLoss(element, energyEV);
        *crossSection += opacity;
        if (kinetic < HighElectronEnergyEV) {
            *lowWeight += opacity * kinetic;
        } else {
            *highWeight += opacity * kinetic;
        }
    }
}

static void AddBoundary(double *energiesEV, size_t *count, double energyEV)
{
    if (energiesEV) {
        energiesEV[*count] = energyEV;
    }
    ++*count;
}

size_t CEQMetalXrayBoundaries(double *energiesEV)
{
    size_t count = 0;
    AddBoundary(energiesEV, &count, MinimumEnergyEV);
    AddBoundary(energiesEV, &count, MaximumEnergyEV);
    for (size_t i = 0; i < CEQ_METAL_XRAY_ELEMENTS; ++i) {
        const MetalXrayElement *element = &Elements[i];
        double lower = MinimumEnergyEV;
        for (size_t j = 0; j < element->segmentCount; ++j) {
            const double upper = element->segments[j].upperEV;
            const double electronBoundary = HighElectronEnergyEV + EnergyLoss(element, lower);
            if (electronBoundary > lower && electronBoundary < upper) {
                AddBoundary(energiesEV, &count, electronBoundary);
            }
            if (upper < MaximumEnergyEV) {
                AddBoundary(energiesEV, &count, upper);
            }
            lower = upper;
        }
    }
    return count;
}
