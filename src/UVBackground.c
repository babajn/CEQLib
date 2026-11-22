#include <gsl/gsl_sf_expint.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Constants.h"
#include "MetalPhotoionization.h"
#include "MetalXray.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "UVBackground.h"

enum { ColumnIntervals = 240, ColumnCount = ColumnIntervals + 2 };

static const double MinimumLogColumn = 12.0;
static const double LogColumnStep = 0.05;
static const double MaximumColumn = 1.e24;
static const double MaxLogEnergyStep = 0.01;
static const double PhotoelectronSplitEV = 100.0;
static const double HeliumAbundance = 0.1; /* N_HeI = 0.1 N_HI, N_HeII = 0. */
static const double Thresholds[CEQ_UVB_SPECIES] = {13.60, 24.59, 54.42};
static const char *const SpeciesNames[CEQ_UVB_SPECIES] = {"HI", "HeI", "HeII"};
static const char DefaultDirectory[] = "../data/FG20";

typedef struct {
    double ionization[CEQ_UVB_SPECIES];
    double heating[CEQ_UVB_SPECIES];
} UVBPrimaryMoments;

typedef struct {
    size_t left;
    size_t absorptionIndex;
    double fraction;
    double opacity;
    UVBPrimaryMoments weight;
    double metalOpacity, metalWeightLow, metalWeightHigh;
    double metalIonization[CEQ_UVB_METAL_SPECIES];
    unsigned int fastElectrons;
} QuadraturePoint;

typedef enum {AtomicQuadrature, ShieldedQuadrature, MetalIonizationQuadrature} UVBQuadratureKind;

typedef struct {
    size_t count;
    QuadraturePoint *points;
    double *transmission;
    UVBPrimaryMoments *integrands;
} QuadratureGrid;

typedef struct {
    UVBPrimaryMoments logRates;
    double highFraction[CEQ_UVB_SPECIES];
} UVBTableRow;

typedef struct {
    size_t index;
    double fraction;
} InterpolationPoint;

/* Resolve the metallicity-dependent shielding between the principal grid
   values without coarsening either the original redshift or column grid. */
static const double MetalRatios[] = {
    0, .015, .03, .05, .075, .1, .125, .15, .175, .2, .225, .25, .275, .3,
    .35, .4, .45, .5, .56, .62, .68, .74, .8, .87, .94, 1, 1.1, 1.2, 1.3,
    1.4, 1.5, 1.65, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.25, 3.5, 3.75, 4
};
enum {
    MetalRatioCount = sizeof MetalRatios/sizeof MetalRatios[0],
    MetalGamma = 0, MetalAtomLow = CEQ_UVB_SPECIES, MetalAtomHigh = 2*CEQ_UVB_SPECIES,
    MetalLow = 3*CEQ_UVB_SPECIES, MetalHigh, MetalMoments
};

typedef struct {
    double value[MetalMoments];
    double ionization[CEQ_UVB_METAL_SPECIES];
} MetalRateRow;

typedef struct {
    size_t redshift;
    MetalRateRow *rows;
} MetalSlice;

typedef struct {
    CEQUVBSpectrum spectrum;
    QuadraturePoint *points, *ionizationPoints;
    size_t pointCount, ionizationPointCount;
    MetalSlice slice[2];
    MetalRateRow *blend;
    const MetalRateRow *current;
} MetalRateTable;

static bool Enabled;
static bool Initialized;
static char *SpectrumDirectory;
static double MaximumEnergyEV;
static UVBTableRow *RateTable; /* [redshift][column]; -INFINITY encodes zero rates. */
static InterpolationPoint CurrentRedshift;
static MetalRateTable MetalTable;

bool CEQUVBackgroundIsEnabled(void)
{
    return Enabled;
}

static void FreeMetalTable(MetalRateTable *table)
{
    free(table->spectrum.redshifts);
    free(table->spectrum.energyEV);
    free(table->spectrum.intensity);
    free(table->points);
    free(table->ionizationPoints);
    free(table->slice[0].rows);
    free(table->slice[1].rows);
    free(table->blend);
    *table = (MetalRateTable){0};
}

EQCLibStatus EQCLibSetRunParamUVBackground(bool use, const char *directory)
{
    if(directory != NULL && *directory == '\0') return EQC_INVALID_INPUT;
    if(!use){
        free(SpectrumDirectory);
        free(RateTable);
        FreeMetalTable(&MetalTable);
        SpectrumDirectory = NULL;
        RateTable = NULL;
        Initialized = false;
        Enabled = false;
        return EQC_SUCCESS;
    }

    const char *path = directory == NULL ? DefaultDirectory : directory;
    if(Enabled && strcmp(SpectrumDirectory,path) == 0) return EQC_SUCCESS;
    const size_t length = strlen(path) + 1;
    char *copy = malloc(length);
    if(copy == NULL) return EQC_INITIALIZATION_FAILED;
    memcpy(copy,path,length);
    free(SpectrumDirectory);
    SpectrumDirectory = copy;
    Initialized = false;
    Enabled = true;
    return EQC_SUCCESS;
}

EQCLibStatus EQCLibSetRunParamUVBackgroundMaximumEnergy(double energyEV)
{
    if(!isfinite(energyEV) || energyEV < 0.0) return EQC_INVALID_INPUT;
    if(MaximumEnergyEV != energyEV){
        MaximumEnergyEV = energyEV;
        Initialized = false;
    }
    return EQC_SUCCESS;
}

static double ColumnAtIndex(size_t index)
{
    return index == 0 ? 0.0 : pow(10.0,MinimumLogColumn + (index-1)*LogColumnStep);
}

static int CompareDouble(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static bool SetQuadraturePoint(QuadraturePoint *point, const double *sourceLogEnergy,
                               size_t left, double logEnergy, double halfWidth, bool highOnly,
                               UVBQuadratureKind kind)
{
    const double energy = exp(logEnergy);
    const double crossSection[CEQ_UVB_SPECIES] = {
        PhotoionizationCrossSectionHI(energy),
        PhotoionizationCrossSectionHeI(energy),
        PhotoionizationCrossSectionHeII(energy)
    };
    const double photonWeight = halfWidth * (2.0*M_PI/PLANK_H_CGS);
    point->left = left;
    point->fraction = (logEnergy - sourceLogEnergy[left]) /
                      (sourceLogEnergy[left+1] - sourceLogEnergy[left]);
    point->opacity = crossSection[CEQ_UVB_HI] + HeliumAbundance*crossSection[CEQ_UVB_HEI];
    if(kind == MetalIonizationQuadrature){
        double low, high;
        CEQMetalXrayMixture(energy,&point->metalOpacity,&low,&high);
        CEQMetalPhotoionizationCrossSections(energy,point->metalIonization,NULL);
        for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++){
            point->metalIonization[i] *= photonWeight;
            if(!isfinite(point->metalIonization[i])) return false;
        }
        return true;
    }
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
        point->weight.ionization[i] = photonWeight*crossSection[i];
        point->weight.heating[i] = point->weight.ionization[i] *
                                   fmax(0.0,energy - Thresholds[i])*EV_TO_ERG;
        if(highOnly && energy < Thresholds[i]+PhotoelectronSplitEV)
            point->weight.heating[i] = 0.0;
        if(!isfinite(point->weight.ionization[i]) || !isfinite(point->weight.heating[i]))
            return false;
    }
    if(kind == ShieldedQuadrature){
        double low, high;
        CEQMetalXrayMixture(energy,&point->metalOpacity,&low,&high);
        point->metalWeightLow = photonWeight*low*EV_TO_ERG;
        point->metalWeightHigh = photonWeight*high*EV_TO_ERG;
        point->fastElectrons = 0;
        for(size_t i = 0; i < CEQ_UVB_SPECIES; i++)
            if(energy >= Thresholds[i]+PhotoelectronSplitEV) point->fastElectrons |= 1U << i;
    }
    return true;
}

static bool BuildQuadrature(const CEQUVBSpectrum *spectrum, bool highOnly, UVBQuadratureKind kind,
                            QuadraturePoint **points, size_t *pointCount)
{
    const size_t count = spectrum->energyCount;
    enum { ExtraEdges = 2 + CEQ_UVB_SPECIES + 1 };
    const size_t metalEdges = kind == AtomicQuadrature ? 0 : CEQMetalXrayBoundaries(NULL);
    const size_t ionizationEdges = kind == MetalIonizationQuadrature ? CEQMetalPhotoionizationBoundaries(NULL) : 0;
    const size_t extra = ExtraEdges + CEQ_UVB_SPECIES + metalEdges + ionizationEdges;
    if(count < 2 || count > SIZE_MAX/sizeof(double)-extra) return false;
    const double kineticMinimum = highOnly ? PhotoelectronSplitEV : 0.0;
    const double lower = fmax(Thresholds[CEQ_UVB_HI]+kineticMinimum,spectrum->energyEV[0]);
    double upper = spectrum->energyEV[count-1];
    if(MaximumEnergyEV > 0.0) upper = fmin(upper,MaximumEnergyEV);
    if(lower >= upper){
        *points = NULL;
        *pointCount = 0;
        return true;
    }

    double *sourceLogEnergy = malloc(count*sizeof(*sourceLogEnergy));
    double *edges = malloc((count+extra)*sizeof(*edges));
    if(sourceLogEnergy == NULL || edges == NULL){
        free(sourceLogEnergy);
        free(edges);
        return false;
    }
    size_t edgeCount = 0;
    edges[edgeCount++] = log(lower);
    edges[edgeCount++] = log(upper);
    for(size_t i = 0; i < count; i++){
        sourceLogEnergy[i] = log(spectrum->energyEV[i]);
        if(spectrum->energyEV[i] > lower && spectrum->energyEV[i] < upper)
            edges[edgeCount++] = sourceLogEnergy[i];
    }
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
        const double boundary = Thresholds[i]+kineticMinimum;
        if(boundary > lower && boundary < upper)
            edges[edgeCount++] = log(boundary);
    }
    /* Preserve the original quadrature grid for the total rates. */
    if(!highOnly && PhotoelectronSplitEV > lower && PhotoelectronSplitEV < upper)
        edges[edgeCount++] = log(PhotoelectronSplitEV);
    if(metalEdges != 0){
        double boundaries[metalEdges+ionizationEdges];
        CEQMetalXrayBoundaries(boundaries);
        if(ionizationEdges != 0) CEQMetalPhotoionizationBoundaries(boundaries+metalEdges);
        for(size_t i = 0; i < metalEdges+ionizationEdges; i++)
            if(boundaries[i] > lower && boundaries[i] < upper) edges[edgeCount++] = log(boundaries[i]);
    }
    if(kind == ShieldedQuadrature){
        for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
            const double boundary = Thresholds[i]+PhotoelectronSplitEV;
            if(boundary > lower && boundary < upper) edges[edgeCount++] = log(boundary);
        }
    }
    qsort(edges,edgeCount,sizeof(*edges),CompareDouble);

    size_t samples = 0;
    for(size_t i = 1; i < edgeCount; i++)
        samples += 2*(size_t)ceil((edges[i]-edges[i-1])/MaxLogEnergyStep);
    QuadraturePoint *mesh = samples <= SIZE_MAX/sizeof(*mesh) ? malloc(samples*sizeof(*mesh)) : NULL;
    if(mesh == NULL){
        free(sourceLogEnergy);
        free(edges);
        return false;
    }

    /* Two interior Gauss–Legendre points per log-energy subinterval. Interior
       points take the correct side of every threshold discontinuity. */
    const double node = 0.57735026918962576451;
    size_t next = 0, left = 0;
    bool ok = true;
    for(size_t i = 1; i < edgeCount && ok; i++){
        const double range = edges[i]-edges[i-1];
        if(range == 0.0) continue;
        const size_t subdivisions = (size_t)ceil(range/MaxLogEnergyStep);
        const double step = range/subdivisions;
        const double halfStep = 0.5*step;
        while(left+2 < count && sourceLogEnergy[left+1] <= edges[i-1]) left++;
        for(size_t j = 0; j < subdivisions && ok; j++){
            const double middle = edges[i-1] + (j+0.5)*step;
            ok = SetQuadraturePoint(&mesh[next++],sourceLogEnergy,left,middle-node*halfStep,halfStep,highOnly,kind) &&
                 SetQuadraturePoint(&mesh[next++],sourceLogEnergy,left,middle+node*halfStep,halfStep,highOnly,kind);
        }
    }
    free(sourceLogEnergy);
    free(edges);
    if(!ok){
        free(mesh);
        return false;
    }
    *points = mesh;
    *pointCount = next;
    return true;
}

static double Attenuation(double opticalDepth)
{
    if(opticalDepth < 1.e-18) return 1.0;
    /* Avoid entering GSL's underflow handler for negligible transmissions. */
    if(opticalDepth >= 690.0) return 0.0;
    return gsl_sf_expint_E2(opticalDepth);
}

static double LogBlend(double first, double second, double fraction)
{
    if(fraction == 0.0 || first == second) return first;
    if(fraction == 1.0) return second;
    /* When either endpoint is zero, interpolate the physical rate instead. */
    if(first == -INFINITY) return second + log(fraction);
    if(second == -INFINITY) return first + log1p(-fraction);
    return first + fraction*(second-first);
}

static double InterpolateIntensity(const double *values, const double *logValues,
                                   const QuadraturePoint *point)
{
    const size_t left = point->left;
    if(values[left] == 0.0 || values[left+1] == 0.0)
        return (1.0-point->fraction)*values[left] + point->fraction*values[left+1];
    return exp(LogBlend(logValues[left],logValues[left+1],point->fraction));
}

static double LogRate(double rate)
{
    return rate > 0.0 ? log(rate) : -INFINITY;
}

static void FreeQuadratureGrid(QuadratureGrid *grid)
{
    free(grid->points);
    free(grid->transmission);
    free(grid->integrands);
}

static bool MakeQuadratureGrid(const CEQUVBSpectrum *spectrum, bool highOnly,
                               const double columns[ColumnCount], QuadratureGrid *grid)
{
    if(!BuildQuadrature(spectrum,highOnly,AtomicQuadrature,&grid->points,&grid->count)) return false;
    const size_t count = grid->count;
    if(count == 0) return true;
    if(count > SIZE_MAX/ColumnCount/sizeof(*grid->transmission) ||
       count > SIZE_MAX/sizeof(*grid->integrands)) return false;
    grid->transmission = malloc(ColumnCount*count*sizeof(*grid->transmission));
    grid->integrands = malloc(count*sizeof(*grid->integrands));
    if(grid->transmission == NULL || grid->integrands == NULL) return false;

    for(size_t n = 0; n < ColumnCount; n++)
        for(size_t j = 0; j < count; j++)
            grid->transmission[n*count+j] = Attenuation(columns[n]*grid->points[j].opacity);
    return true;
}

static void SetIntegrands(QuadratureGrid *grid, const double *intensity, const double *logIntensity)
{
    for(size_t j = 0; j < grid->count; j++){
        const QuadraturePoint *point = &grid->points[j];
        UVBPrimaryMoments *integrand = &grid->integrands[j];
        const double jnu = InterpolateIntensity(intensity,logIntensity,point);
        for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
            integrand->ionization[i] = jnu*point->weight.ionization[i];
            integrand->heating[i] = jnu*point->weight.heating[i];
        }
    }
}

static bool IntegrateRateTable(const CEQUVBSpectrum *spectrum, UVBTableRow *table, FILE *diagnostics)
{
    double columns[ColumnCount];
    for(size_t n = 0; n < ColumnCount; n++) columns[n] = ColumnAtIndex(n);
    QuadratureGrid total = {0};
    QuadratureGrid high = {0};
    double *logIntensity = malloc(spectrum->energyCount*sizeof(*logIntensity));
    /* Separate meshes preserve primary rates at their original quadrature
       points and place the photoelectron cutoffs exactly in the high rates. */
    bool ok = logIntensity != NULL && MakeQuadratureGrid(spectrum,false,columns,&total) &&
              MakeQuadratureGrid(spectrum,true,columns,&high);
    for(size_t z = 0; z < spectrum->redshiftCount && ok; z++){
        const double *intensity = spectrum->intensity + z*spectrum->energyCount;
        for(size_t i = 0; i < spectrum->energyCount; i++) logIntensity[i] = LogRate(intensity[i]);
        SetIntegrands(&total,intensity,logIntensity);
        SetIntegrands(&high,intensity,logIntensity);

        for(size_t n = 0; n < ColumnCount && ok; n++){
            UVBPrimaryMoments rates = {0};
            double highPower[CEQ_UVB_SPECIES] = {0};
            for(size_t j = 0; j < total.count; j++){
                const double attenuated = total.transmission[n*total.count+j];
                for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
                    const double heat = attenuated*total.integrands[j].heating[i];
                    rates.ionization[i] += attenuated*total.integrands[j].ionization[i];
                    rates.heating[i] += heat;
                }
            }
            for(size_t j = 0; j < high.count; j++){
                const double attenuated = high.transmission[n*high.count+j];
                for(size_t i = 0; i < CEQ_UVB_SPECIES; i++)
                    highPower[i] += attenuated*high.integrands[j].heating[i];
            }

            UVBTableRow *row = &table[z*ColumnCount+n];
            for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
                if(!isfinite(rates.ionization[i]) || !isfinite(rates.heating[i]) || !isfinite(highPower[i])){
                    ok = false;
                    break;
                }
                row->logRates.ionization[i] = LogRate(rates.ionization[i]);
                row->logRates.heating[i] = LogRate(rates.heating[i]);
                /* Independent quadrature truncation can place the ratio slightly above one. */
                row->highFraction[i] = rates.heating[i] > 0.0 ? fmin(1.0,highPower[i]/rates.heating[i]) : 0.0;
                if(diagnostics != NULL)
                    fprintf(diagnostics,"%.8g %.8e %s %.12e %.12e %.8g\n",
                            spectrum->redshifts[z],columns[n],SpeciesNames[i],
                            rates.ionization[i],rates.heating[i],row->highFraction[i]);
            }
        }
    }

    FreeQuadratureGrid(&total);
    FreeQuadratureGrid(&high);
    free(logIntensity);
    return ok;
}

static void MatchIonizationNodes(MetalRateTable *table)
{
    size_t k = 0;
    for(size_t j = 0; j < table->ionizationPointCount; j++){
        QuadraturePoint *ion = &table->ionizationPoints[j];
        while(k < table->pointCount){
            const QuadraturePoint *absorption = &table->points[k];
            if(absorption->left > ion->left ||
               (absorption->left == ion->left && absorption->fraction >= ion->fraction)) break;
            k++;
        }
        ion->absorptionIndex = SIZE_MAX;
        if(k < table->pointCount){
            const QuadraturePoint *absorption = &table->points[k];
            if(absorption->left == ion->left && absorption->fraction == ion->fraction &&
               absorption->opacity == ion->opacity && absorption->metalOpacity == ion->metalOpacity)
                ion->absorptionIndex = k;
        }
    }
}

/* Keep each mesh's original summation order. Sharing attenuation at matching
   nodes preserves the previous rates and avoids a second E2 evaluation. */
static bool IntegrateMetalSlices(const MetalRateTable *table, const MetalSlice slices[], size_t count)
{
    const CEQUVBSpectrum *spectrum = &table->spectrum;
    const size_t points = table->pointCount;
    const size_t ionizationPoints = table->ionizationPointCount;
    const size_t maxPoints = SIZE_MAX/count/sizeof(MetalRateRow);
    if(points > maxPoints || ionizationPoints > maxPoints-points) return false;
    MetalRateRow *integrands = malloc(count*(points+ionizationPoints)*sizeof(*integrands));
    double *logIntensity = malloc(spectrum->energyCount*sizeof(*logIntensity));
    double *transmissionAtAbsorption = points == 0 ? NULL : malloc(points*sizeof(*transmissionAtAbsorption));
    if(integrands == NULL || logIntensity == NULL || (points != 0 && transmissionAtAbsorption == NULL)){
        free(integrands);
        free(logIntensity);
        free(transmissionAtAbsorption);
        return false;
    }
    MetalRateRow *ionizationIntegrands = integrands+points*count;
    for(size_t z = 0; z < count; z++){
        const double *intensity = spectrum->intensity+slices[z].redshift*spectrum->energyCount;
        for(size_t i = 0; i < spectrum->energyCount; i++) logIntensity[i] = LogRate(intensity[i]);
        for(size_t j = 0; j < points; j++){
            const QuadraturePoint *point = &table->points[j];
            const double jnu = InterpolateIntensity(intensity,logIntensity,point);
            MetalRateRow terms = {0};
            for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
                const size_t energyBand = point->fastElectrons & (1U << i) ? MetalAtomHigh : MetalAtomLow;
                terms.value[MetalGamma+i] = jnu*point->weight.ionization[i];
                terms.value[energyBand+i] = jnu*point->weight.heating[i];
            }
            terms.value[MetalLow] = jnu*point->metalWeightLow;
            terms.value[MetalHigh] = jnu*point->metalWeightHigh;
            integrands[j*count+z] = terms;
        }
        for(size_t j = 0; j < ionizationPoints; j++){
            const QuadraturePoint *point = &table->ionizationPoints[j];
            const double jnu = InterpolateIntensity(intensity,logIntensity,point);
            for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++)
                ionizationIntegrands[j*count+z].ionization[i] = jnu*point->metalIonization[i];
        }
    }
    free(logIntensity);

    double columns[ColumnCount];
    for(size_t n = 0; n < ColumnCount; n++) columns[n] = ColumnAtIndex(n);
    bool ok = true;
    /* If the spectrum stops before metal absorption, rates depend only on N. */
    const size_t metalCount = points == 0 ? 1 : MetalRatioCount;
    for(size_t m = 0; m < metalCount && ok; m++){
        for(size_t n = 0; n < ColumnCount && ok; n++){
            MetalRateRow totals[2] = {0};
            for(size_t j = 0; j < points; j++){
                const QuadraturePoint *point = &table->points[j];
                const double transmission = Attenuation(columns[n]*(point->opacity+MetalRatios[m]*point->metalOpacity));
                transmissionAtAbsorption[j] = transmission;
                if(transmission == 0.0) continue;
                for(size_t z = 0; z < count; z++)
                    for(size_t k = 0; k < MetalMoments; k++)
                        totals[z].value[k] += transmission*integrands[j*count+z].value[k];
            }
            for(size_t j = 0; j < ionizationPoints; j++){
                const QuadraturePoint *point = &table->ionizationPoints[j];
                const double transmission = point->absorptionIndex != SIZE_MAX ?
                    transmissionAtAbsorption[point->absorptionIndex] :
                    Attenuation(columns[n]*(point->opacity+MetalRatios[m]*point->metalOpacity));
                if(transmission == 0.0) continue;
                for(size_t z = 0; z < count; z++)
                    for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++)
                        totals[z].ionization[i] += transmission*ionizationIntegrands[j*count+z].ionization[i];
            }
            for(size_t z = 0; z < count; z++){
                MetalRateRow *row = &slices[z].rows[m*ColumnCount+n];
                if(points != 0)
                    for(size_t k = 0; k < MetalMoments; k++){
                        if(!isfinite(totals[z].value[k])) ok = false;
                        row->value[k] = LogRate(totals[z].value[k]);
                    }
                if(ionizationPoints != 0)
                    for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++){
                        if(!isfinite(totals[z].ionization[i])) ok = false;
                        row->ionization[i] = LogRate(totals[z].ionization[i]);
                    }
            }
        }
    }
    free(integrands);
    free(transmissionAtAbsorption);
    return ok;
}

static bool SelectMetalRedshift(MetalRateTable *table, InterpolationPoint position)
{
    if(table->pointCount == 0 && table->ionizationPointCount == 0) return true;
    const size_t last = table->spectrum.redshiftCount-1;
    const size_t lower = position.index == last && last > 0 ? last-1 : position.index;
    const size_t count = last == 0 ? 1 : 2;
    const size_t rowCount = table->pointCount == 0 ? ColumnCount : MetalRatioCount*ColumnCount;
    MetalSlice next[2] = {0}, missing[2];
    size_t missingCount = 0;
    bool ok = true;
    for(size_t i = 0; i < count; i++){
        next[i].redshift = lower+i;
        for(size_t k = 0; k < count; k++)
            if(table->slice[k].redshift == next[i].redshift) next[i].rows = table->slice[k].rows;
        if(next[i].rows == NULL){
            next[i].rows = malloc(rowCount*sizeof(*next[i].rows));
            if(next[i].rows == NULL){
                ok = false;
                break;
            }
            missing[missingCount++] = next[i];
        }
    }
    if(ok && missingCount != 0) ok = IntegrateMetalSlices(table,missing,missingCount);
    if(!ok){
        for(size_t i = 0; i < missingCount; i++) free(missing[i].rows);
        return false;
    }

    for(size_t i = 0; i < count; i++){
        if(table->slice[i].rows != next[0].rows && table->slice[i].rows != next[1].rows)
            free(table->slice[i].rows);
        table->slice[i] = next[i];
    }
    if(position.fraction == 0.0){
        table->current = next[position.index-lower].rows;
    }else{
        /* Blend once per global redshift update to keep evaluation read-only
           and avoid interpolating two redshift planes in every cell. */
        for(size_t r = 0; r < rowCount; r++){
            if(table->pointCount != 0)
                for(size_t k = 0; k < MetalMoments; k++)
                    table->blend[r].value[k] = LogBlend(next[0].rows[r].value[k],next[1].rows[r].value[k],position.fraction);
            if(table->ionizationPointCount != 0)
                for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++)
                    table->blend[r].ionization[i] = LogBlend(next[0].rows[r].ionization[i],next[1].rows[r].ionization[i],position.fraction);
        }
        table->current = table->blend;
    }
    return true;
}

EQCLibStatus CEQInitializeUVBackground(int check)
{
    if(!Enabled) return EQC_SUCCESS;

    CEQUVBSpectrum spectrum;
    const EQCLibStatus status = CEQReadUVBackgroundSpectrum(SpectrumDirectory,&spectrum);
    if(status != EQC_SUCCESS) return status;

    UVBTableRow *table = NULL;
    if(spectrum.redshiftCount <= SIZE_MAX/ColumnCount/sizeof(*table))
        table = malloc(spectrum.redshiftCount*ColumnCount*sizeof(*table));
    FILE *diagnostics = NULL;
    if(check == 1 && table != NULL){
        diagnostics = fopen(EQC_DIAGNOSTIC_DIRECTORY "/UVBackgroundRates.dat","w");
        if(diagnostics == NULL) perror("EQCLib could not open UVBackgroundRates.dat");
        else fputs("# z N_HI_cm^-2 absorber Gamma_side_s^-1 photoelectron_power_side_erg_s^-1 fraction_K_ge_100eV\n",
                   diagnostics);
    }
    bool ok = table != NULL && (check != 1 || diagnostics != NULL) &&
              IntegrateRateTable(&spectrum,table,diagnostics);
    if(diagnostics != NULL){
        if(ferror(diagnostics)) ok = false;
        if(fclose(diagnostics) != 0) ok = false;
    }

    MetalRateTable metals = {.spectrum = spectrum};
    const double upper = spectrum.energyEV[spectrum.energyCount-1];
    if(ok && upper > 30.0 && (MaximumEnergyEV == 0.0 || MaximumEnergyEV > 30.0))
        ok = BuildQuadrature(&spectrum,false,ShieldedQuadrature,&metals.points,&metals.pointCount);
    /* A separate mesh preserves the original H/He and BCM92 integration
       points while resolving the additional ionic absorption edges. */
    if(ok)
        ok = BuildQuadrature(&spectrum,false,MetalIonizationQuadrature,&metals.ionizationPoints,&metals.ionizationPointCount);
    if(ok && (metals.pointCount != 0 || metals.ionizationPointCount != 0)){
        MatchIonizationNodes(&metals);
        const size_t rows = metals.pointCount == 0 ? ColumnCount : MetalRatioCount*ColumnCount;
        metals.blend = malloc(rows*sizeof(*metals.blend));
        ok = metals.blend != NULL && SelectMetalRedshift(&metals,(InterpolationPoint){0});
    }
    if(!ok){
        FreeMetalTable(&metals);
        free(table);
        return EQC_INITIALIZATION_FAILED;
    }

    free(RateTable);
    FreeMetalTable(&MetalTable);
    MetalTable = metals;
    RateTable = table;
    CurrentRedshift = (InterpolationPoint){0};
    Initialized = true;
    return EQC_SUCCESS;
}

EQCLibStatus EQCLibTrySetUVBackgroundRedshift(double z)
{
    if(!Initialized) return EQC_INITIALIZATION_FAILED;
    const double *redshifts = MetalTable.spectrum.redshifts;
    const size_t redshiftCount = MetalTable.spectrum.redshiftCount;
    if(!isfinite(z) || z < redshifts[0] || z > redshifts[redshiftCount-1]) return EQC_INVALID_INPUT;

    size_t lower = 0, upper = redshiftCount-1;
    while(upper-lower > 1){
        const size_t middle = lower + (upper-lower)/2;
        if(z < redshifts[middle]) upper = middle;
        else lower = middle;
    }
    const InterpolationPoint requested = z == redshifts[upper] ? (InterpolationPoint){.index = upper} : (InterpolationPoint){
        .index = lower,
        .fraction = (z-redshifts[lower])/(redshifts[upper]-redshifts[lower])
    };
    if(requested.index == CurrentRedshift.index && requested.fraction == CurrentRedshift.fraction)
        return EQC_SUCCESS;
    if(!SelectMetalRedshift(&MetalTable,requested)) return EQC_INITIALIZATION_FAILED;
    CurrentRedshift = requested;
    return EQC_SUCCESS;
}

size_t CEQUVBPersistentBytes(void)
{
    if(RateTable == NULL) return 0;
    const CEQUVBSpectrum *spectrum = &MetalTable.spectrum;
    const size_t grids = (MetalTable.slice[0].rows != NULL) + (MetalTable.slice[1].rows != NULL) +
                         (MetalTable.blend != NULL);
    const size_t rows = MetalTable.pointCount == 0 ? ColumnCount : MetalRatioCount*ColumnCount;
    return ColumnCount*spectrum->redshiftCount*sizeof(*RateTable) +
           (MetalTable.pointCount+MetalTable.ionizationPointCount)*sizeof(*MetalTable.points) +
           (spectrum->redshiftCount + spectrum->energyCount + spectrum->redshiftCount*spectrum->energyCount)*sizeof(double) +
           grids*rows*sizeof(MetalRateRow);
}

static InterpolationPoint ColumnInterpolation(double column)
{
    const double first = 1.e12;
    if(column < first) return (InterpolationPoint){.fraction = column/first};
    if(column >= MaximumColumn) return (InterpolationPoint){.index = ColumnCount-1};
    const double position = (log10(column)-MinimumLogColumn)/LogColumnStep;
    const size_t lower = (size_t)position + 1;
    if(lower >= ColumnCount-1) return (InterpolationPoint){.index = ColumnCount-1};
    return (InterpolationPoint){.index = lower, .fraction = position-(lower-1)};
}

static UVBTableRow InterpolateColumn(size_t redshift, InterpolationPoint column)
{
    const UVBTableRow *a = &RateTable[redshift*ColumnCount+column.index];
    if(column.fraction == 0.0) return *a;
    const UVBTableRow *b = a+1;
    UVBTableRow result;
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
        result.logRates.ionization[i] = LogBlend(a->logRates.ionization[i],b->logRates.ionization[i],column.fraction);
        result.logRates.heating[i] = LogBlend(a->logRates.heating[i],b->logRates.heating[i],column.fraction);
        result.highFraction[i] = a->highFraction[i] + column.fraction*(b->highFraction[i]-a->highFraction[i]);
    }
    return result;
}

static CEQUVBRates EvaluateHemisphere(InterpolationPoint position)
{
    const InterpolationPoint redshift = CurrentRedshift;
    UVBTableRow row = InterpolateColumn(redshift.index,position);
    if(redshift.fraction != 0.0){
        const UVBTableRow upper = InterpolateColumn(redshift.index+1,position);
        for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
            row.logRates.ionization[i] = LogBlend(row.logRates.ionization[i],upper.logRates.ionization[i],redshift.fraction);
            row.logRates.heating[i] = LogBlend(row.logRates.heating[i],upper.logRates.heating[i],redshift.fraction);
            row.highFraction[i] += redshift.fraction*(upper.highFraction[i]-row.highFraction[i]);
        }
    }
    CEQUVBRates rates = {0};
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
        rates.ionization[i] = exp(row.logRates.ionization[i]);
        rates.heating[i] = exp(row.logRates.heating[i]);
        rates.heatingHigh[i] = row.highFraction[i]*rates.heating[i];
    }
    return rates;
}

static bool ValidColumn(double column)
{
    return isfinite(column) && column >= 0.0;
}

static InterpolationPoint MetalInterpolation(double ratio)
{
    size_t lower = 0, upper = MetalRatioCount-1;
    while(upper-lower > 1){
        const size_t middle = lower+(upper-lower)/2;
        if(ratio < MetalRatios[middle]) upper = middle;
        else lower = middle;
    }
    if(ratio == MetalRatios[upper]) return (InterpolationPoint){.index = upper};
    if(ratio == MetalRatios[lower]) return (InterpolationPoint){.index = lower};
    const double base = 1.0+MetalRatios[lower];
    return (InterpolationPoint){
        .index = lower,
        .fraction = log1p((ratio-MetalRatios[lower])/base)/log1p((MetalRatios[upper]-MetalRatios[lower])/base)
    };
}

static MetalRateRow InterpolateMetalColumn(size_t metalIndex, InterpolationPoint column)
{
    const MetalRateRow *a = &MetalTable.current[metalIndex*ColumnCount+column.index];
    if(column.fraction == 0.0) return *a;
    const MetalRateRow *b = a+1;
    MetalRateRow result = {0};
    for(size_t k = 0; k < MetalMoments; k++)
        result.value[k] = LogBlend(a->value[k],b->value[k],column.fraction);
    return result;
}

static CEQUVBRates EvaluateMetalHemisphere(InterpolationPoint position, double ratio, InterpolationPoint metallicity)
{
    MetalRateRow row = InterpolateMetalColumn(metallicity.index,position);
    if(metallicity.fraction != 0.0){
        const MetalRateRow upper = InterpolateMetalColumn(metallicity.index+1,position);
        for(size_t k = 0; k < MetalMoments; k++)
            row.value[k] = LogBlend(row.value[k],upper.value[k],metallicity.fraction);
    }
    CEQUVBRates rates = {
        .metalPowerLow = ratio*exp(row.value[MetalLow]),
        .metalPowerHigh = ratio*exp(row.value[MetalHigh])
    };
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
        rates.ionization[i] = exp(row.value[MetalGamma+i]);
        rates.heatingHigh[i] = exp(row.value[MetalAtomHigh+i]);
        rates.heating[i] = exp(row.value[MetalAtomLow+i])+rates.heatingHigh[i];
    }
    return rates;
}

static double MetalIonizationAtColumn(const MetalRateRow *row, size_t species, double fraction)
{
    return fraction == 0.0 ? row->ionization[species] :
        LogBlend(row->ionization[species],(row+1)->ionization[species],fraction);
}

static CEQUVBRates EvaluateShieldedHemisphere(double column, double ratio, InterpolationPoint metallicity)
{
    const InterpolationPoint position = ColumnInterpolation(column);
    CEQUVBRates rates = ratio == 0.0 || MetalTable.pointCount == 0 ? EvaluateHemisphere(position) :
        EvaluateMetalHemisphere(position,ratio,metallicity);
    if(MetalTable.ionizationPointCount == 0) return rates;

    const MetalRateRow *row = &MetalTable.current[metallicity.index*ColumnCount+position.index];
    for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++){
        double value = MetalIonizationAtColumn(row,i,position.fraction);
        if(metallicity.fraction != 0.0){
            const double upper = MetalIonizationAtColumn(row+ColumnCount,i,position.fraction);
            value = LogBlend(value,upper,metallicity.fraction);
        }
        rates.metalIonization[i] = exp(value);
    }
    return rates;
}

static bool FiniteRates(const CEQUVBRates *rates)
{
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++)
        if(!isfinite(rates->ionization[i]) || !isfinite(rates->heating[i])) return false;
    for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++)
        if(!isfinite(rates->metalIonization[i])) return false;
    return isfinite(rates->metalPowerLow) && isfinite(rates->metalPowerHigh);
}

static double MetalRatio(double Z)
{
    return fmin(MetalRatios[MetalRatioCount-1],Z/ASPLUND_2009_SOLAR_METALLICITY);
}

EQCLibStatus EQCLibTryGetUVBackgroundHemisphereRates(double column, double ionization[3], double heating[3])
{
    if(ionization == NULL || heating == NULL) return EQC_INVALID_INPUT;
    CEQUVBRates rates;
    const EQCLibStatus status = CEQEvaluateUVBackgroundHemisphere(column,&rates);
    if(status != EQC_SUCCESS) return status;
    memcpy(ionization,rates.ionization,sizeof(rates.ionization));
    memcpy(heating,rates.heating,sizeof(rates.heating));
    return EQC_SUCCESS;
}

EQCLibStatus EQCLibTryGetUVBackgroundHemisphereRatesWithMetallicity(
    double column, double Z, double ionization[3], double heating[3], double metalPower[2])
{
    if(ionization == NULL || heating == NULL || metalPower == NULL) return EQC_INVALID_INPUT;
    CEQUVBRates rates;
    const EQCLibStatus status = CEQEvaluateUVBackgroundHemisphereWithMetallicity(column,Z,&rates);
    if(status != EQC_SUCCESS) return status;
    memcpy(ionization,rates.ionization,sizeof(rates.ionization));
    memcpy(heating,rates.heating,sizeof(rates.heating));
    metalPower[0] = rates.metalPowerLow;
    metalPower[1] = rates.metalPowerHigh;
    return EQC_SUCCESS;
}

EQCLibStatus CEQEvaluateUVBackgroundHemisphere(double column, CEQUVBRates *rates)
{
    return CEQEvaluateUVBackgroundHemisphereWithMetallicity(column,0.0,rates);
}

EQCLibStatus CEQEvaluateUVBackgroundHemisphereWithMetallicity(double column, double Z, CEQUVBRates *rates)
{
    if(rates == NULL || !ValidColumn(column) || !isfinite(Z) || Z < 0.0) return EQC_INVALID_INPUT;
    if(!Initialized) return EQC_INITIALIZATION_FAILED;
    const double ratio = MetalRatio(Z);
    const InterpolationPoint metallicity = MetalInterpolation(MetalTable.pointCount == 0 ? 0.0 : ratio);
    const CEQUVBRates result = EvaluateShieldedHemisphere(column,ratio,metallicity);
    if(!FiniteRates(&result)) return EQC_INVALID_INPUT;
    *rates = result;
    return EQC_SUCCESS;
}

EQCLibStatus CEQEvaluateUVBackground(const double columns[2], CEQUVBRates *rates)
{
    return CEQEvaluateUVBackgroundWithMetallicity(columns,0.0,rates);
}

EQCLibStatus CEQEvaluateUVBackgroundWithMetallicity(const double columns[2], double Z, CEQUVBRates *rates)
{
    if(rates == NULL) return EQC_INVALID_INPUT;
    if(!Enabled){
        *rates = (CEQUVBRates){0};
        return EQC_SUCCESS;
    }
    if(columns == NULL || !ValidColumn(columns[0]) || !ValidColumn(columns[1]) || !isfinite(Z) || Z < 0.0)
        return EQC_INVALID_INPUT;
    if(!Initialized) return EQC_INITIALIZATION_FAILED;

    const double ratio = MetalRatio(Z);
    const InterpolationPoint metallicity = MetalInterpolation(MetalTable.pointCount == 0 ? 0.0 : ratio);
    CEQUVBRates total = EvaluateShieldedHemisphere(columns[0],ratio,metallicity);
    const CEQUVBRates other = columns[0] == columns[1] ? total : EvaluateShieldedHemisphere(columns[1],ratio,metallicity);
    for(size_t i = 0; i < CEQ_UVB_SPECIES; i++){
        total.ionization[i] += other.ionization[i];
        total.heating[i] += other.heating[i];
        total.heatingHigh[i] += other.heatingHigh[i];
    }
    total.metalPowerLow += other.metalPowerLow;
    total.metalPowerHigh += other.metalPowerHigh;
    for(size_t i = 0; i < CEQ_UVB_METAL_SPECIES; i++)
        total.metalIonization[i] += other.metalIonization[i];
    if(!FiniteRates(&total)) return EQC_INVALID_INPUT;
    *rates = total;
    return EQC_SUCCESS;
}
