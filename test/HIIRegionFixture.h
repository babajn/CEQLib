#ifndef CEQLIB_TEST_HII_REGION_FIXTURE_H
#define CEQLIB_TEST_HII_REGION_FIXTURE_H

#include <math.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"

typedef struct {
    double nH;
    double metallicity;
    double logIonization;
} HIIRegion;

static inline StructCEQData HIIRegionGasState(HIIRegion region, double temperature)
{
    const double nH = region.nH;
    const double scaledDensity = region.metallicity * nH;
    return (StructCEQData){
        .T = temperature,
        .Z = region.metallicity * SOLAR_METALLICITY,
        .NLyn = pow(10.0, region.logIonization) * nH * LIGHT_C_CGS,
        .Nt = 1.0e19,
        .nel = nH,
        .n = {nH, 0.1*nH, 3.019952e-4*scaledDensity, 4.570882e-4*scaledDensity, 3.548134e-6*scaledDensity},
    };
}

static inline CEQLibStatus InitializeHIIRegionFixture(const char *pegaseDirectory)
{
    /* Dust-free gas with unattenuated stellar ionizing radiation as its only source. */
    CEQLibiSetRunParamHIIRegion(true);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    return CEQLibTryInitialization(0, pegaseDirectory);
}

#endif
