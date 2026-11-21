#ifndef CEQ_MOLECULAR_FIXTURE_H
#define CEQ_MOLECULAR_FIXTURE_H

#include <math.h>
#include "Astro.h"
#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"

typedef struct {
    double nH, T, G0, Nt, logU;   /* logU = -INFINITY for no ionizing flux. */
} Condition;

static StructCEQData GasState(Condition c)
{
    const double nH = c.nH, Z = SOLAR_METALLICITY;
    StructCEQData s = {
        .T = c.T, .Nt = c.Nt, .Z = Z, .Zd = Z, .G0 = c.G0, .X0 = 1.0, .nel = 1.0e-3*nH,
        .NLyn = isfinite(c.logU) ? pow(10.0, c.logU)*nH*LIGHT_C_CGS : 0.0,
        .n = {nH, 0.1*nH, 3.019952e-4*nH, 4.570882e-4*nH, 3.548134e-6*nH},
    };
    for(int i = 0; i < NumberChemicalElement; i++){
        s.natom[i] = s.n[i];
        s.Rho += s.n[i]*ChemicalElementMass[i];
    }
    return s;
}

static inline CEQLibStatus InitializeMolecularFixture(const char *pegaseDirectory)
{
    CEQLibSetRunParamRecombinationCaseB(true);
    CEQLibSetRunParamDustAttenuation(true);
    CEQLibSetRunParamGrainSurfaceProcesses(true);
    CEQLibSetRunParamCosmicRayIonization(true);
    CEQLibSetRunParamSoftXrayIonization(true);
    CEQLibSetRunParamCIEResonanceCooling(true);
    CEQLibSetRunParamHIIRegion(true);
    CEQLibSetRunParamH2Formation(true, 1.0);
    CEQLibSetRunParamCOFormation(true, 1.0);
    CEQLibSetRunParamChargeTransfer(true);
    return CEQLibTryInitialization(0, pegaseDirectory);
}

#endif
