#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "ChemicalEquilibriumCooling.h"
#include "Astro.h"
#include "Constants.h"
#include "CollisionalGasPhaseChemistry.h"
#include "CosmicRay.h"
#include "EquilibriumAbundance.h"
#include "FineStructureLines.h"
#include "HIIregionLineCooling.h"
#include "MoleculeChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"

static StructCEQData GasState(double nH, double temperature, double column)
{
    StructCEQData state = {
        .dt = 10.0 * YEAR_CGS,
        .T = temperature,
        .Nt = column,
        .Z = SOLAR_METALLICITY,
        .Zd = SOLAR_METALLICITY,
        .G0 = 1.7,
        .X0 = 1.0,
        .n = {nH, 0.1*nH, 3.019952e-4*nH, 4.570882e-4*nH, 3.548134e-6*nH},
    };
    state.Td = CEQApproximateDustTemperature(state.Nt,state.G0);
    for(int i = 0; i < NumberChemicalElement; i++){
        state.natom[i] = state.n[i];
        state.Rho += state.n[i] * ChemicalElementMass[i];
    }
    state.MeanMolecularWeight = CEQLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = CEQLibGetMeanDegreeOfFreedom(&state);
    state.U = CEQLibConvertTtoU(&state);
    return state;
}

static void CheckAtomicEquilibrium(void)
{
    const double nH = 0.01584893;
    StructCEQData state = GasState(nH,10083.56,1.0e18);
    state.nel = 0.9*nH;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-9) == CEQ_SUCCESS);
    assert(fabs(state.nel/nH - 0.27455813) < 1.0e-6);
    assert(fabs(CEQLibGetElectronNumberDensity(&state) - state.nel) <= 1.0e-9 * state.nel);
    assert(state.nmol[ipH2] == 0.0 && state.nmol[ipCO] == 0.0);
    assert(isfinite(state.H_tot) && isfinite(state.L_tot));
}

static void CheckHeIIIRecombinationCooling(void)
{
    /* Hydrogenic scaling: Lambda_HeIII(T) = 8 Lambda_HII(T/4) for case A and case B. */
    const bool caseB[] = {true, false};
    for(size_t c = 0; c < sizeof(caseB)/sizeof(caseB[0]); c++){
        CEQLibiSetRunParamRecombinationCaseB(caseB[c]);
        const double temperatures[] = {1.0e4, 3.0e4, 1.0e5, 1.0e6};
        for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
            StructCEQData helium = {.T = temperatures[i], .nel = 1.0, .nHeIII = 1.0};
            StructCEQData hydrogen = {.T = 0.25*temperatures[i], .nel = 1.0, .nion = {[ipHYDROGEN] = 1.0}};
            const double expected = 8.0*CoolingRecombination(&hydrogen);
            /* The He II threshold is 4.0018 times the H threshold, so the scaling holds to ~1e-3. */
            assert(fabs(CoolingHeIIIRecombination(&helium)-expected) <= 2.0e-3*expected);
        }
    }
    CEQLibiSetRunParamRecombinationCaseB(true);
}

static void CheckIonizedGasEquationOfState(void)
{
    const double nH = 1.0;
    StructCEQData state = {
        .T = 1.0e6,
        .Rho = nH * ChemicalElementMass[ipHYDROGEN],
        .n = {[ipHYDROGEN] = nH},
        .nel = nH,
    };
    assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-9) == CEQ_SUCCESS);
    assert(state.nion[ipHYDROGEN] > 0.999 * nH);
    assert(fabs(state.DegreeOfFreedom - 1.5) < 1.0e-12);
    assert(fabs(state.Gamma - 5.0/3.0) < 1.0e-12);

    const double expectedU = 1.5 * (nH + state.nel) * BOLTZMANN_CONSTANT_CGS * state.T / state.Rho;
    assert(fabs(CEQLibConvertTtoU(&state) - expectedU) < 1.0e-12 * expectedU);
    state.U = expectedU;
    assert(fabs(CEQLibConvertUtoT(&state) - state.T) < 1.0e-12 * state.T);
}

static void CheckFailureStatus(void)
{
    StructCEQData state = GasState(0.001,13242.50,1.0e18);
    unsigned char original[sizeof(state)];
    memcpy(original,&state,sizeof(state));
    double rate = 42.0;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == CEQ_NO_EQUILIBRIUM);
    assert(memcmp(original,&state,sizeof(state)) == 0);
    assert(CEQLibTryNetCoolingRateChemcalEquilibrium(&state,&rate) == CEQ_NO_EQUILIBRIUM);
    assert(rate == 42.0 && memcmp(original,&state,sizeof(state)) == 0);
    assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&rate) == CEQ_NO_EQUILIBRIUM);
    assert(rate == 42.0 && memcmp(original,&state,sizeof(state)) == 0);
    assert(isnan(CEQLibNetCoolingRateChemcalEquilibrium(&state)));
    assert(memcmp(original,&state,sizeof(state)) == 0);
}

static void CheckZeroElectronSeed(void)
{
    StructCEQData state = GasState(1.0,100.0,1.0e19);
    state.G0 = state.Z = state.Zd = 0.0;
    for(int i = 1; i < NumberChemicalElement; i++) state.n[i] = state.natom[i] = 0.0;
    CEQLibiSetRunParamCosmicRayIonization(false);
    assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-8) >= 0);
    assert(fabs(state.nel - 0.00240514) < 1.0e-7);
    CEQLibiSetRunParamCosmicRayIonization(true);
}

static void CheckLowIonizationXraySearch(void)
{
    StructCEQData state = {.T = 100.0, .Nt = 1.0e24, .X0 = 1.0, .nel = 1.0, .n = {1.0e3}};
    CEQLibiSetRunParamCosmicRayIonization(false);
    CEQLibiSetRunParamGrainSurfaceProcesses(false);
    /* A seed inside the fitted X-ray range keeps that search domain. */
    assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-8) == CEQ_NO_EQUILIBRIUM);
    state.nel = 0.05;
    assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-8) == CEQ_SUCCESS);
    assert(state.nel > 0.0 && state.nel < 1.0e-4*state.n[ipHYDROGEN]);
    assert(fabs(state.nel-state.nion[ipHYDROGEN]) <= 1.0e-8*state.nel);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
}

static void CheckAtomicThermalNumerics(void)
{
    const double accuracy = 1.0e-8;
    const double samples[][2] = {{870963.6,128.5617}, {2.290868,10004.58}};
    for(size_t i = 0; i < sizeof(samples)/sizeof(samples[0]); i++){
        StructCEQData state = GasState(samples[i][0],samples[i][1],1.0e19);
        state.G0 = 1700.0;
        state.X0 = 1000.0;
        state.nel = 0.9*state.n[ipHYDROGEN];
        state.dt = 1000.0*YEAR_CGS;
        assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,accuracy,&state.Du) == CEQ_SUCCESS);
        const double change = state.dt*state.Du;
        const double source = state.dt*(state.H_tot-state.L_tot)/state.Rho;
        assert(fabs(change-source) <= accuracy*fmax(fabs(change),fabs(source)));
        state.U += change;
        assert(fabs(CEQLibConvertTtoU(&state)-state.U) < 1.0e-12*state.U);
    }
}

/* Isochoric cooling over a fixed interval, split into equal implicit steps. */
static double CoolIsochorically(double nH, double temperature, double duration, int steps)
{
    StructCEQData state = GasState(nH,temperature,1.0e20*sqrt(nH));
    state.nel = 1.0e-3*nH;
    assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-9) == CEQ_SUCCESS);
    state.U = CEQLibConvertTtoU(&state);
    state.dt = duration/steps;
    for(int step = 0; step < steps; step++){
        assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-4,&state.Du) == CEQ_SUCCESS);
        state.U += state.dt*state.Du;
        assert(isfinite(state.U) && state.U > 0.0);
    }
    return state.T;
}

static void CheckHighTemperatureCIEExtrapolation(void)
{
    const double tableMaximum = pow(10.0,8.5);
    const double maximumTemperature = 1.0e9;
    const double nH = 0.1;
    const double justAboveTable = tableMaximum*(1.0 + 1.0e-6);
    StructCEQData endpoint = GasState(nH,tableMaximum,1.0e19);
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&endpoint,1.0e-8) == CEQ_SUCCESS);
    const double endpointCoefficient = endpoint.L_res/(endpoint.nel*nH);
    assert(endpointCoefficient > 0.0);

    const double temperatures[] = {justAboveTable, maximumTemperature};
    for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
        StructCEQData state = endpoint;
        state.T = temperatures[i];
        assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == CEQ_SUCCESS);
        const double expected = endpointCoefficient*sqrt(state.T/tableMaximum);
        const double coefficient = state.L_res/(state.nel*nH);
        assert(fabs(coefficient-expected) <= 1.0e-12*expected);
    }

    const double duration = 1.0e5*YEAR_CGS;
    const double cooled = CoolIsochorically(nH,maximumTemperature,duration,1);
    assert(cooled > tableMaximum && cooled < maximumTemperature);
    assert(CoolIsochorically(nH,justAboveTable,duration,1) < tableMaximum);

    StructCEQData outside = endpoint;
    outside.T = nextafter(maximumTemperature,INFINITY);
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&outside,1.0e-8) == CEQ_INVALID_INPUT);
}

static void CheckLowMetallicityCIECooling(void)
{
    const double minimumTableMetallicity = 1.0e-3*SOLAR_METALLICITY;
    const double metallicityRatio = 1.0e-3;
    const double temperatures[] = {9999.99,1.0e4,1.0e5};
    double coefficients[3];
    for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
        StructCEQData state = {.T = temperatures[i], .Z = minimumTableMetallicity, .n = {1.0,0.1}, .nel = 1.0};
        assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == CEQ_SUCCESS);
        const double coefficient = state.L_res/(state.nel*state.n[ipHYDROGEN]);
        assert(coefficient > 0.0);
        coefficients[i] = coefficient;

        state.Z *= metallicityRatio;
        assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == CEQ_SUCCESS);
        const double traceCoefficient = state.L_res/(state.nel*state.n[ipHYDROGEN]);
        const double expected = metallicityRatio*coefficient;
        assert(fabs(traceCoefficient-expected) <= 1.0e-12*expected);

        state.Z = 0.0;
        assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == CEQ_SUCCESS);
        assert(state.L_res == 0.0);
    }
    assert(fabs(coefficients[0]-coefficients[1]) <= 1.0e-4*coefficients[1]);
}

static void CheckHotGasTimeStepConvergence(void)
{
    const double duration = 2.0e5*YEAR_CGS;

    /* When the step is shorter than the cooling time, one step matches a fine integration. */
    const double diffuse[][2] = {{0.1,1.0e6}, {0.1,1.0e7}};
    for(size_t i = 0; i < sizeof(diffuse)/sizeof(diffuse[0]); i++){
        const double coarse = CoolIsochorically(diffuse[i][0],diffuse[i][1],duration,1);
        const double fine = CoolIsochorically(diffuse[i][0],diffuse[i][1],duration,1000);
        assert(fabs(coarse-fine) <= 0.01*fine);
    }

    /* A step longer than the cooling time succeeds but is not accurate: at nH = 1-10
       a single step leaves the gas 1.6-3 times hotter than the fine integration.
       Refining the step must still converge toward the fine result. */
    const double dense[] = {1.0, 10.0};
    for(size_t i = 0; i < sizeof(dense)/sizeof(dense[0]); i++){
        const double fine = CoolIsochorically(dense[i],1.0e6,duration,1000);
        const double error10 = fabs(CoolIsochorically(dense[i],1.0e6,duration,10)-fine);
        const double error100 = fabs(CoolIsochorically(dense[i],1.0e6,duration,100)-fine);
        assert(fine < 1.0e4);
        assert(error100 < error10);
        assert(error100 <= 0.05*fine);
    }
}

static void CheckMolecularBranches(void)
{
    StructCEQData initial = GasState(1000.0,50.0,1.0e24);
    initial.G0 = 588.5900374411234;
    initial.Zd = 5.0e-4 * SOLAR_METALLICITY;
    initial.n[ipCARBON] = 0.30;
    initial.n[ipOXYGEN] = 0.46;
    initial.Td = CEQApproximateDustTemperature(initial.Nt,initial.G0);
    const double seeds[] = {0.0,0.05,0.20};
    const double roots[] = {0.00218258,0.05198586,0.29849391};
    for(size_t i = 0; i < sizeof(seeds)/sizeof(seeds[0]); i++){
        StructCEQData state = initial;
        state.nmol[ipCO] = seeds[i];
        assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-9) >= 0);
        assert(fabs(state.nmol[ipCO] - roots[i]) < 1.0e-6);
        double hydrogen = state.natom[ipHYDROGEN] + state.nion[ipHYDROGEN] + 2.0*state.nmol[ipH2];
        double carbon = state.natom[ipCARBON] + state.nion[ipCARBON] + state.nCIII + state.nmol[ipCO];
        double oxygen = state.natom[ipOXYGEN] + state.nion[ipOXYGEN] + state.nOIII + state.nmol[ipCO];
        assert(fabs(hydrogen - state.n[ipHYDROGEN]) < 1.0e-9 * state.n[ipHYDROGEN]);
        assert(fabs(carbon - state.n[ipCARBON]) < 1.0e-9 * state.n[ipCARBON]);
        assert(fabs(oxygen - state.n[ipOXYGEN]) < 1.0e-9 * state.n[ipOXYGEN]);
    }
}

static void CheckOxygenLimit(void)
{
    StructCEQData state = GasState(1000.0,50.0,1.0e22);
    state.n[ipCARBON] = 0.30;
    state.n[ipOXYGEN] = 0.10;
    assert(CEQLibChemicalReactionEquilibrium(&state,1.0e-8) == CEQ_NO_EQUILIBRIUM);
    assert(state.nmol[ipCO] == 0.0);
}

static void CheckMolecularXrayLocalSearch(void)
{
    StructCEQData state = GasState(1.0e3,200.0,1.0e21);
    state.nel = 1.0;
    const double accuracy = 1.0e-9;
    assert(CEQLibChemicalReactionEquilibrium(&state,accuracy) == CEQ_SUCCESS);

    /* Reference densities from the full search with the same electron seed. */
    const double ne = 0.297561308014;
    const double nH2 = 499.902127618;
    const double nCO = 0.00386665466694;
    assert(fabs(state.nel - ne) <= 1.0e-8 * ne);
    assert(fabs(state.nmol[ipH2] - nH2) <= 1.0e-8 * nH2);
    assert(fabs(state.nmol[ipCO] - nCO) <= 1.0e-8 * nCO);
    assert(fabs(CEQLibGetElectronNumberDensity(&state) - state.nel) <= accuracy * state.nel);
}

static void CheckEnergyAndEquationOfState(void)
{
    StructCEQData state = GasState(1000.0,6000.0,1.0e22);
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-10) == CEQ_SUCCESS);
    const double oldH2 = state.nmol[ipH2];
    const double oldU = CEQLibConvertTtoU(&state);
    state.U = oldU;
    assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&state.Du) == CEQ_SUCCESS);
    assert(state.U == oldU && state.Du < 0.0);
    assert(state.nmol[ipH2] > oldH2);
    state.U += state.dt * state.Du;
    const double source = state.dt * (state.H_tot - state.L_tot) / state.Rho;
    assert(fabs(state.U - oldU - source) <= 2.0e-8 * fabs(source));
    assert(fabs(CEQLibConvertTtoU(&state) - state.U) <= 1.0e-12 * state.U);

    const double solvedTemperature = state.T;
    double rate;
    assert(CEQLibTryNetCoolingRateChemcalEquilibrium(&state,&rate) == CEQ_SUCCESS);
    assert(fabs(state.T - solvedTemperature) < 1.0e-7 * solvedTemperature);
    assert(fabs(CEQLibConvertTtoU(&state) - state.U) <= 2.0e-8 * state.U);
    assert(isfinite(rate));
}

static void CheckThermalBranchAtCOBoundary(void)
{
    const double nH = pow(10.0,0.5);
    StructCEQData state = GasState(nH,pow(10.0,4.5),1.0e22*sqrt(nH));
    state.nel = 1.0e-3*nH;
    state.dt = 1.0e6*YEAR_CGS;
    assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-4,&state.Du) == CEQ_SUCCESS);
    /* A coarser search can skip the hot root next to the CO regime change. */
    assert(state.T > 1.0e4 && state.T < 1.05e4);
}

static void CheckCIIHydrogenCollisionCooling(void)
{
    StructCEQData line = {.T = 250.0};
    line.nion[ipCARBON] = 0.01;
    line.nmol[ipH2] = 10.0;
    const double atJoin = CoolingCII158um(&line);
    line.T = nextafter(line.T,INFINITY);
    assert(atJoin > 0.0 && fabs(CoolingCII158um(&line)-atJoin) <= 1.0e-12*atJoin);

    const double nH = 1.0e3;
    const double accuracy = 1.0e-4;
    const double initialTemperatures[] = {pow(10.0,5.75),1.0e6,pow(10.0,6.25)};
    for(size_t i = 0; i < sizeof(initialTemperatures)/sizeof(initialTemperatures[0]); i++){
        StructCEQData state = GasState(nH,initialTemperatures[i],1.0e20*pow(nH,0.25));
        state.nel = 1.0e-3*nH;
        state.dt = 1.0e4*YEAR_CGS;
        for(int step = 0; step < 3; step++){
            assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,accuracy,&state.Du) == CEQ_SUCCESS);
            const double change = state.dt*state.Du;
            const double source = state.dt*(state.H_tot-state.L_tot)/state.Rho;
            assert(fabs(change-source) <= accuracy*fmax(fabs(change),fabs(source)));
            state.U += change;
        }
    }
}

static void CheckH2CoolingFitLimits(void)
{
    const double self = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6000.0,0,0,1,0,0);
    const double neutral = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6000.0,0,1,1,1,0);
    const double charged = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(10000.0,1,0,1,0,1);
    const double aboveSelf = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6001.0,0,0,1,0,0);
    const double aboveNeutral = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6001.0,0,1,1,1,0);
    const double aboveCharged = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(10001.0,1,0,1,0,1);
    assert(self > 0.0 && neutral > self && charged > self);
    assert(fabs(aboveSelf-self) <= 1.0e-12*self);
    assert(fabs(aboveNeutral-neutral) <= 1.0e-12*neutral);
    assert(fabs(aboveCharged-charged) <= 1.0e-12*charged);
}

static void CheckDisabledFineStructureLines(void)
{
    double line1 = 1.0, line2 = 2.0, line3 = 3.0;
    assert(CoolingCI(1.0,100.0,0.0,0.0,1.0,0.0,&line1,&line2,&line3) == 0.0);
    assert(line1 == 0.0 && line2 == 0.0 && line3 == 0.0);
    assert(CoolingCI610um(1.0,100.0,0.0,0.0,1.0,0.0) == 0.0);
    assert(CoolingCI230um(1.0,100.0,0.0,0.0,1.0,0.0) == 0.0);
    assert(CoolingCI370um(1.0,100.0,0.0,0.0,1.0,0.0) == 0.0);

    line1 = 1.0;
    line2 = 2.0;
    line3 = 3.0;
    assert(CoolingSiI(1.0,100.0,1.0,0.0,&line1,&line2,&line3) == 0.0);
    assert(line1 == 0.0 && line2 == 0.0 && line3 == 0.0);
    assert(CoolingSiI130um(1.0,100.0,1.0,0.0) == 0.0);
    assert(CoolingSiI45um(1.0,100.0,1.0,0.0) == 0.0);
    assert(CoolingSiI68um(1.0,100.0,1.0,0.0) == 0.0);
}

static void CheckHydrogenFlux(const StructCEQData *state, double accuracy)
{
    const double nH = state->n[ipHYDROGEN];
    const double nHI = state->natom[ipHYDROGEN];
    const double nH2 = state->nmol[ipH2];
    const double dustTemperature = state->Td > 0.0 ?
                                   state->Td : CEQApproximateDustTemperature(state->Nt,state->G0);
    const double formed = nH*H2FormationRate(state->T,state->Zd,dustTemperature)*nHI;
    const double doppler = hypot(12.90*sqrt(1.0e-4*state->T),3.0);
    const double photo = H2PhotoDissociationRate(state->G0,state->Nt,state->Nt*nH2/nH,doppler,state->Zd);
    const double destroyed = nH2*(photo+H2CollisionalDissociationRate(state->T,nH,nHI,nH2,state->nel));
    assert(formed > 0.0 && destroyed > 0.0);
    assert(fabs(formed-destroyed) <= accuracy*fmax(formed,destroyed));
}

static void CheckMolecularFlux(double dustTemperature)
{
    StructCEQData state = GasState(1000.0,75.0,1.0e20);
    state.Td = dustTemperature;
    const double accuracy = 1.0e-9;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == CEQ_SUCCESS);
    CheckHydrogenFlux(&state,accuracy);
    assert(state.H_H2 > 0.0 && isfinite(state.H_H2));

    state.nmol[ipCO] = 1.0e-11;
    assert(CoolingCORovibrationalLine(&state) > 0.0);
}

static void CheckTraceHydrogen(void)
{
    const double accuracy = 1.0e-12;
    StructCEQData initial = GasState(1000.0,50.0,1.0e20);
    initial.X0 = 0.0;
    StructCEQData atomic = initial;
    atomic.G0 = 1.0e100;
    assert(CEQSolveHydrogenEquilibrium(&atomic,accuracy) == CEQ_SUCCESS);
    assert(atomic.nmol[ipH2] < 1.0e-90);
    CheckHydrogenFlux(&atomic,accuracy);

    StructCEQData molecular = initial;
    molecular.G0 = 1.0e-100;
    assert(CEQSolveHydrogenEquilibrium(&molecular,accuracy) == CEQ_SUCCESS);
    assert(molecular.natom[ipHYDROGEN] < 1.0e-90);
    CheckHydrogenFlux(&molecular,accuracy);
}

static void CheckOxygenLineEmissivities(void)
{
    const double tolerance = 2.0e-4;
    /* Independent five-level solve using PyNeb's Z82-WFD96/Kal09 and FFT04-SZ00/SSB14 data. */
    const struct {
        double ne;
        double oxygenII;
        double oxygenIII;
        double opticalDoublet;
        double infraredDoublet;
    } reference[] = {
        {1.0,   3.5638632548e-21, 7.0532399529e-21, 1.4991843, 0.5530553},
        {100.0, 3.4963065702e-21, 6.7333964675e-21, 1.3493476, 0.7989521},
        {1.0e4, 1.7171530434e-21, 4.9690089747e-21, 0.4172276, 7.6530238},
    };
    double oII[OxygenLevelCount][OxygenLevelCount], oIII[OxygenLevelCount][OxygenLevelCount];
    for(size_t i = 0; i < sizeof(reference)/sizeof(reference[0]); i++){
        const double ne = reference[i].ne;
        const double coolingII = CoolingLinesOII(1.0e4,ne,1.0,oII)/ne;
        const double coolingIII = CoolingLinesOIII(1.0e4,ne,1.0,oIII)/ne;
        assert(fabs(coolingII/reference[i].oxygenII-1.0) < tolerance);
        assert(fabs(coolingIII/reference[i].oxygenIII-1.0) < tolerance);
        assert(fabs((oII[1][0]/oII[2][0])/reference[i].opticalDoublet-1.0) < tolerance);    /* 3729/3726 Å */
        assert(fabs((oIII[2][1]/oIII[1][0])/reference[i].infraredDoublet-1.0) < tolerance); /* 52/88 µm */
        assert(fabs(oIII[3][2]/oIII[3][1]-2.983969) < tolerance);                         /* 5007/4959 Å */
    }

    const double lowT = pow(10.0,3.7);
    const double highT = pow(10.0,4.3);
    assert(fabs(CoolingLinesOII(lowT,100.0,1.0,NULL)/9.9799349321e-21-1.0) < tolerance);
    assert(fabs(CoolingLinesOIII(highT,100.0,1.0,NULL)/1.8028688374e-18-1.0) < tolerance);
    assert(CoolingLinesOII(1.0e4,0.0,1.0,oII) == 0.0 && oII[1][0] == 0.0);
    assert(CoolingLinesOIII(1.0e4,100.0,0.0,oIII) == 0.0 && oIII[3][2] == 0.0);
    assert(isnan(CoolingLinesOII(1.0e4,-1.0,1.0,NULL)));
}

static void CheckPhotoionizedOxygenCooling(void)
{
    StructCEQData state = GasState(100.0,1.1e4,1.0e19);
    state.NLyn = 3.0e9;
    state.X0 = state.G0 = state.Zd = 0.0;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == CEQ_SUCCESS);
    assert(state.L_OII > 0.0 && state.L_OIII > 0.0);

    StructCEQData cie = state;
    CEQLibiSetRunParamHIIRegion(false);
    assert(CEQLibCalculateCoolingHeatingRates(&cie) == CEQ_SUCCESS);
    CEQLibiSetRunParamHIIRegion(true);
    assert(cie.L_OII == 0.0 && cie.L_OIII == 0.0 && cie.L_res > 0.0);

    double lines = CoolingLinesHIIregionOxygen(state.T,state.nel,state.nion[ipOXYGEN],state.nOIII);
    const double weight = (state.L_OII+state.L_OIII)/lines;
    assert(weight > 0.99 && weight <= 1.0);
    assert(fabs(state.L_res-(1.0-weight)*cie.L_res) < 1.0e-12*cie.L_res);
    const double expected = cie.L_tot - cie.L_res + state.L_res + state.L_OII + state.L_OIII;
    assert(fabs(state.L_tot-expected) < 1.0e-12*expected);

    StructCEQData transition = state;
    transition.T = 2.0e4;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&transition,1.0e-8) == CEQ_SUCCESS);
    lines = CoolingLinesHIIregionOxygen(transition.T,transition.nel,transition.nion[ipOXYGEN],transition.nOIII);
    const double transitionWeight = (transition.L_OII+transition.L_OIII)/lines;
    assert(transitionWeight > 0.0 && transitionWeight < weight);

    StructCEQData hot = state;
    hot.T = 1.0e6;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&hot,1.0e-8) == CEQ_SUCCESS);
    assert(hot.L_OII == 0.0 && hot.L_OIII == 0.0 && hot.L_res > 0.0);

    StructCEQData unilluminated = state;
    unilluminated.NLyn = 0.0;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&unilluminated,1.0e-8) == CEQ_SUCCESS);
    assert(unilluminated.L_OII == 0.0 && unilluminated.L_OIII == 0.0 && unilluminated.L_res > 0.0);

    state.U = CEQLibConvertTtoU(&state);
    state.dt = 1000.0*YEAR_CGS;
    assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&state.Du) == CEQ_SUCCESS);
    assert(state.Du < 0.0 && state.T < 1.0e4 && state.L_OII > 0.0 && state.L_OIII > 0.0);
    assert(state.nCIII > 0.0); /* Free carbon stays photoionized below the CO temperature threshold. */
}

static void CheckPhotoionizationTurnOn(void)
{
    CEQLibiSetRunParamH2Formation(false,1.0);
    CEQLibiSetRunParamCOFormation(false,1.0);
    const double nH = 100.0;
    StructCEQData initial = GasState(nH,1.0e4,1.0e19);
    initial.G0 = initial.X0 = initial.Zd = 0.0;
    initial.nel = nH;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&initial,1.0e-10) == CEQ_SUCCESS);
    assert(initial.nel < 0.1*nH);
    initial.U = CEQLibConvertTtoU(&initial);
    initial.NLyn = 1.0e-3*nH*LIGHT_C_CGS;

    StructCEQData illuminated = initial, reference = initial;
    reference.nel = nH;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&reference,1.0e-12) == CEQ_SUCCESS);
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&illuminated,1.0e-8) == CEQ_SUCCESS);
    assert(illuminated.nion[ipHYDROGEN] > 0.99*nH);
    assert(fabs(CEQLibGetElectronNumberDensity(&illuminated)-illuminated.nel) <= 1.0e-8*illuminated.nel);
    assert(fabs(illuminated.nel-reference.nel) <= 1.0e-8*reference.nel);
    assert(fabs(illuminated.nOIII-reference.nOIII) <= 1.0e-8*reference.n[ipOXYGEN]);

    initial.dt = 1.0e4*YEAR_CGS;
    StructCEQData updated = initial, tightUpdate = initial;
    tightUpdate.nel = reference.nel;
    const double accuracy = 1.0e-8;
    assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&updated,accuracy,&updated.Du) == CEQ_SUCCESS);
    assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&tightUpdate,1.0e-9,&tightUpdate.Du) == CEQ_SUCCESS);
    const double change = updated.dt*updated.Du;
    const double source = updated.dt*(updated.H_tot-updated.L_tot)/updated.Rho;
    assert(fabs(change-source) <= accuracy*fmax(fabs(change),fabs(source)));
    assert(fabs(updated.T-tightUpdate.T) <= 1.0e-8*tightUpdate.T);
}

static double RelativeDifference(double a, double b)
{
    const double scale = fmax(fabs(a),fabs(b));
    return scale > 0.0 ? fabs(a-b)/scale : 0.0;
}

static void CheckHydrogenOxygenChargeTransfer(void)
{
    CEQLibiSetRunParamH2Formation(false,1.0);
    CEQLibiSetRunParamCOFormation(false,1.0);
    const double nH = 100.0, accuracy = 1.0e-10;
    StructCEQData state = GasState(nH,8000.0,1.0e19);
    state.G0 = state.X0 = state.Zd = 0.0;
    state.nel = nH;
    state.NLyn = 1.0e-4*nH*LIGHT_C_CGS;   /* log U = -4 leaves a few per cent H I. */
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == CEQ_SUCCESS);
    assert(fabs(CEQLibGetElectronNumberDensity(&state)-state.nel) <= 1.0e-9*state.nel);

    const double T = state.T, ne = state.nel;
    const double nHI = state.natom[ipHYDROGEN], nHII = state.nion[ipHYDROGEN];
    const double nOI = state.natom[ipOXYGEN], nOII = state.nion[ipOXYGEN], nOIII = state.nOIII;
    assert(nHI > 0.01*nH && nHI < 0.1*nH);
    const double toOxygenII = nOI*nHII*OIChargeTransferRateToHII(T);
    const double toOxygenI = nOII*nHI*OIIChargeTransferRateToHI(T);
    const double fromOxygenIII = nOIII*nHI*OIIIChargeTransferRateToHI(T);

    /* Each ladder balances with the charge-transfer fluxes included. */
    const double hydrogenIonization = nHI*(HIPhotoIonizationRate(state.NLyn) + ne*HICollisionalIonizationRate(T))
                                    + toOxygenI + fromOxygenIII;
    const double hydrogenRecombination = nHII*ne*HIIRecombinationRate(T) + toOxygenII;
    assert(RelativeDifference(hydrogenIonization,hydrogenRecombination) <= 1.0e-8);
    const double oxygenIIIonization = nOI*(OIPhotoionizationRate(state.NLyn) + ne*OICollisionalIonizationRate(T))
                                    + toOxygenII;
    const double oxygenIIRecombination = nOII*ne*OIIRecombinationRate(T) + toOxygenI;
    assert(RelativeDifference(oxygenIIIonization,oxygenIIRecombination) <= 1.0e-8);
    const double oxygenIIIIonization = nOII*(OIIPhotoionizationRate(state.NLyn) + ne*OIICollisionalIonizationRate(T));
    const double oxygenIIIRecombination = nOIII*ne*OIIIRecombinationRate(T) + fromOxygenIII;
    assert(RelativeDifference(oxygenIIIIonization,oxygenIIIRecombination) <= 1.0e-8);

    /* At this H I fraction, O++ + H dominates the O++ recombination. */
    const double withoutTransfer = (OIIPhotoionizationRate(state.NLyn) + ne*OIICollisionalIonizationRate(T))
                                 / (ne*OIIIRecombinationRate(T));
    assert(nOIII/nOII < 0.25*withoutTransfer);

    /* Dust, cosmic rays and X-rays use the same coupled solve. */
    StructCEQData dusty = GasState(nH,8000.0,1.0e19);
    dusty.nel = nH;
    dusty.NLyn = state.NLyn;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&dusty,accuracy) == CEQ_SUCCESS);
    assert(fabs(CEQLibGetElectronNumberDensity(&dusty)-dusty.nel) <= 1.0e-9*dusty.nel);
}

extern double b_turb;

/* H2 destruction against formation, from the rate functions used by the solver. */
static double MolecularHydrogenImbalance(const StructCEQData *s)
{
    const double nH = s->n[ipHYDROGEN], nHI = s->natom[ipHYDROGEN], nH2 = s->nmol[ipH2];
    const double formation = nH*H2FormationRate(s->T,s->Zd,EffectiveDustTemperature(s))*nHI;
    const double doppler = sqrt(SQ(12.90*sqrt(1.e-4*s->T))+SQ(b_turb));
    const double photo = H2PhotoDissociationRateWithDustFactor(s->G0,nH2*s->Nt/nH,doppler,
                                                               DustShieldingFactorFUV(s->Nt,s->Zd));
    const double collision = H2DissociationFromCoefficients(H2DissociationCoefficients(s->T),s->nel,nHI,nH2);
    return RelativeDifference((photo+collision)*nH2,formation);
}

static void CheckMolecularChargeTransferBalance(void)
{
    /* Partly ionized molecular gas, where coupling O and H used to leave H2
       about 9% out of balance because H I moved after the H2 solve. */
    const struct StructCEQLibParameters_tag saved = CEQLibRunParameters;
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamH2Formation(true,1.0);
    CEQLibiSetRunParamCOFormation(true,1.0);
    const double nH = 100.0, accuracy = 1.0e-8;
    StructCEQData state = GasState(nH,60.0,1.0e21);
    state.nel = 1.0e-3*nH;
    state.NLyn = 1.0e-5*nH*LIGHT_C_CGS;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == CEQ_SUCCESS);
    assert(state.nmol[ipH2] > 0.4*nH && state.nion[ipHYDROGEN] > 1.0e-3*nH);
    assert(MolecularHydrogenImbalance(&state) <= 2.0*accuracy);

    const CEQHydrogenRates rates = CEQHydrogenIonizationRates(&state);
    const double T = state.T, nHI = state.natom[ipHYDROGEN], nHII = state.nion[ipHYDROGEN];
    const double ionization = nHI*rates.ionization
                            + state.nion[ipOXYGEN]*nHI*OIIChargeTransferRateToHI(T)
                            + state.nOIII*nHI*OIIIChargeTransferRateToHI(T);
    const double recombination = nHII*rates.recombination + state.natom[ipOXYGEN]*nHII*OIChargeTransferRateToHII(T);
    assert(RelativeDifference(ionization,recombination) <= 2.0*accuracy);

    /* CO is re-solved after charge transfer moves H2: solving CO again from the
       returned H2 must leave it unchanged. The old order left a 1e-3 offset. */
    const double nCO = 316.22776601683793;
    StructCEQData molecular = GasState(nCO,1000.0,1.0e20*sqrt(nCO));
    molecular.nel = 1.0e-3*nCO;
    molecular.NLyn = 1.0e-5*nCO*LIGHT_C_CGS;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&molecular,accuracy) == CEQ_SUCCESS);
    assert(molecular.nmol[ipCO] > 0.01*molecular.n[ipCARBON]);
    assert(MolecularHydrogenImbalance(&molecular) <= 2.0*accuracy);
    StructCEQData carbon = molecular;
    assert(CEQSolveCarbonEquilibrium(&carbon,molecular.nmol[ipCO],accuracy) == CEQ_SUCCESS);
    assert(RelativeDifference(carbon.nmol[ipCO],molecular.nmol[ipCO]) <= 1.0e-6);
    CEQLibRunParameters = saved;
}

static void CheckH2ShieldingUpperTail(void)
{
    const double edge = 1.016e22, delta = 1.e-6;
    const double atEdge = COShieldingFactor(edge,0.0,3.0);
    const double below = COShieldingFactor(edge/(1.0+delta),0.0,3.0);
    const double above = COShieldingFactor(edge*(1.0+delta),0.0,3.0);
    assert(RelativeDifference(atEdge,4.666e-3) <= 1.e-12);
    assert(below > atEdge && atEdge > above);
    assert(RelativeDifference(below-atEdge,atEdge-above) < 1.e-5);

    /* Match the attenuation per decade in the N_CO=0 column of Visser et al. Table 5. */
    const double nextDecade = COShieldingFactor(10.0*edge,0.0,3.0);
    assert(RelativeDifference(nextDecade/atEdge,3.938e-7/1.108e-2) <= 1.e-12);
    assert(COShieldingFactor(100.0*edge,0.0,3.0) < nextDecade);
}

static void CheckCarbonBalanceAtH2ShieldingBoundary(double nH, double metallicity, double column)
{
    const struct StructCEQLibParameters_tag saved = CEQLibRunParameters;
    CEQLibiSetRunParamHIIRegion(true);
    CEQLibiSetRunParamChargeTransfer(true);
    CEQLibiSetRunParamDustAttenuation(true);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamH2Formation(true,1.0);
    CEQLibiSetRunParamCOFormation(true,1.0);
    /* Both fixtures put the uncoupled H2 above the edge and the coupled H2 below. */
    const double upperH2Column = 1.016e22;
    const double accuracies[] = {1.e-8, 1.e-5};
    for(size_t i=0; i<sizeof(accuracies)/sizeof(accuracies[0]); i++){
        const double accuracy = accuracies[i];
        StructCEQData state = GasState(nH,30.0,column);
        state.nel = 1.e-3*nH;
        state.NLyn = 1.e-5*nH*LIGHT_C_CGS;
        state.Z *= metallicity;
        state.Zd *= metallicity;
        for(int element=ipCARBON; element<NumberChemicalElement; element++){
            state.n[element] *= metallicity;
            state.natom[element] = state.n[element];
        }
        assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == CEQ_SUCCESS);

        StructCEQData uncoupled = state;
        assert(CEQSolveHydrogenEquilibrium(&uncoupled,accuracy) == CEQ_SUCCESS);
        const double h2Column = state.Nt*(state.nmol[ipH2]/nH);
        const double uncoupledColumn = state.Nt*(uncoupled.nmol[ipH2]/nH);
        assert(h2Column < upperH2Column && uncoupledColumn > upperH2Column);

        const double nCO = state.nmol[ipCO];
        const double freeCarbon = state.natom[ipCARBON] + state.nion[ipCARBON] + state.nCIII;
        const double association = 4.e-16*pow(state.T/300.0,-0.2)*state.nmol[ipH2];
        const double photoCI = 3.1e-10*state.G0*DustShieldingFactor(state.Nt,state.Zd,3.0);
        const double formation = (association*photoCI/(association+photoCI))*freeCarbon;
        const double destruction = COPhotoDissociationRate(state.G0,state.Nt,state.Zd)*
                                   COShieldingFactor(h2Column,state.Nt*(nCO/nH),b_turb)*nCO;
        assert(RelativeDifference(formation,destruction) <= accuracy);
        assert(MolecularHydrogenImbalance(&state) <= accuracy);
        assert(RelativeDifference(state.nel,CEQLibGetElectronNumberDensity(&state)) <= accuracy);
    }
    CEQLibRunParameters = saved;
}

static void CheckMolecularHydrogenRoundoffAccuracy(void)
{
    /* At the smallest tolerance the thermal solver uses (16 eps), the local H2
       search started near the root must still meet it when it reports success. */
    const struct StructCEQLibParameters_tag saved = CEQLibRunParameters;
    CEQLibiSetRunParamH2Formation(true,1.0);
    const double accuracy = 16.0*DBL_EPSILON;
    const double temperatures[] = {30.0, 100.0, 1000.0};
    const double offsets[] = {1.0e-9, 1.0e-6, 1.0e-3, -1.0e-3, 0.1};
    for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
        StructCEQData reference = GasState(100.0,temperatures[i],1.0e21);
        reference.nel = 1.0e-4*reference.n[ipHYDROGEN];
        assert(CEQSolveHydrogenEquilibrium(&reference,1.0e-8) == CEQ_SUCCESS);
        for(size_t j = 0; j < sizeof(offsets)/sizeof(offsets[0]); j++){
            StructCEQData state = reference;
            state.nmol[ipH2] = reference.nmol[ipH2]*(1.0+offsets[j]);
            if(CEQSolveHydrogenEquilibrium(&state,accuracy) == CEQ_SUCCESS)
                assert(MolecularHydrogenImbalance(&state) <= 1.0e-12);
        }
    }
    CEQLibRunParameters = saved;
}

static void CheckChargeTransferSwitch(void)
{
    const double nH = 100.0, accuracy = 1.0e-10;
    StructCEQData coupled = GasState(nH,8000.0,1.0e19);
    coupled.G0 = coupled.X0 = coupled.Zd = 0.0;
    coupled.nel = nH;
    coupled.NLyn = 1.0e-4*nH*LIGHT_C_CGS;
    StructCEQData uncoupled = coupled;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&coupled,accuracy) == CEQ_SUCCESS);

    CEQLibiSetRunParamChargeTransfer(false);
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&uncoupled,accuracy) == CEQ_SUCCESS);
    CEQLibiSetRunParamChargeTransfer(true);

    /* Without charge transfer, the O ladder balances on photo- and electron processes alone. */
    const double T = uncoupled.T, ne = uncoupled.nel;
    const double ionization = uncoupled.nion[ipOXYGEN]*(OIIPhotoionizationRate(uncoupled.NLyn)
                                                     + ne*OIICollisionalIonizationRate(T));
    const double recombination = uncoupled.nOIII*ne*OIIIRecombinationRate(T);
    assert(RelativeDifference(ionization,recombination) <= 1.0e-8);
    assert(uncoupled.nOIII/uncoupled.nion[ipOXYGEN] > 4.0*coupled.nOIII/coupled.nion[ipOXYGEN]);
    assert(fabs(CEQLibGetElectronNumberDensity(&uncoupled)-uncoupled.nel) <= 1.0e-9*uncoupled.nel);
}

static void CheckChemistryAccuracySetting(void)
{
    const double invalid[] = {0.0, -1.0e-5, 1.0, NAN, INFINITY};
    for(size_t i = 0; i < sizeof(invalid)/sizeof(invalid[0]); i++)
        assert(CEQLibiSetRunParamChemistryAccuracy(invalid[i]) == CEQ_INVALID_INPUT);
    assert(CEQLibRunParameters.ChemistryAccuracy == 1.0e-8);

    /* A looser chemistry tolerance keeps energy updates within the energy accuracy. */
    const double nH = 100.0, accuracy = 1.0e-4;
    const double fluxes[] = {0.0, 1.0e-3*nH*LIGHT_C_CGS};
    for(size_t i = 0; i < sizeof(fluxes)/sizeof(fluxes[0]); i++){
        StructCEQData state = GasState(nH,1.0e4,1.0e19);
        state.nel = 0.5*nH;
        state.NLyn = fluxes[i];
        state.dt = 1.0e4*YEAR_CGS;
        StructCEQData loose = state;
        assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,accuracy,&state.Du) == CEQ_SUCCESS);
        assert(CEQLibiSetRunParamChemistryAccuracy(1.0e-5) == CEQ_SUCCESS);
        assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&loose,accuracy,&loose.Du) == CEQ_SUCCESS);
        assert(CEQLibiSetRunParamChemistryAccuracy(1.0e-8) == CEQ_SUCCESS);
        assert(fabs(loose.T-state.T) <= accuracy*state.T);
    }
}

static void CheckIonizingFluxNeedsHIIRegion(void)
{
    /* With H II region mode off, NLyn must not reach the chemistry or heating,
       even after the PEGASE2 tables are loaded. Otherwise C would lose its FUV
       photoionization and drop from C II to C I in diffuse gas. */
    const struct StructCEQLibParameters_tag saved = CEQLibRunParameters;
    CEQLibiSetRunParamHIIRegion(false);
    CEQLibiSetRunParamH2Formation(false,1.0);
    CEQLibiSetRunParamCOFormation(false,1.0);
    StructCEQData dark = GasState(1.0,100.0,1.0e19);
    dark.nel = 1.0e-3;
    StructCEQData lit = dark;
    lit.NLyn = 1.0e-3*LIGHT_C_CGS;
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&dark,1.0e-10) == CEQ_SUCCESS);
    assert(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&lit,1.0e-10) == CEQ_SUCCESS);
    assert(lit.nel == dark.nel && lit.H_tot == dark.H_tot && lit.L_tot == dark.L_tot);
    for(int i = 0; i < NumberChemicalElement; i++)
        assert(lit.natom[i] == dark.natom[i] && lit.nion[i] == dark.nion[i]);
    assert(dark.nion[ipCARBON] > 0.99*dark.n[ipCARBON]);
    CEQLibRunParameters = saved;
}

static void CheckPhotoionization(const char *pegaseDirectory)
{
    const double unattenuated = CIFUVPhotoionizationRate(1.7,0.0);
    const double attenuated = CIFUVPhotoionizationRate(1.7,1.0);
    assert(fabs(attenuated/unattenuated-exp(-3.0)) < 1.0e-15);
    CEQLibiSetRunParamDustAttenuation(false);
    assert(CIFUVPhotoionizationRate(1.7,1.0) == unattenuated);
    CEQLibiSetRunParamDustAttenuation(true);

    CEQLibiSetRunParamHIIRegion(true);
    assert(CEQLibTryInitialization(0,pegaseDirectory) == CEQ_SUCCESS);
    const double rate = HIPhotoIonizationRate(1.0);
    const double referenceRate = 0x1.36bcd33264bbp-58;  /* Solar PEGASE2 spectrum. */
    assert(fabs(rate-referenceRate) <= 1.0e-12*referenceRate);
    /* Solar PEGASE2 coefficients normalized by the H-ionizing flux before age averaging. */
    const double referenceHeIIRate = 5.928733248e-21;
    const double referenceHeIIHeating = 1.043341218e-31;
    assert(fabs(HeIIPhotoIonizationRate(1.0)-referenceHeIIRate) <= 1.0e-8*referenceHeIIRate);
    assert(fabs(HeatingHeIIPhotoionization(1.0)-referenceHeIIHeating) <= 1.0e-8*referenceHeIIHeating);

    struct stat info;
    assert(stat("out",&info) != 0 && stat("CEQCooling",&info) != 0);
    assert(CEQLibTryInitialization(0,"missing-pegase") == CEQ_INITIALIZATION_FAILED);
    assert(mkdir("malformed",0700) == 0);
    FILE *input = fopen("malformed/spectra_Kroupa_Z0.0001.dat","w");
    assert(input != NULL);
    for(int line = 0; line < 27; line++) fputc('\n',input);
    fputs("1 1500 0\n",input);
    assert(fclose(input) == 0);
    assert(CEQLibTryInitialization(0,"malformed") == CEQ_INITIALIZATION_FAILED);
    assert(HIPhotoIonizationRate(1.0) == rate);
    CheckPhotoionizedOxygenCooling();
    CheckPhotoionizationTurnOn();
    CheckHydrogenOxygenChargeTransfer();
    CheckMolecularChargeTransferBalance();
    CheckH2ShieldingUpperTail();
    CheckCarbonBalanceAtH2ShieldingBoundary(100.0,0.3,2.0320004172648183e22);
    CheckCarbonBalanceAtH2ShieldingBoundary(1.0,0.03,4.0608174530788829e22);
    CheckMolecularHydrogenRoundoffAccuracy();
    CheckChargeTransferSwitch();
    CheckChemistryAccuracySetting();
    CheckIonizingFluxNeedsHIIRegion();
}

int main(int argc, char **argv)
{
    assert(argc == 2);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(true);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(false);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    CEQLibiSetRunParamH2Formation(false,1.0);
    CEQLibiSetRunParamCOFormation(false,1.0);
    CEQLibInitialization(0);
    CheckDisabledFineStructureLines();

    StructCEQData probe = GasState(1.0,100.0,1.e20);
    probe.nel = 0.01*probe.n[ipHYDROGEN];
    assert(TotalCosmicRayHeatingRate(&probe) == 0.0);
    CEQLibiSetRunParamCosmicRayIonization(true);
    const double cosmicRayHeating = TotalCosmicRayHeatingRate(&probe);
    assert(isfinite(cosmicRayHeating) && cosmicRayHeating > 0.0);
    CheckAtomicEquilibrium();
    CheckHeIIIRecombinationCooling();
    CheckIonizedGasEquationOfState();
    CheckFailureStatus();
    CheckZeroElectronSeed();
    CheckLowIonizationXraySearch();
    CheckAtomicThermalNumerics();
    CheckHighTemperatureCIEExtrapolation();
    CheckLowMetallicityCIECooling();
    CheckHotGasTimeStepConvergence();

    CEQLibiSetRunParamH2Formation(true,1.0);
    CEQLibiSetRunParamCOFormation(true,1.0);
    CEQLibiSetRunParamGrainSurfaceProcesses(false);
    CEQLibiSetRunParamSoftXrayIonization(false);
    CheckMolecularBranches();
    CheckOxygenLimit();

    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CheckMolecularXrayLocalSearch();
    CheckEnergyAndEquationOfState();
    CheckThermalBranchAtCOBoundary();
    CheckCIIHydrogenCollisionCooling();
    CheckH2CoolingFitLimits();
    CheckMolecularFlux(0.0);
    CheckMolecularFlux(100.0);
    CheckTraceHydrogen();
    CheckOxygenLineEmissivities();
    CheckPhotoionization(argv[1]);
    puts("Chemical equilibrium and thermal update checks passed.");
    return 0;
}
