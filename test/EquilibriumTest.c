#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CosmicRay.h"
#include "EquilibriumAbundance.h"
#include "FineStructureLines.h"
#include "MoleculeChemistry.h"
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
    const double samples[][2] = {{870963.6,128.5617}, {2.290868,10004.58}};
    for(size_t i = 0; i < sizeof(samples)/sizeof(samples[0]); i++){
        StructCEQData state = GasState(samples[i][0],samples[i][1],1.0e19);
        state.G0 = 1700.0;
        state.X0 = 1000.0;
        state.nel = 0.9*state.n[ipHYDROGEN];
        state.dt = 1000.0*YEAR_CGS;
        const double oldU = state.U;
        assert(CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&state.Du) == CEQ_SUCCESS);
        state.U += state.dt*state.Du;
        const double source = state.dt*(state.H_tot-state.L_tot)/state.Rho;
        assert(fabs(state.U-oldU-source) < 1.0e-11*oldU);
        assert(fabs(CEQLibConvertTtoU(&state)-state.U) < 1.0e-12*state.U);
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
    CheckFailureStatus();
    CheckZeroElectronSeed();
    CheckLowIonizationXraySearch();
    CheckAtomicThermalNumerics();

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
    CheckPhotoionization(argv[1]);
    puts("Chemical equilibrium and thermal update checks passed.");
    return 0;
}
