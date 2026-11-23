#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"

static const char settingsPath[] = "thermochemistry_mode_test.conf";
static const double accuracy = 1.0e-6;

static void WriteSettings(const char *contents)
{
    FILE *file = fopen(settingsPath,"w");
    assert(file != NULL);
    assert(fputs(contents,file) >= 0);
    assert(fclose(file) == 0);
}

static void CheckModeSettings(void)
{
    EQCLibChemistryMode mode = EQC_CHEMISTRY_NONEQUILIBRIUM;
    assert(EQCLibTryLoadChemistryMode(NULL,&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_EQUILIBRIUM);
    assert(EQCLibTryLoadChemistryMode(NULL,NULL) == EQC_INVALID_INPUT);

    assert(EQCLibTryParseChemistryMode("nonequilibrium",&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);
    assert(EQCLibTryParseChemistryMode("non-equilibrium",&mode) == EQC_INVALID_INPUT);
    assert(EQCLibTryParseChemistryMode(NULL,&mode) == EQC_INVALID_INPUT);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);
    assert(EQCLibTryParseChemistryMode("equilibrium",&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_EQUILIBRIUM);
    assert(EQCLibTryParseChemistryMode("h2_nonequilibrium",&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_H2_NONEQUILIBRIUM);

    WriteSettings("n_H = 10\r\n  chemistry_mode = nonequilibrium  # enable ODEs\r\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);

    fputs("[expected test diagnostics] Duplicate, invalid and malformed chemistry_mode settings.\n",stderr);
    WriteSettings("chemistry_mode = equilibrium\nchemistry_mode = nonequilibrium\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_INVALID_INPUT);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);
    WriteSettings("chemistry_mode = unknown\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_INVALID_INPUT);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);
    WriteSettings("chemistry_mode nonequilibrium\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_INVALID_INPUT);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);

    WriteSettings("# chemistry_mode = nonequilibrium\nn_H = 10\nother_parameter 100\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_EQUILIBRIUM);
    assert(remove(settingsPath) == 0);
    mode = EQC_CHEMISTRY_NONEQUILIBRIUM;
    fputs("[expected test diagnostics] Missing chemistry mode configuration file.\n",stderr);
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_INITIALIZATION_FAILED);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);
}

static void SetCellThermodynamics(StructEQCData *state)
{
    state->nel = EQCLibGetElectronNumberDensity(state);
    state->Rho = EQCLibGetTotalMassDensity(state);
    state->MeanMolecularWeight = EQCLibGetMeanMolecularWeight(state);
    state->DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(state);
    state->U = EQCLibConvertTtoU(state);
    state->invT = 1.0/state->T;
    state->logT = log10(state->T);
}

static StructEQCData InitialCell(void)
{
    StructEQCData state = {
        .GlobalID = 17,
        .t = 2000.0*YEAR_CGS,
        .dt = 1000.0*YEAR_CGS,
        .T = 8000.0,
        .n = {[ipHYDROGEN] = 10.0},
        .natom = {[ipHYDROGEN] = 5.0},
        .nion = {[ipHYDROGEN] = 5.0},
        .Du = 1.0,
    };
    SetCellThermodynamics(&state);
    return state;
}

static StructEQCData InitialSolarCell(void)
{
    StructEQCData state = InitialCell();
    const double nH = state.n[ipHYDROGEN];
    state.n[ipHELIUM] = state.natom[ipHELIUM] = 0.1*nH;
    state.n[ipCARBON] = state.natom[ipCARBON] = 3.019952e-4*nH;
    state.n[ipOXYGEN] = state.natom[ipOXYGEN] = 4.570882e-4*nH;
    const double nSi = 3.548134e-6*nH;
    state.n[ipSILICON] = nSi;
    state.natom[ipSILICON] = 0.25*nSi;
    state.nion[ipSILICON] = 0.5*nSi;
    state.nSiIII = 0.25*nSi;
    SetCellThermodynamics(&state);
    return state;
}

static void CheckAdvance(const StructEQCData *initial, const StructEQCData *expected,
                         StructEQCData *actual, EQCLibChemistryMode mode)
{
    assert(EQCLibTryAdvanceCell(actual,accuracy,mode) == EQC_SUCCESS);
    assert(actual->U == expected->U && actual->U < initial->U);
    assert(actual->T == expected->T);
    assert(actual->natom[ipHYDROGEN] == expected->natom[ipHYDROGEN]);
    assert(actual->nion[ipHYDROGEN] == expected->nion[ipHYDROGEN]);
    assert(actual->nel == expected->nel);
    assert(actual->H_tot == expected->H_tot && actual->L_tot == expected->L_tot);
    assert(actual->Du < 0.0);
    const double energyError = fabs(actual->dt*actual->Du - (actual->U-initial->U));
    assert(energyError <= 4.0*DBL_EPSILON*initial->U);
    assert(actual->GlobalID == initial->GlobalID && actual->t == initial->t);
    assert(actual->dt == initial->dt && actual->Rho == initial->Rho);
}

static void CheckDispatch(void)
{
    const StructEQCData initial = InitialCell();
    StructEQCData expected = initial;
    double du;
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&expected,accuracy,&du) == EQC_SUCCESS);
    expected.U += expected.dt*du;

    EQCLibChemistryMode mode;
    assert(EQCLibTryLoadChemistryMode(NULL,&mode) == EQC_SUCCESS);
    StructEQCData equilibrium = initial;
    CheckAdvance(&initial,&expected,&equilibrium,mode);

    expected = initial;
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&expected,accuracy) == EQC_SUCCESS);
    WriteSettings("chemistry_mode = nonequilibrium\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_SUCCESS);
    assert(remove(settingsPath) == 0);
    StructEQCData nonequilibrium = initial;
    CheckAdvance(&initial,&expected,&nonequilibrium,mode);
    assert(nonequilibrium.nion[ipHYDROGEN] > equilibrium.nion[ipHYDROGEN]);
}

static void CheckSilicon(void)
{
    const StructEQCData initial = InitialSolarCell();
    EQCLibChemistryMode mode;
    WriteSettings("output_interval = 100\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_EQUILIBRIUM);
    StructEQCData equilibrium = initial;
    assert(EQCLibTryAdvanceCell(&equilibrium,accuracy,mode) == EQC_SUCCESS);

    WriteSettings("output_interval = 100\nchemistry_mode = nonequilibrium\n");
    assert(EQCLibTryLoadChemistryMode(settingsPath,&mode) == EQC_SUCCESS);
    assert(mode == EQC_CHEMISTRY_NONEQUILIBRIUM);
    assert(remove(settingsPath) == 0);
    StructEQCData nonequilibrium = initial;
    nonequilibrium.natom[ipSILICON] = 0.0; /* Reconstruct from supplied ions and total. */
    assert(EQCLibTryAdvanceCell(&nonequilibrium,accuracy,mode) == EQC_SUCCESS);
    assert(nonequilibrium.U < initial.U);
    assert(nonequilibrium.L_SiII35 > 0.0);
    assert(nonequilibrium.n[ipSILICON] == initial.n[ipSILICON]);
    assert(nonequilibrium.nion[ipSILICON] == initial.nion[ipSILICON]);
    assert(nonequilibrium.nSiIII == initial.nSiIII);
    const double nSi = initial.n[ipSILICON];
    const double siliconSum = nonequilibrium.natom[ipSILICON] +
                              nonequilibrium.nion[ipSILICON] + nonequilibrium.nSiIII;
    assert(fabs(siliconSum-nSi) <= 4.0*DBL_EPSILON*nSi);
    const double electrons = EQCLibGetElectronNumberDensity(&nonequilibrium);
    assert(fabs(nonequilibrium.nel-electrons) <= 4.0*DBL_EPSILON*electrons);

    /* The initial electrons all come from Si. They seed collisional H ionization. */
    StructEQCData siliconElectrons = initial;
    siliconElectrons.T = 20000.0;
    siliconElectrons.natom[ipHYDROGEN] = siliconElectrons.n[ipHYDROGEN];
    siliconElectrons.nion[ipHYDROGEN] = 0.0;
    SetCellThermodynamics(&siliconElectrons);
    assert(EQCLibTryAdvanceCell(&siliconElectrons,accuracy,mode) == EQC_SUCCESS);
    assert(siliconElectrons.nion[ipHYDROGEN] > 0.0);
}

static void CheckUnchanged(StructEQCData *state, double tolerance, EQCLibChemistryMode mode,
                           EQCLibStatus expected)
{
    unsigned char before[sizeof(*state)];
    memcpy(before,state,sizeof(before));
    assert(EQCLibTryAdvanceCell(state,tolerance,mode) == expected);
    assert(memcmp(before,state,sizeof(before)) == 0);
}

static void CheckCOCompatibility(void)
{
    StructEQCData initial = InitialSolarCell();
    initial.n[ipSILICON] = initial.natom[ipSILICON] = initial.nion[ipSILICON] = initial.nSiIII = 0.0;
    SetCellThermodynamics(&initial);
    StructEQCData state = initial;
    assert(EQCLibTryAdvanceCell(&state,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM) == EQC_SUCCESS);

    state = initial;
    EQCLibSetRunParamCOFormation(true,1.0);
    CheckUnchanged(&state,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM,EQC_INVALID_INPUT);

    unsigned char before[sizeof(state)];
    memcpy(before,&state,sizeof(before));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(before)) == 0);
    assert(EQCLibTryAdvanceCell(&state,accuracy,EQC_CHEMISTRY_EQUILIBRIUM) == EQC_SUCCESS);

    StructEQCData pureHydrogen = InitialCell();
    assert(EQCLibTryAdvanceHydrogenChemistry(&pureHydrogen,accuracy) == EQC_SUCCESS);
    EQCLibSetRunParamCOFormation(false,1.0);

    state = initial;
    state.nmol[ipCO] = 0.5*state.n[ipCARBON];
    state.natom[ipCARBON] -= state.nmol[ipCO];
    state.natom[ipOXYGEN] -= state.nmol[ipCO];
    SetCellThermodynamics(&state);
    CheckUnchanged(&state,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM,EQC_INVALID_INPUT);
}

static void CheckUpdateInputs(void)
{
    assert(EQCLibTryAdvanceCell(NULL,accuracy,EQC_CHEMISTRY_EQUILIBRIUM) == EQC_INVALID_INPUT);
    const StructEQCData initial = InitialCell();
    StructEQCData state = initial;
    CheckUnchanged(&state,accuracy,(EQCLibChemistryMode)42,EQC_INVALID_INPUT);

    const EQCLibChemistryMode modes[] = {
        EQC_CHEMISTRY_EQUILIBRIUM,EQC_CHEMISTRY_NONEQUILIBRIUM,EQC_CHEMISTRY_H2_NONEQUILIBRIUM
    };
    for(size_t i = 0; i < sizeof(modes)/sizeof(modes[0]); i++){
        state = initial;
        CheckUnchanged(&state,NAN,modes[i],EQC_INVALID_INPUT);
        state.dt = -1.0;
        CheckUnchanged(&state,accuracy,modes[i],EQC_INVALID_INPUT);
        state.dt = 0.0;
        CheckUnchanged(&state,accuracy,modes[i],EQC_SUCCESS);
        state.dt = initial.dt;
        state.Rho = 0.0;
        CheckUnchanged(&state,accuracy,modes[i],EQC_INVALID_INPUT);
    }
}

static void CheckMolecularNonequilibrium(void)
{
    const EQCLibChemistryMode mode = EQC_CHEMISTRY_H2_NONEQUILIBRIUM;
    StructEQCData initial = InitialCell();
    const double nH = initial.n[ipHYDROGEN];
    initial.T = 100.0;
    initial.G0 = 1.0;
    initial.dt = 0.5/3.3e-11;
    initial.natom[ipHYDROGEN] = 0.2*nH;
    initial.nion[ipHYDROGEN] = 0.0;
    initial.nmol[ipH2] = 0.4*nH;
    SetCellThermodynamics(&initial);

    EQCLibSetRunParamH2Formation(true,1.0);
    StructEQCData state = initial;
    assert(EQCLibTryCalculateFixedH2CoolingHeatingRates(&state,accuracy) == EQC_SUCCESS);
    assert(state.nmol[ipH2] == initial.nmol[ipH2] && state.U == initial.U);
    assert(EQCLibTryAdvanceCell(&state,accuracy,mode) == EQC_SUCCESS);
    /* At low T and zero dust and column, optically thin photo loss dominates. */
    assert(fabs(state.nmol[ipH2]/initial.nmol[ipH2] - exp(-0.5)) < 1.e-6);
    const double hydrogen = state.natom[ipHYDROGEN] + state.nion[ipHYDROGEN] + 2.0*state.nmol[ipH2];
    assert(fabs(hydrogen/nH-1.0) <= 8.0*DBL_EPSILON);
    assert(fabs(state.nel-EQCLibGetElectronNumberDensity(&state)) <= accuracy*nH);
    assert(fabs(EQCLibConvertTtoU(&state)/state.U-1.0) <= accuracy);
    assert(state.Du == (state.U-initial.U)/initial.dt);

    const double molecular = state.nmol[ipH2];
    state.U *= 0.9;
    const double limitedEnergy = state.U;
    double rate;
    assert(EQCLibTryNetCoolingRateFixedH2(&state,&rate) == EQC_SUCCESS);
    assert(state.nmol[ipH2] == molecular && state.U == limitedEnergy);
    assert(fabs(EQCLibConvertTtoU(&state)/state.U-1.0) <= accuracy);
    assert(rate == (state.H_tot-state.L_tot)/state.Rho);

    state = initial;
    state.nmol[ipH2] = nH;
    CheckUnchanged(&state,accuracy,mode,EQC_INVALID_INPUT);
    EQCLibSetRunParamCOFormation(true,1.0);
    state = initial;
    CheckUnchanged(&state,accuracy,mode,EQC_INVALID_INPUT);
    EQCLibSetRunParamCOFormation(false,1.0);
    EQCLibSetRunParamH2Formation(false,1.0);
}

int main(int argc, char **argv)
{
    assert(argc == 2);
    CheckModeSettings();
    assert(EQCLibTryInitialization(0,argv[1]) == EQC_SUCCESS);
    CheckDispatch();
    CheckUpdateInputs();
    CheckSilicon();
    CheckCOCompatibility();
    CheckMolecularNonequilibrium();
    puts("Common thermochemistry updates and chemistry mode configuration checks passed.");
    return 0;
}
