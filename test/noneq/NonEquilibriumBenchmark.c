#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "PhotochemicalGasPhaseChemistry.h"

enum { Trials = 5, MaximumCalls = 1 << 18 };
static const double Accuracy = 1.0e-6;
static const double BatchSeconds = 0.025;

typedef struct {
    const char *name;
    double density, temperature, column, dust, fuv, photoionizationRate;
    double ionizedHydrogen, molecularHydrogen, singlyIonized, doublyIonized;
    double uvbColumn;
} BenchmarkCase;

static StructEQCData InitialCell(const BenchmarkCase *cell, double years)
{
    const double nH = cell->density;
    StructEQCData state = {
        .T = cell->temperature, .dt = years*YEAR_CGS, .Td = 20.0,
        .Nt = cell->column, .Zd = cell->dust, .G0 = cell->fuv,
        .NLyn = cell->photoionizationRate/HIPhotoIonizationRate(1.0),
        .NHIuvb = {cell->uvbColumn,cell->uvbColumn},
    };
    state.n[ipHYDROGEN] = nH;
    state.n[ipHELIUM] = 0.10*nH;
    state.n[ipCARBON] = 3.0e-4*nH;
    state.n[ipOXYGEN] = 5.0e-4*nH;
    state.natom[ipHYDROGEN] = nH*(1.0-cell->ionizedHydrogen-cell->molecularHydrogen);
    state.nion[ipHYDROGEN] = nH*cell->ionizedHydrogen;
    state.nmol[ipH2] = 0.5*nH*cell->molecularHydrogen;
    for(int i = ipHELIUM; i <= ipOXYGEN; i++){
        state.natom[i] = state.n[i]*(1.0-cell->singlyIonized-cell->doublyIonized);
        state.nion[i] = state.n[i]*cell->singlyIonized;
    }
    state.nHeIII = state.n[ipHELIUM]*cell->doublyIonized;
    state.nCIII = state.n[ipCARBON]*cell->doublyIonized;
    state.nOIII = state.n[ipOXYGEN]*cell->doublyIonized;
    state.nel = EQCLibGetElectronNumberDensity(&state);
    state.Rho = EQCLibGetTotalMassDensity(&state);
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);
    return state;
}

static EQCLibStatus AdvanceEquilibrium(StructEQCData *state, double accuracy)
{
    const EQCLibStatus status =
        EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(state,accuracy,&state->Du);
    if(status == EQC_SUCCESS) state->U += state->dt*state->Du;
    return status;
}

typedef struct {
    const char *name;
    EQCLibStatus (*advance)(StructEQCData *, double);
    bool temperatureTables;
} Method;

static const Method Methods[] = {
    {"equilibrium",AdvanceEquilibrium,false},
    {"nonequilibrium",EQCLibTryAdvanceHydrogenThermochemistry,true},
    {"nonequilibrium_direct",EQCLibTryAdvanceHydrogenThermochemistry,false},
};
enum { MethodCount = sizeof Methods/sizeof Methods[0] };

static bool Measure(const char *scenario, const Method *method, const StructEQCData *initial,
                    clockid_t timer, unsigned int calls, double *elapsed, double *temperature)
{
    StructEQCData state;
    struct timespec start, end;
    EQCLibSetRunParamNonEquilibriumTemperatureTables(method->temperatureTables);
    if(clock_gettime(timer,&start) != 0){
        perror("clock_gettime");
        return false;
    }
    for(unsigned int i = 0; i < calls; i++){
        state = *initial;
        const EQCLibStatus status = method->advance(&state,Accuracy);
        if(status != EQC_SUCCESS){
            fprintf(stderr,"%s %s dt=%g yr: update %u failed (status %d)\n",
                    scenario,method->name,initial->dt/YEAR_CGS,i+1,status);
            return false;
        }
    }
    if(clock_gettime(timer,&end) != 0){
        perror("clock_gettime");
        return false;
    }
    if(!(isfinite(state.U) && state.U > 0 && isfinite(state.T) && state.T > 0 &&
         state.Rho == initial->Rho && state.n[ipHYDROGEN] == initial->n[ipHYDROGEN])){
        fprintf(stderr,"%s %s dt=%g yr: invalid resulting state\n",
                scenario,method->name,initial->dt/YEAR_CGS);
        return false;
    }
    *elapsed = (double)(end.tv_sec-start.tv_sec) + 1.0e-9*(double)(end.tv_nsec-start.tv_nsec);
    *temperature = state.T;
    return true;
}

static bool RunCase(FILE *file, const BenchmarkCase *cell, double years, clockid_t timer, int methodCount)
{
    const StructEQCData initial = InitialCell(cell,years);
    unsigned int calls[MethodCount];
    double elapsed, finalTemperature;
    /* Calibrate and warm up every method before collecting the five batches. */
    for(int method = 0; method < methodCount; method++){
        calls[method] = 1;
        for(;;){
            if(!Measure(cell->name,&Methods[method],&initial,timer,calls[method],&elapsed,&finalTemperature))
                return false;
            if(elapsed >= BatchSeconds || calls[method] == MaximumCalls) break;
            calls[method] *= 2;
        }
    }
    /* Rotate the first method to spread short-term system load. */
    for(int trial = 1; trial <= Trials; trial++){
        for(int i = 0; i < methodCount; i++){
            const int method = (trial-1+i)%methodCount;
            if(!Measure(cell->name,&Methods[method],&initial,timer,calls[method],&elapsed,&finalTemperature))
                return false;
            if(elapsed <= 0.0){
                fprintf(stderr,"%s %s dt=%g yr: batch measured zero elapsed time\n",cell->name,Methods[method].name,years);
                return false;
            }
            fprintf(file,"%s,%g,%g,%g,%s,%d,%u,%.12e,%.12e\n",cell->name,cell->density,cell->temperature,years,
                    Methods[method].name,trial,calls[method],elapsed,finalTemperature);
        }
    }
    return true;
}

int main(int argc, char **argv)
{
    if(argc != 2 && argc != 3){
        fprintf(stderr,"Usage: %s PEGASE2_DIRECTORY [UVB_DIRECTORY|--uvb-off]\n",argv[0]);
        return EXIT_FAILURE;
    }
    const bool compareUVB = argc == 3;
    const bool enableUVB = compareUVB && strcmp(argv[2],"--uvb-off") != 0;
    if(enableUVB && EQCLibSetRunParamUVBackground(true,argv[2]) != EQC_SUCCESS) return EXIT_FAILURE;
    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamChargeTransfer(true);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCOFormation(false,1.0);
    EQCLibSetRunParamCIEResonanceCooling(false);
    EQCLibSetRunParamCosmicRayIonization(false);
    EQCLibSetRunParamSoftXrayIonization(false);
    if(EQCLibTryInitialization(0,argv[1]) != EQC_SUCCESS) return EXIT_FAILURE;
    if(enableUVB && EQCLibTrySetUVBackgroundRedshift(0.0) != EQC_SUCCESS) return EXIT_FAILURE;

    const BenchmarkCase originalCases[] = {
        {.name = "ionizing", .density = 10.0, .temperature = 8000.0, .photoionizationRate = 1.0e-7},
        {.name = "recombining", .density = 10.0, .temperature = 1.35e4,
         .ionizedHydrogen = 0.98, .singlyIonized = 0.20, .doublyIonized = 0.70},
        {.name = "fuv", .density = 100.0, .temperature = 100.0, .column = 2.0e20,
         .dust = SOLAR_METALLICITY, .fuv = 1.0e3, .molecularHydrogen = 0.8},
    };
    const BenchmarkCase uvbCases[] = {
        {.name = "diffuse", .density = 0.04, .temperature = 9000.0, .uvbColumn = 1.0e16,
         .ionizedHydrogen = 0.1, .singlyIonized = 0.1},
        {.name = "shielded", .density = 1.0, .temperature = 8000.0, .uvbColumn = 1.0e20,
         .ionizedHydrogen = 0.01, .singlyIonized = 0.1},
    };
    const double originalSteps[] = {10.0,1.0e3,1.0e5};
    const double uvbSteps[] = {1.0e5};
    const BenchmarkCase *cases = compareUVB ? uvbCases : originalCases;
    const size_t caseCount = compareUVB ? sizeof uvbCases/sizeof uvbCases[0] : sizeof originalCases/sizeof originalCases[0];
    const double *steps = compareUVB ? uvbSteps : originalSteps;
    const size_t stepCount = compareUVB ? sizeof uvbSteps/sizeof uvbSteps[0] : sizeof originalSteps/sizeof originalSteps[0];
    const clockid_t timer = compareUVB ? CLOCK_PROCESS_CPUTIME_ID : CLOCK_MONOTONIC;
    const int methodCount = compareUVB ? 2 : MethodCount;
    FILE *file = fopen("benchmark.csv","w");
    if(file == NULL){
        perror("benchmark.csv");
        return EXIT_FAILURE;
    }
    fputs("scenario,nH,T_initial,dt_yr,method,trial,calls,elapsed_s,T_final\n",file);
    bool success = true;
    for(size_t i = 0; success && i < caseCount; i++){
        EQCLibSetRunParamH2Formation(cases[i].molecularHydrogen > 0.0,1.0);
        for(size_t j = 0; success && j < stepCount; j++)
            success = RunCase(file,&cases[i],steps[j],timer,methodCount);
    }
    if(fclose(file) != 0){
        perror("benchmark.csv");
        return EXIT_FAILURE;
    }
    if(!success) return EXIT_FAILURE;
    if(compareUVB) printf("Equilibrium and non-equilibrium CPU timings measured with UVB %s.\n",enableUVB ? "on" : "off");
    else puts("Equilibrium and both non-equilibrium timings measured.");
    return EXIT_SUCCESS;
}
