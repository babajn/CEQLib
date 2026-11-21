#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"

enum { DensityCount = 13, TemperatureCount = 25, SampleCount = DensityCount*TemperatureCount, StepCount = 3 };
static const double Accuracy = 1.0e-4;

static StructEQCData GasState(double density, double temperature, double timestep)
{
    StructEQCData state = {
        .dt = timestep, .T = temperature, .Nt = 1.0e20*sqrt(density),
        .Z = ASPLUND_2009_SOLAR_METALLICITY, .Zd = SOLAR_METALLICITY, .G0 = 1.7, .X0 = 1.0,
        .nel = 1.0e-3*density,
        .n = {density, 0.1*density, 3.019952e-4*density, 4.570882e-4*density, 3.548134e-6*density},
    };
    for(int i = 0; i < NumberChemicalElement; i++){
        state.natom[i] = state.n[i];
        state.Rho += state.n[i]*ChemicalElementMass[i];
    }
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);
    return state;
}

static bool SameSolution(const StructEQCData *a, const StructEQCData *b)
{
    if(a->T != b->T || a->U != b->U || a->Du != b->Du || a->nel != b->nel ||
       a->H_tot != b->H_tot || a->L_tot != b->L_tot ||
       a->nHeIII != b->nHeIII || a->nCIII != b->nCIII ||
       a->nOIII != b->nOIII || a->nSiIII != b->nSiIII ||
       a->nmol[ipH2] != b->nmol[ipH2] || a->nmol[ipCO] != b->nmol[ipCO]) return false;
    for(int i = 0; i < NumberChemicalElement; i++){
        if(a->natom[i] != b->natom[i] || a->nion[i] != b->nion[i]) return false;
    }
    return true;
}

static EQCLibStatus AdvanceState(StructEQCData *state)
{
    const EQCLibStatus status =
        EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(state,Accuracy,&state->Du);
    if(status == EQC_SUCCESS) state->U += state->dt*state->Du;
    return status;
}

static bool RunBenchmark(const char *mode, double years)
{
    static StructEQCData serial[SampleCount], parallel[SampleCount];
    static EQCLibStatus status[SampleCount];
    for(int i = 0; i < DensityCount; i++){
        const double density = pow(10.0,-2.0+0.5*i);
        for(int j = 0; j < TemperatureCount; j++){
            const int sample = i*TemperatureCount+j;
            serial[sample] = GasState(density,pow(10.0,1.0+0.25*j),years*YEAR_CGS);
            parallel[sample] = serial[sample];
        }
    }

    bool ok = true;
    for(int step = 1; step <= StepCount; step++){
        double serialTime = 0.0, worstTime = 0.0, worstDensity = 0.0, worstTemperature = 0.0;
        int successes = 0;
        for(int i = 0; i < SampleCount; i++){
            const double density = serial[i].n[ipHYDROGEN];
            const double temperature = serial[i].T;
            const double start = omp_get_wtime();
            status[i] = AdvanceState(&serial[i]);
            const double elapsed = omp_get_wtime()-start;
            serialTime += elapsed;
            if(elapsed > worstTime){
                worstTime = elapsed;
                worstDensity = density;
                worstTemperature = temperature;
            }
            if(status[i] == EQC_SUCCESS) successes++;
            else fprintf(stderr,"%s dt=%g yr step=%d sample=%d status=%d nH=%.9g cm^-3 T=%.9g K\n",
                         mode,years,step,i,status[i],density,temperature);
        }

        int mismatches = 0, threads = 0;
        const double start = omp_get_wtime();
        #pragma omp parallel reduction(+:mismatches)
        {
            #pragma omp single
            threads = omp_get_num_threads();
            #pragma omp for schedule(dynamic)
            for(int i = 0; i < SampleCount; i++){
                const EQCLibStatus result = AdvanceState(&parallel[i]);
                if(result != status[i] || !SameSolution(&parallel[i],&serial[i])) mismatches++;
            }
        }
        const double parallelTime = omp_get_wtime()-start;
        printf("%s dt=%g yr step=%d: %d/%d converged, serial mean %.3f ms, worst %.3f ms (nH=%.3g, T=%.3g K)\n",
               mode,years,step,successes,SampleCount,1.0e3*serialTime/SampleCount,1.0e3*worstTime,
               worstDensity,worstTemperature);
        printf("OpenMP: %d threads, mean %.3f ms, %.3f s total, %d mismatches\n",
               threads,1.0e3*parallelTime/SampleCount,parallelTime,mismatches);
        if(successes != SampleCount || mismatches != 0) ok = false;
    }
    return ok;
}

int main(int argc, char **argv)
{
    if(argc != 2 || (strcmp(argv[1],"atomic") != 0 && strcmp(argv[1],"molecular") != 0)){
        fprintf(stderr,"Usage: %s atomic|molecular\n",argv[0]);
        return EXIT_FAILURE;
    }
    const bool molecules = strcmp(argv[1],"molecular") == 0;
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamH2Formation(molecules,1.0);
    EQCLibSetRunParamCOFormation(molecules,1.0);
    if(EQCLibTryInitialization(0,NULL) != EQC_SUCCESS) return EXIT_FAILURE;

    const double timesteps[] = {10.0,1.0e4,1.0e6};
    bool ok = true;
    for(size_t i = 0; i < sizeof(timesteps)/sizeof(timesteps[0]); i++){
        if(!RunBenchmark(argv[1],timesteps[i])) ok = false;
    }
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
