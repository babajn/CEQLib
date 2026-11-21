#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "HIIRegionFixture.h"

enum Illumination {
    NoFlux,
    Continuous,
    NewlyIlluminated,
    IlluminationCount
};

static const char *const illuminationName[IlluminationCount] = {
    "no_flux", "continuous", "newly_illuminated",
};
static const double densities[] = {1.0, 10.0, 100.0, 1000.0};
static const double metallicities[] = {0.0, 0.1, 1.0};
static const double logIonizations[] = {-4.0, -3.0, -2.0};
static const double temperatures[] = {7500.0, 10000.0, 15000.0};
static const double chemistryAccuracy = 1.0e-8;
static const double energyAccuracy = 1.0e-4;
enum {
    DensityCount = sizeof(densities)/sizeof(densities[0]),
    MetallicityCount = sizeof(metallicities)/sizeof(metallicities[0]),
    IonizationCount = sizeof(logIonizations)/sizeof(logIonizations[0]),
    TemperatureCount = sizeof(temperatures)/sizeof(temperatures[0]),
    SampleCount = DensityCount * MetallicityCount * IonizationCount * TemperatureCount,
    BlockCount = 7
};

typedef EQCLibStatus (*UpdateGas)(StructEQCData *);

static EQCLibStatus SolveChemistry(StructEQCData *state)
{
    return EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(state, chemistryAccuracy);
}

static EQCLibStatus AdvanceEnergy(StructEQCData *state)
{
    const EQCLibStatus status =
        EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(state, energyAccuracy, &state->Du);
    if(status == EQC_SUCCESS) state->U += state->dt * state->Du;
    return status;
}

static bool CheckStatus(EQCLibStatus status, const StructEQCData *state, const char *operation,
                        enum Illumination mode)
{
    if(status == EQC_SUCCESS) return true;
    fprintf(stderr, "%s / %s failed: status=%d n_H=%g Z/Z_sun=%g T=%g K NLyn=%g dt=%g yr\n",
            operation, illuminationName[mode], status, state->n[ipHYDROGEN],
            state->Z/ASPLUND_2009_SOLAR_METALLICITY, state->T, state->NLyn, state->dt/YEAR_CGS);
    return false;
}

static bool PrepareSamples(StructEQCData states[IlluminationCount][SampleCount])
{
    int index = 0;
    for(int d = 0; d < DensityCount; d++)
    for(int z = 0; z < MetallicityCount; z++)
    for(int u = 0; u < IonizationCount; u++)
    for(int t = 0; t < TemperatureCount; t++, index++){
        const HIIRegion region = {densities[d], metallicities[z], logIonizations[u]};
        StructEQCData steady = HIIRegionGasState(region, temperatures[t]);
        for(int i = 0; i < NumberChemicalElement; i++) steady.Rho += steady.n[i]*ChemicalElementMass[i];
        StructEQCData dark = steady;
        dark.NLyn = 0.0;
        if(!CheckStatus(SolveChemistry(&dark), &dark, "initial chemistry", NoFlux) ||
           !CheckStatus(SolveChemistry(&steady), &steady, "initial chemistry", Continuous))
            return false;
        dark.U = EQCLibConvertTtoU(&dark);
        steady.U = EQCLibConvertTtoU(&steady);

        states[NoFlux][index] = dark;
        states[Continuous][index] = steady;
        /* Switch on the same flux as the continuous case, keeping the preceding state. */
        states[NewlyIlluminated][index] = dark;
        states[NewlyIlluminated][index].NLyn = steady.NLyn;
    }
    return true;
}

static double ThreadCpuSeconds(void)
{
    struct timespec now;
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now) != 0){
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return now.tv_sec + 1.0e-9 * now.tv_nsec;
}

static bool TimeUpdates(FILE *csv, StructEQCData states[IlluminationCount][SampleCount],
                        const char *operation, UpdateGas update, double years, int batches)
{
    const int calls = batches * SampleCount;
    for(int block = 0; block < BlockCount; block++){
        /* Rotate execution order to spread short-term CPU timing variation across modes. */
        for(int k = 0; k < IlluminationCount; k++){
            const enum Illumination mode = (block + k) % IlluminationCount;
            const double start = ThreadCpuSeconds();
            for(int batch = 0; batch < batches; batch++){
                for(int i = 0; i < SampleCount; i++){
                    StructEQCData state = states[mode][i];
                    state.dt = years * YEAR_CGS;
                    if(!CheckStatus(update(&state), &state, operation, mode)) return false;
                }
            }
            const double microseconds = 1.0e6 * (ThreadCpuSeconds() - start) / calls;
            if(fprintf(csv, "%s,%s,%.0f,%d,%d,%d,%.8g\n", operation,
                       illuminationName[mode], years, block + 1, SampleCount, calls, microseconds) < 0)
                return false;
        }
    }
    return true;
}

int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr, "Usage: %s PEGASE2_DIRECTORY > hii_benchmark.csv\n", argv[0]);
        return EXIT_FAILURE;
    }
    if(InitializeHIIRegionFixture(argv[1]) != EQC_SUCCESS) return EXIT_FAILURE;
    StructEQCData states[IlluminationCount][SampleCount];
    if(!PrepareSamples(states)) return EXIT_FAILURE;
    if(puts("operation,illumination,dt_years,block,cases,calls,cpu_us_per_call") < 0)
        return EXIT_FAILURE;
    const bool success =
        TimeUpdates(stdout, states, "chemistry", SolveChemistry, 0.0, 200) &&
        TimeUpdates(stdout, states, "energy", AdvanceEnergy, 10.0, 60) &&
        TimeUpdates(stdout, states, "energy", AdvanceEnergy, 1.0e4, 40) &&
        TimeUpdates(stdout, states, "energy", AdvanceEnergy, 1.0e6, 16);
    return success && fflush(stdout) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
