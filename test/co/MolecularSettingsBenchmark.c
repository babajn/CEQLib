#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "MolecularFixture.h"

enum { DensityCount = 13, BlockCount = 5 };
static const double temperatures[] = {30., 100., 300., 1000., 3000., 8000., 15000.};
static const double logIonization[] = {-INFINITY, -5., -3.};
static const double years[] = {10., 1.e4, 1.e6};
static const double chemistryAccuracy = 1.e-5;
static const double energyAccuracy = 1.e-4;

typedef struct {
    const char *name;
    bool h2, co, chargeTransfer;
} Settings;

static const Settings modes[] = {
    {"h2_co_ct", true, true, true}, {"h2_co", true, true, false},
    {"h2_ct", true, false, true},   {"h2", true, false, false},
    {"atomic_ct", false, false, true}, {"atomic", false, false, false},
};

enum {
    TemperatureCount = sizeof(temperatures)/sizeof(temperatures[0]),
    FluxCount = sizeof(logIonization)/sizeof(logIonization[0]),
    ModeCount = sizeof(modes)/sizeof(modes[0]),
    CaseCount = DensityCount*TemperatureCount*FluxCount,
};

typedef struct {
    Condition input;
    StructEQCData gas;
    EQCLibStatus status;
} Sample;

static void SetMode(Settings mode)
{
    EQCLibSetRunParamH2Formation(mode.h2, 1.0);
    EQCLibSetRunParamCOFormation(mode.co, 1.0);
    EQCLibSetRunParamChargeTransfer(mode.chargeTransfer);
}

static void PrepareSamples(Sample samples[ModeCount][CaseCount])
{
    for(int mode=0; mode<ModeCount; mode++){
        SetMode(modes[mode]);
        int index=0;
        for(int u=0; u<FluxCount; u++)
        for(int t=0; t<TemperatureCount; t++)
        for(int d=0; d<DensityCount; d++, index++){
            const double density = pow(10., -2.+0.5*d);
            Sample *sample = &samples[mode][index];
            sample->input = (Condition){density, temperatures[t], 1.7, 1.e20*sqrt(density), logIonization[u]};
            sample->gas = GasState(sample->input);
            /* Equilibrate each model at the same T and convert with its own EOS. */
            sample->status = EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&sample->gas, 1.e-9);
            if(sample->status == EQC_SUCCESS) sample->gas.U = EQCLibConvertTtoU(&sample->gas);
        }
    }
}

static double ThreadCpuSeconds(void)
{
    struct timespec now;
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now) != 0){
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return now.tv_sec + 1.e-9*now.tv_nsec;
}

static void RunSample(const Sample *sample, int mode, int block, int index, double dt)
{
    StructEQCData gas = sample->gas;
    EQCLibStatus status = sample->status;
    double elapsed = NAN;
    double energyFraction = NAN;
    if(status == EQC_SUCCESS){
        gas.dt = dt*YEAR_CGS;
        double du;
        const double start = ThreadCpuSeconds();
        status = EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&gas, energyAccuracy, &du);
        elapsed = 1.e6*(ThreadCpuSeconds()-start);
        if(status == EQC_SUCCESS) energyFraction = gas.dt*du/sample->gas.U;
    }
    const bool ok = status == EQC_SUCCESS;
    printf("%s,%d,%d,%.0f,%g,%g,%g,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
           modes[mode].name, block+1, index, dt, sample->input.nH, sample->input.T, sample->input.logU,
           sample->status, status, elapsed, ok ? gas.T : NAN, energyFraction,
           ok ? 2.0*gas.nmol[ipH2]/gas.n[ipHYDROGEN] : NAN,
           ok ? gas.nmol[ipCO]/gas.n[ipCARBON] : NAN,
           ok ? gas.nion[ipHYDROGEN]/gas.n[ipHYDROGEN] : NAN,
           ok ? gas.nion[ipOXYGEN]/gas.n[ipOXYGEN] : NAN,
           ok ? gas.nOIII/gas.n[ipOXYGEN] : NAN);
}

int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr, "Usage: %s PEGASE2_DIRECTORY > molecular_settings.csv\n", argv[0]);
        return EXIT_FAILURE;
    }
    if(EQCLibSetRunParamChemistryAccuracy(chemistryAccuracy) != EQC_SUCCESS ||
       InitializeMolecularFixture(argv[1]) != EQC_SUCCESS)
        return EXIT_FAILURE;

    static Sample samples[ModeCount][CaseCount];
    PrepareSamples(samples);
    puts("mode,block,case,dt_years,nH,T0,log_U,seed_status,status,cpu_us,T,delta_U_over_U0,H2,CO,HII,OII,OIII");
    for(size_t t=0; t<sizeof(years)/sizeof(years[0]); t++)
    for(int block=0; block<BlockCount; block++)
    for(int i=0; i<ModeCount; i++){
        const int mode = (i+block)%ModeCount;
        SetMode(modes[mode]);
        for(int index=0; index<CaseCount; index++)
            RunSample(&samples[mode][index], mode, block, index, years[t]);
    }
    return fflush(stdout) == 0 && !ferror(stdout) ? EXIT_SUCCESS : EXIT_FAILURE;
}
