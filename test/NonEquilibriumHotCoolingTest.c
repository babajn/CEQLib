#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"
#include "PhotochemicalGasPhaseChemistry.h"

enum { HotSamples = 601, OffSamples = 121, Trials = 5, MaximumCalls = 1 << 18 };
static const double HistoryAccuracy = 1.e-8;
static const double BenchmarkAccuracy = 1.e-6;
static const double ColdBoundary = 1.5e4;
static const double HotBoundary = 2.5e4;
static const char HistoryHeader[] = "scenario,method,phase,time_yr,nH,T,U,H_PI,L_tot,L_OII,L_OIII,L_res,OII,OIII\n";

typedef struct {
    const char *name;
    double temperature, hydrogenIonized, singlyIonized, doublyIonized;
} Scenario;

static StructEQCData InitialCell(const Scenario *scenario)
{
    const double nH = 10.0;
    const double neutral = 1.0-scenario->singlyIonized-scenario->doublyIonized;
    StructEQCData state = {.T = scenario->temperature, .Z = SOLAR_METALLICITY};
    state.n[ipHYDROGEN] = nH;
    state.n[ipHELIUM] = 0.10*nH;
    state.n[ipCARBON] = 3.e-4*nH;
    state.n[ipOXYGEN] = 5.e-4*nH;
    state.natom[ipHYDROGEN] = nH*(1.0-scenario->hydrogenIonized);
    state.nion[ipHYDROGEN] = nH*scenario->hydrogenIonized;
    for(int i = ipHELIUM; i <= ipOXYGEN; i++){
        state.natom[i] = state.n[i]*neutral;
        state.nion[i] = state.n[i]*scenario->singlyIonized;
    }
    state.nHeIII = state.n[ipHELIUM]*scenario->doublyIonized;
    state.nCIII = state.n[ipCARBON]*scenario->doublyIonized;
    state.nOIII = state.n[ipOXYGEN]*scenario->doublyIonized;
    state.nel = EQCLibGetElectronNumberDensity(&state);
    state.Rho = EQCLibGetTotalMassDensity(&state);
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);
    state.invT = 1.0/state.T;
    state.logT = log10(state.T);
    return state;
}

static bool Advance(StructEQCData *state, double years, double accuracy, const char *scenario)
{
    state->dt = years*YEAR_CGS;
    const EQCLibStatus status = EQCLibTryAdvanceCell(state,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM);
    if(status == EQC_SUCCESS) return true;
    fprintf(stderr,"%s: update failed at T=%g K, dt=%g yr (status %d)\n",
            scenario,state->T,years,status);
    return false;
}

static bool WriteState(FILE *file, const char *scenario, const char *method, const char *phase,
                       double years, const StructEQCData *state)
{
    const double nO = state->n[ipOXYGEN];
    return fprintf(file,"%s,%s,%s,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                   scenario,method,phase,years,state->n[ipHYDROGEN],state->T,state->U,state->H_PI,
                   state->L_tot,state->L_OII,state->L_OIII,state->L_res,state->nion[ipOXYGEN]/nO,
                   state->nOIII/nO) >= 0;
}

static bool RunHotHistory(FILE *file, const Scenario *scenario, const char *method, bool legacy)
{
    const double firstYear = 1.e-2, lastYear = 2.e6;
    StructEQCData state = InitialCell(scenario);
    if(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) != EQC_SUCCESS) return false;
    double priorYear = 0.0;
    for(int i = -1; i < HotSamples; i++){
        double year = 0.0;
        if(i >= 0){
            year = firstYear*pow(lastYear/firstYear,(double)i/(HotSamples-1));
            if(!Advance(&state,year-priorYear,HistoryAccuracy,scenario->name)) return false;
        }
        if((legacy && state.L_res != 0.0) ||
           (!legacy && state.T >= HotBoundary &&
            (state.L_res <= 0.0 || state.L_OII != 0.0 || state.L_OIII != 0.0)) ||
           (state.T <= ColdBoundary && state.L_res != 0.0)){
            fprintf(stderr,"%s %s: unexpected oxygen cooling at T=%g K\n",scenario->name,method,state.T);
            return false;
        }
        if(!WriteState(file,scenario->name,method,"cooling",year,&state)) return false;
        priorYear = year;
    }
    if(state.T < ColdBoundary) return true;
    fprintf(stderr,"%s %s: history did not cool through the oxygen-line range\n",scenario->name,method);
    return false;
}

static bool RunSwitch(FILE *file, const char *method)
{
    const Scenario pulse = {"ionizing", 8000.0, 0.0, 0.0, 0.0};
    StructEQCData state = InitialCell(&pulse);
    state.NLyn = 1.e-7/HIPhotoIonizationRate(1.0);
    if(!Advance(&state,5.e8/YEAR_CGS,HistoryAccuracy,pulse.name)) return false;
    if(!(state.T < ColdBoundary && state.H_PI > 0.0 && state.L_OII > 0.0 && state.L_OIII > 0.0 &&
         state.L_res == 0.0)){
        fprintf(stderr,"%s: photoionized gas did not enter the warm oxygen-line regime\n",method);
        return false;
    }
    if(!WriteState(file,pulse.name,method,"on",0.0,&state)) return false;
    const StructEQCData illuminated = state;
    state.NLyn = 0.0;
    if(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) != EQC_SUCCESS ||
       !(state.H_PI == 0.0 && state.L_res == 0.0 &&
         state.L_OII == illuminated.L_OII && state.L_OIII == illuminated.L_OIII)){
        fprintf(stderr,"%s: oxygen lines changed instantaneously at light-off\n",method);
        return false;
    }
    if(!WriteState(file,pulse.name,method,"off",0.0,&state)) return false;

    double priorYear = 0.0;
    const double firstYear = 1.e-3, lastYear = 1.e6;
    for(int i = 0; i < OffSamples; i++){
        const double year = firstYear*pow(lastYear/firstYear,(double)i/(OffSamples-1));
        if(!Advance(&state,year-priorYear,HistoryAccuracy,pulse.name) ||
           !WriteState(file,pulse.name,method,"off",year,&state)) return false;
        priorYear = year;
    }
    if(state.T < illuminated.T && state.nion[ipOXYGEN] < illuminated.nion[ipOXYGEN] &&
       state.nOIII < illuminated.nOIII && state.L_OII < illuminated.L_OII && state.L_OIII < illuminated.L_OIII)
        return true;
    fprintf(stderr,"%s: the dark oxygen lines did not fade as oxygen recombined\n",method);
    return false;
}

static bool Measure(const StructEQCData *initial, unsigned calls, double *seconds, double *temperature)
{
    struct timespec start, end;
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&start) != 0) return false;
    StructEQCData state;
    for(unsigned int i = 0; i < calls; i++){
        state = *initial;
        const EQCLibStatus status = EQCLibTryAdvanceCell(&state,BenchmarkAccuracy,EQC_CHEMISTRY_NONEQUILIBRIUM);
        if(status != EQC_SUCCESS) return false;
    }
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&end) != 0) return false;
    *seconds = (double)(end.tv_sec-start.tv_sec)+1.e-9*(double)(end.tv_nsec-start.tv_nsec);
    *temperature = state.T;
    return *seconds > 0.0 && isfinite(*temperature);
}

static bool Benchmark(FILE *file, const Scenario *scenario, const char *method)
{
    const double steps[] = {10.0, 1.e3, 1.e5};
    for(size_t i = 0; i < sizeof steps/sizeof steps[0]; i++){
        StructEQCData initial = InitialCell(scenario);
        initial.dt = steps[i]*YEAR_CGS;
        double elapsed, finalTemperature;
        unsigned calls = 1;
        /* Calibrate each batch and warm the solver before measuring. */
        for(;;){
            if(!Measure(&initial,calls,&elapsed,&finalTemperature)) return false;
            if(elapsed >= 0.025 || calls == MaximumCalls) break;
            calls *= 2;
        }
        for(int trial = 1; trial <= Trials; trial++){
            if(!Measure(&initial,calls,&elapsed,&finalTemperature) ||
               fprintf(file,"%s,%s,%.0f,%d,%u,%.12e,%.12e\n",
                       scenario->name,method,steps[i],trial,calls,elapsed,finalTemperature) < 0) return false;
        }
    }
    return true;
}

int main(int argc, char **argv)
{
    if(argc != 3 || (strcmp(argv[2],"legacy") != 0 && strcmp(argv[2],"cie_linear") != 0 &&
                     strcmp(argv[2],"cie_bridge") != 0)){
        fprintf(stderr,"Usage: %s PEGASE2_DIRECTORY legacy|cie_linear|cie_bridge\n",argv[0]);
        return EXIT_FAILURE;
    }
    const char *method = argv[2];
    const bool legacy = strcmp(method,"legacy") == 0;
    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    if(EQCLibTryInitialization(0,argv[1]) != EQC_SUCCESS) return EXIT_FAILURE;

    const Scenario cases[] = {
        {"hot_1e6", 1.e6, 1.0, 0.0, 1.0},
        {"hot_1e7", 1.e7, 1.0, 0.0, 1.0},
        {"transition", 2.e4, 0.98, 0.2, 0.7},
    };
    FILE *hot = fopen("hot_cooling.csv","w");
    FILE *pulse = fopen("oxygen_after_switch.csv","w");
    FILE *benchmark = fopen("hot_benchmark.csv","w");
    bool success = hot && pulse && benchmark;
    if(!success) perror("Opening hot cooling output");
    if(success){
        success = fputs(HistoryHeader,hot) >= 0 && fputs(HistoryHeader,pulse) >= 0 &&
                  fputs("scenario,method,dt_yr,trial,calls,elapsed_s,T_final\n",benchmark) >= 0 &&
                  RunHotHistory(hot,&cases[0],method,legacy) && RunHotHistory(hot,&cases[1],method,legacy) &&
                  RunSwitch(pulse,method);
        for(size_t i = 0; success && i < sizeof cases/sizeof cases[0]; i++){
            success = Benchmark(benchmark,&cases[i],method);
            if(!success) fprintf(stderr,"%s: benchmark failed for %s\n",method,cases[i].name);
        }
    }
    if(hot && fclose(hot) != 0) success = false;
    if(pulse && fclose(pulse) != 0) success = false;
    if(benchmark && fclose(benchmark) != 0) success = false;
    if(!success) return EXIT_FAILURE;
    printf("%s: hot cooling, oxygen turn-off and update timing recorded.\n",method);
    return EXIT_SUCCESS;
}
