#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"
#include "chimes_proto.h"

typedef struct {
    const char *name;
    double density, temperature, history_end;
    double h_ionized, h_molecular, heco_singly, heco_doubly;
    bool thermal, molecules;
} Case;

static const Case cases[] = {
    {"hydrogen", 10., 1000., 1.e5, 1., 0., 0., 0., false, false},
    {"recombining", 10., 13500., 3.e5, .98, 0., .20, .70, true, false},
    {"molecular", 30., 100., 1.e5, 1.e-6, .80, 0., 0., true, true},
    {"hot", 10., 1.e6, 2.e5, 1., 0., 0., 1., true, false},
};

static const double steps[] = {10., 1.e3, 1.e5};
static const double nominal_tolerance = 1.e-6;
static const double tight_tolerance = 1.e-8;
/* CHIMES d8f0fdb leaves CR rates uninitialized when cr_rate is exactly zero. */
static const double negligible_cr_rate = 1.e-50;
static const double chimes_temperature_floor = 10.1;
enum { EQCLib, CHIMES, SolverCount, Batches = 5, HistoryPoints = 160 };
static const char *const solver_names[] = {"EQCLib", "CHIMES"};

typedef struct {
    StructEQCData eqc;
    struct gasVariables chimes;
    ChimesFloat abundances[CHIMES_TOTSIZE];
} State;

static struct globalVariables chimes;

static void SetAbundance(State *state, int species, double value)
{
    const int index = chimes.speciesIndices[species];
    if (index >= 0) state->abundances[index] = value;
}

static double Abundance(const State *state, int species)
{
    const int index = chimes.speciesIndices[species];
    return index >= 0 ? state->abundances[index] : 0.;
}

static State InitialState(const Case *test)
{
    const double nH = test->density;
    const double He = test->thermal ? .1 : 0.;
    const double C = test->thermal ? 3.e-4 : 0.;
    const double O = test->thermal ? 5.e-4 : 0.;
    const double metal_mass = 12.*C + 16.*O;
    const double metal_fraction = metal_mass / (1. + 4.*He + metal_mass);
    const double neutral = 1. - test->heco_singly - test->heco_doubly;
    State state = {
        .eqc = {.T = test->temperature, .Z = metal_fraction},
        .chimes = {
            .element_abundances = {He, C, 0., O}, .nH_tot = nH,
            .temperature = test->temperature, .ThermEvolOn = test->thermal,
            .metallicity = metal_fraction / .0129, .TempFloor = chimes_temperature_floor,
            .cr_rate = negligible_cr_rate,
            .light_speed_reduction_factor = 1.,
        },
    };
    StructEQCData *eqc = &state.eqc;
    eqc->n[ipHYDROGEN] = nH;
    eqc->n[ipHELIUM] = He*nH;
    eqc->n[ipCARBON] = C*nH;
    eqc->n[ipOXYGEN] = O*nH;
    eqc->natom[ipHYDROGEN] = nH*(1. - test->h_ionized - test->h_molecular);
    eqc->nion[ipHYDROGEN] = nH*test->h_ionized;
    eqc->nmol[ipH2] = .5*nH*test->h_molecular;
    for (int i = ipHELIUM; i <= ipOXYGEN; i++) {
        eqc->natom[i] = eqc->n[i]*neutral;
        eqc->nion[i] = eqc->n[i]*test->heco_singly;
    }
    eqc->nHeIII = He*nH*test->heco_doubly;
    eqc->nCIII = C*nH*test->heco_doubly;
    eqc->nOIII = O*nH*test->heco_doubly;
    eqc->nel = EQCLibGetElectronNumberDensity(eqc);
    eqc->Rho = EQCLibGetTotalMassDensity(eqc);
    eqc->MeanMolecularWeight = EQCLibGetMeanMolecularWeight(eqc);
    eqc->DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(eqc);
    eqc->U = EQCLibConvertTtoU(eqc);

    SetAbundance(&state, sp_HI, eqc->natom[ipHYDROGEN]/nH);
    SetAbundance(&state, sp_HII, test->h_ionized);
    SetAbundance(&state, sp_H2, .5*test->h_molecular);
    SetAbundance(&state, sp_elec, eqc->nel/nH);
    const int stages[][3] = {{sp_HeI, sp_HeII, sp_HeIII},
                             {sp_CI, sp_CII, sp_CIII},
                             {sp_OI, sp_OII, sp_OIII}};
    const double abundances[] = {He, C, O};
    for (int i = 0; i < 3; i++) {
        SetAbundance(&state, stages[i][0], abundances[i]*neutral);
        SetAbundance(&state, stages[i][1], abundances[i]*test->heco_singly);
        SetAbundance(&state, stages[i][2], abundances[i]*test->heco_doubly);
    }
    return state;
}

static void SetTolerance(double tolerance)
{
    chimes.relativeTolerance = tolerance;
    chimes.absoluteTolerance = 1.e-3*tolerance;
    chimes.explicitTolerance = tolerance;
}

static void Advance(State *state, const Case *test, int solver, double years, double tolerance)
{
    int status;
    if (solver == EQCLib) {
        state->eqc.dt = years*YEAR_CGS;
        status = test->thermal
            ? EQCLibTryAdvanceCell(&state->eqc, tolerance, EQC_CHEMISTRY_NONEQUILIBRIUM)
            : EQCLibTryAdvanceHydrogenChemistry(&state->eqc, tolerance);
    } else {
        state->chimes.abundances = state->abundances;
        state->chimes.hydro_timestep = years*YEAR_CGS;
        status = chimes_network(&state->chimes, &chimes);
    }
    if (status != 0) {
        fprintf(stderr, "%s %s dt=%.9g yr tol=%g failed (status %d)\n",
                test->name, solver_names[solver], years, tolerance, status);
        exit(EXIT_FAILURE);
    }
}

static void WriteState(FILE *output, State *state, const Case *test,
                       int solver, double years, double tolerance)
{
    double t, h1, h2, el, o2, o3, co, cooling, heating;
    const double nH = test->density;
    if (solver == EQCLib) {
        StructEQCData *eqc = &state->eqc;
        if (CEQCalculateNonEquilibriumCoolingHeatingRates(eqc) != EQC_SUCCESS) {
            fprintf(stderr, "%s EQCLib cooling-rate evaluation failed\n", test->name);
            exit(EXIT_FAILURE);
        }
        t = eqc->T;
        h1 = eqc->nion[ipHYDROGEN]/nH;
        h2 = 2.*eqc->nmol[ipH2]/nH;
        el = eqc->nel/nH;
        o2 = eqc->nion[ipOXYGEN]/nH;
        o3 = eqc->nOIII/nH;
        co = eqc->nmol[ipCO]/nH;
        cooling = eqc->L_tot;
        heating = eqc->H_tot;
    } else {
        state->chimes.abundances = state->abundances;
        t = state->chimes.temperature;
        h1 = Abundance(state, sp_HII);
        h2 = 2.*Abundance(state, sp_H2);
        el = Abundance(state, sp_elec);
        o2 = Abundance(state, sp_OII);
        o3 = Abundance(state, sp_OIII);
        co = Abundance(state, sp_CO);
        cooling = chimes_cooling_rate_external(&state->chimes, &chimes, 1);
        heating = chimes_cooling_rate_external(&state->chimes, &chimes, 2);
    }
    if (!(t > chimes_temperature_floor && isfinite(t) && h1 >= 0. && h2 >= 0. && h1+h2 <= 1.000001 &&
          el >= 0. && isfinite(el) && isfinite(o2) && isfinite(o3) && isfinite(co) &&
          isfinite(cooling) && isfinite(heating) && fabs(heating) < 1.e-20)) {
        fprintf(stderr, "%s %s invalid state at %.9g yr: T=%g cooling=%g heating=%g\n",
                test->name, solver_names[solver], years, t, cooling, heating);
        exit(EXIT_FAILURE);
    }
    fprintf(output, "%s,%s,%.12g,%.0e,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
            test->name, solver_names[solver], years, tolerance, t, h1, h2, el, o2, o3, co, cooling, heating);
}

static double ThreadSeconds(void)
{
    struct timespec time;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time) != 0) {
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return time.tv_sec + 1.e-9*time.tv_nsec;
}

static double Measure(const State *initial, const Case *test, int solver, double years, unsigned calls)
{
    const double start = ThreadSeconds();
    for (unsigned i = 0; i < calls; i++) {
        State state = *initial;
        Advance(&state, test, solver, years, nominal_tolerance);
    }
    return ThreadSeconds() - start;
}

static void Benchmark(FILE *timings, const State *initial, const Case *test, double years)
{
    unsigned calls[SolverCount] = {1, 1};
    SetTolerance(nominal_tolerance);
    for (int solver = 0; solver < SolverCount; solver++) {
        while (Measure(initial, test, solver, years, calls[solver]) < .025 && calls[solver] < (1u << 18))
            calls[solver] *= 2;
    }
    for (int batch = 0; batch < Batches; batch++) {
        for (int i = 0; i < SolverCount; i++) {
            const int solver = (i + batch)%SolverCount;
            const double elapsed = Measure(initial, test, solver, years, calls[solver]);
            fprintf(timings, "%s,%s,%.12g,%d,%u,%.12g\n",
                    test->name, solver_names[solver], years, batch+1, calls[solver], elapsed);
        }
    }
}

static void RunCase(FILE *history, FILE *endpoints, FILE *timings, const Case *test)
{
    chimes.Tmol_K = test->molecules ? 1.e5 : 0.;
    EQCLibSetRunParamH2Formation(test->molecules, 1.);
    const State initial = InitialState(test);
    SetTolerance(tight_tolerance);
    for (int solver = 0; solver < SolverCount; solver++) {
        State state = initial;
        WriteState(history, &state, test, solver, 0., tight_tolerance);
        double previous = 0.;
        for (int i = 0; i < HistoryPoints; i++) {
            const double years = .1*pow(test->history_end/.1, (double)i/(HistoryPoints-1));
            Advance(&state, test, solver, years-previous, tight_tolerance);
            WriteState(history, &state, test, solver, years, tight_tolerance);
            previous = years;
        }
    }
    for (size_t i = 0; i < sizeof steps/sizeof steps[0]; i++) {
        const double tolerances[] = {nominal_tolerance, tight_tolerance};
        for (size_t j = 0; j < sizeof tolerances/sizeof tolerances[0]; j++) {
            SetTolerance(tolerances[j]);
            for (int solver = 0; solver < SolverCount; solver++) {
                State state = initial;
                Advance(&state, test, solver, steps[i], tolerances[j]);
                WriteState(endpoints, &state, test, solver, steps[i], tolerances[j]);
            }
        }
        Benchmark(timings, &initial, test, steps[i]);
    }
    printf("Completed %s (CHIMES allocated species: %d)\n", test->name, chimes.totalNumberOfSpecies);
}

static FILE *OpenOutput(const char *directory, const char *case_group, const char *filename, const char *header)
{
    char path[4096];
    if (snprintf(path, sizeof path, "%s/%s-%s.csv", directory, case_group, filename) >= (int)sizeof path) {
        fputs("Output path is too long\n", stderr);
        exit(EXIT_FAILURE);
    }
    FILE *file = fopen(path, "w");
    if (file == NULL) {
        perror(path);
        exit(EXIT_FAILURE);
    }
    fputs(header, file);
    return file;
}

int main(int argc, char **argv)
{
    if (argc != 5 || (strcmp(argv[4], "hydrogen") != 0 && strcmp(argv[4], "mixed") != 0)) {
        fprintf(stderr, "Usage: %s PEGASE_DIR CHIMES_DATA_FILE OUTPUT_DIR hydrogen|mixed\n", argv[0]);
        return EXIT_FAILURE;
    }
    const bool hydrogen = strcmp(argv[4], "hydrogen") == 0;
    chimes = (struct globalVariables){
        .redshift_dependent_UVB_index = -1, .StaticMolCooling = 1, .grain_temperature = 20.,
        .cmb_temperature = 2.725, .scale_metal_tolerances = 1, .chimes_debug = 1,
        .element_included = {!hydrogen, 0, !hydrogen, 0, 0, 0, 0, 0, 0},
    };
    if (snprintf(chimes.MainDataTablePath, sizeof chimes.MainDataTablePath, "%s", argv[2]) >=
        (int)sizeof chimes.MainDataTablePath) {
        fputs("CHIMES data path is too long\n", stderr);
        return EXIT_FAILURE;
    }
    strcpy(chimes.EqAbundanceTablePath, "None");
    init_chimes(&chimes);

    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamRecombinationCaseB(false);
    if (EQCLibTryInitialization(0, argv[1]) != EQC_SUCCESS) return EXIT_FAILURE;

    const char state_columns[] = "case,solver,time_yr,rtol,T_K,x_HII,f_H2,x_e,x_OII,x_OIII,x_CO,cooling,heating\n";
    FILE *history = OpenOutput(argv[3], argv[4], "history", state_columns);
    FILE *endpoints = OpenOutput(argv[3], argv[4], "steps", state_columns);
    FILE *timings = OpenOutput(argv[3], argv[4], "timings", "case,solver,dt_yr,batch,calls,cpu_s\n");
    for (size_t i = 0; i < sizeof cases/sizeof cases[0]; i++) {
        if (hydrogen != cases[i].thermal) RunCase(history, endpoints, timings, &cases[i]);
    }
    const int history_status = fclose(history);
    const int endpoint_status = fclose(endpoints);
    const int timing_status = fclose(timings);
    return history_status || endpoint_status || timing_status ? EXIT_FAILURE : EXIT_SUCCESS;
}
