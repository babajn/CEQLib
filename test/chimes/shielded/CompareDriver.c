/* Shared driver for the EQCLib-CHIMES comparison.
 *
 * Usage: compare_<backend> MODE VARIANT DATA_DIRECTORY OUTPUT.csv
 *   MODE      equilibrium | timing | courant | cooling
 *   VARIANT   EQCLib: default | co      CHIMES: matched | full
 *   DATA_DIRECTORY  EQCLib: PEGASE2 directory, CHIMES: chimes-data directory
 */
#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "CompareBackend.h"

enum { Repeats = 5 };

static double CpuSeconds(void)
{
    struct timespec now;
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now) != 0){
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return (double)now.tv_sec + 1.0e-9*(double)now.tv_nsec;
}

static int CompareDouble(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void PrintSnapshot(FILE *f, const Snapshot *s)
{
    fprintf(f, ",%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g",
            s->T, s->xe, s->xHII, s->fH2, s->xHeII, s->xHeIII, s->xCII, s->xCIII, s->xCO,
            s->xOII, s->xOIII, s->xSiII, s->heat, s->cool);
}
static const char *SnapshotHeader =
    "T,xe,xHII,fH2,xHeII,xHeIII,xCII,xCIII,xCO,xOII,xOIII,xSiII,heat,cool";

static void Equilibrium(FILE *f)
{
    Cell *cell = BackendNewCell();
    fprintf(f, "code,nH,T0,status,iterations,cpu_us,%s\n", SnapshotHeader);
    for(int d = 0; d <= 24; d++){
        const double nH = pow(10.0, -2.0+0.25*d);
        for(int t = 0; t <= 55; t++){
            const double T = pow(10.0, 1.5+0.1*t);
            int iterations = 0;
            const double start = CpuSeconds();
            const int status = BackendEquilibrium(cell, nH, T, &iterations);
            const double elapsed = 1.0e6*(CpuSeconds()-start);
            Snapshot s;
            BackendSnapshot(cell, &s, 1);
            fprintf(f, "%s,%.9g,%.9g,%d,%d,%.6g", BackendName(), nH, T, status, iterations, elapsed);
            PrintSnapshot(f, &s);
            fputc('\n', f);
        }
        fflush(f);
        fprintf(stderr, "%s equilibrium: nH = %g done\n", BackendName(), nH);
    }
}

/* Courant step of an SPH particle of the given mass: dt = 0.3 h / c_s with
 * h = (m/rho)^(1/3), gamma = 5/3, mu = 1.27 below 1e4 K and 0.6 above, capped
 * at 1 Myr.  Both codes receive the same dt for the same (nH, T0, mass). */
static double CourantYears(double nH, double T, double massSolar)
{
    const double mH = 1.6726e-24, k = 1.380649e-16, solar = 1.989e33;
    const double rho = 1.4*mH*nH, mu = T < 1.0e4 ? 1.27 : 0.6;
    const double h = cbrt(massSolar*solar/rho), cs = sqrt(5.0/3.0*k*T/(mu*mH));
    return fmin(0.3*h/cs/CMP_YEAR, 1.0e6);
}

/* courant == 0: fixed steps of 10, 1e4 and 1e6 yr.
 * courant == 1: Courant steps for particle masses of 1e2 and 1e4 Msun. */
static void Timing(FILE *f, int courant)
{
    static const double temperatures[] = {30., 100., 300., 1.e3, 3.e3, 8.e3, 2.e4, 1.e5, 1.e6, 1.e7};
    static const double years[] = {10., 1.e4, 1.e6};
    static const double masses[] = {1.e2, 1.e4};
    enum { DensityCount = 13, TemperatureCount = sizeof temperatures/sizeof temperatures[0] };
    Cell *seed[DensityCount][TemperatureCount];
    int seedStatus[DensityCount][TemperatureCount];
    Cell *work = BackendNewCell();
    for(int d = 0; d < DensityCount; d++)
        for(int t = 0; t < TemperatureCount; t++){
            int iterations;
            seed[d][t] = BackendNewCell();
            seedStatus[d][t] = BackendEquilibrium(seed[d][t], pow(10.0, -2.0+0.5*d), temperatures[t], &iterations);
        }
    fprintf(stderr, "%s timing: seeds ready\n", BackendName());
    fprintf(f, "code,method,mass_msun,nH,T0,dt_yr,seed_status,status,repeats,cpu_us_median,cpu_us_min,cpu_us_max,%s\n",
            SnapshotHeader);
    const size_t outer = courant ? sizeof masses/sizeof masses[0] : sizeof years/sizeof years[0];
    for(size_t y = 0; y < outer; y++){
        for(int d = 0; d < DensityCount; d++)
            for(int t = 0; t < TemperatureCount; t++){
                const double nH = pow(10.0, -2.0+0.5*d), mass = courant ? masses[y] : 0.0;
                const double dtYears = courant ? CourantYears(nH, temperatures[t], mass) : years[y];
                for(int m = 0; m < BackendMethodCount(); m++){
                    double cpu[Repeats];
                    int status = seedStatus[d][t], done = 0;
                    Snapshot s = {.T = NAN};
                    if(status == 0){
                        for(int r = 0; r < Repeats; r++){
                            BackendCopyCell(work, seed[d][t]);
                            const double start = CpuSeconds();
                            status = BackendAdvance(work, dtYears*CMP_YEAR, m);
                            cpu[r] = 1.0e6*(CpuSeconds()-start);
                            if(status != 0) break;
                            done++;
                        }
                        if(status == 0) BackendSnapshot(work, &s, 0);
                    }
                    double median = NAN, lo = NAN, hi = NAN;
                    if(done == Repeats){
                        qsort(cpu, Repeats, sizeof cpu[0], CompareDouble);
                        median = cpu[Repeats/2]; lo = cpu[0]; hi = cpu[Repeats-1];
                    }
                    fprintf(f, "%s,%s,%g,%.9g,%.9g,%.9g,%d,%d,%d,%.6g,%.6g,%.6g", BackendName(), BackendMethodName(m),
                            mass, nH, temperatures[t], dtYears, seedStatus[d][t], status, done, median, lo, hi);
                    if(status == 0) PrintSnapshot(f, &s);
                    else fputs(",nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan", f);
                    fputc('\n', f);
                }
            }
        fflush(f);
        fprintf(stderr, "%s timing: %s %g done\n", BackendName(), courant ? "mass [Msun]" : "dt [yr]",
                courant ? masses[y] : years[y]);
    }
}

static void Cooling(FILE *f)
{
    static const double densities[] = {0.1, 1.0, 10.0};
    static const double temperatures[] = {1.0e7, 1.0e4};
    const int perDex = 20, steps = 9*perDex;   /* 1 yr to 1 Gyr */
    Cell *seed = BackendNewCell(), *cell = BackendNewCell();
    fprintf(f, "code,method,nH,T0,step,t_yr,status,cpu_us,%s\n", SnapshotHeader);
    for(size_t d = 0; d < sizeof densities/sizeof densities[0]; d++)
        for(size_t t = 0; t < sizeof temperatures/sizeof temperatures[0]; t++){
            int iterations;
            const int seedStatus = BackendEquilibrium(seed, densities[d], temperatures[t], &iterations);
            for(int m = 0; m < BackendMethodCount(); m++){
                Snapshot s;
                BackendCopyCell(cell, seed);
                BackendSnapshot(cell, &s, 1);
                fprintf(f, "%s,%s,%g,%g,0,0,%d,0", BackendName(), BackendMethodName(m), densities[d],
                        temperatures[t], seedStatus);
                PrintSnapshot(f, &s);
                fputc('\n', f);
                if(seedStatus != 0) continue;
                double time = 0.0;
                for(int k = 0; k <= steps; k++){
                    const double next = pow(10.0, (double)k/perDex);
                    const double start = CpuSeconds();
                    const int status = BackendAdvance(cell, (next-time)*CMP_YEAR, m);
                    const double elapsed = 1.0e6*(CpuSeconds()-start);
                    time = next;
                    BackendSnapshot(cell, &s, status == 0);
                    fprintf(f, "%s,%s,%g,%g,%d,%.9g,%d,%.6g", BackendName(), BackendMethodName(m), densities[d],
                            temperatures[t], k+1, time, status, elapsed);
                    PrintSnapshot(f, &s);
                    fputc('\n', f);
                    if(status != 0) break;
                }
                fflush(f);
            }
            fprintf(stderr, "%s cooling: nH = %g, T0 = %g done\n", BackendName(), densities[d], temperatures[t]);
        }
}

int main(int argc, char **argv)
{
    if(argc != 5){
        fprintf(stderr, "Usage: %s equilibrium|timing|courant|cooling VARIANT DATA_DIRECTORY OUTPUT.csv\n", argv[0]);
        return EXIT_FAILURE;
    }
    if(BackendInit(argv[2], argv[3]) != 0) return EXIT_FAILURE;
    FILE *f = fopen(argv[4], "w");
    if(f == NULL){
        perror(argv[4]);
        return EXIT_FAILURE;
    }
    if(strcmp(argv[1], "equilibrium") == 0) Equilibrium(f);
    else if(strcmp(argv[1], "timing") == 0) Timing(f, 0);
    else if(strcmp(argv[1], "courant") == 0) Timing(f, 1);
    else if(strcmp(argv[1], "cooling") == 0) Cooling(f);
    else{
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    return fclose(f) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
