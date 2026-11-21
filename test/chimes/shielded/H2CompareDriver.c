#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "H2Diagnostics.h"

static const double Densities[] = {3.0, 10.0, 30.0, 100.0};
static const double InitialTemperature = 5000.0;
static const char RateColumns[] =
    "fH2,xe,xHII,form_grain,form_other,destroy_photo,destroy_coll,destroy_other,"
    "k_form,k_photo0,dust_shield,self_shield,shield_other,"
    "cool_h2,cool_h2_diss,heat_h2,heat,cool,energy,net_oh,net_chplus";

static void WriteRates(FILE *out, const H2Rates *r)
{
    const double values[] = {
        r->fH2,r->xe,r->xHII,r->form_grain,r->form_other,
        r->destroy_photo,r->destroy_coll,r->destroy_other,
        r->k_form,r->k_photo0,r->dust_shield,r->self_shield,r->shield_other,
        r->cool_h2,r->cool_h2_diss,r->heat_h2,r->heat,r->cool,r->energy,
        r->net_oh,r->net_chplus,
    };
    for(size_t i=0;i<sizeof values/sizeof values[0];i++)
        fprintf(out,",%.12g",values[i]);
    fputc('\n',out);
}

static void WriteFixedTemperature(FILE *out, Cell *cell, double nH, double T, double fH2, int status)
{
    H2Rates rates = {0};
    if(status == 0) status = BackendH2Rates(cell,nH,T,fH2,&rates);
    fprintf(out,"%s,%s,%.9g,%.9g,%d,%.9g",BackendName(),fH2 < 0.0 ? "native" : "shared",nH,T,status,T);
    WriteRates(out,&rates);
}

static void FixedTemperature(FILE *out, int useSharedGrid)
{
    Cell *cell = BackendNewCell();
    static const double extraTemperatures[] = {30.0,300.0,3000.0,5000.0};
    static const double molecularFractions[] = {
        0.4,0.5,0.6,0.7,0.8,0.9,0.93,0.95,0.97,0.99,0.995,0.999,
        0.9999,0.99999,0.999999,0.9999999,1.0,
    };
    fprintf(out,"code,sample,nH,T0,status,T,%s\n",RateColumns);
    for(size_t i=0;i<sizeof Densities/sizeof Densities[0];i++){
        const double nH = Densities[i];
        const int extraCount = sizeof extraTemperatures/sizeof extraTemperatures[0];
        for(int point=-extraCount;point<=100;point++){
            const double T = point < 0 ? extraTemperatures[point+extraCount] : pow(10.0,1.5+0.025*point);
            int iterations;
            const int status = BackendEquilibrium(cell,nH,T,&iterations);
            WriteFixedTemperature(out,cell,nH,T,-1.0,status);
            if(!useSharedGrid) continue;

            WriteFixedTemperature(out,cell,nH,T,0.0,0);
            for(int molecular=0;molecular<=92;molecular++)
                WriteFixedTemperature(out,cell,nH,T,pow(10.0,-12.0+0.125*molecular),0);
            for(size_t j=0;j<sizeof molecularFractions/sizeof molecularFractions[0];j++)
                WriteFixedTemperature(out,cell,nH,T,molecularFractions[j],0);
        }
        fprintf(stderr,"%s fixed-temperature grid: nH = %g done\n",BackendName(),nH);
    }
}

static void CoolingHistory(FILE *out)
{
    Cell *cell = BackendNewCell();
    Cell *initial = BackendNewCell();
    fprintf(out,"code,method,nH,T0,step,t_yr,status,T,%s\n",RateColumns);
    for(size_t i=0;i<sizeof Densities/sizeof Densities[0];i++){
        const double nH = Densities[i];
        int iterations;
        const int initialStatus = BackendEquilibrium(initial,nH,InitialTemperature,&iterations);
        for(int method=0;method<BackendMethodCount();method++){
            BackendCopyCell(cell,initial);
            double previousYear = 0.0;
            for(int step=0;step<=181;step++){
                const double year = step == 0 ? 0.0 : pow(10.0,(step-1)*0.05);
                int status = step == 0 ? initialStatus : BackendAdvance(cell,(year-previousYear)*CMP_YEAR,method);
                Snapshot snapshot = {0};
                H2Rates rates = {0};
                if(status == 0){
                    BackendSnapshot(cell,&snapshot,0);
                    status = BackendH2Rates(cell,nH,snapshot.T,-1.0,&rates);
                }
                fprintf(out,"%s,%s,%.9g,%.9g,%d,%.9g,%d,%.9g",
                        BackendName(),BackendMethodName(method),nH,InitialTemperature,
                        step,year,status,snapshot.T);
                WriteRates(out,&rates);
                if(status != 0) break;
                previousYear = year;
            }
        }
        fprintf(stderr,"%s cooling histories: nH = %g done\n",BackendName(),nH);
    }
}

static double CpuSeconds(void)
{
    struct timespec now;
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&now) != 0){
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

/* Use the same Courant step prescription as the shielded comparison. */
static double CourantYears(double nH, double T, double mass)
{
    const double mH = 1.6726e-24, k = 1.380649e-16, solar = 1.989e33;
    const double h = cbrt(mass*solar/(1.4*mH*nH));
    const double cs = sqrt((5.0/3.0)*k*T/(1.27*mH));
    return fmin(0.3*h/cs/CMP_YEAR,1.0e6);
}

static void CourantTiming(FILE *out)
{
    static const double temperatures[] = {30.0,100.0,300.0,1000.0,3000.0,5000.0,8000.0};
    static const double masses[] = {100.0,1.0e4};
    enum { Repeats = 7 };
    Cell *seed = BackendNewCell(), *work = BackendNewCell();
    fprintf(out,"code,method,mass_msun,nH,T0,dt_yr,seed_status,status,repeats,"
                "cpu_us_median,cpu_us_min,cpu_us_max\n");
    for(size_t d=0;d<sizeof Densities/sizeof Densities[0];d++){
        const double nH = Densities[d];
        for(size_t t=0;t<sizeof temperatures/sizeof temperatures[0];t++){
            const double T = temperatures[t];
            int iterations;
            const int seedStatus = BackendEquilibrium(seed,nH,T,&iterations);
            for(size_t m=0;m<sizeof masses/sizeof masses[0];m++){
                const double dtYears = CourantYears(nH,T,masses[m]);
                int status = seedStatus, completed = 0;
                double elapsed[Repeats];
                if(status == 0){
                    BackendCopyCell(work,seed);
                    status = BackendAdvance(work,dtYears*CMP_YEAR,0);
                }
                for(int repeat=0;status==0 && repeat<Repeats;repeat++){
                    BackendCopyCell(work,seed);
                    const double start = CpuSeconds();
                    status = BackendAdvance(work,dtYears*CMP_YEAR,0);
                    elapsed[repeat] = 1.0e6*(CpuSeconds()-start);
                    if(status == 0) completed++;
                }
                double median = NAN, minimum = NAN, maximum = NAN;
                if(completed == Repeats){
                    qsort(elapsed,Repeats,sizeof elapsed[0],CompareDouble);
                    median = elapsed[Repeats/2];
                    minimum = elapsed[0];
                    maximum = elapsed[Repeats-1];
                }
                fprintf(out,"%s,%s,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%.9g,%.9g,%.9g\n",
                        BackendName(),BackendMethodName(0),masses[m],nH,T,dtYears,seedStatus,
                        status,completed,median,minimum,maximum);
            }
        }
    }
}

int main(int argc, char **argv)
{
    if(argc != 5){
        fprintf(stderr,"Usage: %s equilibrium|native|cooling|courant VARIANT DATA_DIRECTORY OUTPUT.csv\n",argv[0]);
        return 2;
    }
    if(strcmp(argv[1],"equilibrium") != 0 && strcmp(argv[1],"native") != 0 &&
       strcmp(argv[1],"cooling") != 0 && strcmp(argv[1],"courant") != 0){
        fprintf(stderr,"Unknown calculation %s\n",argv[1]);
        return 2;
    }
    if(BackendInit(argv[2],argv[3]) != 0) return 1;
    FILE *out = fopen(argv[4],"w");
    if(out == NULL){
        perror(argv[4]);
        return 1;
    }
    if(strcmp(argv[1],"equilibrium") == 0){
        FixedTemperature(out,strcmp(argv[2],"form025") != 0);
    }else if(strcmp(argv[1],"native") == 0){
        FixedTemperature(out,0);
    }else if(strcmp(argv[1],"cooling") == 0){
        CoolingHistory(out);
    }else{
        CourantTiming(out);
    }
    if(fclose(out) != 0){
        perror(argv[4]);
        return 1;
    }
    return 0;
}
