#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CoolingTableGnatFerland.h"
#include "CoolingCIE.h"

/* Tracked elements use the cell number densities. The last table holds
   N, Ne, Mg, S, Ca and Fe for the Asplund et al. (2009) solar mixture per H nucleus
   and is scaled by Z / ASPLUND_2009_SOLAR_METALLICITY. */
enum { TrackedElementCount = 3, CoolingTableCount = TrackedElementCount+1, UntrackedTable = TrackedElementCount };
static const int CoolingElement[TrackedElementCount] = {ipCARBON,ipOXYGEN,ipSILICON};
static const double *const TableLogLambda[CoolingTableCount] = {
    CoolingTableCIELogLambdaCarbon,
    CoolingTableCIELogLambdaOxygen,
    CoolingTableCIELogLambdaSilicon,
    CoolingTableCIELogLambdaUntrackedSolar
};

static bool CoolingTableInitialized;
static gsl_spline *LinearLogT[CoolingTableCount];
static gsl_spline *SteffenLogT[CoolingTableCount];

bool InitializeCoolingTableCIE(int Check)
{
    if(CoolingTableInitialized) return true;

    for(int i=0;i<CoolingTableCount;i++){
        LinearLogT[i] = gsl_spline_alloc(gsl_interp_linear,COOLING_TABLE_CIE_SIZE_T);
        SteffenLogT[i] = gsl_spline_alloc(gsl_interp_steffen,COOLING_TABLE_CIE_SIZE_T);
        if(!LinearLogT[i] || !SteffenLogT[i] ||
           gsl_spline_init(LinearLogT[i],CoolingTableCIELogT,TableLogLambda[i],COOLING_TABLE_CIE_SIZE_T) != 0 ||
           gsl_spline_init(SteffenLogT[i],CoolingTableCIELogT,TableLogLambda[i],COOLING_TABLE_CIE_SIZE_T) != 0){
            for(int j=0;j<=i;j++){
                if(LinearLogT[j]) gsl_spline_free(LinearLogT[j]);
                if(SteffenLogT[j]) gsl_spline_free(SteffenLogT[j]);
                LinearLogT[j] = SteffenLogT[j] = NULL;
            }
            return false;
        }
    }
    if(Check == 1){
        FILE *output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CIECooling.dat","w");
        if(output == NULL){
            perror(EQC_DIAGNOSTIC_DIRECTORY "/CIECooling.dat");
        }else{
            fputs("# log10(T/K) Lambda_C Lambda_O Lambda_Si [erg cm^3 s^-1]"
                  " Lambda_untracked_solar [erg cm^3 s^-1 per H nucleus, before the ramp]\n",output);
            for(int it=0;it<COOLING_TABLE_CIE_SIZE_T;it++){
                fprintf(output,"%.8f %.8e %.8e %.8e %.8e\n",CoolingTableCIELogT[it],
                        pow(10.0,TableLogLambda[0][it]),
                        pow(10.0,TableLogLambda[1][it]),
                        pow(10.0,TableLogLambda[2][it]),
                        pow(10.0,TableLogLambda[UntrackedTable][it]));
            }
            fclose(output);
        }
    }
    CoolingTableInitialized = true;
    return true;
}

#ifndef EQC_CIE_TRACKED_ONLY
/* The supplementary term ramps from 10^4 to 10^5 K in log temperature. */
static double UntrackedCoolingWeight(double logT)
{
    if(logT <= 4.0) return 0.0;
    if(logT >= 5.0) return 1.0;
    const double s = logT-4.0;
    return s*s*(3.0-2.0*s);
}
#endif

/* Gnat & Ferland (2012), Table 1: each coefficient multiplies n_e n_element. */
static double CoolingRateWithSpline(const StructEQCData *state, gsl_spline *const spline[])
{
    if(!EQCLibRunParameters.UseCIEResonanceCooling) return 0.0;
    if(!CoolingTableInitialized || !isfinite(state->T) || state->T <= 0.0 || state->T > 1.0e9 ||
       !isfinite(state->Z) || state->Z < 0.0) return NAN;
    if(state->Z == 0.0) return 0.0;

    const double inputLogT = log10(state->T);
    const double minLogT = CoolingTableCIELogT[0];
    const double maxLogT = CoolingTableCIELogT[COOLING_TABLE_CIE_SIZE_T-1];
    const double tableLogT = fmin(fmax(inputLogT,minLogT),maxLogT);
    /* Reuse the bracket on the shared grid without mutating shared state. */
    gsl_interp_accel accel = {0};
    double cooling = 0.0;
    for(int i=0;i<TrackedElementCount;i++){
        const double density = state->n[CoolingElement[i]];
        if(!isfinite(density) || density < 0.0) return NAN;
        if(density > 0.0){
            const double logLambda = gsl_spline_eval(spline[i],tableLogT,&accel);
            cooling += density*pow(10.0,logLambda);
        }
    }

#ifndef EQC_CIE_TRACKED_ONLY
    const double weight = UntrackedCoolingWeight(inputLogT);
    const double nH = state->n[ipHYDROGEN];
    if(!isfinite(nH) || nH < 0.0) return NAN;
    if(weight > 0.0 && nH > 0.0){
        const double logLambda = gsl_spline_eval(spline[UntrackedTable],tableLogT,&accel);
        cooling += weight*(state->Z/ASPLUND_2009_SOLAR_METALLICITY)*nH*pow(10.0,logLambda);
    }
#endif

    /* Match the existing low-temperature taper and high-temperature free-free extension. */
    if(inputLogT < minLogT) cooling *= exp(100.0*(inputLogT-minLogT));
    if(inputLogT > maxLogT) cooling *= pow(10.0,0.5*(inputLogT-maxLogT));
    return state->nel*cooling;
}

double CoolingRateCIE(StructEQCData *state)
{
    return CoolingRateWithSpline(state,LinearLogT);
}

double CEQNonEquilibriumCoolingRateCIE(const StructEQCData *state)
{
    return CoolingRateWithSpline(state,SteffenLogT);
}
