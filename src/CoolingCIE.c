#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "Astro.h"
#include "Constants.h"
#include "ChemicalEquilibriumCooling.h"
#include "CoolingTableSutherlandDopita.h"
#include "CoolingCIE.h"

static bool CoolingTableInitialized;
static gsl_spline *LinearLogT[COOLING_TABLE_CIE_SIZE_Z];
static gsl_spline *SteffenLogT[COOLING_TABLE_CIE_SIZE_Z];

bool InitializeCoolingTableCIE(int Check)
{
    if(CoolingTableInitialized) return true;

    FILE *output = NULL;
    if(Check == 1){
        output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CIECooling.dat","w");
        if(output == NULL) perror(EQC_DIAGNOSTIC_DIRECTORY "/CIECooling.dat");
    }
    bool ready = true;
    for(int im=1;im<COOLING_TABLE_CIE_SIZE_Z;im++){
        double metalLogLambda[COOLING_TABLE_CIE_SIZE_T];
        for(int it=0;it<COOLING_TABLE_CIE_SIZE_T;it++){
            const double LambZ = pow(10.0,CoolingTableCIELogLambda[im][it]);
            const double Lamb0 = pow(10.0,CoolingTableCIELogLambda[0][it]);
            metalLogLambda[it] = log10(MAX(1.e-100,LambZ-Lamb0));
            if(output != NULL){
                fprintf(output,"%f %f %e %e %e\n",
                        CoolingTableCIELogZ[im],CoolingTableCIELogT[it],
                        LambZ,Lamb0,pow(10.0,metalLogLambda[it]));
            }
        }
        LinearLogT[im] = gsl_spline_alloc(gsl_interp_linear,COOLING_TABLE_CIE_SIZE_T);
        SteffenLogT[im] = gsl_spline_alloc(gsl_interp_steffen,COOLING_TABLE_CIE_SIZE_T);
        if(!LinearLogT[im] || !SteffenLogT[im] ||
           gsl_spline_init(LinearLogT[im],CoolingTableCIELogT,metalLogLambda,COOLING_TABLE_CIE_SIZE_T) != 0 ||
           gsl_spline_init(SteffenLogT[im],CoolingTableCIELogT,metalLogLambda,COOLING_TABLE_CIE_SIZE_T) != 0){
            ready = false;
            break;
        }
    }
    if(output != NULL) fclose(output);
    if(!ready){
        for(int im=1;im<COOLING_TABLE_CIE_SIZE_Z;im++){
            if(LinearLogT[im]) gsl_spline_free(LinearLogT[im]);
            if(SteffenLogT[im]) gsl_spline_free(SteffenLogT[im]);
            LinearLogT[im] = SteffenLogT[im] = NULL;
        }
        return false;
    }
    CoolingTableInitialized = true;
    return true;
}

/* Sutherland & Dopita (1993), with the primordial contribution subtracted. */
static double CoolingRateWithSpline(const StructEQCData *state, gsl_spline *const spline[])
{
    if(!EQCLibRunParameters.UseCIEResonanceCooling) return 0.0;
    if(!CoolingTableInitialized || !isfinite(state->T) || state->T <= 0.0 || state->T > 1.0e9 ||
       !isfinite(state->Z) || state->Z < 0.0) return NAN;
    if(state->Z == 0.0) return 0.0;

    const int lastMetallicity = COOLING_TABLE_CIE_SIZE_Z-1;
    const double inputLogT = log10(state->T);
    const double logT = fmin(inputLogT,CoolingTableCIELogT[COOLING_TABLE_CIE_SIZE_T-1]);
    const double tableLogT = fmax(logT,CoolingTableCIELogT[0]);
    const double logZ = log10(state->Z/SOLAR_METALLICITY);
    double val;
    if(logZ <= CoolingTableCIELogZ[1]){
        /* Scale the lowest-metallicity contribution linearly toward zero. */
        val = gsl_spline_eval(spline[1],tableLogT,NULL) + (logZ-CoolingTableCIELogZ[1]);
    }else if(logZ >= CoolingTableCIELogZ[lastMetallicity]){
        val = gsl_spline_eval(spline[lastMetallicity],tableLogT,NULL);
    }else{
        int i=1;
        while(CoolingTableCIELogZ[i]<logZ) i++;
        const double val1 = gsl_spline_eval(spline[i-1],tableLogT,NULL);
        const double val2 = gsl_spline_eval(spline[i],tableLogT,NULL);
        val = ((CoolingTableCIELogZ[i]-logZ)*val1+(logZ-CoolingTableCIELogZ[i-1])*val2)
                /(CoolingTableCIELogZ[i]-CoolingTableCIELogZ[i-1]);
    }
    const double reduce = MIN(1.0,exp(100*(logT-4.0)));
    /* Continue the metal contribution above the table with the free-free T^0.5 scaling. */
    return pow(10.0,val + 0.5*(inputLogT-logT))*reduce*state->nel*state->n[ipHYDROGEN];
}

double CoolingRateCIE(StructEQCData *state)
{
    return CoolingRateWithSpline(state,LinearLogT);
}

double CEQNonEquilibriumCoolingRateCIE(const StructEQCData *state)
{
    return CoolingRateWithSpline(state,SteffenLogT);
}
