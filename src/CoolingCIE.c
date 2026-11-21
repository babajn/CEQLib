#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "Astro.h"
#include "Constants.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "CoolingTableSutherlandDopita.h"
#include "CoolingCIE.h"

#define Zsol (SOLAR_METALLICITY)

static int FlagCoolingTable = OFF;
static gsl_spline       *SplineAccessLogT[COOLING_TABLE_CIE_SIZE_Z];

void InitializeCoolingTableCIE(int Check){

    //InitializeCoolingRateLookUpTableCIE(Check);
    if(FlagCoolingTable == ON)  return;

    FILE *output = NULL;
    if(Check == 1){
        output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CIECooling.dat","w");
        if(output == NULL) perror(EQC_DIAGNOSTIC_DIRECTORY "/CIECooling.dat");
    }
    for(int im=1;im<COOLING_TABLE_CIE_SIZE_Z;im++){
        for(int it=0; it<COOLING_TABLE_CIE_SIZE_T; it++){
            double LambZ = pow(10.0,CoolingTableCIELogLambda[im][it]);
            double Lamb0 = pow(10.0,CoolingTableCIELogLambda[0][it]);
            CoolingTableCIELogLambda[im][it] = log10(MAX(1.e-100,LambZ-Lamb0));
            if(output != NULL){
                //fprintf(stderr,"logZ=%f logT=%f LambZ=%e Lamb0=%e LambZ-Lamb0=%e\n",
                //        CoolingTableCIELogZ[im],CoolingTableCIELogT[it],
                //        LambZ,Lamb0,
                //        pow(10.0,CoolingTableCIELogLambda[im][it])
                //        );
                fprintf(output,"%f %f %e %e %e\n",
                        CoolingTableCIELogZ[im],CoolingTableCIELogT[it],
                        LambZ,Lamb0,pow(10.0,CoolingTableCIELogLambda[im][it])
                        );
            }
        }
        //SplineAccessLogT[im] = gsl_spline_alloc(gsl_interp_cspline,COOLING_TABLE_CIE_SIZE_T);
        SplineAccessLogT[im] = gsl_spline_alloc(gsl_interp_linear,COOLING_TABLE_CIE_SIZE_T);
        gsl_spline_init(SplineAccessLogT[im],CoolingTableCIELogT,
                CoolingTableCIELogLambda[im],COOLING_TABLE_CIE_SIZE_T);
    }
    if(output != NULL) fclose(output);
    //exit(1);

    //if(Check==1){
    //    double metal[] = {1.e-3*Zsol,3.e-3*Zsol,1.e-2*Zsol,3.e-2*Zsol,
    //                      1.e-1*Zsol,3.e-1*Zsol,Zsol,3.e0*Zsol,1.e+1*Zsol};
    //    char fname[MaxCharactersInLine];
    //    for(int i=0;i<9;i++){
    //        FILE *fp;
    //        sprintf(fname,"out/CoolingCIE.Z%02d",i);
    //        FileOpen(fp,fname,"w");
    //        fprintf(stderr,"im=%d: [Fe/H]=%f\n",i,log10(metal[i]/Zsol));
    //        for(double T=5.e3;T<3.e+8;T*=1.1){
    //            fprintf(fp,"%e %e\n",
    //                    log10(T),log10(CoolingRateCIE(log10(T),metal[i])));
    //        }
    //        fclose(fp);
    //    }
    //}

    FlagCoolingTable = ON;
    return;
}


//// normalized cooling rate for CIE gas [erg cm^3/s] 
// ref. Sutherland & Dopita, ApJS, (1993)
//      http://www.mso.anu.edu.au/~ralph/data/cool/
double CoolingRateCIE(StructEQCData *PhydroPP){
    
    if(EQCLibRunParameters.UseCIEResonanceCooling==false)   return 0.e0;
    if(FlagCoolingTable == OFF || !isfinite(PhydroPP->T) || PhydroPP->T <= 0.0 || PhydroPP->T > 1.0e9
       || !isfinite(PhydroPP->Z) || PhydroPP->Z < 0.0) return NAN;
    if(PhydroPP->Z == 0.0) return 0.0;

    const double inputLogT = log10(PhydroPP->T);
    const double logT = fmin(inputLogT,CoolingTableCIELogT[COOLING_TABLE_CIE_SIZE_T-1]);
    double Z    = PhydroPP->Z; 
    double ne   = PhydroPP->nel;
    double nH   = PhydroPP->n[ipHYDROGEN];
    //double reduce = MAX(0.0,1.e0-exp(-20.0*(logT-3.9)));
    //double reduce = MAX(0.0,1.e0-exp(-20.0*(logT-3.95)));
    double FeH = log10(Z/Zsol);
    //fprintf(stderr,"FeH=%f\n",FeH);
    //exit(1);
    double val = 1.e-100;
    if(FeH <= CoolingTableCIELogZ[1]){
        /* Scale the lowest-metallicity contribution linearly toward zero. */
        const double tableLogT = fmax(logT,CoolingTableCIELogT[0]);
        val = gsl_spline_eval(SplineAccessLogT[1],tableLogT,NULL) + (FeH-CoolingTableCIELogZ[1]);
    } else if(FeH >= CoolingTableCIELogZ[COOLING_TABLE_CIE_SIZE_Z-1]){
        if(logT<CoolingTableCIELogT[0])
            val = CoolingTableCIELogLambda[COOLING_TABLE_CIE_SIZE_Z-1][0];
        else{
            val = gsl_spline_eval(SplineAccessLogT[COOLING_TABLE_CIE_SIZE_Z-1],
                    logT,NULL);
        }
    } else {
        int i=0;
        double val1 = 0.0, val2 = 0.0;
        while(CoolingTableCIELogZ[i]<FeH)   i++;
        if(logT<CoolingTableCIELogT[0]){
            val1 = CoolingTableCIELogLambda[i-1][0];
            val2 = CoolingTableCIELogLambda[i  ][0];
        }
        else{
            val1 = gsl_spline_eval(SplineAccessLogT[i-1],logT,NULL);
            val2 = gsl_spline_eval(SplineAccessLogT[i  ],logT,NULL);
        }
        val = ((CoolingTableCIELogZ[i]-FeH)*val1+(FeH-CoolingTableCIELogZ[i-1])*val2)
                /(CoolingTableCIELogZ[i]-CoolingTableCIELogZ[i-1]);
    }
    double reduce = MIN(1.0,exp(100*(logT-4.0)));
    //double reduce = 1.e0;
    //if(logT<4.0)    reduce = MIN(1.0,exp(100*(logT-4.01)));
    /* Continue the metal contribution above the table with the free-free T^0.5 scaling. */
    return pow(10.0,val + 0.5*(inputLogT-logT))*reduce*ne*nH;
}
