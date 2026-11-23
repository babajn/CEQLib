#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalElements.h"
#include "ChemicalEquilibriumCooling.h"
#include "Xray.h"

static double SoftXrayIonizationRate(double xe, double Ntot, CEQSoftXrayColumn *column);
static double SoftXrayHeatingRateBody(double xe, double Ntot);

void InitializeSoftXrayIonization(int Check){
    
    if(Check==1){
        FILE *output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/SoftXrayIonization.dat","w");
        for(double logNt=17.0; logNt<=22.0; logNt+=0.01){
            double Ntot  = pow(10.0,logNt);
            fprintf(output,"%e ",Ntot);
            for(double logxe = -4.0; logxe<=0.0; logxe+=1.0){
                double xe    = pow(10.0,logxe);
                double xi_XR = SoftXrayIonizationRate(xe,Ntot,NULL);
                double Gm_XR = SoftXrayHeatingRateBody(xe,Ntot);
                fprintf(output,"%e ",xi_XR);
                fprintf(output,"%e ",Gm_XR);
            }
            fprintf(output,"\n");
        }
        fclose(output);
    }
    //exit(1);

    return ;
}

/// X-ray ionization per H [1/s/H]
static double SoftXrayIonizationRate(double xe, double Ntot, CEQSoftXrayColumn *column){
    // ref. (A6)(A7)(A8) of Wolfire et al. (1995)
    if(xe<=0.0) return 0.0;
#ifdef EQC_LEGACY_XRAY_CUTOFF
    if(xe > EQC_XRAY_IONIZATION_MAX_XE) return 0.0;
#else
    if(xe >= EQC_XRAY_IONIZATION_MAX_XE) return 0.0;
#endif
    CEQSoftXrayColumn local = {0};
    if(column == NULL) column = &local;
    if(!column->ready){
        const double p1 = MAX(-0.3,log10(1.e-18*Ntot));
        column->base = -15.6-1.10*p1+9.13e-2*SQ(p1);
        column->attenuation = exp(-SQ((p1-0.41)/0.84));
        column->ready = true;
    }
    double p2 = log10(xe);
    double f4 = 1.06 + 4.08e-2*p2 + 6.51e-3*SQ(p2);
    double f5 = 1.90 + 0.678*p2 + 0.113*SQ(p2);
    double rate = f4*column->base + f5*0.87*column->attenuation;
    double ionization = pow(10.0,rate);
#ifndef EQC_LEGACY_XRAY_CUTOFF
    /* Appendix A fits the rate through xe=0.1. Join the fit smoothly to zero
       at the previous cutoff, using the remaining fraction to avoid
       cancellation close to the cutoff. */
    const double taperStart = 0.1;
    if(xe > taperStart){
        const double remaining = (EQC_XRAY_IONIZATION_MAX_XE - xe) /
                                 (EQC_XRAY_IONIZATION_MAX_XE - taperStart);
        ionization *= SQ(remaining) * (3.0 - 2.0*remaining);
    }
#endif
    return ionization;
}

double CEQSoftXrayIonizationRate(double xe, double Ntot, CEQSoftXrayColumn *column){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return 0.5*SoftXrayIonizationRate(xe,Ntot,column);
}

double SoftXrayIonizationRateHI(double xe, double Ntot){
    return CEQSoftXrayIonizationRate(xe,Ntot,NULL);
}

double SoftXrayIonizationRateHeI(double xe, double Ntot){
    return CEQSoftXrayIonizationRate(xe,Ntot,NULL);
}

double SoftXrayIonizationRateHeII(double xe, double Ntot){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return 0.e0;
}

double SoftXrayIonizationRateCI(double xe, double Ntot){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return 0.e0;
}

double SoftXrayIonizationRateOI(double xe, double Ntot){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return 0.e0;
}

double SoftXrayIonizationRateSiI(double xe, double Ntot){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return 0.e0;
}

double SoftXrayIonizationRateSiII(double xe, double Ntot){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return 0.e0;
}

/// Soft X-ray ionization heating per H [erg/s/H]
static double SoftXrayHeatingRateBody(double xe, double Ntot){
    // eq(A9)(A10) of Wolfire et al. (1995)
    //// note: this fitting formula might be obtained in the equilibrium state (Wolfire et al. 1995)
    if(xe<=0.0) return 0.0;
    //if(xe>0.5)  return 0.e0;
    double p1 = MAX(-0.3,log10(1.e-18*Ntot));
    double p2 = log10(xe);
    double f6 = 0.990 + (-2.74e-3 + 1.13e-3*p2)*p2;
    return pow(10.0,f6*(-26.5+(-0.920+5.89e-2*p1)*p1) + f6*0.96*exp(-SQ((p1-0.38)/0.87))); // [erg/s]
}

/// Soft X-ray heating [erg/cc/s]
//double TotalSoftXrayHeatingRate(double xe, double nH, double Ntot){
double TotalSoftXrayHeatingRate(StructEQCData *PhydroPP){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    double nH   = PhydroPP->n[ipHYDROGEN];
    double xe   = PhydroPP->nel/nH;
    double Ntot = PhydroPP->Nt;
    double X0   = PhydroPP->X0;
    //// note: this fitting formula might be obtained in the equilibrium state (Wolfire et al. 1995)
    return nH*X0*SoftXrayHeatingRateBody(xe,Ntot); // [erg/s/cc]
}

double TotalSoftXrayHeatingRateFixed(double nH, double X0){
    if(EQCLibRunParameters.UseSoftXrayIonization==false)    return 0.e0;
    return nH*X0*(1.6e-17); // [erg/s/cc]
}

#if 0
// eq.(5) and Table 1 of Wolfire et al. (1995)
static double T1 = 1.e+6.16, EM1 = 5.3e-3, N1 = 0.e0;
static double T2 = 1.e+6.33, EM2 = 2.5e-3, N2 = 3.6e+20; 
static double N3 = 3.9e+20;
double SoftXrayIntensity_4piJnyu(double nyu){
    double hnyu = PLANK_H_CGS*nyu;
    double Lamb1;
    double Lamb2;
    double C3;
    return 0.5*(Lamb1*EM1 + Lamb2*EM2*exp(-sig_nyu*N2) + C3*pow(hnyu,-0.4)*exp(-sig_nyu)*N3);
}
#endif
