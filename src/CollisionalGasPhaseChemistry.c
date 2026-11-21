#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Astro.h"
#include "Constants.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "CollisionalGasPhaseChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "CosmicRay.h"
#include "Xray.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

//#define RATES_GLOVER_JAPPSEN_2007

#define T_HI    (157807.0) // [K]
#define T_HeI   (285335.0) // [K]
#define T_HeII  (631515.0) // [K]

void InitializeCollisionalIonization(int Check){
        
    if(Check==1){
        FILE *output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CollisionalIonization.dat","w");
        for(double logT=1.0; logT<=6.0; logT+=0.1){
            double T = pow(10.0,logT);
            fprintf(output,"%e ",T);
            fprintf(output,"%e ",HICollisionalIonizationRate(T));
            fprintf(output,"%e ",HeICollisionalIonizationRate(T));
            fprintf(output,"%e ",HeIICollisionalIonizationRate(T));
            fprintf(output,"%e ",CICollisionalIonizationRate(T));
            fprintf(output,"%e ",CIICollisionalIonizationRate(T));
            fprintf(output,"%e ",OICollisionalIonizationRate(T));
            fprintf(output,"%e ",OIICollisionalIonizationRate(T));
            fprintf(output,"\n");
        }
        fclose(output);

        FILE *output2 = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CollisionalRecombination.dat","w");
        for(double logT=1.0; logT<=6.0; logT+=0.1){
            double T = pow(10.0,logT);
            fprintf(output2,"%e ",T);
            fprintf(output2,"%e ",HIIRecombinationRate(T));
            fprintf(output2,"%e ",HIIRecombinationRateCaseB(T));
            fprintf(output2,"%e ",HeIIRecombinationRate(T));
            fprintf(output2,"%e ",HeIIRecombinationRateCaseB(T));
            fprintf(output2,"%e ",HeIIIRecombinationRate(T));
            fprintf(output2,"%e ",HeIIIRecombinationRateCaseB(T));
            fprintf(output2,"%e ",CIIRecombinationRate(T));
            fprintf(output2,"%e ",CIIIRecombinationRate(T));
            fprintf(output2,"%e ",OIIRecombinationRate(T));
            fprintf(output2,"%e ",OIIIRecombinationRate(T));
            fprintf(output2,"\n");
        }
        fclose(output2);
    }
    //exit(1);

    return ;
}


// eq.(1) and Table 1 of Voronov (1997)
static double CollisionalIonizationRateVoronov1997(double T, double dE, double P, double A, double X, double K){
    double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    double U  = dE/Te;
    return A*(1+P*sqrt(U))*pow(U,K)*exp(-U)/(X+U);
}

//// Ionization rate of collision with electrons [cm^3/s]
double HICollisionalIonizationRate(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    // R11 in Table 1 of Glover & Jappsen (2007)
    double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    double lnT  = log(Te),  lnT2 = SQ(lnT),  lnT3 = lnT*lnT2, lnT4 = lnT*lnT3,
           lnT5 = lnT*lnT5, lnT6 = lnT*lnT5, lnT7 = lnT*lnT6, lnT8 = lnT*lnT7;
    return exp(- 3.271396786e+1      + 1.35365560e+1*lnT
               - 5.739328750e+0*lnT2 + 1.56315498e+0*lnT3
               - 2.877056000e-1*lnT4 + 3.48255977e-2*lnT5
               - 2.631976170e-3*lnT6 + 1.11954395e-4*lnT7
               - 2.039149850e-6*lnT8
               );
#else
    //// ref. Cen (1992); see also Abel et al. (1997)
    return 5.85e-11*sqrt(T)*exp(-T_HI/T)/(1.0+sqrt(1.e-5*T));
    //return CollisionalIonizationRateVoronov1997(T,13.6,0.0,0.291e-7,0.232,0.39);
#endif
}
double HeICollisionalIonizationRate(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    // R24 in Table 1 of Glover & Jappsen (2007)
    double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    double lnT  = log(Te),  lnT2 = SQ(lnT),  lnT3 = lnT*lnT2, lnT4 = lnT*lnT3,
           lnT5 = lnT*lnT5, lnT6 = lnT*lnT5, lnT7 = lnT*lnT6, lnT8 = lnT*lnT7;
    return exp(- 4.409864886e+1      + 2.391596563e+1*lnT
               - 1.075323020e+1*lnT2 + 3.058038750e+0*lnT3
               - 5.685118900e-1*lnT4 + 6.795391230e-2*lnT5
               - 5.009056100e-3*lnT6 + 2.067236160e-4*lnT7
               - 3.649161410e-6*lnT8
               );
#else
    //// ref. eq.(26b) of Cen (1992); Abel et al. (1997)
    return 2.38e-11*sqrt(T)*exp(-T_HeI/T)/(1.0+sqrt(1.e-5*T));    
    //return CollisionalIonizationRateVoronov1997(T,24.6,0.0,0.175e-7,0.180,0.35);
#endif
}
double HeIICollisionalIonizationRate(double T){
    //// ref. eq.(26c) of Cen (1992)
    return 5.68e-12*sqrt(T)*exp(-T_HeII/T)/(1.0+sqrt(1.e-5*T));
    //return CollisionalIonizationRateVoronov1997(T,54.4,1.0,0.205e-8,0.265,0.25);
}

double CICollisionalIonizationRate(double T){
    // ref. Voronov (1997); R33 in Table 1 of Glover & Jappsen (2007) 
    return CollisionalIonizationRateVoronov1997(T,11.3,0.0,0.685e-7,0.193,0.25);
    //double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    //double u  = 11.26/Te;
    //return 6.85e-8*pow(u,0.25)*exp(-u)/(0.193 + u);
}
double CIICollisionalIonizationRate(double T){
    // ref. Voronov (1997)
    return CollisionalIonizationRateVoronov1997(T,24.4,1.0,0.186e-7,0.286,0.24);
}
double CIIICollisionalIonizationRate(double T){
    // ref. Voronov (1997)
    return CollisionalIonizationRateVoronov1997(T,47.9,1.0,0.635e-8,0.427,0.21);
}
double CIVCollisionalIonizationRate(double T){
    // ref. Voronov (1997)
    return CollisionalIonizationRateVoronov1997(T,64.5,1.0,0.150e-8,0.416,0.13);
}
double CVCollisionalIonizationRate(double T){
    // ref. Voronov (1997)
    return CollisionalIonizationRateVoronov1997(T,392.1,1.0,0.299e-9,0.666,0.02);
}
double CVICollisionalIonizationRate(double T){
    // ref. Voronov (1997)
    return CollisionalIonizationRateVoronov1997(T,490.0,1.0,0.123e-9,0.620,0.16);
}

double NICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,14.5,0.0,0.482e-7,0.0652,0.42);
}
double NIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,29.6,0.0,0.298e-7,0.310,0.30);
}
double NIIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,47.5,1.0,0.810e-8,0.350,0.24);
}
double NIVCollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,77.5,1.0,0.371e-8,0.549,0.18);
}
double NVCollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,97.9,0.0,0.151e-8,0.0167,0.74);
}
double NVICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,552.1,0.0,0.371e-9,0.546,0.29);
}
double NVIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,667.0,1.0,0.777e-10,0.624,0.16);
}

double OICollisionalIonizationRate(double T){
    // ref. Voronov (1997); R35 in Table 1 of Glover & Jappsen (2007) 
    return CollisionalIonizationRateVoronov1997(T,13.6,0.0,0.359e-7,0.073,0.34);
    //double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    //double u  = 13.6/Te;
    //return 3.59e-8*pow(u,0.34)*exp(-u)/(0.073 + u);
}
double OIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,35.1,1.0,0.139e-7,0.212,0.22);
}
double OIIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,54.9,1.0,0.931e-8,0.270,0.27);
}
double OIVCollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,77.4,0.0,0.102e-7,0.614,0.27);
}
double OVCollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,113.9,1.0,0.219e-8,0.630,0.17);
}
double OVICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,138.1,0.0,0.195e-8,0.360,0.54);
}
double OVIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,739.3,0.0,0.212e-9,0.396,0.35);
}
double OVIIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,871.4,1.0,0.521e-10,0.629,0.16);
}

double SiICollisionalIonizationRate(double T){
    // ref. Voronov (1997); R34 in Table 1 of Glover & Jappsen (2007) 
    return CollisionalIonizationRateVoronov1997(T,8.2,1.0,0.188e-6,0.376,0.25);
    //double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    //double u  = 8.2/Te;
    //return 1.88e-7*(1.e0+sqrt(u))*pow(u,0.25)*exp(-u)/(0.376 + u);
}
double SiIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,16.4,1.0,0.643e-7,0.632,0.20);
}
double SiIIICollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,33.5,1.0,0.201e-7,0.473,0.22);
}
double SiIVCollisionalIonizationRate(double T){
    return CollisionalIonizationRateVoronov1997(T,54.0,1.0,0.494e-8,0.172,0.23);
}
#if 0
// ref. eq.(6)(7) and Tables 2,3 in Lennon et al. (1988)
//      see also p.134 of Draine's textbook
static double Tion_CI = 11.3*EV_TO_ERG/BOLTZMANN_CONSTANT_CGS;  // [K]
static double aCI[] = {5.9849e-8,1.1903e-8,-3.0141e-8,-1.3693e-8,8.3749e-9,4.0150e-9};
static double alphaCI = 3.7442e-7;
static double betaCI[] = {-6.5826e-7,2.0520e-6,-4.4688e-6};
double CICollisionalIonizationRate(double T){
    double kci = 0.0;
    double kToverI = T/Tion_CI;
    double IoverkT = Tion_CI/T;
    if(kToverI<0.1)
        return 0.0;
    else if(kToverI<=10.0){
        for(int i=0; i<=5; i++)
            kci += aCI[i]*pow(log10(kToverI),i);
        return kci*exp(-IoverkT)*sqrt(kToverI);
    }
    else{
        kci = alphaCI*log(kToverI);
        for(int i=0; i<=2; i++)
            kci += betaCI[i]*pow(IoverkT,i);
        return kci*sqrt(IoverkT);
    }
}
// ref. eq.(6)(7) and Tables 2,3 in Lennon et al.,JPCRD,17,1285,(1988)
//      see also p.134 of Draine's textbook
static double Tion_OI = 13.6*EV_TO_ERG/BOLTZMANN_CONSTANT_CGS;  // [K]
static double aOI[] = {3.3523e-8,1.3434e-8,-6.7039e-9,-1.9954e-8,1.6196e-9,6.5780e-9};
static double alphaOI = 3.2684e-7;
static double betaOI[] = {-6.7409e-7,2.3910e-6,-6.0366e-6};
double OICollisionalIonizationRate(double T){
    double kci = 0.0;
    double kToverI = T/Tion_OI;
    double IoverkT = Tion_OI/T;
    if(kToverI<0.1)
        return 0.0;
    else if(kToverI<=10.0){
        for(int i=0; i<=5; i++)
            kci += aOI[i]*pow(log10(kToverI),i);
        return kci*exp(-IoverkT)*sqrt(kToverI);
    }
    else{
        kci = alphaOI*log(kToverI);
        for(int i=0; i<=2; i++)
            kci += betaOI[i]*pow(IoverkT,i);
        return kci*sqrt(IoverkT);
    }
}
#endif

//// Recombination Rate [cm^3/s]
double HIIRecombinationRate(double T){
    if(EQCLibRunParameters.UseRecombinationCaseB==true)
        return HIIRecombinationRateCaseB(T);
    else
        return HIIRecombinationRateCaseA(T);
}
double HIIRecombinationRateCaseA(double T){
    // ref. Hui & Gnedin (1997); Ferland et al. (1992)
    // R13 in Table 1 of Glover & Jappsen (2007)
    double lamb = 2.0*(T_HI)/T;
    return 1.269e-13*pow(lamb,1.503)*pow(1.0+pow(1.91571*lamb,0.470),-1.923);
}
double HIIRecombinationRateCaseB(double T){
    // ref. Hui & Gnedin (1997); Ferland et al. (1992)
    // R13 in Table 1 of Glover & Jappsen (2007)
    double lamb = 2.0*(T_HI)/T;
    return 2.753e-14*pow(lamb,1.500)*pow(1.0+pow(0.36496*lamb,0.407),-2.242);
}
double HeIIRecombinationRate(double T){
    if(EQCLibRunParameters.UseRecombinationCaseB==true)
        return HeIIRecombinationRateCaseB(T);
    else
        return HeIIRecombinationRateCaseA(T);
}
double HeIIRecombinationRateCaseA(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    // R25 in Table 1 of Glover & Jappsen (2007)
    double lnT = log(T), lnT2 = SQ(lnT), lnT3 = lnT*lnT2;
    return 1.e-11*pow(T,-0.5)*(12.72-1.615*lnT-0.3162*lnT2+0.0493*lnT3);
#else
    // ref. Hui & Gnedin (1997); Burgess & Seaton (1960)
    double lamb = 2.0*(T_HeI)/T;
    return (3.0e-14)*pow(lamb,0.654);
#endif
}
double HeIIRecombinationRateCaseB(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    // R25 in Table 1 of Glover & Jappsen (2007)
    double lnT = log(T), lnT2 = SQ(lnT), lnT3 = lnT*lnT2;
    return 1.e-11*pow(T,-0.5)*(11.19-1.676*lnT-0.2852*lnT2+0.04433*lnT3);
#else
    // ref. Hui & Gnedin (1997); Burgess & Seaton (1960)
    double lamb = 2.0*(T_HeI)/T;
    return 1.26e-14*pow(lamb,0.750);
#endif
}
double HeIIDielectronicRecombinationRate(double T){
    // ref. eq.(28) in Cen (1992)
    // R25 in Table 1 of Glover & Jappsen (2007)
    double invT = T_HeII/T;
    return 1.90e-3*pow(T,-1.5)*exp(-0.75*invT)*(1.e0+0.3*exp(-0.15*invT));
}
double HeIIIRecombinationRate(double T){
    if(EQCLibRunParameters.UseRecombinationCaseB==true)
        return HeIIIRecombinationRateCaseB(T);
    else
        return HeIIIRecombinationRateCaseA(T);
}
double HeIIIRecombinationRateCaseA(double T){
    // ref. Hui & Gnedin (1997)
    // The recombination rate for He++ in Cen (1992) has an incorrect temperature dependence. 
    // See Fukugita & Kawasaki (1994) for the correct scaling.
    double lamb = 2.0*(T_HeII)/T;
    return 2.538e-13*pow(lamb,1.503)*pow(1.0+pow(1.91571*lamb,0.470),-1.923);
}
double HeIIIRecombinationRateCaseB(double T){
    // ref. Hui & Gnedin (1997)
    double lamb = 2.0*(T_HeII)/T;
    return 5.506e-14*pow(lamb,1.500)*pow(1.0+pow(0.36496*lamb,0.407),-2.242);
}
double CIIRecombinationRate(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    // radiative + dielectronic recombination 
    // ref. Nahar & Pradhan (1997); R30 in Table 1 of Glover & Jappsen (2007)
    if(T<=7950)         return 4.67e-12*pow(T/300,-0.6 );
    else if(T<=21140)   return 1.23e-17*pow(T/300, 2.49)*exp(  21845.6/T); 
    else                return 9.62e-8 *pow(T/300,-1.37)*exp(-115786.2/T);
#else
    // ref. Millar et al., A&AS, 121, 139, (1997)
    //      eq.(1) and Table 4; Reaction 3208
    return 4.40e-12*pow(3.3333e-3*T,-0.61);
#endif
}
double CIIIRecombinationRate(double T){
    // original fitting to Nahar (1999)
    // radiative + dielectronic recombination 
    double invT = 1.e0/T;
    return 2.0e-11*pow(T/300,-0.61) + 2.6e-6*pow(T,-1.5)*exp(-5150*invT)*(1+0.5*exp(-1.e+3*invT))
         + 8.5e-3*pow(T,-1.7)*exp(-125000.0*invT)*(1.+4.5*exp(-55000.0*invT)) + 1.9e-12*exp(-1.e-6*T);
    //fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
    //exit(1);
    //// radiative + dielectronic recombination 
    //// ref. Verner (1999); Table 14.7 of Draine's text
    //return 0.e0;
}
double OIIRecombinationRate(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    // radiative + dielectronic recombination 
    // ref. Nahar (1999); R32 in Table 1 of Glover & Jappsen (2007)
    double invT = 1.e0/T;
    return 1.41e-10*pow(T,-0.66) + 7.4e-4*pow(T,-1.5)*exp(-175000.0*invT)*(1.0+0.062*exp(-145000.0*invT));
    //if(T<=400)  return 1.30e-10*pow(T,-0.64);
    //else        return 1.41e-10*pow(T,-0.66) 
    //                 + 7.4e-4*pow(T,-1.5)*exp(-175000/T)*(1.0+0.062*exp(-145000/T));
#else
    // ref. Millar et al., A&AS, 121, 139, (1997)
    //      eq.(1) and Table 4; Reaction 3218
    return 3.40e-12*pow(3.3333e-3*T,-0.63);
#endif
}
double OIIIRecombinationRate(double T){
    // original fitting to Nahar (1999)
    // radiative + dielectronic recombination 
    double invT = 1.e0/T;
    return 5.5e-10*pow(T,-0.585) + 4.8e-3*pow(T,-1.55)*exp(-205000.0*invT)*(1.0+0.55*exp(-150000.0*invT));
    //fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
    //exit(1);
    //// radiative + dielectronic recombination 
    //// ref. Verner (1999); Table 14.7 of Draine's text
    //return 0.0;
}
double SiIIRecombinationRate(double T){
    // radiative + dielectronic recombination 
    // ref. Nahar (2000); R31 in Table 1 of Glover & Jappsen (2007)
    if(T<=2000)         return 7.50e-12*pow(T/300,-0.55 );
    else if(T<=1.e+4)   return 4.86e-12*pow(T/300,-0.32 );
    else                return 9.08e-14*pow(T/300, 0.818);
}
double SiIIIRecombinationRate(double T){
    // radiative recombination
    // ref. Nahar (1995,1996); R47 in Table 1 of Glover & Jappsen (2007)
    double rr = 1.75e-12*pow(1.e-4*T,-0.6346);
   
    // dielectronic recombination
    // ref. Mazzotta (1998); R47 in Table 1 of Glover & Jappsen (2007)
    double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    double invTe = 1.e0/Te;
    double di = pow(Te,-1.5)*(2.2552e-11*exp(-2.76*invTe) + 5.6058e-9*exp(-10.13*invTe)); 
    
    return rr + di;
}

///////// Charge Transfer Rate [cm^3/s]
//// H+ + He -> He+ + H
double HeIChargeTransferRateToHII(double T){
    /// R27 in Table 1 of Glover & Jappsen (2007)
    if(T<=1.e+4)    return 1.26e-9*pow(T,-0.75)*exp(-127500/T); 
    else            return 4.0e-37*pow(T, 4.74);
}
//// He+ + H -> He + H+
double HeIIChargeTransferRateToHI(double T){
    /// R26 in Table 1 of Glover & Jappsen (2007)
    return 1.25e-15*pow(T/300.0,0.25); 
}

//// H+ + O -> O+ + H
double OIChargeTransferRateToHII(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    /// ref. Stancil et al. (1998); R37 in Table 1 of  & Jappsen (2007)
    return (1.08e-11*pow(T,0.517) + 4.00e-10*pow(T,0.00669))*exp(-227.0/T);
#else
    // ref. Millar et al., A&AS, 121, 139, (1997)
    return 7.0e-10*exp(-232.0/T);
#endif
}
//// O+ + H -> H+ + O
double OIIChargeTransferRateToHI(double T){
#ifdef RATES_GLOVER_JAPPSEN_2007
    /// ref. Stancil et al. (1999); R36 in Table 1 of Glover & Jappsen (2007)
    return 4.99e-11*pow(T,0.405) + 7.45e-10*pow(T,-0.458);
#else
    // ref. Millar et al., A&AS, 121, 139, (1997)
    //return 6.80e-10;
    return 7.7e-10;
#endif
}
//// O++ + H -> O+ + H+
double OIIIChargeTransferRateToHI(double T){
    /// ref. Table 1 of Kingdon & Ferland (1996), fitted for 1e2-1e5 K
    const double T4 = fmin(10.0,fmax(1.e-2,1.e-4*T));
    return 1.04e-9*pow(T4,0.27)*(1.0+2.02*exp(-5.92*T4));
}
//// He+ + O -> O+ + He
double OIChargeTransferRateToHeII(double T){
    /// ref. Zhao et al. (2004); Glover & Jappsen (2007)
    return 4.991e-15*pow(1.e-4*T, 0.3794)*exp(-T/1121000.0) 
         + 2.780e-15*pow(1.e-4*T,-0.2163)*exp( T/815800.0);
}
//// H+ + C -> C+ + H
double CIChargeTransferRateToHII(double T){
    /// ref. Stancil et al. (1998); Glover & Jappsen (2007)
    return 3.9e-16*pow(T,0.213); 
}
//// H + C+ -> C + H+
double CIIChargeTransferRateToHI(double T){
    /// ref. Stancil et al. (1998); Glover & Jappsen (2007)
    return 6.08e-14*pow(1.e-4*T,1.96)*exp(-170000.0/T);
}
//// He+ + C -> C+ + He
double CIChargeTransferRateToHeII(double T){
    /// ref. Kimura et al. (1993); Glover & Jappsen (2007)
    if(T<=200)          return 8.58e-17*pow(T,0.757);
    else if(T<=2000)    return 3.25e-17*pow(T,0.968);
    else                return 2.77e-19*pow(T,1.597);
}
//// H+ + Si -> Si+ + H
double SiIChargeTransferRateToHII(double T){
    /// ref. Kingdon & Ferland (1996); Glover & Jappsen (2007)
    if(T<=1.e+4)    return 5.88e-13*pow(T,0.848);
    else            return 1.45e-13*T;
}
//// He+ + Si -> Si+ + He
double SiIChargeTransferRateToHeII(double T){
    /// ref. Glover & Jappsen (2007)
    return 3.3e-9;
}
//// C+ + Si -> Si+ + C
double SiIChargeTransferRateToCII(double T){
    /// ref. Glover & Jappsen (2007)
    return 2.1e-9;
}
//// H+ + Si+ -> Si++ + H
double SiIIChargeTransferRateToHII(double T){
    /// ref. Glover & Jappsen (2007)
    double Te = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T; 
    return 4.10e-10*pow(1.e-4*T,0.24)*(1.0+3.17*exp(1.e-6*T/2.39))*exp(-3.178/Te);
}
//// H + Si++ -> Si+ + H+
double SiIIIChargeTransferRateToHI(double T){
    /// ref. Glover & Jappsen (2007)
    return 1.23e-10*pow(1.e-4*T,0.24)*(1.0+3.17*exp(1.e-6*T/2.39));
}


//////////////////////////
//// Collisional Ionization Cooling [erg/cc/s]
// ref. Cen (1992)
static double HydrogenIonizationCooling(double T){
    return 1.27e-21*sqrt(T)*exp(-T_HI/T)/(1.0+sqrt(1.e-5*T));
}
// ref. eq.(12b) of Cen (1992)
static double NeutralHeliumIonizationCooling(double T){
    return 9.38e-22*sqrt(T)*exp(-T_HeI/T)/(1.0+sqrt(1.e-5*T));
}
// ref. eq.(12c) of Cen (1992)
static double SinglyIonizedHeliumIonizationCooling(double T){
    return 4.95e-22*sqrt(T)*exp(-T_HeII/T)/(1.0+sqrt(1.e-5*T));
}

//// Collisional Excitation Cooling [erg/cc/s]
// ref.) eq(6-12) in p. 140 of Spitzer (1978)
//       modified by Cen (1992)
// electron excitation pf the excited levels in neutral H can be 
// an important source of  energy loss. The Ly-a photons
static double HydrogenExcitationCooling(double T){
    return 7.5e-19*exp(-118348.0/T)/(1.0+sqrt(1.e-5*T));
}
// Collisional Excitation Cooling of HeI (n=2,3,4 triplets) [erg/cc/s]
// ref. eq.(15c) of Cen (1992) 
static double HeliumIExcitationCooling(double T){
    return 9.10e-27*pow(T,-0.1687)*exp(-13179.0/T)/(1.0+sqrt(1.e-5*T));
}
// Collisional Excitation Cooling of HeII (n=2) [erg/cc/s]
// ref. eq.(15c) of Cen (1992) 
static double HeliumIIExcitationCooling(double T){
    return 5.54e-17*pow(T,-0.397)*exp(-473638.0/T)/(1.0+sqrt(1.e-5*T));
}

// ref. Hui & Gnedin (1997)
static double HydrogenRecombinationCoolingCaseA(double T){
    const double lambda = 2.0*T_HI/T;
    return 1.778e-29*T*pow(lambda,1.965)*pow(1.0+pow(1.84843*lambda,0.502),-2.697);
}

static double HydrogenRecombinationCoolingCaseB(double T){
    const double lambda = 2.0*T_HI/T;
    return 3.435e-30*T*pow(lambda,1.970)*pow(1.0+pow(0.44444*lambda,0.376),-3.720);
}
// ref. Hui & Gnedin (1997)
static double HeliumIIRecombinationCoolingCaseA(double T){
    return BOLTZMANN_CONSTANT_CGS*T*HeIIRecombinationRateCaseA(T);
}

static double HeliumIIRecombinationCoolingCaseB(double T){
    return BOLTZMANN_CONSTANT_CGS*T*HeIIRecombinationRateCaseB(T);
}

/* Hui & Gnedin (1997). Lambda_HeIII(T) = 8 Lambda_HII(T/4) in the
   hydrogenic scaling. The fits carry a factor of T. */
static double HeliumIIIRecombinationCoolingCaseA(double T){
    const double lambda = 2.0*T_HeII/T;
    return 3.556e-29*T*pow(lambda,1.965)*pow(1.0+pow(1.84843*lambda,0.502),-2.697);
}

static double HeliumIIIRecombinationCoolingCaseB(double T){
    const double lambda = 2.0*T_HeII/T;
    return 6.870e-30*T*pow(lambda,1.970)*pow(1.0+pow(0.44444*lambda,0.376),-3.720);
}

// ref. eq.(14) of Cen (1992)
static double HeliumIIDielectronicRecombinationCooling(double T){
    return 0.75*BOLTZMANN_CONSTANT_CGS*T_HeII*HeIIDielectronicRecombinationRate(T);
}

double CoolingCollisionalIonization(StructEQCData *state){
    return state->nel*state->natom[ipHYDROGEN]*HydrogenIonizationCooling(state->T);
}

double CoolingHeICollisionalIonization(StructEQCData *state){
    return state->nel*state->natom[ipHELIUM]*NeutralHeliumIonizationCooling(state->T);
}

double CoolingHeIICollisionalIonization(StructEQCData *state){
    return state->nel*state->nion[ipHELIUM]*SinglyIonizedHeliumIonizationCooling(state->T);
}

double CoolingCollisionalExcitation(StructEQCData *state){
    return state->nel*state->natom[ipHYDROGEN]*HydrogenExcitationCooling(state->T);
}

double CoolingHeICollisionalExcitation(StructEQCData *state){
    return SQ(state->nel)*state->nion[ipHELIUM]*HeliumIExcitationCooling(state->T);
}

double CoolingHeIICollisionalExcitation(StructEQCData *state){
    return state->nel*state->nion[ipHELIUM]*HeliumIIExcitationCooling(state->T);
}

double CoolingRecombination(StructEQCData *state){
    const double eta = EQCLibRunParameters.UseRecombinationCaseB ? HydrogenRecombinationCoolingCaseB(state->T)
                                                               : HydrogenRecombinationCoolingCaseA(state->T);
    return state->nion[ipHYDROGEN]*state->nel*eta;
}

double CoolingHeIIRecombination(StructEQCData *state){
    const double eta = EQCLibRunParameters.UseRecombinationCaseB ? HeliumIIRecombinationCoolingCaseB(state->T)
                                                               : HeliumIIRecombinationCoolingCaseA(state->T);
    return state->nion[ipHELIUM]*state->nel*eta;
}

double CoolingHeIIIRecombination(StructEQCData *state){
    const double eta = EQCLibRunParameters.UseRecombinationCaseB ? HeliumIIIRecombinationCoolingCaseB(state->T)
                                                               : HeliumIIIRecombinationCoolingCaseA(state->T);
    return state->nHeIII*state->nel*eta;
}

double CoolingHeIIDielectronicRecombination(StructEQCData *state){
    return state->nel*state->nion[ipHELIUM]*HeliumIIDielectronicRecombinationCooling(state->T);
}

/* Combine coefficients only when they multiply the same collider and ion densities. */
static double NeutralHydrogenCooling(double T){
    return HydrogenExcitationCooling(T)+HydrogenIonizationCooling(T);
}

static double SinglyIonizedHeliumCooling(double T){
    return HeliumIIExcitationCooling(T)+SinglyIonizedHeliumIonizationCooling(T);
}

static double TotalHeliumIIRecombinationCoolingCaseA(double T){
    return HeliumIIRecombinationCoolingCaseA(T)+HeliumIIDielectronicRecombinationCooling(T);
}

static double TotalHeliumIIRecombinationCoolingCaseB(double T){
    return HeliumIIRecombinationCoolingCaseB(T)+HeliumIIDielectronicRecombinationCooling(T);
}

static CEQTemperatureTable neutralHydrogenCoolingTable, hydrogenRecombinationCoolingTable[2];
static CEQTemperatureTable neutralHeliumIonizationCoolingTable, neutralHeliumExcitationCoolingTable;
static CEQTemperatureTable singlyIonizedHeliumCoolingTable, heliumIIRecombinationCoolingTable[2];
static CEQTemperatureTable heliumIIIRecombinationCoolingTable[2];

void CEQInitializeAtomicCoolingTemperatureTables(void){
    CEQFillTemperatureTable(&neutralHydrogenCoolingTable,NeutralHydrogenCooling,NULL,0);
    CEQFillTemperatureTable(&hydrogenRecombinationCoolingTable[0],HydrogenRecombinationCoolingCaseA,NULL,0);
    CEQFillTemperatureTable(&hydrogenRecombinationCoolingTable[1],HydrogenRecombinationCoolingCaseB,NULL,0);
    CEQFillTemperatureTable(&neutralHeliumIonizationCoolingTable,NeutralHeliumIonizationCooling,NULL,0);
    CEQFillTemperatureTable(&neutralHeliumExcitationCoolingTable,HeliumIExcitationCooling,NULL,0);
    CEQFillTemperatureTable(&singlyIonizedHeliumCoolingTable,SinglyIonizedHeliumCooling,NULL,0);
    CEQFillTemperatureTable(&heliumIIRecombinationCoolingTable[0],TotalHeliumIIRecombinationCoolingCaseA,NULL,0);
    CEQFillTemperatureTable(&heliumIIRecombinationCoolingTable[1],TotalHeliumIIRecombinationCoolingCaseB,NULL,0);
    CEQFillTemperatureTable(&heliumIIIRecombinationCoolingTable[0],HeliumIIIRecombinationCoolingCaseA,NULL,0);
    CEQFillTemperatureTable(&heliumIIIRecombinationCoolingTable[1],HeliumIIIRecombinationCoolingCaseB,NULL,0);
}

void CEQAtomicCoolingWithTables(StructEQCData *state, const CEQTemperaturePoint *point){
    const double T = state->T, ne = state->nel;
    const bool caseB = EQCLibRunParameters.UseRecombinationCaseB;
    const double nHI = state->natom[ipHYDROGEN], nHII = state->nion[ipHYDROGEN];
    const double nHeI = state->natom[ipHELIUM], nHeII = state->nion[ipHELIUM], nHeIII = state->nHeIII;
    state->L_HI = state->L_HeI = state->L_HeII = 0.0;
    if(ne == 0.0) return;
    if(nHI != 0.0)
        state->L_HI = ne*nHI*CEQReadTemperatureTable(&neutralHydrogenCoolingTable,point,T,NeutralHydrogenCooling);
    if(nHII != 0.0)
        state->L_HI += nHII*ne*CEQReadTemperatureTable(&hydrogenRecombinationCoolingTable[caseB],point,T,
                                                      caseB ? HydrogenRecombinationCoolingCaseB : HydrogenRecombinationCoolingCaseA);
    if(nHeI != 0.0)
        state->L_HeI = ne*nHeI*CEQReadTemperatureTable(&neutralHeliumIonizationCoolingTable,point,T,NeutralHeliumIonizationCooling);
    if(nHeII != 0.0){
        state->L_HeI += SQ(ne)*nHeII*CEQReadTemperatureTable(&neutralHeliumExcitationCoolingTable,point,T,HeliumIExcitationCooling)
                       + nHeII*ne*CEQReadTemperatureTable(&heliumIIRecombinationCoolingTable[caseB],point,T,
                                                         caseB ? TotalHeliumIIRecombinationCoolingCaseB : TotalHeliumIIRecombinationCoolingCaseA);
        state->L_HeII = ne*nHeII*CEQReadTemperatureTable(&singlyIonizedHeliumCoolingTable,point,T,SinglyIonizedHeliumCooling);
    }
    if(nHeIII != 0.0)
        state->L_HeII += nHeIII*ne*CEQReadTemperatureTable(&heliumIIIRecombinationCoolingTable[caseB],point,T,
                                                          caseB ? HeliumIIIRecombinationCoolingCaseB : HeliumIIIRecombinationCoolingCaseA);
}

