#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalElements.h"
#include "ChemicalEquilibriumCooling.h"
#include "HIIregionLineCooling.h"
#include "StatisticalEquilibrium.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


// Principal Collisionally Excited Cooling Lines in HII regions
//  lines : Table 27.2 in Draine's textbook
//  A_ul  : Table 3.8, 3.9 in Osterbrock's textbook (1989)
//  n_crit: Table 18.1, 18.2 in Draine'x textbook
//  electreon collision strength for OII: Tayal (2007), p.499 of Draine's textbook
// ion      lins                    A_ul [1/s]  n_cirt_u [1/cc] (T=1.e+4 K)
//--------------------------------------------------------------------------
// NII      6585A (1D2->3P2)         3.0e-3     7.68e+4
//          6550A (1D2->3P1)         1.0e-3     7.68e+4
// OII      3730A (2D5/2->4S3/2)     3.6e-5     3.31e+3
//          3727A (2D3/2->4S3/2)     1.8e-4     4.49e+3
// OIII     88.36um (3P1->3P0)       2.6e-5     1.74e+3
//          51.81um (3P2->3P1)       9.8e-5     3.79e+3
//          5008A (1D2->3P2)         2.0e-2     6.40e+5
//          4960A (1D2->3P1)         6.7e-3     6.40e+5
// NeII     12.81um         
// SII      6733A (2D3/2->4S3/2)     8.8e-4     1.49e+4
//          6718A (2D5/2->4S3/2)     2.6e-4     1.57e+3
// SIII     33.48um (3P1->3P0)       4.7e-4     4.22e+3
//          18.71um (3P2->3P1)       2.1e-3     1.31e+4
//          9071A (1D2->3P1)         2.2e-2     7.33e+5
//          9533A (1D2->3P2)         5.8e-2     7.33e+5
//--------------------------------------------------------------------------

//static double LambdaOIII[4];
//static double E_OIII[4];
//static double Aul_OIII[4];
//
//static double LambdaOII[2];
//static double E_OII[2];
//static double Aul_OII[2];
//
//// electreon collision strength for OII: Tayal (2007), p.499 of Draine's textbook
//// coefficients for collision strength
//static double aOII[2], bOII[2], cOII[2];     
//
//// statistical weight of the upper and lower levels
//static double WeightOII1[2];
//static double WeightOII2[2];

//#define OIII_OII_RATIO 2.0  // assumed a typical factor (see Fig.27.1 in Draine's textbook)

// Spontaneous Emission Rate [1/sec]
#define AulOIII32 (2.0e-2)    // n = 3 -> 2 5008 
#define AulOIII31 (6.7e-3)    // n = 3 -> 1 4960
#define AulOIII21 (9.8e-5)    // n = 2 -> 1 52um
#define AulOIII10 (2.6e-5)    // n = 1 -> 0 88um

#define AulOII20  (1.8e-2)    // n = 2 -> 0
#define AulOII10  (3.6e-2)    // n = 1 -> 0

// Excitation Temperature [K]
// ref. Fig. 18.1, 18.3 of Draine's text
#define TOIII32   (28728.0)
#define TOIII31   (29006.0)
#define TOIII21   (278.0)
#define TOIII10   (163.0)

#define TOII20    (38604.0)
#define TOII10    (38575.0)

// Statistical Weight Ratio: g1/g0 (upper/lower)
// ref. Fig. 18.1, 18.3 of Draine's text
#define gulOIII32 (1.0000)
#define gulOIII31 (1.6667)
#define gulOIII21 (1.6667)
#define gulOIII10 (3.0000)

#define gulOII20  (1.0000)
#define gulOII10  (1.5000)

//static double Energy_;     
//static double IntensityCMB_;
//static double EinsteinA10_; 
//static double EinsteinB10_;
//static double EinsteinB01_;

static double Lambda_OIII32;     
static double Lambda_OIII31;     
static double Lambda_OIII21;     
static double Lambda_OIII10;     
static double Lambda_OII20;     
static double Lambda_OII10;     

static double Energy_OIII32;     
static double Energy_OIII31;     
static double Energy_OIII21;     
static double Energy_OIII10;     
static double Energy_OII20;     
static double Energy_OII10;     

static double IntensityCMB_OIII32;
static double IntensityCMB_OIII31;
static double IntensityCMB_OIII21;
static double IntensityCMB_OIII10;
static double IntensityCMB_OII20;
static double IntensityCMB_OII10;

static double EinsteinAul_OIII32; 
static double EinsteinAul_OIII31; 
static double EinsteinAul_OIII21; 
static double EinsteinAul_OIII10; 
static double EinsteinAul_OII20; 
static double EinsteinAul_OII10; 

static double EinsteinBul_OIII32;
static double EinsteinBul_OIII31;
static double EinsteinBul_OIII21;
static double EinsteinBul_OIII10;
static double EinsteinBul_OII20;
static double EinsteinBul_OII10;

static double EinsteinBlu_OIII32;
static double EinsteinBlu_OIII31;
static double EinsteinBlu_OIII21;
static double EinsteinBlu_OIII10;
static double EinsteinBlu_OII20;
static double EinsteinBlu_OII10;


void InitializeCoolingLinesHIIregion(int Check){

    if(CEQLibRunParameters.UseHIIRegion==false) return;

    //Energy_CI610um        = TCI610um*BOLTZMANN_CONSTANT_CGS;     
    //IntensityCMB_CI610um  = MeanSpecificIntensityCMB(Energy_CI610um);
    //EinsteinA10_CI610um   = AulCI610um;
    //EinsteinB10_CI610um   = EinsteinCoefficientB10(EinsteinA10_CI610um,Energy_CI610um);         
    //EinsteinB01_CI610um   = EinsteinCoefficientB01(EinsteinB10_CI610um,gulCI610um);
    //if(Check==1)
    //    fprintf(stderr,"[CI]610um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
    //            EinsteinA10_CI610um,EinsteinB10_CI610um,EinsteinB01_CI610um,gulCI610um,Energy_CI610um,IntensityCMB_CI610um); 
    
    // [OIII]
    //  5008A  (1D2->3P2; n = 3 -> 2)
    //  4960A  (1D2->3P1; n = 3 -> 1)
    //  51.8um (3P2->3P1; n = 2 -> 1)
    //  88.4um (3P1->3P0; n = 1 -> 0)
    Energy_OIII32      = TOIII32*BOLTZMANN_CONSTANT_CGS; 
    Lambda_OIII32      = PLANK_H_CGS*LIGHT_C_CGS/Energy_OIII32;
    EinsteinAul_OIII32 = AulOIII32;
    EinsteinBul_OIII32 = EinsteinCoefficientB10(EinsteinAul_OIII32,Energy_OIII32);         
    EinsteinBlu_OIII32 = EinsteinCoefficientB01(EinsteinBul_OIII32,gulOIII32);
    
    Energy_OIII31      = TOIII31*BOLTZMANN_CONSTANT_CGS; 
    Lambda_OIII31      = PLANK_H_CGS*LIGHT_C_CGS/Energy_OIII31;
    EinsteinAul_OIII31 = AulOIII31;
    EinsteinBul_OIII31 = EinsteinCoefficientB10(EinsteinAul_OIII31,Energy_OIII31);         
    EinsteinBlu_OIII31 = EinsteinCoefficientB01(EinsteinBul_OIII31,gulOIII31);
    
    Energy_OIII21      = TOIII21*BOLTZMANN_CONSTANT_CGS; 
    Lambda_OIII21      = PLANK_H_CGS*LIGHT_C_CGS/Energy_OIII21;
    EinsteinAul_OIII21 = AulOIII21;
    EinsteinBul_OIII21 = EinsteinCoefficientB10(EinsteinAul_OIII21,Energy_OIII21);         
    EinsteinBlu_OIII21 = EinsteinCoefficientB01(EinsteinBul_OIII21,gulOIII21);
    
    Energy_OIII10      = TOIII10*BOLTZMANN_CONSTANT_CGS; 
    Lambda_OIII10      = PLANK_H_CGS*LIGHT_C_CGS/Energy_OIII10;
    EinsteinAul_OIII10 = AulOIII10;
    EinsteinBul_OIII10 = EinsteinCoefficientB10(EinsteinAul_OIII10,Energy_OIII10);         
    EinsteinBlu_OIII10 = EinsteinCoefficientB01(EinsteinBul_OIII10,gulOIII10);

    // [OII]
    //  3727A (2D3/2->4S3/2; n = 2 -> 0)     1.8e-4     4.49e+3
    //  3730A (2D5/2->4S3/2; n = 1 -> 0)     3.6e-5     3.31e+3
    Energy_OII20      = TOII20*BOLTZMANN_CONSTANT_CGS; 
    Lambda_OII20      = PLANK_H_CGS*LIGHT_C_CGS/Energy_OII20;
    EinsteinAul_OII20 = AulOII20;
    EinsteinBul_OII20 = EinsteinCoefficientB10(EinsteinAul_OII20,Energy_OII20);         
    EinsteinBlu_OII20 = EinsteinCoefficientB01(EinsteinBul_OII20,gulOII20);
    
    Energy_OII10      = TOII10*BOLTZMANN_CONSTANT_CGS; 
    Lambda_OII10      = PLANK_H_CGS*LIGHT_C_CGS/Energy_OII10;
    EinsteinAul_OII10 = AulOII10;
    EinsteinBul_OII10 = EinsteinCoefficientB10(EinsteinAul_OII10,Energy_OII10);         
    EinsteinBlu_OII10 = EinsteinCoefficientB01(EinsteinBul_OII10,gulOII10);

    if(Check==1){
        //fprintf(stderr,"OII : E(%.1fA)=%feV\n",
        //        LambdaOII[0]/ANGSTROM_CGS,E_OII[0]*ERG_TO_EV);
        //fprintf(stderr,"OII : E(%.1fA)=%feV\n",
        //        LambdaOII[1]/ANGSTROM_CGS,E_OII[1]*ERG_TO_EV);
        //fprintf(stderr,"OIII : E(%.1fum)=%feV\n",
        //        LambdaOIII[0]/MICRON_CGS,E_OIII[0]*ERG_TO_EV);
        //fprintf(stderr,"OIII : E(%.1fum)=%feV\n",
        //        LambdaOIII[1]/MICRON_CGS,E_OIII[1]*ERG_TO_EV);
        //fprintf(stderr,"OIII : E(%.1fA)=%feV\n",
        //        LambdaOIII[2]/ANGSTROM_CGS,E_OIII[2]*ERG_TO_EV);
        //fprintf(stderr,"OIII : E(%.1fA)=%feV\n",
        //        LambdaOIII[3]/ANGSTROM_CGS,E_OIII[3]*ERG_TO_EV);
    }

    return ;
}

double CoolingLinesHIIregion(double T, double ne, double nX, int ipX){
    fprintf(stderr," #%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
    exit(1);
    return 0.0;
}

double CoolingLinesHIIregionOxygen(double T, double ne, double nOII, double nOIII){
    fprintf(stderr," #%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
    exit(1);
    return 0.0;
#if 0
   // eq.(27.25) in Draine's textbook
   // eq.(3.21)(3.26) in Osterbrock (1989)
   double T4   = 1.e-4*T;
   double lnT4 = log(T4);
   double invT = 1.e0/T;
   double invkT = invT/(BOLTZMANN_CONSTANT_CGS);
   double L_OII = 0.e0;
   double gammOII12;
   double q21, q12;
   for(int i=0; i<2; i++){
       // electreon collision strength for OII
       // ref. Tayal (2007), p.499 of Draine's textbook
       gammOII12 = aOII[i]*pow(T4,bOII[i]-cOII[i]*lnT4);

       // collision rate coefficients [cm^3/s]
       q21 = 8.629e-6*sqrt(invT)*gammOII12/WeightOII2[i];
       q12 = (WeightOII2[i]/WeightOII1[i])*q21*exp(-E_OII[i]*invkT);
       
       // cooling rate [erg*cm^3/s]
       L_OII += E_OII[i]*Aul_OII[i]*q12/(Aul_OII[i]+ne*q21);
   }
   // [erg/cc/s]
   return ne*nOII*(1.e0+OIII_OII_RATIO)*L_OII;
#endif
}

// 4*pi* line-integrated emissivity [erg/cc/s]
// ref. eq(4)(5) of Peters et al. MNRAS, 466, 3293, (2017)
static double G13(double T){
    //if(T<4.3e+3)    return 0.0;
    if(T<=2.5e+4)
        return 0.350 + T*(-2.62e-7 + T*(-8.15e-11 + 6.19e-15*T));
    else if(T<5.e+5)
        return 0.276 + T*(4.99e-6 - T*(8.85e-12+ 7.18e-18*T)); 
    else
        return 0.e0;
}
double CoolingHalphaEmission(StructCEQData *PhydroPP){
    double T    = PhydroPP->T;
    double T4   = 1.e-4*T;
    double lnT4 = log(T4);
    double ne   = PhydroPP->nel;
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nHII = PhydroPP->nion[ipHYDROGEN];
    // Two processes contribute to the emis- sion of Hα radiation. 
    // The first one is the recombination of ionised hydrogen with a free electron. 
    double j_a_R = 4.0*PI*2.82e-26*ne*nHII*pow(T4,-0.942-0.031*lnT4); 

    // The second contribution comes from collisional excitation of neutral hydrogen by free electrons. 
    double kT    = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T;
    double j_a_C = 1.30e-17*ne*nHI*exp(-12.1/kT)*G13(T)/sqrt(T); 

    return j_a_R + j_a_C;   // [erg/cc/s]
}


