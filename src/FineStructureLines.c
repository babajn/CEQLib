#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "FineStructureLines.h"
#include "StatisticalEquilibrium.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

//#define FINE_STRUCUTRE_COOLING_CI
#define FINE_STRUCUTRE_COOLING_OI
#define FINE_STRUCUTRE_COOLING_CII
//#define FINE_STRUCUTRE_COOLING_SiI
#define FINE_STRUCUTRE_COOLING_SiII

#define OrthH2Frac  (0.75)
#define ParaH2Frac  (0.25) 

static const double Zsol = SOLAR_METALLICITY;
//static void InitializeCIIH2CollisionCoolingCoefficient(void);

//////// De-excitation Rates [cm^3/s]
// ref.) Table 6 of Glover & Jappsen (2007)
static double q10_CI_oH2(double T);
static double q20_CI_oH2(double T);
static double q21_CI_oH2(double T);
static double q10_CI_pH2(double T);
static double q20_CI_pH2(double T);
static double q21_CI_pH2(double T);
static double q10_CI_HI(double T);
static double q20_CI_HI(double T);
static double q21_CI_HI(double T);
static double q10_CI_HII(double T);
static double q20_CI_HII(double T);
static double q21_CI_HII(double T);
static double q10_CI_e(double T);
static double q20_CI_e(double T);
static double q21_CI_e(double T);

static double q10_CII_oH2(double T);
static double q10_CII_pH2(double T);
static double q10_CII_HI(double T);
static double q10_CII_e(double T);

static double q10_OI_oH2(double T);
static double q20_OI_oH2(double T);
static double q21_OI_oH2(double T);
static double q10_OI_pH2(double T);
static double q20_OI_pH2(double T);
static double q21_OI_pH2(double T);
static double q10_OI_HI(double T);
static double q20_OI_HI(double T);
static double q21_OI_HI(double T);
static double q10_OI_HII(double T);
static double q20_OI_HII(double T);
static double q21_OI_HII(double T);
static double q10_OI_e(double T);
static double q20_OI_e(double T);
static double q21_OI_e(double T);

static double q10_SiI_HI(double T);
static double q20_SiI_HI(double T);
static double q21_SiI_HI(double T);
static double q10_SiI_HII(double T);
static double q20_SiI_HII(double T);
static double q21_SiI_HII(double T);

static double q10_SiII_HI(double T);
static double q10_SiII_e(double T);

// Spontaneous Emission Rate [1/sec]
#define AulCI610um   (7.9e-8)    // [1/s]
#define AulCI230um   (2.1e-14)   // [1/s]
#define AulCI370um   (2.7e-7)    // [1/s]
#define AulCII158um  (2.3e-6)    // [1/s]
#define AulOI063um   (8.9e-5)    // [1/s]
#define AulOI044um   (1.3e-10)   // [1/s]
#define AulOI146um   (1.8e-5)    // [1/s]
#define AulSiI139um  (8.4e-6)
#define AulSiI045um  (2.4e-10)
#define AulSiI068um  (4.2e-5)
#define AulSiII035um (2.2e-4)    // [1/s]

// Excitation Temperature [K]
#define TCI610um     ( 24.0)     // [K]
#define TCI230um     ( 63.0)     // [K]
#define TCI370um     ( 39.0)     // [K]
#define TCII158um    ( 92.0)     // [K]
#define TOI063um     (230.0)     // [K]
#define TOI044um     (330.0)     // [K]
#define TOI146um     ( 98.0)     // [K]
#define TSiI130um    (110.0)     // [K]
#define TSiI045um    (320.0)     // [K]
#define TSiI068um    (210.0)     // [K]
#define TSiII035um   (410.0)     // [K]

// Statistical Weight Ratio: g1/g0 (upper/lower)
#define gulCI610um     (3.0000)
#define gulCI230um     (5.0000)
#define gulCI370um     (1.6667)
#define gulCII158um    (2.0000)    // g(^2P3/2)g(^2P1/2) 
#define gulOI063um     (0.6000)
#define gulOI044um     (0.2000)    
#define gulOI146um     (0.3333)
#define gulSiI130um    (3.0000)
#define gulSiI045um    (5.0000)
#define gulSiI068um    (1.6667)
#define gulSiII035um   (2.0000)

//// Statistical Weight Ratio: g0/g1 (lower/upper)
//#define gluCI610um     (0.3333)    // g(^3P1)/g(^3P0)
//#define gluCI230um     (0.2000)
//#define gluCI370um     (0.6000)    // g(^3P2)/g(^3P1)
//#define gluCII158um    (0.5000)    // g(^2P1/2)/g(^2P3/2)
//#define gluOI063um     (1.6667)    // g(^3Ps)/g(^3P1)
//#define gluOI044um     (5.0000)    
//#define gluOI146um     (3.0000)    // g(^3P1)/g(^3P0)
//#define gluSiI130um    (0.3333)
//#define gluSiI045um    (0.2000)
//#define gluSiI068um    (0.6000)
//#define gluSiII035um   (0.5000)    // g(^2P1/2)/g(^2P3/2)


static double Energy_CI610um;     
static double IntensityCMB_CI610um;
static double EinsteinA10_CI610um; 
static double EinsteinB10_CI610um;
static double EinsteinB01_CI610um;

static double Energy_CI230um;     
static double IntensityCMB_CI230um;
static double EinsteinA10_CI230um; 
static double EinsteinB10_CI230um;
static double EinsteinB01_CI230um;

static double Energy_CI370um;     
static double IntensityCMB_CI370um;
static double EinsteinA10_CI370um; 
static double EinsteinB10_CI370um;
static double EinsteinB01_CI370um;

static double Energy_CII158um;     
static double IntensityCMB_CII158um;
static double EinsteinA10_CII158um; 
static double EinsteinB10_CII158um;
static double EinsteinB01_CII158um;

static double Energy_OI063um;     
static double IntensityCMB_OI063um;
static double EinsteinA10_OI063um; 
static double EinsteinB10_OI063um;
static double EinsteinB01_OI063um;

static double Energy_OI044um;     
static double IntensityCMB_OI044um;
static double EinsteinA10_OI044um; 
static double EinsteinB10_OI044um;
static double EinsteinB01_OI044um;

static double Energy_OI146um;     
static double IntensityCMB_OI146um;
static double EinsteinA10_OI146um; 
static double EinsteinB10_OI146um;
static double EinsteinB01_OI146um;

static double Energy_SiII035um;     
static double IntensityCMB_SiII035um;
static double EinsteinA10_SiII035um; 
static double EinsteinB10_SiII035um;
static double EinsteinB01_SiII035um;


static double CoolingCIDetail(double nCI, double T, double nel, double nHII, double nHI, double nH2, 
        double *nCI0, double *nCI1, double *nCI2, double *CoolingCI610um, double *CoolingCI230um, double *CoolingCI370um);
static double CoolingCII158umDetail(double nCII, double T, double nel, double nHI, double nH2, double *nCII0, double *nCII1, double *ncrit);
static double CoolingOIDetail(double nOI, double T, double nel, double nHII, double nHI, double nH2,
        double *nOI0, double *nOI1, double *nOI2, double *CoolingOI063um, double *CoolingOI044um, double *CoolingOI146um);
static double CoolingSiII35umDetail(double nSiII, double T, double nel, double nHI, double *nSiII0, double *nSiII1, double *ncrit);

///// Einstein Coefficient (Stimulated Emission)
//static double EinsteinCoefficientB10(double A10, double E10){
//    double nyu10 = E10/PLANK_H_CGS;  // [Hz]
//    return A10*CUBE(LIGHT_C_CGS/nyu10)/(8.0*PI*PLANK_H_CGS); 
//
//}
///// Einstein Coefficient (Absorption)
//static double EinsteinCoefficientB01(double B10, double g1g0){
//    return g1g0*B10;
//}

static double CriticalDensityFSLines(double A10, double B10, double g1g0, double I10, double n1){
    return (A10 + B10*I10)*n1/(g1g0*B10*I10);
}

/// [erg/s/cm2/Hz/str]
static double SpecificIntensityBlackBody(double T, double nyu){
    double kT = BOLTZMANN_CONSTANT_CGS*T; // [erg]
    double fact = 2.0*PLANK_H_CGS*nyu*SQ(nyu/LIGHT_C_CGS); 
    return fact/(exp(PLANK_H_CGS*nyu/kT)-1.0);
}
static double PhotonOccupationNumberBlackBody(double T, double nyu){
    double kT = BOLTZMANN_CONSTANT_CGS*T; // [erg]
    return 1.e0/(exp(PLANK_H_CGS*nyu/kT)-1.0);
}

/// [erg/s/cm2/Hz/str]
static double MeanSpecificIntensityCMB(double E10){
    double nyu10 = E10/PLANK_H_CGS;  // [Hz]
    return SpecificIntensityBlackBody(CMB_TEMPERATURE,nyu10);
}
static double PhotonOccupationNumberCMB(double E10){
    double nyu10 = E10/PLANK_H_CGS;  // [Hz]
    return PhotonOccupationNumberBlackBody(CMB_TEMPERATURE,nyu10);
}

#if 0
// ref) eq(56) of Glover & Jappsen (2007)
static void SolveTwoLevelSystemStatisticalEquilibrium(
        double E10, double g1g0,
        double I10,double A10, double B10, double B01, 
        int NcolSpecies, double ncol[], double q10[], 
        double n, double kT, double *n0, double *n1){

    double I01 = I10;
    double C10 = 0.e0;
    for(int k=0;k<NcolSpecies;k++)   C10 += q10[k]*ncol[k];
    double C01 = C10*g1g0*exp(-E10/kT);
    double e1 = B01*I01 + C01;
    double e2 = A10 + B10*I10 + C10;

    if(e1+e2==0.0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"Error: Zero e1+e2\n"); 
        fprintf(stderr," e1=%e e2=%e\n",e1,e2); 
        fprintf(stderr," kT=%e \n"); 
        exit(1);
    }
    *n1 = n*e1/(e1 + e2);
    *n0 = n - *n1;
    //*n0 = n*e2/(e1 + e2);
    //*n1 = n - *n0;

    //double f10 = (B01*I01 + C01)/(A10 + B10*I10 + C10);
    //*n1 = n*f10/(1.0+f10);
    //*n0 = n - *n1;

    return ;
}

// ref) eqs(61)(62) of Glover & Jappsen (2007)
static void SolveThreeLevelSystemStatisticalEquilibrium(
        double E10, double E20, double E21,
        double g1g0,double g2g0,double g2g1,
        double I10, double I20, double I21,
        double A10, double A20, double A21,
        double B10, double B20, double B21,
        double B01, double B02, double B12,
        int NcolSpecies, double ncol[], double q10[], double q20[], double q21[],
        double n, double kT, double *n0, double *n1, double *n2){
    
    double I01 = I10, I02 = I20, I12 = I21;
    double C10 = 0.e0, C20 = 0.e0, C21 = 0.e0;
    for(int k=0;k<NcolSpecies;k++){
        C10 += q10[k]*ncol[k];
        C20 += q20[k]*ncol[k];
        C21 += q21[k]*ncol[k];
    }
    double C01 = C10*g1g0*exp(-E10/kT);
    double C02 = C20*g2g0*exp(-E20/kT);
    double C12 = C21*g2g1*exp(-E21/kT);

    double e1 = B01*I01 + C01 + B02*I02 + C02;
    double e2 = A10 + B10*I10 + C10;
    double e3 = A20 + B20*I20 + C20;
    double e4 = B10*I10 + C10 + B12*I12 + C12;
    double e5 = B01*I01 + C01;
    double e6 = A21 + B21*I21 + C21;

    //*n0 = n*(e3*e4 + e2*e6)/(e1*e4 + e1*e6 + e3*e4 - e2*e5 + e3*e5 + e2*e6);
    //*n1 = n*(e3*e5 + e1*e6)/(e1*e4 + e1*e6 + e3*e4 - e2*e5 + e3*e5 + e2*e6);
    //*n2 = n - *n0 - *n1;

    double f10 = (e1*e6+e3*e5)/(e2*e6+e3*e4);
    double f20 = (e1*e4-e2*e5)/(e3*e4+e2*e6);
    *n0 = n*1.0/(1.0+f10+f20);
    *n1 = n*f10/(1.0+f10+f20);
    *n2 = n - *n0 - *n1; 

    return ;
}
#endif

void InitializeFineStructureLines(int Check){

    double n1 = 1.0;
    //Pall.TCMB = CMB_TEMPERATURE;

#ifdef FINE_STRUCUTRE_COOLING_CI
    Energy_CI610um        = TCI610um*BOLTZMANN_CONSTANT_CGS;     
    IntensityCMB_CI610um  = MeanSpecificIntensityCMB(Energy_CI610um);
    EinsteinA10_CI610um   = AulCI610um;
    EinsteinB10_CI610um   = EinsteinCoefficientB10(EinsteinA10_CI610um,Energy_CI610um);         
    EinsteinB01_CI610um   = EinsteinCoefficientB01(EinsteinB10_CI610um,gulCI610um);
    if(Check==1)
        fprintf(stderr,"[CI]610um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
                EinsteinA10_CI610um,EinsteinB10_CI610um,EinsteinB01_CI610um,gulCI610um,Energy_CI610um,IntensityCMB_CI610um); 

    Energy_CI230um        = TCI230um*BOLTZMANN_CONSTANT_CGS;     
    IntensityCMB_CI230um  = MeanSpecificIntensityCMB(Energy_CI230um);
    EinsteinA10_CI230um   = AulCI230um;
    EinsteinB10_CI230um   = EinsteinCoefficientB10(EinsteinA10_CI230um,Energy_CI230um);         
    EinsteinB01_CI230um   = EinsteinCoefficientB01(EinsteinB10_CI230um,gulCI230um);
    if(Check==1)
        fprintf(stderr,"[CI]230um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
                EinsteinA10_CI230um,EinsteinB10_CI230um,EinsteinB01_CI230um,gulCI230um,Energy_CI230um,IntensityCMB_CI230um); 

    Energy_CI370um        = TCI370um*BOLTZMANN_CONSTANT_CGS;     
    IntensityCMB_CI370um  = MeanSpecificIntensityCMB(Energy_CI370um);
    EinsteinA10_CI370um   = AulCI370um;
    EinsteinB10_CI370um   = EinsteinCoefficientB10(EinsteinA10_CI370um,Energy_CI370um);         
    EinsteinB01_CI370um   = EinsteinCoefficientB01(EinsteinB10_CI370um,gulCI370um);
    if(Check==1)
        fprintf(stderr,"[CI]370um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
                EinsteinA10_CI370um,EinsteinB10_CI370um,EinsteinB01_CI370um,gulCI370um,Energy_CI370um,IntensityCMB_CI370um); 
#endif

    Energy_CII158um       = TCII158um*BOLTZMANN_CONSTANT_CGS;    // [erg]
    IntensityCMB_CII158um = MeanSpecificIntensityCMB(Energy_CII158um);
    EinsteinA10_CII158um  = AulCII158um;
    EinsteinB10_CII158um  = EinsteinCoefficientB10(EinsteinA10_CII158um,Energy_CII158um);         
    EinsteinB01_CII158um  = EinsteinCoefficientB01(EinsteinB10_CII158um,gulCII158um);
    if(Check==1)
        fprintf(stderr,"[CII]158um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e n_CMB = %e n_crit = %e\n",
                EinsteinA10_CII158um,EinsteinB10_CII158um,EinsteinB01_CII158um,gulCII158um,Energy_CII158um,
                IntensityCMB_CII158um,PhotonOccupationNumberCMB(Energy_CII158um),
                CriticalDensityFSLines(EinsteinA10_CII158um,EinsteinB10_CII158um,gulCII158um,IntensityCMB_CII158um,n1));
    
    Energy_OI063um        = TOI063um*BOLTZMANN_CONSTANT_CGS;    // [erg]
    IntensityCMB_OI063um  = MeanSpecificIntensityCMB(Energy_OI063um);
    EinsteinA10_OI063um   = AulOI063um;
    EinsteinB10_OI063um   = EinsteinCoefficientB10(EinsteinA10_OI063um,Energy_OI063um);         
    EinsteinB01_OI063um   = EinsteinCoefficientB01(EinsteinB10_OI063um,gulOI063um);
    if(Check==1)
        fprintf(stderr,"[OI]63um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
                EinsteinA10_OI063um,EinsteinB10_OI063um,EinsteinB01_OI063um,gulOI063um,Energy_OI063um,IntensityCMB_OI063um); 

    Energy_OI044um        = TOI044um*BOLTZMANN_CONSTANT_CGS;    // [erg]
    IntensityCMB_OI044um  = MeanSpecificIntensityCMB(Energy_OI044um);
    EinsteinA10_OI044um   = AulOI044um;
    EinsteinB10_OI044um   = EinsteinCoefficientB10(EinsteinA10_OI044um,Energy_OI044um);         
    EinsteinB01_OI044um   = EinsteinCoefficientB01(EinsteinB10_OI044um,gulOI044um);
    if(Check==1)
        fprintf(stderr,"[OI]44um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
                EinsteinA10_OI044um,EinsteinB10_OI044um,EinsteinB01_OI044um,gulOI044um,Energy_OI044um,IntensityCMB_OI044um); 

    Energy_OI146um        = TOI146um*BOLTZMANN_CONSTANT_CGS;    // [erg]
    IntensityCMB_OI146um  = MeanSpecificIntensityCMB(Energy_OI146um);
    EinsteinA10_OI146um   = AulOI146um;
    EinsteinB10_OI146um   = EinsteinCoefficientB10(EinsteinA10_OI146um,Energy_OI146um);         
    EinsteinB01_OI146um   = EinsteinCoefficientB01(EinsteinB10_OI146um,gulOI146um);
    if(Check==1)
        fprintf(stderr,"[OI]146um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e\n",
                EinsteinA10_OI146um,EinsteinB10_OI146um,EinsteinB01_OI146um,gulOI146um,Energy_OI146um,IntensityCMB_OI146um); 
   
    Energy_SiII035um       = TSiII035um*BOLTZMANN_CONSTANT_CGS;    // [erg]
    IntensityCMB_SiII035um = MeanSpecificIntensityCMB(Energy_SiII035um);
    EinsteinA10_SiII035um  = AulSiII035um;
    EinsteinB10_SiII035um  = EinsteinCoefficientB10(EinsteinA10_SiII035um,Energy_SiII035um);         
    EinsteinB01_SiII035um  = EinsteinCoefficientB01(EinsteinB10_SiII035um,gulSiII035um);
    if(Check==1)
        fprintf(stderr,"[SiII]035um: (A10,B10,B01)=(%e,%e,%e) g1/g0 = %f E10=%e I10=%e n_CMB = %e n_crit = %e\n",
                EinsteinA10_SiII035um,EinsteinB10_SiII035um,EinsteinB01_SiII035um,gulSiII035um,Energy_SiII035um,
                IntensityCMB_SiII035um,PhotonOccupationNumberCMB(Energy_SiII035um),
                CriticalDensityFSLines(EinsteinA10_SiII035um,EinsteinB10_SiII035um,gulSiII035um,IntensityCMB_SiII035um,n1));


    if(Check==1){
        int nRho = 80; 
        int nTmp =  9;
        // cf. Fig. 17.4, 17.7 of Drain's text  
        double CalcTmp[] = {20,30,50,100,200,500,1000,2000,5000};
        double xe   = 1.e-3;
        double nH2  = 0.e0, nHII = 0.e0, nHI;
        double nCII = 1.0, nCII0, nCII1, nCIIcrit;
        double nOI = 1.0, nOI0, nOI1, nOI2, nOIcrit;
        double nSiII = 1.0, nSiII0, nSiII1, nSiIIcrit;
        double CoolingOI063um, CoolingOI044um, CoolingOI146um;
        FILE *outCII  = fopen(CEQ_DIAGNOSTIC_DIRECTORY "/FineStructureLineTestCII.dat","w");
        FILE *outOI   = fopen(CEQ_DIAGNOSTIC_DIRECTORY "/FineStructureLineTestOI.dat","w");
        FILE *outSiII = fopen(CEQ_DIAGNOSTIC_DIRECTORY "/FineStructureLineTestSiII.dat","w");
        for(int i=0;i<nRho;i++){
            double nH  = pow(10.0,-1+i*0.1);
            double nel = xe*nH;
            nHII = nel;
            nHI  = nH - nHII;
            fprintf(outCII,"%e ",nHI);
            fprintf(outOI ,"%e ",nHI);
            fprintf(outSiII,"%e ",nHI);
            for(int it=0;it<nTmp;it++){
                double T      = CalcTmp[it]; //20.0*(1+it);
                double L_CII  = CoolingCII158umDetail(nCII,T,nel,nHI,nH2,&nCII0,&nCII1,&nCIIcrit);
                double L_OI   = CoolingOIDetail(nOI,T,nel,nHII,nHI,nH2,&nOI0,&nOI1,&nOI2,&CoolingOI063um,&CoolingOI044um,&CoolingOI146um);
                double L_SiII = CoolingSiII35umDetail(nSiII,T,nel,nHI,&nSiII0,&nSiII1,&nSiIIcrit);
                //fprintf(outCII ,"%e ",L_CII );
                //fprintf(outOI  ,"%e ",L_OI  );
                //fprintf(outSiII,"%e ",L_SiII);
                fprintf(outCII ,"%e ",nCII1/nCII0);
                fprintf(outOI  ,"%e ",nOI1/nOI0);
                fprintf(outSiII,"%e ",nSiII1/nSiII0);
                //fprintf(outOI ,"%e %e ",nOI1/nOI0,nOI2/nOI1);
                if(i==0)    fprintf(stderr," T=%f nCIIcrit=%e nOIcrit=%e nSiIIcrit=%e\n",T,nCIIcrit,nOIcrit,nSiIIcrit);
            }
            fprintf(outCII,"\n");
            fprintf(outOI ,"\n");
            fprintf(outSiII,"\n");
        }
        fclose(outCII);
        fclose(outOI);
        fclose(outSiII);
       
        /*
#define CarbonAbundance     (pow(10.0,-3.52))
#define OxygenAbundance     (pow(10.0,-3.34))
#define SiliconAbundance    (pow(10.0,-5.45))
        FILE *outCool = fopen("out/FineStructureLineCoolingTest.dat","w");
        //double nH = 0.6, nel = 1.7e-2*nH;
        double nH = 20.0, nel = 4.0e-4*nH;
        nHII  = nel;
        nHI   = nH - nHII;
        nH2   = 0.0;
        nCII  = (CarbonAbundance )*nH;
        nOI   = (OxygenAbundance )*nH * nHI/nH;
        nSiII = (SiliconAbundance)*nH;
        for(int it=0;it<300;it++){
            double T = pow(10.0,1.0+it*0.01);
            double L_CII  = CoolingCII158um(nCII,T,nel,nHI,nH2);
            double L_OI   = CoolingOI63um(nOI,T,nel,nHII,nHI,nH2);
            double L_SiII = CoolingSiII35um(nSiII,T,nel,nHI);
            fprintf(outCool,"%e ",T            );
            fprintf(outCool,"%e ",L_CII /SQ(nH));
            fprintf(outCool,"%e ",L_OI  /SQ(nH));
            fprintf(outCool,"%e ",L_SiII/SQ(nH));
            fprintf(outCool,"\n");
        }
        fclose(outCool);
        */
    }

    return ;
}

// total cooling rate of fine structure lines
//double TotalFineStructureLinesCoolingRate(double T, double nel, double nHII, double nHI, double nH2, 
//        double nCI, double nCII, double nOI, double nSiI, double nSiII){
//    
//    double L_CI610 = 0.e0, L_CI230 = 0.e0, L_CI370 = 0.e0;
//    double L_OI063 = 0.e0, L_OI044 = 0.e0, L_OI146 = 0.e0;
//    double L_SiI130 = 0.e0, L_SiI045 = 0.e0, L_SiI068 = 0.e0;
//
//    double L_CI     = CoolingCI(nCI,T,nel,nHII,nHI,nH2,&L_CI610,&L_CI230,&L_CI370);
//    double L_CII158 = CoolingCII158um(nCII,T,nel,nHI,nH2);
//    double L_OI     = CoolingOI(nOI,T,nel,nHII,nHI,nH2,&L_OI063,&L_OI044,&L_OI146);
//    double L_SiI    = CoolingSiI(nSiI,T,nHI,nHII,&L_SiI130,&L_SiI045,&L_SiI068);
//    double L_SiII35 = CoolingSiII35um(nSiII,T,nel,nHI);
//    if(isnan(L_CI)){
//        fprintf(stderr,"L_CI=%e nCI=%e T=%e nel=%e nHII=%e nHI=%e nH2=%e\n",
//                L_CI,nCI,T,nel,nHII,nHI,nH2);
//        exit(1);
//    }
//    if(isnan(L_CII158)){
//        fprintf(stderr,"L_CII158=%e\n",L_CII158);
//        exit(1);
//    }
//    if(isnan(L_OI)){
//        fprintf(stderr,"L_OI=%e\n",L_OI);
//        exit(1);
//    }
//    if(isnan(L_SiI)){
//        fprintf(stderr,"L_SiI=%e\n",L_SiI);
//        exit(1);
//    }
//    if(isnan(L_SiII35)){
//        fprintf(stderr,"L_SiII35=%e\n",L_SiII35);
//        exit(1);
//    }
//    return L_CI + L_CII158 + L_OI + L_SiI + L_SiII35;
//}

// Fine structure transition between CI
double CoolingCI(double nCI, double T, double nel, double nHII, double nHI, double nH2, 
        double *CoolingCI610um, double *CoolingCI230um, double *CoolingCI370um){
#ifdef FINE_STRUCUTRE_COOLING_CI
    double n0, n1, n2; //, ncrit;
    double L_CI610um, L_CI230um, L_CI370um;
    double L_CI = CoolingCIDetail(nCI,T,nel,nHII,nHI,nH2,&n0,&n1,&n2,&L_CI610um,&L_CI230um,&L_CI370um);
    *CoolingCI610um = L_CI610um;
    *CoolingCI230um = L_CI230um;
    *CoolingCI370um = L_CI370um;
    return L_CI;
#else
    *CoolingCI610um = 0.0;
    *CoolingCI230um = 0.0;
    *CoolingCI370um = 0.0;
    return 0.e0;
#endif
}
static double CoolingCIDetail(double nCI, double T, double nel, double nHII, double nHI, double nH2, 
        double *nCI0, double *nCI1, double *nCI2, double *CoolingCI610um, double *CoolingCI230um, double *CoolingCI370um){

    double E10 = Energy_CI610um;     
    double E20 = Energy_CI230um;     
    double E21 = Energy_CI370um;     

    double I10 = IntensityCMB_CI610um;
    double I20 = IntensityCMB_CI230um;
    double I21 = IntensityCMB_CI370um;

    double A10 = EinsteinA10_CI610um; 
    double A20 = EinsteinA10_CI230um; 
    double A21 = EinsteinA10_CI370um; 

    double B10 = EinsteinB10_CI610um;
    double B20 = EinsteinB10_CI230um;
    double B21 = EinsteinB10_CI370um;

    double B01 = EinsteinB01_CI610um;
    double B02 = EinsteinB01_CI230um;
    double B12 = EinsteinB01_CI370um;
    
    double g1g0 = gulCI610um; 
    double g2g0 = gulCI230um;
    double g2g1 = gulCI370um;

    int NcolSpecies = 5;
    double ncol[] = {OrthH2Frac*nH2,ParaH2Frac*nH2,nHI          ,nHII          ,nel        };
    double q10[]  = {q10_CI_oH2(T) ,q10_CI_pH2(T) ,q10_CI_HI(T) ,q10_CI_HII(T) ,q10_CI_e(T)};
    double q20[]  = {q20_CI_oH2(T) ,q20_CI_pH2(T) ,q20_CI_HI(T) ,q20_CI_HII(T) ,q20_CI_e(T)};
    double q21[]  = {q21_CI_oH2(T) ,q21_CI_pH2(T) ,q21_CI_HI(T) ,q21_CI_HII(T) ,q21_CI_e(T)};

    double kT = BOLTZMANN_CONSTANT_CGS*T;
    double n2 = 0.e0;    // uppest level [1/cc]
    double n1 = 0.5*nCI; // middle level [1/cc]
    double n0 = 0.5*nCI; // lowest level [1/cc]
    SolveThreeLevelSystemStatisticalEquilibrium(
            E10,E20,E21,g1g0,g2g0,g2g1,I10,I20,I21,
            A10,A20,A21,B10,B20,B21,B01,B02,B12,
            NcolSpecies,ncol,q10,q20,q21,
            nCI,kT,&n0,&n1,&n2);
    
    *nCI0 = n0;
    *nCI1 = n1;
    *nCI2 = n2;

    /// spontaneous + stimulated emissions [erg/s/cm^3]
    double lamb610um = (A10 + B10*I10)*E10*n1;
    double lamb230um = (A20 + B20*I20)*E20*n2;
    double lamb370um = (A21 + B21*I21)*E21*n2;
    /// absorption [erg/s/cm^3]
    double gamm610um = B01*I10*E10*n0;
    double gamm230um = B02*I20*E20*n0;
    double gamm370um = B12*I21*E21*n1;
    /// net energy loss [erg/s/cm^3]
    double net610um  = lamb610um - gamm610um;
    double net230um  = lamb230um - gamm230um;
    double net370um  = lamb370um - gamm370um;

    *CoolingCI610um = net610um;
    *CoolingCI230um = net230um;
    *CoolingCI370um = net370um;
    return net610um + net230um + net370um;
}

double CoolingCI610um(double nCI, double T, double nel, double nHII, double nHI, double nH2){
    double CoolingCI610um, CoolingCI230um, CoolingCI370um;
    CoolingCI(nCI,T,nel,nHII,nHI,nH2,&CoolingCI610um,&CoolingCI230um,&CoolingCI370um);
    return CoolingCI610um;
}

double CoolingCI230um(double nCI, double T, double nel, double nHII, double nHI, double nH2){
    double CoolingCI610um, CoolingCI230um, CoolingCI370um;
    CoolingCI(nCI,T,nel,nHII,nHI,nH2,&CoolingCI610um,&CoolingCI230um,&CoolingCI370um);
    return CoolingCI230um;
}

double CoolingCI370um(double nCI, double T, double nel, double nHII, double nHI, double nH2){
    double CoolingCI610um, CoolingCI230um, CoolingCI370um;
    CoolingCI(nCI,T,nel,nHII,nHI,nH2,&CoolingCI610um,&CoolingCI230um,&CoolingCI370um);
    return CoolingCI370um;
}

// Fine structure transition between CII(^2P1/2-^2P3/2)
double CoolingCII158um(StructCEQData *PhydroPP){
    double nCII = PhydroPP->nion [ipCARBON  ];
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    double nel  = PhydroPP->nel; 
    double T    = PhydroPP->T;
#ifdef FINE_STRUCUTRE_COOLING_CII
    double n0, n1, ncrit;
    return CoolingCII158umDetail(nCII,T,nel,nHI,nH2,&n0,&n1,&ncrit);
#else
    return 0.e0;
#endif
}
static double CoolingCII158umDetail(double nCII, double T, double nel, double nHI, double nH2, double *nCII0, double *nCII1, double *ncrit){
    double I10   = IntensityCMB_CII158um;
    double E10   = Energy_CII158um;
    double A10   = EinsteinA10_CII158um;
    double B10   = EinsteinB10_CII158um;         
    double B01   = EinsteinB01_CII158um;
    double g1g0  = gulCII158um; 

    int NcolSpecies = 4;
    double ncol[] = {OrthH2Frac*nH2,ParaH2Frac*nH2,nHI          ,nel         };
    double q10[]  = {q10_CII_oH2(T),q10_CII_pH2(T),q10_CII_HI(T),q10_CII_e(T)};

    double kT = BOLTZMANN_CONSTANT_CGS*T;
    double n1 = 0.5*nCII;    // upper level [1/cc]
    double n0 = 0.5*nCII;    // lower level [1/cc]
    SolveTwoLevelSystemStatisticalEquilibrium(
        E10,g1g0,
        I10,A10,B10,B01,
        NcolSpecies,ncol,q10, 
        nCII,kT,&n0,&n1);

    double C10 = 0.e0;
    for(int k=0;k<NcolSpecies;k++)   C10 += q10[k]*ncol[k];
    *nCII0 = n0;
    *nCII1 = n1;
    *ncrit = A10/C10;

    /// spontaneous + stimulated emissions [erg/s/cm^3]
    double lamb = (A10 + B10*I10)*E10*n1;
    /// absorption [erg/s/cm^3]
    double gamm = B01*I10*E10*n0;
    /// net energy loss [erg/s/cm^3]
    // if n0 > (A10+B10*I10)*n1/{(g1/g0)*B10*I10}, 
    // then the ions will absorb more energy than they emit (= heating)
    return lamb - gamm;
}

// Fine structure transition between ^3P2-^3P1
double CoolingOI(double nOI, double T, double nel, double nHII, double nHI, double nH2, 
        double *CoolingOI063um, double *CoolingOI044um, double *CoolingOI146um){
#ifdef FINE_STRUCUTRE_COOLING_OI
    double n0, n1, n2; //, ncrit;
    double L_OI063um, L_OI044um, L_OI146um;
    double L_OI = CoolingOIDetail(nOI,T,nel,nHII,nHI,nH2,&n0,&n1,&n2,&L_OI063um,&L_OI044um,&L_OI146um);
    *CoolingOI063um = L_OI063um;
    *CoolingOI044um = L_OI044um;
    *CoolingOI146um = L_OI146um;
    return L_OI;
#else
    *CoolingOI063um = 0.0;
    *CoolingOI044um = 0.0;
    *CoolingOI146um = 0.0;
    return 0.e0;
#endif
}
static double CoolingOIDetail(double nOI, double T, double nel, double nHII, double nHI, double nH2,
        double *nOI0, double *nOI1, double *nOI2, double *CoolingOI063um, double *CoolingOI044um, double *CoolingOI146um){
    
    double E10 = Energy_OI063um;     
    double E20 = Energy_OI044um;     
    double E21 = Energy_OI146um;     

    double I10 = IntensityCMB_OI063um;
    double I20 = IntensityCMB_OI044um;
    double I21 = IntensityCMB_OI146um;

    double A10 = EinsteinA10_OI063um; 
    double A20 = EinsteinA10_OI044um; 
    double A21 = EinsteinA10_OI146um; 

    double B10 = EinsteinB10_OI063um;
    double B20 = EinsteinB10_OI044um;
    double B21 = EinsteinB10_OI146um;

    double B01 = EinsteinB01_OI063um;
    double B02 = EinsteinB01_OI044um;
    double B12 = EinsteinB01_OI146um;
    
    double g1g0 = gulOI063um; 
    double g2g0 = gulOI044um;
    double g2g1 = gulOI146um;

    int NcolSpecies = 4;
    double ncol[] = {OrthH2Frac*nH2,ParaH2Frac*nH2,nHI          ,nHII         };
    double q10[]  = {q10_OI_oH2(T) ,q10_OI_pH2(T) ,q10_OI_HI(T) ,q10_OI_HII(T)};
    double q20[]  = {q20_OI_oH2(T) ,q20_OI_pH2(T) ,q20_OI_HI(T) ,q20_OI_HII(T)};
    double q21[]  = {q21_OI_oH2(T) ,q21_OI_pH2(T) ,q21_OI_HI(T) ,q21_OI_HII(T)};

    double kT = BOLTZMANN_CONSTANT_CGS*T;
    double n2 = 0.e0;    // uppest level [1/cc]
    double n1 = 0.5*nOI; // middle level [1/cc]
    double n0 = 0.5*nOI; // lowest level [1/cc]
    SolveThreeLevelSystemStatisticalEquilibrium(
            E10,E20,E21,g1g0,g2g0,g2g1,I10,I20,I21,
            A10,A20,A21,B10,B20,B21,B01,B02,B12,
            NcolSpecies,ncol,q10,q20,q21,
            nOI,kT,&n0,&n1,&n2);

    *nOI0 = n0;
    *nOI1 = n1;
    *nOI2 = n2;
    //*ncrit = A10;

    /// spontaneous + stimulated emissions [erg/s/cm^3]
    double lamb063um = (A10 + B10*I10)*E10*n1;
    double lamb044um = (A20 + B20*I20)*E20*n2;
    double lamb146um = (A21 + B21*I21)*E21*n2;
    /// absorption [erg/s/cm^3]
    double gamm063um = B01*I10*E10*n0;
    double gamm044um = B02*I20*E20*n0;
    double gamm146um = B12*I21*E21*n1;
    /// net energy loss [erg/s/cm^3]
    double net063um  = lamb063um - gamm063um;
    double net044um  = lamb044um - gamm044um;
    double net146um  = lamb146um - gamm146um;

    *CoolingOI063um = net063um;
    *CoolingOI044um = net044um;
    *CoolingOI146um = net146um;
    return net063um + net044um + net146um;
}

double CoolingOI63um(StructCEQData *PhydroPP){
    double nHII = PhydroPP->nion [ipHYDROGEN];
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    double nOI  = PhydroPP->natom[ipOXYGEN  ];
    double nel  = PhydroPP->nel; 
    double T    = PhydroPP->T;
    double CoolingOI063um, CoolingOI044um, CoolingOI146um;
    CoolingOI(nOI,T,nel,nHII,nHI,nH2,&CoolingOI063um,&CoolingOI044um,&CoolingOI146um);
    return CoolingOI063um;
}

double CoolingOI44um(StructCEQData *PhydroPP){
    double nHII = PhydroPP->nion [ipHYDROGEN];
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    double nOI  = PhydroPP->natom[ipOXYGEN  ];
    double nel  = PhydroPP->nel; 
    double T    = PhydroPP->T;
    double CoolingOI063um, CoolingOI044um, CoolingOI146um;
    CoolingOI(nOI,T,nel,nHII,nHI,nH2,&CoolingOI063um,&CoolingOI044um,&CoolingOI146um);
    return CoolingOI044um;
}

double CoolingOI146um(StructCEQData *PhydroPP){
    double nHII = PhydroPP->nion [ipHYDROGEN];
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    double nOI  = PhydroPP->natom[ipOXYGEN  ];
    double nel  = PhydroPP->nel; 
    double T    = PhydroPP->T;
    double CoolingOI063um, CoolingOI044um, CoolingOI146um;
    CoolingOI(nOI,T,nel,nHII,nHI,nH2,&CoolingOI063um,&CoolingOI044um,&CoolingOI146um);
    return CoolingOI146um;
}

double CoolingSiI(double nSiI, double T, double nHI, double nHII,
        double *CoolingSiI130um, double *CoolingSiI45um, double *CoolingSiI68um){
#ifdef FINE_STRUCUTRE_COOLING_SiI
    const double rate = NAN; // Si I cooling is unimplemented.
#else
    const double rate = 0.0;
#endif
    *CoolingSiI130um = rate;
    *CoolingSiI45um = rate;
    *CoolingSiI68um = rate;
    return rate;
}

double CoolingSiI130um(double nSiI, double T, double nHI, double nHII){
    double CoolingSiI130um, CoolingSiI45um, CoolingSiI68um;
    CoolingSiI(nSiI,T,nHI,nHII,&CoolingSiI130um,&CoolingSiI45um,&CoolingSiI68um);
    return CoolingSiI130um;
}

double CoolingSiI45um(double nSiI, double T, double nHI, double nHII){
    double CoolingSiI130um, CoolingSiI45um, CoolingSiI68um;
    CoolingSiI(nSiI,T,nHI,nHII,&CoolingSiI130um,&CoolingSiI45um,&CoolingSiI68um);
    return CoolingSiI45um;
}

double CoolingSiI68um(double nSiI, double T, double nHI, double nHII){
    double CoolingSiI130um, CoolingSiI45um, CoolingSiI68um;
    CoolingSiI(nSiI,T,nHI,nHII,&CoolingSiI130um,&CoolingSiI45um,&CoolingSiI68um);
    return CoolingSiI68um;
}

// Fine structure transition between ^2P1/2-^2P3/2
double CoolingSiII35um(StructCEQData *PhydroPP){
    double nSiII = PhydroPP->nion [ipSILICON ];
    double nHI   = PhydroPP->natom[ipHYDROGEN];
    double nel   = PhydroPP->nel; 
    double T     = PhydroPP->T;
#ifdef FINE_STRUCUTRE_COOLING_SiII
    double n0, n1, ncrit;
    return CoolingSiII35umDetail(nSiII,T,nel,nHI,&n0,&n1,&ncrit);
#else
    return 0.0;
#endif
}
static double CoolingSiII35umDetail(double nSiII, double T, double nel, double nHI, double *nSiII0, double *nSiII1, double *ncrit){
    double I10   = IntensityCMB_SiII035um;
    double E10   = Energy_SiII035um;
    double A10   = EinsteinA10_SiII035um;
    double B10   = EinsteinB10_SiII035um;         
    double B01   = EinsteinB01_SiII035um;
    double g1g0  = gulSiII035um; 

    int NcolSpecies = 2;
    double ncol[] = {nHI           ,nel         };
    double q10[]  = {q10_SiII_HI(T),q10_SiII_e(T)};

    double kT = BOLTZMANN_CONSTANT_CGS*T;
    double n1 = 0.5*nSiII;    // upper level [1/cc]
    double n0 = 0.5*nSiII;    // lower level [1/cc]
    SolveTwoLevelSystemStatisticalEquilibrium(
        E10,g1g0,
        I10,A10,B10,B01,
        NcolSpecies,ncol,q10, 
        nSiII,kT,&n0,&n1);

    double C10 = 0.e0;
    for(int k=0;k<NcolSpecies;k++)   C10 += q10[k]*ncol[k];
    *nSiII0 = n0;
    *nSiII1 = n1;
    *ncrit = A10/C10;

    /// spontaneous + stimulated emissions [erg/s/cm^3]
    double lamb = (A10 + B10*I10)*E10*n1;
    /// absorption [erg/s/cm^3]
    double gamm = B01*I10*E10*n0;
    /// net energy loss [erg/s/cm^3]
    // if n0 > (A10+B10*I10)*n1/{(g1/g0)*B10*I10}, 
    // then the ions will absorb more energy than they emit (= heating)
    return lamb - gamm;
}


//////// De-excitation Rates [cm^3/s] for CI
// ref.) Table 6 of Glover & Jappsen (2007)
static double q10_CI_oH2(double T){
    return 8.7e-11 - 6.6e-11*exp(-T/218.3) + 6.6e-11*exp(-2*T/218.3);
}
static double q20_CI_oH2(double T){
    return 1.2e-10 - 6.1e-11*exp(-T/387.3);
}
static double q21_CI_oH2(double T){
    return 2.9e-10 - 1.9e-10*exp(-T/348.9);
}

static double q10_CI_pH2(double T){
    return 7.9e-11 - 8.7e-11*exp(-T/126.4) + 1.3e-10*exp(-2*T/126.4);
}
static double q20_CI_pH2(double T){
    return 1.1e-10 - 8.6e-11*exp(-T/223.0) + 8.7e-11*exp(-2*T/223.0);
}
static double q21_CI_pH2(double T){
    return 2.7e-10 - 2.6e-10*exp(-T/250.7) + 1.8e-10*exp(-2*T/250.7);
}

static double q10_CI_HI(double T){
    return 1.6e-10*pow(0.01*T,0.14);
}
static double q20_CI_HI(double T){
    return 9.2e-11*pow(0.01*T,0.26);
}
static double q21_CI_HI(double T){
    return 2.9e-10*pow(0.01*T,0.26);
}

static double q10_CI_HII(double T){
    if(T<=5000.0)   return (9.6e-11 - 1.8e-14*T + 1.9e-18*SQ(T))*pow(T,0.45);
    else            return 8.9e-10*pow(T,0.117);
}
static double q20_CI_HII(double T){
    if(T<=5000.0)   return (3.1e-12 - 6.0e-16*T + 3.9e-20*SQ(T))*T;
    else            return 2.3e-9*pow(T,0.0965);
}
static double q21_CI_HII(double T){
    if(T<=5000.0)   return (1.0e-10 - 2.2e-14*T + 1.7e-18*SQ(T))*pow(T,0.70);
    else            return 9.2e-9*pow(T,0.0535);
}

static double q10_CI_e(double T){
    double lnT  = log(T);
    double lnT2 = SQ(lnT), lnT3 = lnT*lnT2, lnT4 = lnT*lnT3;
    if(T<=1.e+3)
        return 2.88e-6*pow(T,-0.5)*exp(-9.25141 - 7.73782e-1*lnT + 3.61184e-1*lnT2 
                - 1.50892e-2*lnT3 - 6.56325e-4*lnT4);
    else
        return 2.88e-6*pow(T,-0.5)*exp(-4.44600e+2 - 2.27913e+2*lnT + 4.25950e+1*lnT2 
                - 3.47620e+0*lnT3 + 1.0508e-1*lnT4);
}
static double q20_CI_e(double T){
    double lnT  = log(T);
    double lnT2 = SQ(lnT), lnT3 = lnT*lnT2, lnT4 = lnT*lnT3;
    if(T<=1.e+3)
        return 1.73e-6*pow(T,-0.5)*exp(-7.69735 - 1.30743*lnT + 0.697638*lnT2 
                - 0.111338*lnT3 + 0.705277e-2*lnT4);
    else
        return 1.73e-6*pow(T,-0.5)*exp(3.50609e+2 - 1.87474e+2*lnT + 3.61803e+1*lnT2 
                - 3.03283e+0*lnT3 + 9.38138e-2*lnT4);
}
static double q21_CI_e(double T){
    double lnT  = log(T);
    double lnT2 = SQ(lnT), lnT3 = lnT*lnT2, lnT4 = lnT*lnT3;
    if(T<=1.e+3)
        return 1.73e-6*pow(T,-0.5)*exp(-7.4387 - 0.57443*lnT + 0.358264*lnT2
                - 4.18166e-2*lnT3 + 2.35272e-3*lnT4);
    else
        return 1.73e-6*pow(T,-0.5)*exp(3.86186e+2 - 2.02192e+2*lnT + 3.85049e+1*lnT2
                - 3.19268e+0*lnT3 + 9.78573e-2*lnT4);
}


//////// De-excitation Rates [cm^3/s] for [CII] 157.7um
// Glover & Jappsen (2007), Table 6 and references 8–9.
// Anchor the published T^0.07 extrapolation to the low-temperature fit
// at 250 K to remove the jump in the tabulated formula.
static const double CIIH2FitUpperTemperature = 250.0;

static double q10_CII_oH2(double T){
    const double fitTemperature = fmin(T,CIIH2FitUpperTemperature);
    const double rate = 4.70e-10 + 4.6e-13*fitTemperature;
    return T <= CIIH2FitUpperTemperature ? rate : rate*pow(T/CIIH2FitUpperTemperature,0.07);
}
static double q10_CII_pH2(double T){
    const double fitTemperature = fmin(T,CIIH2FitUpperTemperature);
    const double rate = 2.50e-10*pow(fitTemperature,0.12);
    return T <= CIIH2FitUpperTemperature ? rate : rate*pow(T/CIIH2FitUpperTemperature,0.07);
}
static double q10_CII_HI(double T){
    if(T<=2000.0)   return 8.0e-10*pow(0.01*T,0.07);
    else            return 3.1e-10*pow(0.01*T,0.385);
}
static double q10_CII_e(double T){
    if(T<=2000.0)   return 3.86e-7*pow(0.01*T,-0.5);
    else            return 2.43e-7*pow(0.01*T,-0.345);
}


//////// De-excitation Rates [cm^3/s] for [OI] 63.1um, 44.2um, 145.6um
// ref.) Table 6 of Glover & Jappsen (2007)
static double q10_OI_oH2(double T){
    return 2.7e-11*pow(T,0.362);
}
static double q20_OI_oH2(double T){
    return 5.49e-11*pow(T,0.317);
}
static double q21_OI_oH2(double T){
    return 2.74e-14*pow(T,1.060);
}

static double q10_OI_pH2(double T){
    return 3.46e-11*pow(T,0.316);
}
static double q20_OI_pH2(double T){
    return 7.07e-11*pow(T,0.268);
}
static double q21_OI_pH2(double T){
    return 3.33e-15*pow(T,1.360);
}

static double q10_OI_HI(double T){
    return 9.2e-11*pow(0.01*T,0.67);
}
static double q20_OI_HI(double T){
    return 4.3e-11*pow(0.01*T,0.80);
}
static double q21_OI_HI(double T){
    return 1.1e-10*pow(0.01*T,0.44);
}

static double q10_OI_HII(double T){
    if(T<=194.0)        return 6.38e-11*pow(T,0.40);
    else if(T<=3686.0)  return 7.75e-12*pow(T,0.80);
    else                return 2.65e-10*pow(T,0.37);
}
static double q20_OI_HII(double T){
    if(T<=511.0)        return 6.10e-13*pow(T,1.10);
    else if(T<=7510.0)  return 2.12e-12*pow(T,0.90);
    else                return 4.49e-10*pow(T,0.30);
}
static double q21_OI_HII(double T){
    if(T<=2090.0)   return 2.03e-11*pow(T,0.56);
    else            return 3.43e-10*pow(T,0.19);
}

static double q10_OI_e(double T){
    return 5.12e-10*pow(T,-0.075);
}
static double q20_OI_e(double T){
    return 4.86e-10*pow(T,-0.026);
}
static double q21_OI_e(double T){
    return 1.08e-14*pow(T,0.926);
}


//////// De-excitation Rates [cm^3/s] for [SiI] 130um, 45um, 68um
// ref.) Table 6 of Glover & Jappsen (2007)
static double q10_SiI_HI(double T){
    return 3.5e-10*pow(0.01*T,-0.03);
}
static double q20_SiI_HI(double T){
    return 1.7e-11*pow(0.01*T,0.17);
}
static double q21_SiI_HI(double T){
    return 5.0e-10*pow(0.01*T,0.17);
}

static double q10_SiI_HII(double T){
    return 7.2e-9;
}
static double q20_SiI_HII(double T){
    return 7.2e-9;
}
static double q21_SiI_HII(double T){
    return 2.2e-8;
}

//////// De-excitation Rates [cm^3/s] for [SiII] 35um
// ref.) Table 6 of Glover & Jappsen (2007)
static double q10_SiII_HI(double T){
    return 4.95e-10*pow(0.01*T,0.24);
}
static double q10_SiII_e(double T){
    return 1.2e-6*pow(0.01*T,-0.5);
}

#if 0
////////////////
static bool CIIcolfirst = true;
double TmpCIIH2col[] = { // Table 4 in Flower & Launay (1977)    
    10,20,30,40,50,60,70,80,90,100,
    120,140,160,180,200,220,250
};
double gamCIIparaH2col[] = { // Table 4 in Flower & Launay (1977)
    8.2e-4,8.9e-2,0.43,0.97,1.6,2.2,2.8,3.3,3.9,4.3,
    5.1,5.8,6.4,6.9,7.3,7.7,8.2
};
double gamCIIorthoH2col[] = { // Table 4 in Flower & Launay (1977)
    1.2e-3,0.12,0.57,1.2,2.0,2.7,3.4,4.1,4.7,5.2,
    6.2,7.0,7.8,8.4,9.0,9.5,10.0
};
static int NumberofCIIH2col = 0;
static gsl_interp_accel *TmpCIIH2colAccess;
static gsl_spline       *gamCIIparaH2colAccess;
static gsl_spline       *gamCIIorthoH2colAccess;
static double CIIH2collisionCoolingCoefficient(double logT){   // [cm3/s]
    double val1 = gsl_spline_eval(gamCIIparaH2colAccess ,logT,TmpCIIH2colAccess);
    double val2 = gsl_spline_eval(gamCIIorthoH2colAccess,logT,TmpCIIH2colAccess);
    return (val1+val2)*1.e+24;
    //return pow(10.0,val1+val2);
}
static void InitializeCIIH2CollisionCoolingCoefficient(void){

    NumberofCIIH2col = sizeof(TmpCIIH2col)/sizeof(double);
    if(MPIGetMyID() == MPI_ROOT_RANK)
        dprintlmpi(NumberofCIIH2col);

    for(int i=0; i<NumberofCIIH2col; i++){
        TmpCIIH2col[i] = log10(TmpCIIH2col[i]);
        //gamCIIparaH2col[i]  = log10(gamCIIparaH2col[i]);
        //gamCIIorthoH2col[i] = log10(gamCIIorthoH2col[i]);
    }
     
    TmpCIIH2colAccess = gsl_interp_accel_alloc();
    gamCIIparaH2colAccess  = gsl_spline_alloc(gsl_interp_linear,NumberofCIIH2col);
    gamCIIorthoH2colAccess = gsl_spline_alloc(gsl_interp_linear,NumberofCIIH2col);
    gsl_spline_init(gamCIIparaH2colAccess ,TmpCIIH2col,gamCIIparaH2col ,NumberofCIIH2col);
    gsl_spline_init(gamCIIorthoH2colAccess,TmpCIIH2col,gamCIIorthoH2col,NumberofCIIH2col);
    
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *output = fopen("CII-H2CollisionDestructionCoefficiency.dat","w");
        for(double logT=1.0; logT<=5.0; logT+=0.1)
            fprintf(output,"%e %e\n",logT,CIIH2collisionCoolingCoefficient(logT)); 
        fclose(output);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    CIIcolfirst = false;
         
    return ;
}
#endif
