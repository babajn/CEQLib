#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "GrainSurfaceChemistry.h"
#include "ChemicalEquilibriumCooling.h"

static double Zsol = SOLAR_METALLICITY;

extern double DustShieldingFactor(double Ntot, double Zd, double gamma);
static double PhotoelectricHeatingEfficiency(double T, double k);

void InitializeGrainSurfaceChemistry(int Check){

    if(Check==1){
        FILE *output = fopen(CEQ_DIAGNOSTIC_DIRECTORY "/GrainSurfaceChemistry.dat","w");
        for(double logk=0.0; logk<=7.0; logk+=0.01){
            double k = pow(10.0,logk);
            fprintf(output,"%e ",k);
            for(double logT = 1.0; logT<=5.0; logT+=1.0){
                double T   = pow(10.0,logT);
                double eff = PhotoelectricHeatingEfficiency(T,k);
                fprintf(output,"%e ",eff);
            }
            fprintf(output,"\n");
        }
        fclose(output);
    }

    return ;
}

//// Recombination on the Dust [cm^3/s]
//// ref. Weingartner & Draine, ApJ, 563, 842, (2001) 
static double C_HII[]  = {12.25,8.074e-6,1.378,5.087e+2,1.586e-2,0.4723,1.102e-5};
static double C_HeII[] = {5.572,3.185e-7,1.512,5.115e+3,3.903e-7,0.4956,5.494e-7};
static double C_CII[]  = {45.58,6.089e-3,1.128,4.331e+2,4.845e-2,0.8120,1.333e-4};
static double C_SiII[] = {2.166,5.678e-8,1.874,4.375e+4,1.635e-6,0.8964,7.538e-5};
static double RecombinationRateOnDust(double T, double G0, double ne, double Z, double C[7]){
    if(CEQLibRunParameters.UseGrainSurfaceProcesses==false || Z==0.0) return 0.0;
    double numerator = 1.e-14*(Z/Zsol)*C[0];
    if(G0==0.0) return numerator;
    if(ne==0.0) return 0.0;
    double psi = G0*sqrt(T)/ne;
    double secondExponent = C[2] - C[5] - C[6]*log(T);
    return numerator/(1.0 + C[1]*(pow(psi,C[2]) + C[3]*pow(T,C[4])*pow(psi,secondExponent)));
}
double HIIRecombinationRateOnDust(double T, double G0, double ne, double Z){
    return RecombinationRateOnDust(T,G0,ne,Z,C_HII);
}
double DIIRecombinationRateOnDust(double T, double G0, double ne, double Z){
    // see R62 in Table 3 of Glover & Jappsen (2007)
    return 0.71*HIIRecombinationRateOnDust(T,G0,ne,Z);  // 1/sqrt(2)
}
double HeIIRecombinationRateOnDust(double T, double G0, double ne, double Z){
    return RecombinationRateOnDust(T,G0,ne,Z,C_HeII);
}
double CIIRecombinationRateOnDust(double T, double G0, double ne, double Z){
    return RecombinationRateOnDust(T,G0,ne,Z,C_CII);
}
double OIIRecombinationRateOnDust(double T, double G0, double ne, double Z){
    // see R65 in Table 3 of Glover & Jappsen (2007)
    return 0.25*HIIRecombinationRateOnDust(T,G0,ne,Z);
}
double SiIIRecombinationRateOnDust(double T, double G0, double ne, double Z){
    return RecombinationRateOnDust(T,G0,ne,Z,C_SiII);
}


//#define DUST_TEMPERATURE    (8.0)  // [K]
//#define DUST_TEMPERATURE    (10.0)  // [K]
#define DUST_TEMPERATURE    (20.0)  // [K]
//double TotalDustCoolingRate(double T, double nH, double ne, double Z, double G0){
//    double Tdust = DUST_TEMPERATURE;
//    double Lrec  = CoolingDustRecombination(T,nH,ne,Z,G0);
//    double Lgd   = CoolingGasDustEnergyTransfer(T,Tdust,nH,Z);
//    double Hpho  = HeatingDustPhotoelectic(T,nH,ne,Z,G0);
//    if(isnan(Lrec)){
//        fprintf(stderr,"Ldustrec=%e ne=%e nH=%e T=%e G0=%e\n",Lrec,ne,nH,T,G0);
//        exit(1);
//    }
//    if(isnan(Lgd)){
//        fprintf(stderr,"Lgd=%e\n",Lgd);
//        exit(1);
//    }
//    if(isnan(Hpho)){
//        fprintf(stderr,"Hphoel=%e\n",Hpho);
//        exit(1);
//    }
//    return Lrec+Lgd-Hpho;
//}
//
//double TotalDustCoolingRateHIIregion(double T, double nH, double ne, double Z, double G0){
//    double Tdust= DUST_TEMPERATURE;
//    double Lrec = CoolingDustRecombination(T,nH,ne,Z,G0);
//    double Hpho = HeatingDustPhotoelectic(T,nH,ne,Z,G0);
//    if(isnan(Lrec)){
//        fprintf(stderr,"Ldustrec=%e ne=%e nH=%e T=%e G0=%e\n",
//                Lrec,ne,nH,T,G0);
//        exit(1);
//    }
//    if(isnan(Hpho)){
//        fprintf(stderr,"Hphoel=%e\n",Hpho);
//        exit(1);
//    }
//    return Lrec-Hpho;
//}


//// Dust recombination cooling [erg/cc/s]
// eq.(44) of Bakes & Tielens (1994)
// eq.(9) of Wolfire et al. (1995)
#define al (0.944)
#define bt (0.735)
#define phiPAH (1.0)   
// adjustable parameter introduced by Wolfire et al. (2003) to make the heating 
// rate consistent with the values of the electron recombination rates that are 
// inferred observationally for PAHs; the original Bakes & Tielens (1994) 
// treatment corresponds to phiPAH=1
double CoolingDustRecombination(StructCEQData *PhydroPP){
    double T    = PhydroPP->T;
    double nH   = PhydroPP->n[ipHYDROGEN]; 
    double ne   = PhydroPP->nel; 
    double Zd   = PhydroPP->Zd;
    double G0   = PhydroPP->G0;
    double Ntot = PhydroPP->Nt;
    if(CEQLibRunParameters.UseGrainSurfaceProcesses==false || ne==0.0 || Zd==0.0) return 0.0;
    if(CEQLibRunParameters.UseDustAttenuation==true)
        G0 *= DustShieldingFactor(Ntot,Zd,2.5);
    double k = G0*sqrt(T)/(ne*phiPAH);
    double beta = bt*pow(T,-0.068); 
    //return (4.65e-30)*phiPAH*pow(T,al)*pow(k,beta)*(Z/Zsol)*ne*nH; // eq.(9) of W+95
    return (3.49e-30)*phiPAH*pow(T,al)*pow(k,beta)*(Zd/Zsol)*ne*nH;  // eq.(44) of BT94
}

//// Dust Photoelectric Heating [erg/cc/s]
// eq.(42)(43) of Bakes & Tielens (1994)
static double PhotoelectricHeatingEfficiency(double T, double k){
    return (4.87e-2)/(1.0+(4.0e-3)*pow(k,0.73)) + (3.65e-2)*pow(1.e-4*T,0.7)/(1.0+(2.0e-4)*k);
}
double HeatingDustPhotoelectic(StructCEQData *PhydroPP){
    double nH    = PhydroPP->n[ipHYDROGEN]; 
    double ne    = PhydroPP->nel; 
    double T     = PhydroPP->T;
    double Zd    = PhydroPP->Zd;
    double G0    = PhydroPP->G0;
    double Ntot  = PhydroPP->Nt;
    if(ne==0.0 || G0==0.0 || Zd==0.0) return 0.0;
    if(CEQLibRunParameters.UseDustAttenuation==true)
        G0 *= DustShieldingFactor(Ntot,Zd,2.5);
    double k = G0*sqrt(T)/(ne*phiPAH);
    double eff = PhotoelectricHeatingEfficiency(T,k);
    return 1.0e-24*eff*G0*nH*(Zd/Zsol); // eq.(11) of Wolfire et al. (1995) 
    //return (1.3e-24)*eff*G0*nH*(Z/Zsol); // eq.(42) of Glover & MacLow (2007)
}

//// Gas-Dust Energy Transfer Cooling/Heating [erg/cc/s]
// ref. eq.(2.15) of Hollenbach & McKee (1989)
#define MinimunDustSize (100.0) // [A]
double CoolingGasDustEnergyTransfer(double T, double Tdust, double nH, double Z){
#ifdef GAS_DUST_ENERGY_EXCHANG
    return (1.2e-31)*SQ(nH)*(Z/Zsol)*sqrt(1.e-3*T)*sqrt(100.0/MinimunDustSize)*(1.e0-0.8*exp(-75.0/T))*(T-Tdust);
#endif
    return 0.e0;
}
double HeatingGasDustEnergyTransfer(double T, double Tdust, double nH, double Z){
    return CoolingGasDustEnergyTransfer(T,Tdust,nH,Z);
}
//// Planck mean opacity due to the grains
// ref. Pollack et al. (1994)
double PlanckMeanGrainOpacity(double Tdust){
    return 4.0e-4*SQ(Tdust);    // cm^2/g
}

//// Approximate Dust Temperature
// ref. Hollenbach et al., ApJ, 377, 192, (1991)
//      see also Hosokawa & Inutsuka, ApJ, 646, 240, (2006)
// to determine the dust temperature in the cloud, we assume that the
// incident radiation on a grain consists of the attenuated incident UV
// fieldm the 2.7 K background, and the infrared radiation field
// represented by dust emission from an effective slab thickness of Av~1
// corresponding to 100um emission optical depth of tau100=0.001 and a
// constant temperature given by T0=12.2*Go^0.2 [K].
#define tau100  (1.e-3)     // the emission optical depth at 100um
#define nu0     (3.e+15)    // [1/s]: characteristic frequency
#define sigmaV  (4.9e-22)   // dust cross-section @ V-band [cm^2] 
double ApproximateDustTemperature(double Ntot, double G0){
#ifdef SOLVE_DUST_TEMPERATURE
    double Av = 1.086*sigmaV*Ntot;
    //double Av = sigmaV*Ntot;   // eq.(21.7) in Draine's textbook (2011)
    double T0 = 12.2*pow(G0,0.2);   // eq.(6) [K]
    return pow(8.9e-11*nu0*G0*exp(-1.8*Av)+143.48907
            +3.4e-2*(0.42-log(3.5e-2*tau100*T0))*tau100*pow(T0,6.0),0.2);
#endif
    return DUST_TEMPERATURE;
}

double EffectiveDustTemperature(const StructCEQData *state){
    return isfinite(state->Td) && state->Td > 0.0 ?
           state->Td : ApproximateDustTemperature(state->Nt,state->G0);
}

double ApproximateDustToGasRatio(double T){
    //return 1.e0;
    double logT = log10(T);
    return (0.5*PI-atan(logT-5.0))/PI;
}

