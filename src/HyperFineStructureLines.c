#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "HyperFineStructureLines.h"

static const double Zsol = SOLAR_METALLICITY;

#define A10HI21cm   (2.8843e-15)    // [1/s]
#define THI21cm     (6.81e-2)       // [K]
#define nyuHI21cm   (1.4204e+9)     // [Hz]

static double Energy_HI21cm;
static double EinsteinA10_HI21cm;

void InitializeHyperFineStructureLines(int Check){

    Energy_HI21cm      = PLANK_H_CGS*nyuHI21cm;
    EinsteinA10_HI21cm = A10HI21cm;

    return ;
}

static double NormalizedLineProfileHI21cm(double v_ch, double T, double vlos){
    return 1.0;
    //double dv = sqrt(2.0*BOLTZMANN_CONSTANT_CGS*T/PROTON_MASS_CGS);
}

static double k10H(double T){
    // collisional de-excitation rate coefficient with HI (eq.17.8 of Draine 2011)
    double T2   = 1.e-2*T;
    double lnT2 = log(T2); 
    if(T<20.0)          return 0.0;
    else if(T<300.0)    return 1.19e-10*pow(T2,0.74-0.20*lnT2);
    else                return 2.24e-10*pow(T2,0.207)*exp(-0.876/T2);
}

// Spin Temperature [K]
double SpinTemperatureHI21cm(double nHI, double T_k){
    double T_R  = 3.77;  // [K]; inc. CMBR + synchrotron (Draine 2011, Chap.7)
    double T_a  = T_k;   // Field (1959)
    double n_a  = 1.e-6; // [/cc] (Liszt 2001)
    double Rc10 = nHI*k10H(T_k);
    double y_c  = THI21cm*Rc10/(T_k*EinsteinA10_HI21cm);
    //double y_a = THI21cm*Ra10/(T_a*EinsteinA10_HI21cm);
    double y_a  = 5.9e+11*n_a/(T_a*sqrt(T_k));  // Field (1958)
    return (T_R + y_c*T_k + y_a*T_a)/(1 + y_c + y_a);
}

// Emission coefficient [erg/s/cc/str]
double EmisivityCoefficientHI21cm(double nHI, double T){
    return (3.0/16.0/PI)*PLANK_H_CGS*nyuHI21cm*EinsteinA10_HI21cm*nHI;
}

// Absorpsion coefficient [/cm]
double AbsorptionCoefficientHI21cm(double nHI, double T){
    double Ts = SpinTemperatureHI21cm(nHI,T);
    return (3.0/32.0/PI)*(PLANK_H_CGS*SQ(LIGHT_C_CGS)/BOLTZMANN_CONSTANT_CGS/nyuHI21cm)*EinsteinA10_HI21cm*(nHI/Ts);
}

// [erg/s/cc]
double CoolingHI21cm(double nHI, double T){
    return 4.0*PI*EmisivityCoefficientHI21cm(nHI,T);
}

double SourceFunctionHI21cm(double nHI, double T, double vlos, double v_ch, double *j_nyu, double *k_nyu){
    double j = EmisivityCoefficientHI21cm(nHI,T);
    double k = AbsorptionCoefficientHI21cm(nHI,T);
    //double phi = NormalizedLineProfileHI21cm(v_ch,T,vlos);
    *j_nyu = j;
    *k_nyu = k;
    return j/k;
}
