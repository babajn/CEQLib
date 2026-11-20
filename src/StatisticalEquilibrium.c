#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalElements.h"
#include "StatisticalEquilibrium.h"

/// Einstein Coefficient (Stimulated Emission)
double EinsteinCoefficientB10(double A10, double E10){
    double nyu10 = E10/PLANK_H_CGS;  // [Hz]
    return A10*CUBE(LIGHT_C_CGS/nyu10)/(8.0*PI*PLANK_H_CGS); 

}
/// Einstein Coefficient (Absorption)
double EinsteinCoefficientB01(double B10, double g1g0){
    return g1g0*B10;
}

void SolveTwoLevelSystemStatisticalEquilibrium(
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
        fprintf(stderr," kT=%e \n",kT); 
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
void SolveThreeLevelSystemStatisticalEquilibrium(
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
