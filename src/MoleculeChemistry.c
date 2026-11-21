#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalElements.h"
#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#include "ShieldingFactorCO.h"
#include "CosmicRay.h"
#include "Xray.h"
#include "FineStructureLines.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

static double Zsol = SOLAR_METALLICITY;
#define FUVPumpingHeatingFactorDB96 (8.5)   // Draine & Bertoldi (1996)

extern double b_turb;

//static double RateEquation(int index, double fH2);
//static double Sh(double Tk);
//static double SelfShieldingFactor(double N2);
//static void InitializeCollisionDissociationCoefficient(void);
//static double HIcollisionCoefficientMartin1998(double logT);
//static double HIcollisionCoefficient(double T, double nH, double nHI, double nH2);
//static double H2collisionCoefficient(double T, double nH, double nHI, double nH2);
static double CriticalDensity(double T, double yHI, double yH2);
static double f_co_self(double logNco);
static double f_co_H2(double logNh2);

static int NumberTableShieldingCO;
static int NumberTableShieldingH2;
static gsl_spline *logShieldCOAccess;
static gsl_spline *logShieldH2Access;
static double h2ShieldingTailPower;
static double h2ShieldingTailOpticalDepth;

void InitializeMoleculeChemistory(int Check){

    if(EQCLibRunParameters.UseH2Formation==true&&Check==1){
        FILE *output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/H2ShieldingFactor.dat","w");
        for(double logNh2=12.0; logNh2<=23.0; logNh2+=0.01){
            double b5   = 3.0;
            double f_H2 = H2SelfShieldingFactor(pow(10,logNh2),b5);
            fprintf(output,"%e %e\n",logNh2,f_H2);
        }
        fclose(output);
        
        output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CoolingRateH2.dat","w");
        double nHI = 1.0, nHII = 0.e0, ne = nHII, nHe = 0.0; 
        //for(double lognH2=-1.0; lognH2<=4.0; lognH2+=1.0){
        for(double lognH2=2.0; lognH2<=2.0; lognH2+=1.0){
            double nH2 = pow(10.0,lognH2);
            for(double logT=1.0; logT<=4.0; logT+=0.01){
                double T = pow(10.0,logT);
                double L_rv_NK93 = CoolingH2RovibrationalLineNeufeldKaufman1993(T,nHI,nH2);
                double L_rv_HM79 = CoolingH2RovibrationalLineHollenbachMcKee1979(T,nHI,nH2);
                double L_rv_GA08 = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(T,nHII,nHI,nH2,nHe,ne);
                double L_rv_GP98 = CoolingH2RovibrationalLineLowDensityLimitGalliPalla1998(T,nHI,nH2);
                fprintf(output,"%e ",T);
                fprintf(output,"%e ",nH2);
                fprintf(output,"%e ",L_rv_NK93);
                fprintf(output,"%e ",L_rv_HM79);
                fprintf(output,"%e ",L_rv_GA08);
                fprintf(output,"%e ",L_rv_GP98);
                fprintf(output,"\n");
            }
            fprintf(output,"\n");
        }
        fclose(output);
    }

    if(NumberTableShieldingCO == 0){
        NumberTableShieldingCO = sizeof(TableNco)/sizeof(double);
        NumberTableShieldingH2 = sizeof(TableNh2)/sizeof(double);

        for(int i=0; i<NumberTableShieldingCO; i++){
            TableNco[i]         = log10(TableNco[i]);
            TableShieldingCO[i] = log10(TableShieldingCO[i]);
            //fprintf(stderr,"%f %f\n",TableNco[i],TableShieldingCO[i]);
        }
        logShieldCOAccess = gsl_spline_alloc(gsl_interp_linear,NumberTableShieldingCO);
        gsl_spline_init(logShieldCOAccess,TableNco,TableShieldingCO,NumberTableShieldingCO);

        for(int i=0; i<NumberTableShieldingH2; i++){
            TableNh2[i]         = log10(TableNh2[i]);
            TableShieldingH2[i] = log10(TableShieldingH2[i]);
            //fprintf(stderr,"%f %f\n",TableNh2[i],TableShieldingH2[i]);
        }
        logShieldH2Access = gsl_spline_alloc(gsl_interp_linear,NumberTableShieldingH2);
        gsl_spline_init(logShieldH2Access,TableNh2,TableShieldingH2,NumberTableShieldingH2);

        /* Continue the Lee et al. H2 shielding as S_edge*x^-p*exp(-tau*(x-1)),
           where x=N_H2/N_edge. Match the slope at the edge and the decrease
           over the next decade from Visser et al. (2009), Table 5, 12CO with
           N_CO=0. This keeps the tail normalized to the Lee et al. table. */
        const int last = NumberTableShieldingH2-1;
        const double edgeSlope = (TableShieldingH2[last-1]-TableShieldingH2[last]) /
                                 (TableNh2[last]-TableNh2[last-1]);
        const double logDecadeAttenuation = log(1.108e-2/3.938e-7);
        h2ShieldingTailOpticalDepth = (logDecadeAttenuation-edgeSlope*log(10.0))/(9.0-log(10.0));
        h2ShieldingTailPower = edgeSlope-h2ShieldingTailOpticalDepth;
    }

    if(EQCLibRunParameters.UseCOFormation==true && Check==1){
        FILE *output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/COShieldingFactor.dat","w");
        for(double logNco=12.0; logNco<=23.0; logNco+=0.01){
            double f_self = f_co_self(logNco);
            fprintf(output,"%e %e\n",logNco,pow(10,f_self));
        }
        for(double logNh2=12.0; logNh2<=23.0; logNh2+=0.01){
            double f_H2   = f_co_H2(logNh2);
            fprintf(output,"%e %e\n",logNh2,pow(10,f_H2));
        }
        fclose(output);
    }

    return ;
}


//static bool HIcolfirst = true;

//double TotalH2FormationRate(double T, double Metal, double Tdust){
//    // [cm3/s]
//    //return MAX(CHEM_TINY,H2formationRate(T,Metal));
//    return H2FormationRate(T,Metal,Tdust);
//}
//
//double TotalH2DissociationRate(double nH, double nHI, double nH2, double ne, double T, double Metal, double G0, double Nt){
//    // [/s]
//    double N2 = Nt*nH2/nH;
//    double b5 = 12.90*sqrt(1.e-4*T);    // eq(6.34) of Draine's textbook
//    double dn2_diss = H2PhotoDissociationRate(G0,Nt,N2,b5,Metal);
//    double dn2_dest = H2CollisionalDissociationRate(T,nH,nHI,nH2,ne);
//    return MAX(CHEM_TINY,dn2_diss+dn2_dest);
//    //double dn2_dest1= nHI*HIcollisionCoefficient(T,nH,nHI,nH2);
//    //double dn2_dest2= nH2*H2collisionCoefficient(T,nH,nHI,nH2);
//    //return MAX(CHEM_TINY,dn2_diss + dn2_dest1 + dn2_dest2);
//}

/// see eq(1) and Table B2 of Glover et al. (2010)
double DustShieldingFactor(double Ntot, double Zd, double gamma){
    if(EQCLibRunParameters.UseDustAttenuation==false)   return 1.e0;

    // eq(1) of Glover et al. (2010)
    double Av = VisualExtinction(Ntot,Zd);
    return exp(-gamma*Av);
}
double VisualExtinction(double Ntot, double Zd){
    if(EQCLibRunParameters.UseDustAttenuation==false)   return 0.e0;

    // eq(2) of Glover et al. (2010)
    return 0.535e-21*(Zd/Zsol)*Ntot;
}

#define sigma1000MW     (2.0e-21)   // [cm^-2]
double DustShieldingFactorFUV(double Ntot, double Zd){
    if(EQCLibRunParameters.UseDustAttenuation==false)   return 1.e0;

    double tau1000 = (double)sigma1000MW*(Zd/Zsol)*Ntot;
    return exp(-tau1000);
}


// Tcr: cirtical temperatur of dusts
//#define TcrH2   (100.0) // [K]
//      Not considered in these equations is the thermal evapolation of 
//      H atoms from grains. Above a certain critical grain temperature Tcr,
//      the rate at which H2 is formed on a grain surface drops rapidly
//      because of this evapolation.
//  ref. Tcr=75 K (Appendix C of ielens & Hollenbach 1985)
//       25 K < Tcr < 50 K (Hollenbach & Salpeter 1971)
//       Tcr~100K (Cazaux & Tielens 2002,2004)

///// H2 formation on the grain surface 
//  ref. Buch & Zhang (1991)
//       Pelpessy et al. (2006)
/// Uncertainties of Rf: 
// H2 formation rate parameter (Pelpessy et al. 2006)
// The canonical formation rate Rf=3.e-17 cm^3/s (Jura 1974,1975) for
// typical CNM HI gas (T~100K), corresponds to mu = 3.5, but values of
// mu~17-18 are not excluded.
//static double mu = 3.5;
/// H2 formation Probability
// probability that H atoms eventually form an H2 molecule that detaches
// itself from the grain
//static double gamH2 = 1.0;
/// Sticking Probability
//  probability that H atoms stick on the grain surface
//  ref. Buch & Zhang (1991)
//static double T0 = 100.0;
//static double Sh(double T){
//    return 1.e0/SQ(1.0+T/T0);
//}
//double H2formationRate(double T, double Z){
//    // [cm3/s]
//    double sigma_d = Sh(T)*gamH2;
//    return 3.5e-17*mu*(Z/Zsol)*sqrt(1.e-2*T)*sigma_d*(CLUMPING_FACTOR_H2);
//}

// ref. Hollenbach & McKee (1979)
//      eq.(4,5) of Tielens & Hollenbach (1985)
//      eq.(A13,14) of Koyama & Inutsuka (2000)
// Cp: clumping facort
//  ref. Gnedin et al.,ApJ, 697,55,(2009)
//       Krumholz et al.,ApJ,699,850,(2009)
//  the gas is clustered on scales unresolved in the simulation,
//  so that the effective H2 formation rate is higher in proportion to
//  the gas clumping factor Cp. Numerical somulations of turbulent
//  molecular clouds with the width that depends on the average Much 
//  number of the flows: sigma_ln^2 = 1-2 (Padoan et al. 1997; Ostriker
//  et al. 2001; McKee & Ostriker 2007). Since for a lognormal
//  distribution, Cp = exp(sigma_ln^2), most likely values for the
//  clumping factor are in the range Cp=3-10. Empirical evidence also
//  indicates that the gas in molecular clouds is hightly clumped with
//  the ratio of typical molecualr gas density to the mean cloud density
//  of ~30-100 ($3.1.2 in McKee & Ostriker 2007).
//  Note that the clumping factor is expected to approach 1 as the
//  resosulution does ~100pc.
//   e.g. Cp = 30 (Krumholz & Gnedin, 2009)
double H2FormationRate(double T, double Z, double Tdust){
    double Cp = EQCLibRunParameters.CpForH2;
    // R60 in Table 3 of Glover & Jappsen (2007) 
    double St = 1.e0/(1.0 + 4.e-2*sqrt(T+Tdust) + 2.e-3*T + 8.e-6*SQ(T));
    //return 3.0e-18*Cp*sqrt(T)*(Z/Zsol)*St/(1.0 + 1.e+4*exp(-600.0/Tdust));
    return 3.0e-17*Cp*sqrt(T)*(Z/Zsol)*St/(1.0 + 1.e+4*exp(-600.0/Tdust));
    // (A13) of Koyama & Inutsuka 2000
    //return 3.e-17*sqrt(T/300.0)*St;
}

///// H2 dissociation by FUV 
// H2 Self-Sheilding factor: Draine & Bertoldi (1996)
double H2SelfShieldingFactor(double N2, double b5){
    // simplified formula
    // eq.(36) of Draine & Bertoldi (1996) 
    if(N2<1.e+14)   return 1.e0;
    else            return pow(1.e-14*N2,-0.75);

    //// more accurate formula
    // eq.(37) of Draine & Bertoldi (1996) 
    if(N2==0)   return 1.e0;
    double x  = 2.e-15*N2;
    double sqrt1x = sqrt(1.0+x);
    return 0.965/SQ(1.0+x/b5) + 0.035*exp(-8.5e-4*sqrt1x)/sqrt1x;
}

double H2PhotoDissociationRateWithDustFactor(double G0, double N2, double b5, double dustFactor){
    // eq.(29,31) of Glover & MacLow (2007)
    return G0*(3.3e-11*H2SelfShieldingFactor(N2,b5)*dustFactor);
}
double H2PhotoDissociationRate(double G0, double Nt, double N2, double b5, double Zd){
    return H2PhotoDissociationRateWithDustFactor(G0,N2,b5,DustShieldingFactorFUV(Nt,Zd));
}

H2CollisionalDissociationCoefficients H2DissociationCoefficients(double T){
    // R8, R9, R10 in Table 1 of Glover & Jappsen (2007)
    const double invT = 1.e0/T;
    return (H2CollisionalDissociationCoefficients) {
        .electron = 3.73e-9*pow(T,0.1121)*exp(-99430*invT),
        .atomicHydrogen = 6.67e-12*sqrt(T)*exp(-(1+63590*invT)),
        .molecularHydrogen = 5.996e-30*pow(T,4.1881)*exp(-54657.4*invT)*pow(1.0+6.761e-6*T,-5.6881)
    };
}

double H2CollisionalDissociationRate(double T, double nH, double nHI, double nH2, double ne){
    return H2DissociationFromCoefficients(H2DissociationCoefficients(T),ne,nHI,nH2);
}


/// Critical Density for H2-H collisions
// ref.) eq.(6) in Lepp & Shull (1983)
// ref.) eq.(B6a) in Shapiro & Kang (1987).
static double ncrHI(double logT4){
    return pow(10.0,4.00-0.416*logT4-0.327*SQ(logT4));  // [/cc] 
}
/// Critical Density for H2-H2 collisions
// ref.) eq.(B6b) in Shapiro & Kang (1987).
static double ncrH2(double logT4){
    return pow(10.0,4.845-1.3*logT4+1.62*SQ(logT4));  // [/cc] 
}
///// Harmonic Mean of Critical Densities: 
// ref.) Glover & MacLow, ApJ, 169, 239, (2007)
// a weighted harmonic mean of the value for H2-H collisions given by Lepp & Shull (1983)
// and the value for H2-H2 collisions given by Shapiro & Kang (1987).
static double CriticalDensity(double T, double yHI, double yH2){
    double logT4 = log10(T) - 4.0;
    return 2.0*ncrHI(logT4)*ncrH2(logT4)/(ncrH2(logT4)*yHI+ncrHI(logT4)*yH2);

    // see eq(A5)(A6) in Koyama & Inustuka (2000) 
    //return 1.e+6*sqrt(1/T)*(1.6*yHI*exp(-SQ(400.0/T))+1.4*yH2*exp(-1.2e+4/(T+1.2e+3)));
}


////////////////////////////////////////////////
//
// Cooling and Heating Rates associating with H2 chemistory
//
////////////////////////////////////////////////
//double TotalH2CoolingRate(double T, double nH, double nHI, double nH2, double Z, double G0, 
//        double Nt, double N2, double Tdust){
//    double b5 = 12.90*sqrt(1.e-4*T);    // eq(6.34) of Draine's textbook
//    double Lvibrot   = CoolingH2RovibrationalLine(T,nHI,nH2);
//    double Lh2coll   = CoolingH2CollisionalDissociation(T,nH,nHI,nH2);
//    double Hphotdiss = HeatingH2Photodissociate(nH,nHI,nH2,Nt,N2,G0,Z,b5);
//    double Hphotpump = HeatingH2PhotomPumping(T,nH,nHI,nH2,Nt,N2,G0,Z,b5);
//    double Hh2form   = HeatingH2Formation(T,nH,nHI,nH2,Z,Tdust);
//    if(isnan(Lvibrot)||isnan(Lh2coll)||isnan(Hphotdiss)
//            ||isnan(Hphotpump)||isnan(Hh2form)){
//        fprintf(stderr,"Lvibrot=%e\n",Lvibrot);
//        fprintf(stderr,"Lh2coll=%e\n",Lh2coll);
//        fprintf(stderr,"Hphotdiss=%e\n",Hphotdiss);
//        fprintf(stderr,"Hphotpump=%e\n",Hphotpump);
//        fprintf(stderr,"Hh2form=%e\n",Hh2form);
//        fprintf(stderr,"nH=%e nHI=%e nH2=%e\n",nH,nHI,nH2);
//        exit(1);
//    }
//    // [erg/cc/s]
//    return Lvibrot + Lh2coll - Hphotdiss - Hphotpump - Hh2form;
//}

//// Cooling associating with H2 rovibrational lines [erg/cc/s]
// eq.(A7) of Galli & Palla, A&A, 335, 403, (1998)
double CoolingH2RovibrationalLineLowDensityLimitGalliPalla1998(double T, double nHI, double nH2){
    double logT = log10(T);
    // [erg/s/cc]
    return nHI*nH2*pow(10.0,-103.0 + logT*(97.59 + logT*(-48.05 + logT*(10.80 - logT*0.9032))));
}

// Table 8 of Glover & Abel, MNRAS, 388, 1627, (2008)
static double a_wH_10_T_100[]     = {-16.818342, 37.383713,  58.145166, 48.656103,  20.159831,   3.8479610};
static double a_wH_100_T_1000[]   = {-24.311209, 3.5692468, -11.332860,-27.850082, -21.328264,  -4.2519023};
static double a_wH_1000_T_6000[]  = {-24.311209, 4.6450521, -3.7209846, 5.9369081, -5.5108047,   1.5538288};
static double a_wH2_100_T_6000[]  = {-23.962112, 2.0943374,-0.77151436,0.43693353,-0.14913216,-0.033638326};
static double a_wHe_10_T_6000[]   = {-23.689237, 2.1892372,-0.81520438,0.29036281,-0.16596184,  0.19191375};
static double a_wHp_10_T_10000[]  = {-21.716699, 1.3865783,-0.37915285,0.11453688,-0.23214154, 0.058538864};
static double a_wel_10_T_200[]    = {-34.286155,-48.537163, -77.121176,-51.352459, -15.169160, -0.98120322};
static double a_wel_200_T_10000[] = {-22.190316, 1.5728955,-0.21335100,0.96149759,-0.91023195,  0.13749749};
static double CoolingH2collisionLowDensityLimitBody(double logT3, double a[]){
    double logL = 0.e0;
    for(int k=0;k<=5;k++)    logL += a[k]*pow(logT3,k);
    return pow(10.0,logL);  // [erg*cc/s]
}
//// Glover & Abel, MNRAS, 388, 1627, (2008)
double CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(double T, double nHII, double nHI, double nH2,
        double nHe, double ne){
    const double logT3 = log10(T)-3.0;
    // Extend each Table 8 fit continuously with its upper endpoint value.
    const double neutralLogT3 = T > 6000.0 ? log10(6000.0)-3.0 : logT3;
    const double chargedLogT3 = T > 10000.0 ? 1.0 : logT3;

    double L_H = 0.0, L_He = 0.0, L_Hp = 0.0, L_el = 0.0;
    const double L_H2 = CoolingH2collisionLowDensityLimitBody(neutralLogT3,a_wH2_100_T_6000);
    if(T > 10.0){
        if(T <= 100.0)        L_H = CoolingH2collisionLowDensityLimitBody(logT3,a_wH_10_T_100);
        else if(T <= 1000.0)  L_H = CoolingH2collisionLowDensityLimitBody(logT3,a_wH_100_T_1000);
        else                  L_H = CoolingH2collisionLowDensityLimitBody(neutralLogT3,a_wH_1000_T_6000);
        L_He = CoolingH2collisionLowDensityLimitBody(neutralLogT3,a_wHe_10_T_6000);
        L_Hp = CoolingH2collisionLowDensityLimitBody(chargedLogT3,a_wHp_10_T_10000);
        L_el = T <= 200.0 ?
               CoolingH2collisionLowDensityLimitBody(logT3,a_wel_10_T_200) :
               CoolingH2collisionLowDensityLimitBody(chargedLogT3,a_wel_200_T_10000);
    }
    // total H2 cooling [erg/s/cc]
    return nH2*(nHI*L_H + nH2*L_H2 + nHe*L_He + nHII*L_Hp + ne*L_el);
}

// collisional de-excitation rates for HI and H2 exciting the J to J-2 transition
// ref. eq(6.31)(6.32) in Hollenbach & McKee (1979)
//      eq(2.12) in Hollenbach & McKee (1989)
// ref. A 1.5 of Koyama & Inutsuka (2000)
//      -Hollenbach & McKee (1979)
//      -Galli & Palla (1998)
static double Ev10_H2 =  5860.0; // [K]
static double Ev20_H2 = 11720.0; // [K] 
static double Er20_H2 =   510.0; // [K]: 0-0S(0), 28.2um 
static double Er31_H2 =   850.0; // [K]: 0-0S(1), 17.0um
//static double Bk      =   85.25; // [K]   
static double fH2para  = 0.25;
static double fH2ortho = 0.75;
static double gamJ(double T, int J, int id){   // [cm^3/s]
    if(id==0){  // HI
        //// (6.31) in Hollenbach & McKee (1979)
        //double T3 = 1.e-3*T;
        //return ((1.e-11*sqrt(T3))/(1+60.0*pow(T3,-4))+1.e-12*T3)
        //    *(0.33+0.9*exp(-SQ((J-3.5)/0.9)));

        //// eq(2.12) in Hollenbach & McKee (1989)
        // Draine et al. (1983); Flower et al. (1986)
        //double x = 2*Bk*(2*J-1)/T;
        //return 4.6e-12*(2*J-3)*sqrt(T)*sqrt(1+x)
        //    *exp(-(10.0*Bk*(2*J-1)/(T+Bk*J*(J+1)))+0.1187*(4*J-2));

        //// eq(A4) in Galli & Palla (1998)
        if(J==2){    // J=2->0 transition
            double a[] = {2.93e-14,1.21e-15,2.16e-19,1.32e-21};
            return a[0]+a[1]*T+a[2]*SQ(T)+a[3]*CUBE(T);
        }
        else if(J==3){
            double a[] = {8.34e-14,5.97e-16,7.76e-19,1.72e-21};
            return a[0]+a[1]*T+a[2]*SQ(T)+a[3]*CUBE(T);
        }
        else if(J==4){
            double a[] = {7.54e-14,2.65e-16,2.14e-19,2.65e-21};
            return a[0]+a[1]*T+a[2]*SQ(T)+a[3]*CUBE(T);
        }
        else if(J==5){
            double a[] = {2.95e-14,1.21e-16,5.53e-19,2.51e-21};
            return a[0]+a[1]*T+a[2]*SQ(T)+a[3]*CUBE(T);
        }
        else
            exit(1);
    }
    else if(id==1){ // H2
        //// (6.32) in Hollenbach & McKee (1979)
        return (3.3e-12+6.6e-15*T)*(0.276*SQ(J)*exp(-pow(J/3.18,1.7)));
    }
    else{
        exit(1);
    }
}
double CoolingH2RovibrationalLineHollenbachMcKee1979(double T, double nHI, double nH2){
   
    double T3 = 1.e-3*T;
    //if(T<=10.0)  return 0.0;
    //if(T3>5.0)  return TINY;
    double invT3 = 1.e0/T3;
    double sqrtT = sqrt(T);
    double invT  = 1.e0/T;
    
    //////////
    // vibratio-rotation transition of H2
    // v_u->v_l S(J_l) if J_l = J_u - 2
    // v_u->v_l Q(J_l) if J_l = J_u
    // v_u->v_l O(J_l) if J_l = J_u + 2
    //////////

    // cooling rate per H,H2 in rotational LTE and vibrational LTE
    //  n*L_r,v_LTE = A_J*(g_J/g_J')*exp(-DeltaE/kT)*DeltaE
    //  where g_J = 2*J+1 is the statistical weight, J denotes the total
    //  angular momentum quantum number of the rotational energy level,
    //  A_J is the Einstein coefficient for the J->J' transition at the
    //  same vibrational energy level or between two vibrational energy
    //  levels, and DeltaE is the energy difference between E_J and E_J'. 
    //  eq.(6.37),(6.38): [erg/s]
    double nLr_LTE = (9.5e-22*pow(T3,3.76)/(1.e0+0.12*pow(T3,2.1)))*exp(-CUBE(0.13*invT3)) + 3.0e-24*exp(-0.51*invT3);
    double nLv_LTE = 6.7e-19*exp(-5.86*invT3) + 1.6e-18*exp(-11.7*invT3); 
    //nLr_LTE = nLr_LTE;
    //nLv_LTE = nLv_LTE;

    // collisional de-excitaiion rate for vibrational levels,
    // gamma_vv'(T) [cm^3/s]
    //  eq.(6.29) in Hollenbach & McKee (1979)
    //  eq.(2.8) in Hollenbach & McKee (1989)
    double gam10[2], gam20[2];
    // w/ H
    gam10[0] = 1.0e-12*sqrtT*exp(-1000.0*invT); 
    gam20[0] = 1.6e-12*sqrtT*exp(-SQ(400.0*invT));
    // w/ H2
    gam10[1] = 1.4e-12*sqrtT*exp(-18100.0/(T+1200.0));
    gam20[1] = 0.e0;

    // sum
    double n[] = {nHI,nH2};
    double Lvr[] = {0.e0,0.e0};
    for(int i=0; i<2; i++){
        // cooling coefficients in rotational LTE and vibrational LTE
        //  eq.(6.37),(6.38): [erg cm^3/s] of HM79
        //double Lr_LTE = nLr_LTE/n[i]; // [erg cm^3/s] 
        //double Lv_LTE = nLv_LTE/n[i]; // [erg cm^3/s]

        // low-density limit (n->0) cooling coefficients [erg*cc/s]:
        //  L_r,v_n->0 = C_J*(g_J/g_J')*exp(-DeltaE/kT)*DeltaE
        //  where C_J is the collisional de-excitation coefficients with
        //  H or H2 for the J->J' transition in the same vibrational
        //  level or between v=i and v=j, which are estimated by the
        //  fitting formula given in HM79 and HM89.
        //  eq.(6.40),(6.41): [erg cm^3/s] of HM79
        //  Lr_n0: J=2->0, J=3->1 at v=0, i.e., 0-0S(0), 0-0S(1)
        //  Lv_n0: v=1->0, v=2->0         i.e., 1-0 
        double Lr_n0 = fH2para *5.00*gamJ(T,2,i)*exp(-Er20_H2*invT)*Er20_H2
                     + fH2ortho*2.33*gamJ(T,3,i)*exp(-Er31_H2*invT)*Er31_H2;    
        double Lv_n0 = gam10[i]*exp(-Ev10_H2*invT)*Ev10_H2 
                     + gam20[i]*exp(-Ev20_H2*invT)*Ev20_H2;
        Lr_n0 *= BOLTZMANN_CONSTANT_CGS;
        Lv_n0 *= BOLTZMANN_CONSTANT_CGS;

        // eq.(6.39)
        //double ncr_rot_over_ni = Lr_LTE/Lr_n0;
        //double ncr_vir_over_ni = Lv_LTE/Lv_n0;

        // eq.(6.36) in Hollenbach & McKee (1979):[erg cm^3/s]
        //Lvr[i] = Lr_LTE/(1.e0+ncr_rot_over_ni) + Lv_LTE/(1.e0+ncr_vir_over_ni); 
        Lvr[i] = nLr_LTE*Lr_n0/(n[i]*Lr_n0+nLr_LTE) + nLv_LTE*Lv_n0/(n[i]*Lv_n0+nLv_LTE);
        if(isnan(Lvr[i])){
            fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
            fprintf(stderr,"Lvr[%d]=%e n=%e T=%e\n",i,Lvr[i],n[i],T);
            fprintf(stderr,"nLr_LTE=%e nLv_LTE=%e Lr_n0=%e Lv_n0=%e gam10=%e gam20=%e\n",
                    nLr_LTE,nLv_LTE,Lr_n0,Lv_n0,gam10[i],gam20[i]);
            exit(1);
        }
    }

    // eq.(6.35) in Hollenbach & McKee (1979):[erg/cc/s]
    return nH2*(nHI*Lvr[0] + nH2*Lvr[1]);
    //return nH2*(nHI*Lvr[0] + nH2*Lvr[1])*MAX(0.0,1.e0-exp(T3-3.0));
}

#if 1
//// Neufeld & Kaufman (1993) and Neufeld et al. (1995)
// ref. C.4 of Meijerink & Spaans (2005)
static double L0_vib_NK93(double T){
    // Table 5 of Neufeld & Kaufman 1993
    // [erg cm3/s]
    return 1.19e-24*sqrt(T)*exp(-18100/(T+1190))*exp(-5987/T);
}
static double LLTE_vid_NK93(double T){
    // Table 5 of Neufeld & Kaufman 1993
    // [erg/s]
    return 1.10e-18*exp(-6744/T);
}
double CoolingH2RovibrationalLineNeufeldKaufman1993(double T, double nHI, double nH2){
    // 100 < T < 4000 K; ~N <= 10^19 /cm2/(km/s)

    // vibrational line cooling (Neufeld & Kaufman 1993)
    // only first two terms of eq.(5)
    double invL0_vir   = 1.e0/L0_vib_NK93(T);
    double invLLTE_vir = 1.e0/LLTE_vid_NK93(T);
    double invL_vir    = invL0_vir + nH2*invLLTE_vir;
   
    // [erg/cc/s]
    // TODO: Replace the placeholder 1.e0 with a rotational cooling coefficient in matching units.
    return SQ(nH2)*(1.e0/invL_vir + 1.e0);
}
#endif

//// Le Bourlot et al. (1999)
//   http://ccp7.dur.ac.uk/cooling_by_h2/
double CoolingH2RovibrationalLineLeBourlot1999(double T, double nHI, double nH2){
    fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
    exit(1);

    // cooling rate per H2 [erg/s]
    //  1 < nH < 10^3 /cc
    //  100 < T < 10^4 K
    //  10^-8 < nH/nH2 < 10^6
    //  0.1 < orttho/para < 3
    double W_H2 = 0.e0;  

    return W_H2;
}

//// Cooling via Rotational & Vibrational Lines of H2 [erg/cc/s]
double CoolingH2RovibrationalLine(StructEQCData *PhydroPP){
    double nHII = PhydroPP->nion[ipHYDROGEN];
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nH2  = PhydroPP->nmol [ipH2];
    if(nH2==0.e0)   return 0.e0;
    double nHe  = PhydroPP->n[ipHELIUM];
    double ne   = PhydroPP->nel;
    double T    = PhydroPP->T;
    //return CoolingH2RovibrationalLineLowDensityLimitGalliPalla1998(T,nHI,nH2);
    return CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(T,nHII,nHI,nH2,nHe,ne);    // default
    //return CoolingH2RovibrationalLineHollenbachMcKee1979(T,nHI,nH2);
    //return CoolingH2RovibrationalLineNeufeldKaufman1993(T,nHI,nH2);
    //return CoolingH2RovibrationalLineLeBourlot1999(T,nHI,nH2);
}

//// Cooling associating with H2 dissociation by collisions [erg/cc/s]
// Following Hollenbach & McKee (1979), collisional dissociation absorb
// the same energy as the biding energy, that is, 4.48 eV per H2.
// see p.587 of Hollenbach & McKee (1979)
static double CoolingH2CollisionalDissociationBody(double nHI, double nH2, double nH, double ne, double T){
    if(nH2==0.e0)   return 0.e0;
    return 4.48*(EV_TO_ERG)*nH2*H2CollisionalDissociationRate(T,nH,nHI,nH2,ne);
    //double dissHI = nHI*HIcollisionCoefficient(T,nH,nHI,nH2);
    //double dissH2 = nH2*H2collisionCoefficient(T,nH,nHI,nH2);
    //return 4.48*(EV_TO_ERG)*nH2*(dissHI + dissH2);
}
double CoolingH2CollisionalDissociation(StructEQCData *PhydroPP){
    double nHI = PhydroPP->natom[ipHYDROGEN];
    double nH2 = PhydroPP->nmol [ipH2];
    if(nH2==0.e0)   return 0.e0;
    double nH  = PhydroPP->n[ipHYDROGEN];
    double ne  = PhydroPP->nel;
    double T   = PhydroPP->T;
    return CoolingH2CollisionalDissociationBody(nHI,nH2,nH,ne,T);
    //double diss = H2CollisionalDissociationRate(T,nH,nHI,nH2,0.0);
    //return 4.48*(EV_TO_ERG)*nH2*diss;
}

//// Heating associating with photodissociation of H2 by FUV [erg/cc/s]
// Following Black & Dalgarno (1977) and Glover & MacLow (2007), we assume that
// each photodissociation deposits 0.4 eV (= 6.4 x 10^-13 erg) of heat into the gas.
// see also  eq(3.18)(3.19) of Teilens' textbook (2005)
double HeatingH2Photodissociate(StructEQCData *PhydroPP){
    double nH2 = PhydroPP->nmol[ipH2];
    double Nt  = PhydroPP->Nt;
    double N2  = Nt*nH2/PhydroPP->n[ipHYDROGEN]; 
    double G0  = PhydroPP->G0;
    double Zd  = PhydroPP->Zd;
    //double b5  = 12.90*sqrt(1.e-4*PhydroPP->T);    // eq(6.34) of Draine's textbook
    double b5  = sqrt(SQ(12.90*sqrt(1.e-4*PhydroPP->T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
    // see eq.(44) in Glover & MacLow, ApJ, 169, 239, (2007)
    return (6.4e-13)*H2PhotoDissociationRate(G0,Nt,N2,b5,Zd)*nH2;
}

//// Heating associating with excitaion of H2 by FUV [erg/cc/s]
// Radiative pupmping rate that is 8.5 times larger than the photodissociation rate (Draine & Bertoldi 1996).
// Each excitation transfers an average of 2*n/(n+ncr) [eV] = 3.2e-12*n/(n+ncr) [erg] to the gas (Appendix A of Burton et al. 1990).
// eq.(6.46) of Hollenbach & McKee (1979)
double HeatingH2PhotomPumping(StructEQCData *PhydroPP){
    double nH  = PhydroPP->n[ipHYDROGEN];
    double nHI = PhydroPP->natom[ipHYDROGEN];
    double nH2 = PhydroPP->nmol[ipH2];
    double Nt  = PhydroPP->Nt;
    double N2  = Nt*nH2/nH;
    double G0  = PhydroPP->G0;
    double Zd  = PhydroPP->Zd;
    double T   = PhydroPP->T;
    //double b5  = 12.90*sqrt(1.e-4*T);    // eq(6.34) of Draine's textbook
    double b5  = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
    double yHI = nHI/nH;
    double yH2 = nH2/nH;
    // see pp.247 in Glover & MacLow, ApJ, 169, 239, (2007)
    return (FUVPumpingHeatingFactorDB96)*(3.2e-12)*H2PhotoDissociationRate(G0,Nt,N2,b5,Zd)*nH2*nH/(nH+CriticalDensity(T,yHI,yH2));
    //return (3.2e-12)*H2dissociationRate(G0,Nt,N2,b5,Z)*nH2*nH/(nH+CriticalDensity(T,yHI,yH2));
    
    // eq(A5) of Koyama & Inustuka (2000)
    //return (31.7e-11)*H2dissociationRate(G0,Nt,N2,Z)*nH2*nH/(nH+CriticalDensity(T,yHI,yH2));
}

//// Heating associating with H2 formation on grains [erg/cc/s]
// Hollenbach & McKee (1989)
// eq.(6.43) of Hollenbach & McKee (1979)
double HeatingH2Formation(StructEQCData *PhydroPP){
    double nH  = PhydroPP->n[ipHYDROGEN];
    double nHI = PhydroPP->natom[ipHYDROGEN];
    double nH2 = PhydroPP->nmol[ipH2];
    //double G0  = PhydroPP->G0;
    double Zd  = PhydroPP->Zd;
    double T   = PhydroPP->T;
    double Td  = EffectiveDustTemperature(PhydroPP);
    double yHI  = nHI/nH;
    double yH2  = nH2/nH;
    // see eq.(45) in Glover & MacLow, ApJ, 169, 239, (2007)
    return (7.2e-12)*H2FormationRate(T,Zd,Td)*nHI*SQ(nH)/(nH+CriticalDensity(T,yHI,yH2));
    
    // eq(A5) of Koyama & Inustuka (2000)
    //return EV_TO_ERG*H2formationRate(T,Z,Tdust)*nHI*nH*(0.2+4.2*nH/(nH+CriticalDensity(T,yHI,yH2)));
}


//extern double DustShieldingFactor(double Ntot, double Z, double gamma);

////// Self-shielding factor for CO photodissociation
//// see eq.(C3) of Tielens & Hollenbach (1985)
//double COSelfShieldingFactor(double Nco, double b5){
//    double tau   = 1.1e-14*Nco/b5;
//    double alpha = 0.0;
//    double beta  = 0.0; // eq.(A10)(A11)
//    return (1.e0+alpha*beta)/(1.e0+alpha);
//}
    
//// CO self-shielding and shielding of CO by H2 L-W lines.
// ref. Table 10, 11 in Lee et al. A&A, 311, 690, (1996)
static double f_co_self(double logNco){
    if(logNco<TableNco[0])                          return 0.e0;
    if(logNco>TableNco[NumberTableShieldingCO-1])   return -7.0;
    return gsl_spline_eval(logShieldCOAccess,logNco,NULL);
}
static double f_co_H2(double logNh2){
    if(logNh2<TableNh2[0])                          return 0.e0;
    const int last = NumberTableShieldingH2-1;
    if(logNh2>TableNh2[last]){
        const double logRatio = logNh2-TableNh2[last];
        return TableShieldingH2[last] - h2ShieldingTailPower*logRatio
               - (h2ShieldingTailOpticalDepth/log(10.0))*expm1(logRatio*log(10.0));
    }
    return gsl_spline_eval(logShieldH2Access,logNh2,NULL);
}

int CEQCOShieldingPointCount(void){
    return NumberTableShieldingCO;
}

void CEQCOShieldingPoint(int index, double *logColumn, double *logFactor){
    *logColumn = TableNco[index];
    *logFactor = TableShieldingCO[index];
}
double COShieldingFactor(double Nh2, double Nco, double b5){
    //double f_self = COSelfShieldingFactor(Nco,b5);
    //double f_H2   = H2SelfShieldingFactor(N2,b5);
    //return f_self*f_H2;
    double f_self = f_co_self(log10(MAX(Nco,1.e0)));
    double f_H2   = f_co_H2(log10(MAX(Nh2,1.e0)));
    return pow(10.0,f_self+f_H2);
}

///// CO dissociation by FUV 
double COPhotoDissociationRate(double G0, double Ntot, double Zd){
    //// eq.(3) of  Glover & Clark, MNRAS, 421, 116, (2012)
    return 2.e-10*G0*DustShieldingFactor(Ntot,Zd,2.5); // [1/s]
    // Visser et al. (2009) 
    //return 2.6e-10*G0*DustShieldingFactor(Ntot,Z,3.57); // [1/s]
}

////// photodiossociation rate of CO [cm^3/s]
//// see eq(C2) and Tables 5 of Tielens & Hollenbach (1985)
//double COPhotodissociationRate(double G0dust){
//    return MAX(CHEM_TINY,(1.4e-11)*pow(G0dust,1.23));
//    //return MAX(CHEM_TINY,(1.4e-11)*G0dust);
//    //return (1.4e-11)*G0*exp(-3.2*Av);
//}


#ifdef CHEMISTRY_CO_NELSON_LANGER_1997
// ref. Nelson & Langer, ApJ, 482, 796, (1997)
//      see footnote 1 of Glover & Clark, MNRAS, 421, 116, (2012)
///// CO formation rate 
#define k0_CO (5.e-16)  // [cm^3/s]: rate coefficient for the formation of the intermediate CHx or dacial 
#define k1    (5.e-10)  // [cm^3/s]: rate coefficient for the formation of CO from O+CHx
static double CHxPhotodissociationRate(double G0, double Ntot, double Zd){
    // eq.(5) of Glover & Clark, MNRAS, 421, 116, (2012)
    return 5.e-10*G0*DustShieldingFactor(Ntot,Zd,2.5); // [1/s]
}
double COFormationRate(double nH, double xO, double G0, double Ntot, double Zd){
    // eq.(4) of Glover & Clark, MNRAS, 421, 116, (2012)
    double beta = k1*xO/(k1*xO + CHxPhotodissociationRate(G0,Ntot,Zd)/nH);
    return k0_CO*beta; // [cm^3/s]
}
#endif

#if 0
//#ifdef CHEMISTRY_CO_KETO_CASELLI_2008
// equilibrium carbon chemistry
// ref. Keto & Caselli, ApJ, 683, 238, (2008)
//      see footnote 2 of Glover & Clark, MNRAS, 421, 116, (2012)
static double COFormationRateCH2(double nH2, double T){
    // CH2 formation rate coefficient from CII and H2 (via radiative association to form CH2+)
    //double kRA = 6.e-16; // ref.  Table 12 of Tielens & Hollenbach (1985)
    double Cp = EQCLibRunParameters.CpForCO;
    double kRA = 4.e-16*pow(T/300.0,-0.2);  // k148 in Table B1 of Glover et al. (2010)
    return kRA*Cp*nH2; // [/s] 
}

static double COtoCIIratio(double nH2, double T, double G0, double Ntot, double N2, double Nco, double Z){
    // eq.(14) in Glover & Clark, MNRAS, 421, 116, (2012)
    double b_turb = 0.e0;    // km/s
    double b5   = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));    // eq(6.34) of Draine's textbook
    double r_CO = COPhotoDissociationRate(G0,Ntot,Z)*COShieldingFactor(N2,Nco,b5);
    return COFormationRateCH2(nH2,T)/r_CO;
}

static double CItoCIIratio(double nH2, double T, double G0, double Ntot, double Z){
    // eq.(15) in Glover & Clark, MNRAS, 421, 116, (2012)
    double r_C = 3.1e-10*G0*DustShieldingFactor(Ntot,Z,3.0); // R171 in Table B2 of Glover et al. 2010
    return COFormationRateCH2(nH2,T)/r_C;
}

// eq(11)-(13) in Glover & Clark, MNRAS, 421, 116, (2012)
void EquilibriumCarbonFractionalAbundance(double xC, double nH, double ne, double nH2, double T, double G0, double X0, double Ntot, double N2, 
        double Z, double *fCII, double *fCI, double *fCO){
    if(EQCLibRunParameters.UseCOFormation==true){
        double Nco = Ntot*(*fCO)*xC;
        double COtoCII = COtoCIIratio(nH2,T,G0,Ntot,N2,Nco,Z);  // eq.(14)
        double CItoCII = CItoCIIratio(nH2,T,G0,Ntot,Z);         // eq.(15)
        double CIItoCtotal = 1.e0/(1.0 + COtoCII + CItoCII);    // eq.(11)
        *fCII = CIItoCtotal;
        *fCI  = CItoCII*CIItoCtotal;
        *fCO  = MAX(0.0,1.e0 - *fCII - *fCI);
    } else{
        double xe         = ne/nH;
        double Av         = 0.0; //VisualExtinction(Ntot,Z);    for Wolfire et al. 1995
        double dnCI_phti  = CIFUVPhotoionizationRate(G0,Av);
        double dnCI_coli  = ne*CICollisionalIonizationRate(T);
        double dnCI_cr    = X0*CosmicRayIonizationRateCI(xe);
        double dnCI_xr    = X0*SoftXrayIonizationRateCI(xe,Ntot);
        double dnCII_rec  = ne*CIIRecombinationRate(T);
        double dnCII_dust = nH*CIIRecombinationRateOnDust(T,G0,ne,Z); 
        double CIItoCI   = (dnCI_phti + dnCI_coli + dnCI_cr + dnCI_xr)/(dnCII_rec + dnCII_dust);
        *fCII = CIItoCI/(1.0+CIItoCI);
        *fCI  = MAX(0.0,1.0 - (*fCII));
        *fCO  = 0.e0;
    }
    if(*fCII<0||*fCI<0||*fCO<0){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"  fCII=%e fCI=%e fCO=%e\n",*fCII,*fCI,*fCO);
        exit(1);
    }
    return ;
}
//void EquilibriumCarbonFractionalAbundanceIterative(int IterMax, double xC, double nH2, double T, double G0, double Ntot, double N2, 
//        double Z, double *fCII, double *fCI, double *fCO){
//    double FracCII = 0.9, FracCI = 0.01, FracCO = 1.0-FracCII-FracCI;
//    for(int itr=0; itr<IterMax; itr++){
//        double Ncoi = Ntot*(FracCO)*xC;
//        EquilibriumCarbonFractionalAbundance(xC,nH,ne,nH2,T,G0,Ntot,N2,Z,&FracCII,&FracCI,&FracCO);
//        double Ncof = Ntot*(FracCO)*xC;
//        if(itr>100&&fabs(Ncof-Ncoi)<1.e-3*Ncoi)  break;
//    }
//    *fCII = FracCII;
//    *fCI  = FracCI ;
//    *fCO  = FracCO ;
//    return ;
//}
//#endif
#endif

// see A 1.6 of Koyama & Inustuka, ApJ, 532, 980, (2000)
// ref. McKee et al. (1982)
//      Hollenbach & McKee (1989)
#define A0      (9.7e-8)    // [1/s]
#define E0      (2.76*BOLTZMANN_CONSTANT_CGS)      // [erg]
#define dE10    (3080.0*BOLTZMANN_CONSTANT_CGS)    // [erg]
// 12CO rotational cooling (McKee et al. 1982)
// eq.(6.23) in Hollenbach & McKee (1979)
// eq.(5.2a) in McKee et al. (1982)
static double CoolingRotationCOLinesHollenbachMcKee(double T, double nCO, double nH, double nHI, double nH2){
    double ncr = (3.3e+6)*pow(1.e-3*T,0.75);    // [/cc]
    double ncrOvernH = ncr/nH;
    double rot = 4.0*SQ(BOLTZMANN_CONSTANT_CGS*T)*A0/(E0*(1.0+ncrOvernH+1.5*sqrt(ncrOvernH)));  //[erg/s]
    return nCO*rot;
}
// CO vibrotational cooling (Hollenbach & McKee 1989)
// eq.(2.13a,b) in Hollenbach & McKee (1989)  
static double CoolingVibrationCOLinesHollenbachMcKee(double T, double nCO, double nH, double nHI, double nH2){
    double invT  = 1.e0/T;
    double fact  = exp(-3080.0*invT); 
    double gamHI = (3.0e-12)*sqrt(T)*exp(-pow(2000.0*invT,3.43))*fact;  // [cm3/s]
    double gamH2 = (4.3e-14)*T*exp(-pow(3.14e+5*invT,0.333))*fact;      // [cm3/s]
    double vib   = (gamHI*nHI+gamH2*nH2)*dE10;  // [erg/s]
    return nCO*vib;
}
static double CoolingCORovibrationalLineHollenbachMcKee(double T, double nCO, double nH, double nHI, double nH2){
    double Lrot = CoolingRotationCOLinesHollenbachMcKee(T,nCO,nH,nHI,nH2);
    double Lvib = CoolingVibrationCOLinesHollenbachMcKee(T,nCO,nH,nHI,nH2);
    return Lrot + Lvib;
}

//// CO cooling [erg/cc/s]
double CoolingCORovibrationalLine(StructEQCData *PhydroPP){
    double nH  = PhydroPP->n[ipHYDROGEN];
    double nHI = PhydroPP->natom[ipHYDROGEN];
    double nH2 = PhydroPP->nmol [ipH2];
    double nCO = PhydroPP->nmol [ipCO];
    double T   = PhydroPP->T;
    if(nCO<=0.0) return 0.e0;
    return CoolingCORovibrationalLineHollenbachMcKee(T,nCO,nH,nHI,nH2);
}


#if 0
//// Total Cooling Rate For Oxygen Chemistry
double TotalOxygenCoolingRate(double T, double ne, double nHI, double nH2, double nH2O, double Nh2o){ 
#ifdef CHEMISTRY_H2O
    //////////
    // under construction....
    //////////
    exit(1);
    double L_H2O = CoolingH2ORovibrationalLines(T,ne,nHI,nH2,nH2O,Nh2o);
    double L_OH  = CoolingOHRotationalLines();

    return L_H2O + L_OH;
#else
    return 0.0;
#endif
}

//// H2O vibrational cooling
// ref. Table 2 and 5 of Neufeld & Kaufman, ApJ, 418, 263, (1993)
//      rotational cooling : eq.(11,13,14) of Glover et al. (2010)
//      vibrational cooling: eq.(16,20,21,22) of Glover et al. (2010)
//      see Appendix C.4 of Meijerink & Spaans (2005)
#define NumberofTmpLlte  (6) 
#define NumberoflogNLlte (8) 
static double TableTmp[NumberofTmpLlte] = { // Table 5 in Neufeld & Kaufman (1993)   
    1.e+2,2.e+2,4.e+2,1.e+3,2.e+3,4.e+3
};
static double TableN[NumberoflogNLlte] = {
    13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0
};
static double TableLlteVib[NumberoflogNLlte][NumberofTmpLlte] = { 
    {10.98,11.05,11.07,10.88,10.67,10.67},  // logN=13
    {10.99,11.05,11.07,10.88,10.67,10.67},  // logN=14
    {11.06,11.09,11.09,10.88,10.67,10.67},  // logN=15
    {11.36,11.31,11.23,10.93,10.69,10.68},  // logN=16
    {11.91,11.82,11.65,11.15,10.79,10.72},  // logN=17
    {12.59,12.50,12.28,11.61,11.11,10.90},  // logN=18
    {13.29,13.20,12.99,12.25,11.65,11.30},  // logN=19
    {13.90,13.90,13.76,13.01,12.34,11.88}   // logN=20
};
gsl_interp_accel *TmpAccess;
gsl_spline       *vibLlteAccess;
static void InitializeH2ORoVibrationalCoolingCoeffient(void){

    //NumberofHIcol = sizeof(gamHIcol)/sizeof(double);
    //if(MPIGetMyID() == MPI_ROOT_RANK)
    //    dprintlmpi(NumberofHIcol);

    //for(int i=0; i<NumberofHIcol; i++){
    //    TmpHIcol[i] = log10(TmpHIcol[i]);
    //    gamHIcol[i] = log10(gamHIcol[i]);
    //}
    // 
    //TmpHIcolAccess = gsl_interp_accel_alloc();
    //gamHIcolAccess = gsl_spline_alloc(gsl_interp_linear,NumberofHIcol);
    //gsl_spline_init(gamHIcolAccess,TmpHIcol,gamHIcol,NumberofHIcol);
    //
    //if(MPIGetMyID() == MPI_ROOT_RANK){
    //    FILE *output = fopen("HI-H2CollisionDestructionCoefficiency.dat","w");
    //    for(double logT=1.0; logT<=5.0; logT+=0.1)
    //        fprintf(output,"%e %e\n",logT,HIcollisionCoefficientMartin1998(logT)); 
    //    fclose(output);
    //}
    //MPI_Barrier(MPI_COMM_WORLD);
}
static double LlteForH2Ovib(double T, double logN){
    double val = 0.0;   //gsl_spline_eval(gamHIcolAccess,logT,TmpHIcolAccess);
    return pow(10.0,-val)*exp(-2325.0/T);
}
double CoolingH2ORovibrationalLines(double T, double ne, double nHI, double nH2, double nH2O, double N_H2O){
   
    //// under construction....
    
    double invT = 1.e0/T;

    /// rotational cooling
    // eq(14,15) of Glover et al. (2010)
    double ke  = pow(10.0,-8.020+15.749*pow(T,-1.6667)-47.137*pow(T,-0.33333)+76.648*sqrt(invT)-60.191*pow(T,-0.66667));
    double kH2 = 7.4e-12*sqrt(T);   // [cm^3/s]
    double neffH2Orot = nH2 + 10.0*nHI + (ke/kH2)*ne;
    // Table 5 of Neufeld & Kaufman (1993)
    double L0rot = 1.0;
    double invL0rot = 1.e0/L0rot;
    double Llterot = 1.0;
    double invLlterot = 1.e0/Llterot;
    double n12 = 1.0;
    double a = 1.0;
    // eq (11) of Glover et al. (2010)
    double invLrot = invL0rot + neffH2Orot*invLlterot +invL0rot*pow(neffH2Orot/n12,a)*(1.e0-n12*L0rot/Llterot);

    /// vibrational cooling
    // eq(20,21,22) of Glover et al. (2010)
    double Le = 2.6e-6*sqrt(invT)*exp(-2325.0*invT);
    double LH = 6.4e-15*exp(-47.5*pow(invT,0.33333))*exp(-2325.0*invT); 
    double neffH2Ovib = nH2 + 10.0*nHI + (Le/LH)*ne;
    // Table 5 of Neufeld & Kaufman (1993)
    double L0vib = 1.03e-26*exp(-47.5*pow(invT,0.33333))*exp(-2325.0*invT);    // [erg cm^3/s] 
    double invL0vib   = 1.e0/L0vib;
    double invLltevib = 1.e0/LlteForH2Ovib(T,log10(N_H2O));
    // eq.(16) of Glover et al. (2010)
    double invLvib = invL0vib + neffH2Ovib*invLltevib;
   
    // eq.(C.2) of Meijerink & Spaans (2005)
    return (1.e0/invLrot)*nH2O*neffH2Orot + (1.e0/invLvib)*nH2O*neffH2Ovib; // [erg/cc/s]
}


//// OH cooling
// ref. Pavlovski et al. (2002)
//      Hollenbach & McKee (1979)
//      see $2.3.3 of Glover et al. (2010)
double CoolingOHRotationalLines(void){
    return 0.0;
}
#endif
