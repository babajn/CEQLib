#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "FineStructureLines.h"
#include "CosmicRay.h"
#include "PhotoelectronDeposition.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
//#include <gsl/gsl_roots.h>

// ref. Stahler & Palla (2004) "The Formation of Stars"
#define E_ELECTRON  (35.0)              // average energy of primary electron (Spitzer 1978) [eV] 
#define zetaHI      (1.8e-17)           // ionization rate [1/s] (McKee 1994; Wolfire et al. 1995)
#define zetaDI      (zetaHI)            // ionization rate [1/s] (Le Teuff et al. 2000; Glover et al. 2010)
#define zetaHeI     (1.1*zetaHI)        // ionization rate [1/s] (Le Teuff et al. 2000; Glover et al. 2010) 
#define zetaH2      (2.0*zetaHI)        // ionization rate [1/s] (Le Teuff et al. 2000; Glover et al. 2010)
#define zetaCI      (3.8*zetaHI)        // ionization rate [1/s] (Le Teuff et al. 2000; Glover et al. 2010)
#define zetaOI      (5.7*zetaHI)        // ionization rate [1/s] (Le Teuff et al. 2000; Glover et al. 2010)
#define zetaSiI     (6.5*zetaHI)        // ionization rate [1/s] (Table 4 in Glover & Jappsen 2007)
#define zetaSiII    (2.5*zetaHI)        // ionization rate [1/s] (Table 4 in Glover & Jappsen 2007)
#define zetaCO      (6.5*zetaHI)        // ionization rate [1/s] (Le Teuff et al. 2000; Glover et al. 2010)
#define dE_HI       (6.0*EV_TO_ERG)     // deposite energy [erg]
#define dE_H2       (7.0*EV_TO_ERG)     // deposite energy [erg]

static double HeatDepositedForEachPrimaryElectron(double E_in_eV, double xe);

static double x_e[] = {1.e-4,1.e-3,1.e-2,1.e-1};
static double f_3[] = {3.25,7.20,16.2,31.0};
enum { ElectronFractionTableSize = sizeof(x_e)/sizeof(x_e[0]) };
static gsl_spline *f3Access;
void InitializeCosmicRayIonization(int Check){
    if(f3Access == NULL){
        f3Access = gsl_spline_alloc(gsl_interp_linear,ElectronFractionTableSize);
        gsl_spline_init(f3Access,x_e,f_3,ElectronFractionTableSize);
    }

    if(Check == 1 && EQCLibRunParameters.UseCosmicRayIonization){
        FILE *output = fopen(EQC_DIAGNOSTIC_DIRECTORY "/CosmicRayIonization.dat","w");
        StructEQCData state = {.n = {1.0}, .X0 = 1.0};
        for(double xe = 1.e-6; xe <= 1.0; xe *= 2.0){
            state.nel = xe*state.n[ipHYDROGEN];
            double deposited = HeatDepositedForEachPrimaryElectron((double)E_ELECTRON,xe);
            fprintf(output,"%e %e %e\n",xe,deposited,TotalCosmicRayHeatingRate(&state));
        }
        fclose(output);
    }
}

// CR ionization
double CosmicRayIonizationRateHI(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // see Table B3 in Glover et al. (2010)
    return zetaHI;  // [s^-1]

    // ref. Shull & van Steenberg (1985)
    //      (A1) of Wolfire et al. (1995)
    if(xe<0.0)  return 0.e0;
    //xe = MIN(xe,1.0);
    double phiHI = 0.3908*pow(1.0-pow(xe,0.4092),1.7592)*(E_ELECTRON/13.6-1.0)*pow(1.e-3*E_ELECTRON,2.313*xe);
    return zetaHI*phiHI;  // [s^-1]
}

/// Cosmica ray ionization rate of D
//// D + cr -> D+ + e
double CosmicRayIonizationRateDI(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;
    
    // Le Teuff et al. (2000)
    // see Table B3 in Glover et al. (2010)
    return zetaDI;  // [s^-1]
}

/// Cosmica ray ionization rate of H2
// H2 + cr -> H2+ + e
double CosmicRayIonizationRateH2(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table B3 in Glover et al. (2010)
    return zetaH2;  // [s^-1]
}

double CosmicRayIonizationRateHeI(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table B3 in Glover et al. (2010)
    return zetaHeI;  // [s^-1]

    //// ref. Shull & van Steenberg (1985)
    ////      (A2) of Wolfire et al. (1995)
    //if(xe<0.0)  return 0.e0;
    //xe = MIN(xe,1.0);
    //double phiHeI = 0.0554*pow(1.0-pow(xe,0.4614),1.6660)*(E_ELECTRON/24.6-1.0)*pow(1.e-3*E_ELECTRON,2.313*xe);
    //return zetaHeI*phiHeI;  // [s^-1]
}

/// Cosmica ray ionization rate of CI
// C + cr -> C+ + e
double CosmicRayIonizationRateCI(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table B3 in Glover et al. (2010)
    return zetaCI;
}

/// Cosmica ray ionization rate of OI
// O + cr -> O+ + e
double CosmicRayIonizationRateOI(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table B3 in Glover et al. (2010)
    return zetaOI;
}

/// Cosmica ray ionization rate of SiI
// Si + cr -> Si+ + e
double CosmicRayIonizationRateSiI(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table 4 in Glover & Jappsen (2007)
    return zetaSiI;
}

/// Cosmica ray ionization rate of SiI
// Si+ + cr -> Si++ + e
double CosmicRayIonizationRateSiII(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table 4 in Glover & Jappsen (2007)
    return zetaSiII;
}

/// Cosmica ray ionization rate of CO
// CO + cr -> CO+ + e
double CosmicRayIonizationRateCO(double xe){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;

    // Le Teuff et al. (2000)
    // see Table B3 in Glover et al. (2010)
    return zetaCO;
}


// CR ionization heating
static double HeatDepositedForEachPrimaryElectron(double E_in_eV, double xe){
    // ref. Shull & van Steenberg (1985)
    //      eq(A3)(A4)(A5) of Wolfire et al. (1995)
    if(xe<0.0||xe>=1.0)  return 0.e0;
    //if(xe>0.5)  return 0.e0;
    double f1  = CEQPhotoelectronHeatingFraction(xe);
    double f2  = 2.0*E_in_eV-13.6;
    const int last = ElectronFractionTableSize - 1;
    double f3;
    if(xe < x_e[0] || xe > x_e[last]){
        // Continue the endpoint segment linearly, as in the reference diagnostic table.
        const int lo = xe < x_e[0] ? 0 : last - 1;
        const int hi = lo + 1;
        f3 = f_3[lo] + (xe - x_e[lo]) * (f_3[hi] - f_3[lo]) / (x_e[hi] - x_e[lo]);
    }else{
        f3 = gsl_spline_eval(f3Access,xe,NULL);
    }
    if(isnan(f1)||isnan(f2)||isnan(f3)){
        fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
        fprintf(stderr,"  f1=%e f2=%e f3=%e xe=%e\n",f1,f2,f3,xe);
        exit(1);
    }
    return (E_in_eV*EV_TO_ERG)*(f1 + 1.0/(f2/f3+1.0/(1.0-f1))); // [erg]
}
double TotalCosmicRayHeatingRate(StructEQCData *PhydroPP){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;
    double nH = PhydroPP->n[ipHYDROGEN];
    double xe = PhydroPP->nel/nH;
    double X0 = PhydroPP->X0;
    if(xe<0.0||xe>=1.0)  return 0.e0;

    //return MAX(0.0,1-xe)*nH*zetaHI*(double)E_ELECTRON*EV_TO_ERG;

    // ref. Shull & van Steenberg (1985)
    //      eq(4) and eq(A3)(A4)(A5) of Wolfire et al. (1995)
    return nH*(1-xe)*X0*zetaHI*HeatDepositedForEachPrimaryElectron((double)E_ELECTRON,xe);
}

double TotalCosmicRayHeatingRateFixed(double nH, double X0){
    if(EQCLibRunParameters.UseCosmicRayIonization==false)   return 0.0;
    // ref. Shull & van Steenberg (1985)
    //      eq(4) and eq(A3)(A4)(A5) of Wolfire et al. (1995)
    double xe   = 0.001;
    double Frac = 0.3;  // fro xe=0.001
    return nH*(1.-xe)*X0*zetaHI*(E_ELECTRON*EV_TO_ERG)*Frac;
}

