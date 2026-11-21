#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "Constants.h"
#include "Astro.h"
#include "ChemicalElements.h"
#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

static enum{
    ipHI,
    ipHeI,
    ipHeII,
    ipCI,
    ipCII,
    ipOI,
    ipOII,
    ipSiI,
    ipSiII,
    ipNpi
};

#define Eth_HI_eV       (13.60)
#define Eth_HeI_eV      (24.59)
#define Eth_HeII_eV     (54.42)
#define Eth_CI_eV       (11.26)
#define Eth_CII_eV      (24.38)
#define Eth_OI_eV       (13.62)
#define Eth_OII_eV      (35.12)
#define Eth_SiI_eV      (8.152)
#define Eth_SiII_eV     (16.35)

enum { nZ_PI = 7 };
const int iZs = 4;  // for solar metallicity
static const double Z_PI[nZ_PI] = {0.0001,0.0004,0.004,0.008,0.02,0.05,0.1};
static const char *const SpectrumFiles[nZ_PI] = {
    "spectra_Kroupa_Z0.0001.dat", "spectra_Kroupa_Z0.0004.dat",
    "spectra_Kroupa_Z0.0040.dat", "spectra_Kroupa_Z0.0080.dat",
    "spectra_Kroupa_Z0.0200.dat", "spectra_Kroupa_Z0.0500.dat",
    "spectra_Kroupa_Z0.1000.dat"
};

typedef struct {
    double photonRatio[ipNpi];
    double ionization[ipNpi];
    double heating[ipNpi];
} PhotoionizationTable;

static PhotoionizationTable Photoionization[nZ_PI];
static bool ReadPEGASESpectrum(const char *directory, const char *name, bool check,
                               PhotoionizationTable *table);

static bool CloseSpectrumFile(FILE *stream){
    if(stream == NULL) return true;
    bool ok = !ferror(stream);
    return fclose(stream) == 0 && ok;
}

CEQLibStatus CEQInitializePhotochemistry(int check, const char *pegaseDirectory){
    if(!CEQLibRunParameters.UseHIIRegion) return CEQ_SUCCESS;
    if(pegaseDirectory == NULL) pegaseDirectory = "../data/PEGASE2";
    if(*pegaseDirectory == '\0') return CEQ_INVALID_INPUT;

    PhotoionizationTable tables[nZ_PI] = {0};
    for(int i = 0; i < nZ_PI; i++){
        if(!ReadPEGASESpectrum(pegaseDirectory,SpectrumFiles[i],check == 1,&tables[i]))
            return CEQ_INITIALIZATION_FAILED;
    }

    if(check == 1){
        FILE *output = fopen(CEQ_DIAGNOSTIC_DIRECTORY "/PhotoionizationNormalizedNLyn.dat","w");
        if(output == NULL){
            perror("CEQLib could not open PhotoionizationNormalizedNLyn.dat");
            return CEQ_INITIALIZATION_FAILED;
        }
        for(int i = 0; i < nZ_PI; i++){
            fprintf(output,"%f ",Z_PI[i]);
            for(int ip = 0; ip < ipNpi; ip++)
                fprintf(output,"%e %e %e ",tables[i].photonRatio[ip],
                        tables[i].ionization[ip],tables[i].heating[ip]);
            fputc('\n',output);
        }
        if(!CloseSpectrumFile(output)){
            fprintf(stderr,"CEQLib could not write PhotoionizationNormalizedNLyn.dat\n");
            return CEQ_INITIALIZATION_FAILED;
        }
    }

    memcpy(Photoionization,tables,sizeof(tables));
    return CEQ_SUCCESS;
}

void InitializePhotochemistory(int Check){
    if(CEQInitializePhotochemistry(Check,NULL) != CEQ_SUCCESS) exit(EXIT_FAILURE);
}

//// Photo-ionization rate [1/s]
static double SimplifiedPhotoionizationRate(double NLyn, int iZ, int ip){
    if(NLyn==0.0)   return 0.e0;
    return Photoionization[iZ].ionization[ip]*NLyn;   // [1/s]
}

double HIPhotoIonizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipHI);
}
double HeIPhotoIonizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipHeI);
}
double HeIIPhotoIonizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipHeII);
}
double CIPhotoionizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipCI);
}
double CIIPhotoionizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipCII);
}
double OIPhotoionizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipOI);
}
double OIIPhotoionizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipOII);
}
double SiIPhotoionizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipSiI);
}
double SiIIPhotoionizationRate(double NLyn){
    return SimplifiedPhotoionizationRate(NLyn,iZs,ipSiII);
}

//// FUV photo-ionization [1/s]
double CIFUVPhotoionizationRate(double G0, double Av){
    // see eq(C2) and Tables 5 of Tielens & Hollenbach (1985)
    //return (2.1e-10)*G0*exp(-2.6*Av);
    
    // R171 in Table B2 of Glover et al. 2010
    if(CEQLibRunParameters.UseDustAttenuation==true)
        G0 *= exp(-3.0*Av);
    return 3.1e-10*G0; 
}
double SiIFUVPhotoionizationRate(double G0, double Av){
    // see eq(C2) and Tables 5 of Tielens & Hollenbach (1985)
    if(CEQLibRunParameters.UseDustAttenuation==true)
        G0 *= exp(-1.6*Av);
    return 1.2e-9*G0;
}


////// Photo-ionization Heating [erg/s]
static double SimplifiedPhotoionizationHeating(double NLyn, int iZ, int ip){
    if(NLyn==0.0)   return 0.e0;
    return Photoionization[iZ].heating[ip]*NLyn; // [erg/s]
}
double HeatingPhotoionization(StructCEQData *PhydroPP){
    double NLyn = PhydroPP->NLyn;
    if(NLyn==0.e0)  return 0.e0;
    return PhydroPP->natom[ipHYDROGEN]*HeatingHIPhotoionization(NLyn)
          +PhydroPP->natom[ipHELIUM  ]*HeatingHeIPhotoionization(NLyn)
          +PhydroPP->natom[ipCARBON  ]*HeatingCIPhotoionization(NLyn)
          +PhydroPP->natom[ipOXYGEN  ]*HeatingOIPhotoionization(NLyn)
          +PhydroPP->natom[ipSILICON ]*HeatingSiIPhotoionization(NLyn)
          +PhydroPP->nion [ipHELIUM  ]*HeatingHeIIPhotoionization(NLyn)
          +PhydroPP->nion [ipCARBON  ]*HeatingCIIPhotoionization(NLyn)
          +PhydroPP->nion [ipOXYGEN  ]*HeatingOIIPhotoionization(NLyn)
          +PhydroPP->nion [ipSILICON ]*HeatingSiIIPhotoionization(NLyn);
}
double HeatingHIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipHI);
}
double HeatingHeIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipHeI);
}
double HeatingHeIIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipHeII);
}
double HeatingCIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipCI);
}
double HeatingCIIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipCII);
}
double HeatingOIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipOI);
}
double HeatingOIIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipOII);
}
double HeatingSiIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipSiI);
}
double HeatingSiIIPhotoionization(double NLyn){
    return SimplifiedPhotoionizationHeating(NLyn,iZs,ipSiII);
}

////// Photo-ionization Heating [erg/cc/s]
//// ref. Osterblock (1989)
////      Draine's text (2011) p.316-317
//#define Tphi    (1.e+4) // [K] typical energy of photoelectron
//#define T_c     (1.e+4) // (0.8*Tphi)
//static double Phi0(double Tcolor){
//    // Table 6.1 in Spitzer (1978)
//    // Table 27.1 in Draine's textbook
//    // Tcolor   4000    8000    16000   32000   64000 [K]
//    // phi0     0.977   0.959   0.922   0.864   0.775
//    // <phi>    1.051   1.101   1.199   1.380   1.655
//    return 1.101;
//} 
//double HeatingHIPhotoionization(double nHI, double Epi, double NLyn){
//    //double Tcolor = Tphi;
//    // eq.(27.6) of Draine's textbook
//    //double Eel = Phi0(Tcolor)*(BOLTZMANN_CONSTANT_CGS)*Tcolor;
//    return nHI*Epi*HIPhotoIonizationRate(T_c,NLyn);
//}
//double HeatingHeIPhotoionization(double nHeI, double Epi, double NLyn){
//    return nHeI*Epi*HeIPhotoIonizationRate(0.0,NLyn);
//    //return 1.5*(BOLTZMANN_CONSTANT_CGS*Tphi)*nHeI*HeIPhotoIonizationRate(T,Fion);
//}
//
//double HeatingHIPhotoionizationForIonizationBalance(double T, double ne, double nHII, double NLyn){
//    // assuming photoionization balance
//    // ref. Osterblock's text eq.(3.2)
//    if(NLyn==0.e0)   return 0.0;
//    if(CEQLibRunParameters.UseRecombinationCaseB==true)
//        return ne*nHII*1.5*(BOLTZMANN_CONSTANT_CGS*T_c)*HIIRecombinationRateCaseB(T);
//    else
//        return ne*nHII*1.5*(BOLTZMANN_CONSTANT_CGS*T_c)*HIIRecombinationRateCaseA(T);
//}
//double HeatingHeIPhotoionizationForIonizationBalance(double T, double ne, double nHeII, double NLyn){
//    // assuming photoionization balance
//    // ref. Osterblock's text eq.(3.2)
//    if(NLyn==0.e0)   return 0.0;
//    if(CEQLibRunParameters.UseRecombinationCaseB==true)
//        return ne*nHeII*1.5*(BOLTZMANN_CONSTANT_CGS)*T_c*HeIIRecombinationRateCaseB(T);
//    else
//        return ne*nHeII*1.5*(BOLTZMANN_CONSTANT_CGS)*T_c*HeIIRecombinationRateCaseA(T);
//}

        
//// CI ionization heating [erg/cc/s]
// ref. eq(A6) of Tielens & Hollenbach (1985)
//      eq(7.19) of Stahler & Palla (2004) "The Formation of Stars"
#define dE_CI   (1.0*EV_TO_ERG) // typical energy deposited to electrons
double HeatingCarbonPhotoinozation(double T, double nCI, double G0, double Av, double Nc, double Nh2){
    double tauC  = 1.0e-17*Nc;  // attenuation by self-absorptio (Werner 1970)
    double tauH2 = 1.4e-22*Nh2; // attenuation by H2 (de Jong et al. 1980)
    return 2.2e-22*nCI*G0*exp(-2.4*Av-tauC-tauH2)/(1.e0+tauH2);
    //return dE_CI*CarbonPhotoionizationRate(G0dust)*nCI;
}


// Verner et al., ApJ, 465, 487, (1996) [cm^2]
static double PhotoionizationCrossSectionVerer1996(double E, double E0, double y0, double y1, double yw, double ya, double P){
    // eq.(1) of Verner+1996
    double x = E/E0 - y0;
    double y = sqrt(SQ(x)+SQ(y1));
    return (SQ(x-1)+SQ(yw))*pow(y,0.5*P-5.5)*pow(1+sqrt(y/ya),-P);
}

double PhotoionizationCrossSectionHI(double E){
    if(E<Eth_HI_eV)  return 0.e0;
    return 5.475e-14*PhotoionizationCrossSectionVerer1996(E,4.298e-1,0.e0,0.e0,0.e0,3.288e+1,2.963);
}

double PhotoionizationCrossSectionHeI(double E){
    if(E<Eth_HeI_eV)    return 0.e0;
    return 9.492e-16*PhotoionizationCrossSectionVerer1996(E,1.361e+1,4.434e-1,2.136,2.039,1.469,3.188);
}

double PhotoionizationCrossSectionHeII(double E){
    if(E<Eth_HeII_eV) return 0.e0;
    return 1.369e-14*PhotoionizationCrossSectionVerer1996(E,1.720,0.e0,0.e0,0.e0,3.288e+1,2.963);
}

double PhotoionizationCrossSectionCI(double E){
    if(E<Eth_CI_eV) return 0.e0;
    return 5.027e-16*PhotoionizationCrossSectionVerer1996(E,2.144,1.133,1.607,9.157e-2,6.216e+1,5.101);
}

double PhotoionizationCrossSectionCII(double E){
    if(E<Eth_CII_eV) return 0.e0;
    return 8.709e-18*PhotoionizationCrossSectionVerer1996(E,0.4058,49.29,3.234,2.093,1.261e+2,8.578);
}

double PhotoionizationCrossSectionOI(double E){
    if(E<Eth_OI_eV) return 0.e0;
    return 1.745e-15*PhotoionizationCrossSectionVerer1996(E,1.240,8.698,0.1271,0.07589,3.784,17.64);
}

double PhotoionizationCrossSectionOII(double E){
    if(E<Eth_OII_eV) return 0.e0;
    return 5.967e-17*PhotoionizationCrossSectionVerer1996(E,1.386,21.31,1.503e-2,1.934e-2,3.175e+1,8.943);
}

double PhotoionizationCrossSectionSiI(double E){
    if(E<Eth_SiI_eV) return 0.e0;
    return 2.506e-17*PhotoionizationCrossSectionVerer1996(E,23.17,1.672e-5,0.4207,0.2837,20.57,3.546);
}

double PhotoionizationCrossSectionSiII(double E){
    if(E<Eth_SiII_eV) return 0.e0;
    return 4.140e-18*PhotoionizationCrossSectionVerer1996(E,2.556,6.634,0.1272,1.570,13.37,11.91);
}

static FILE *OpenSpectrumFile(const char *directory, const char *prefix,
                              const char *name, const char *mode){
    const int length = snprintf(NULL,0,"%s/%s%s",directory,prefix,name);
    if(length < 0) return NULL;
    char *path = malloc((size_t)length+1);
    if(path == NULL) return NULL;
    snprintf(path,(size_t)length+1,"%s/%s%s",directory,prefix,name);
    FILE *stream = fopen(path,mode);
    if(stream == NULL) perror(path);
    free(path);
    return stream;
}

static bool ReadSpectrumNumber(FILE *input, double *value){
    return fscanf(input,"%lf",value) == 1 && isfinite(*value);
}

static bool ReadPEGASESpectrum(const char *directory, const char *name, bool check,
                               PhotoionizationTable *table){
    enum { ContinuumCapacity = 1500, LineCapacity = 100, HeaderLines = 27, Properties = 19 };
    FILE *input = OpenSpectrumFile(directory,"",name,"r");
    if(input == NULL) return false;

    FILE *spectrum = NULL, *rates = NULL, *scenario = NULL;
    bool ok = false;
    if(check){
        spectrum = OpenSpectrumFile("out","",name,"w");
        rates = OpenSpectrumFile("out","Photoionization_",name,"w");
        scenario = OpenSpectrumFile("out","","PEGASE_evolutionary_scenario.txt","w");
        if(spectrum == NULL || rates == NULL || scenario == NULL) goto done;
    }

    for(int j = 0; j < HeaderLines; j++){
        int ch;
        do {
            ch = fgetc(input);
            if(ch == EOF) goto done;
            if(scenario != NULL && fputc(ch,scenario) == EOF) goto done;
        } while(ch != '\n');
    }
    double timestepCount, continuumCount, lineCount;
    if(!ReadSpectrumNumber(input,&timestepCount) || !ReadSpectrumNumber(input,&continuumCount) ||
       !ReadSpectrumNumber(input,&lineCount) ||
       timestepCount < 1 || timestepCount > INT_MAX || timestepCount != floor(timestepCount) ||
       continuumCount < 2 || continuumCount >= ContinuumCapacity || continuumCount != floor(continuumCount) ||
       lineCount < 0 || lineCount >= LineCapacity || lineCount != floor(lineCount)) goto done;

    const int nbtimesteps = (int)timestepCount;
    const int nbcontinuum = (int)continuumCount;
    const int nblines = (int)lineCount;
    double wavelength[ContinuumCapacity] = {0.0};
    double energy[ContinuumCapacity], flux[ContinuumCapacity], photonFlux[ContinuumCapacity];
    for(int k = 1; k <= nbcontinuum; k++){
        if(!ReadSpectrumNumber(input,&wavelength[k]) || wavelength[k] <= wavelength[k-1]) goto done;
    }
    double ignored;
    for(int k = 0; k < nblines; k++){
        if(!ReadSpectrumNumber(input,&ignored)) goto done;
    }

    int youngSpectra = 0;
    for(int j = 0; j < nbtimesteps; j++){
        double age;
        if(!ReadSpectrumNumber(input,&age) || age < 0.0) goto done;
        for(int k = 1; k < Properties; k++){
            if(!ReadSpectrumNumber(input,&ignored)) goto done;
        }
        for(int k = 1; k <= nbcontinuum; k++){
            if(!ReadSpectrumNumber(input,&flux[k])) goto done;
        }
        for(int k = 0; k < nblines; k++){
            if(!ReadSpectrumNumber(input,&ignored)) goto done;
        }

        for(int k = 1; k <= nbcontinuum; k++){
            const double energyErg = PLANK_H_CGS*LIGHT_C_CGS/(wavelength[k]*ANGSTROM_CGS);
            photonFlux[k] = flux[k]/energyErg;
            energy[k] = energyErg*ERG_TO_EV;
            if(!isfinite(energy[k]) || energy[k] <= 0.0 || !isfinite(photonFlux[k])) goto done;
            if(spectrum != NULL)
                fprintf(spectrum,"%lf %lf %le %le \n",wavelength[k],energy[k],flux[k],photonFlux[k]);
        }
        if(spectrum != NULL) fputc('\n',spectrum);

        // photoionization rate [cm^2/s/Msun]
        double N_LynHI = 0.e0;
        double N_Lyn[ipNpi], xi[ipNpi], gamma[ipNpi];
        for(int ip=0; ip<ipNpi; ip++)
            N_Lyn[ip] = xi[ip] = gamma[ip] = 0.e0;
        for(int k=1; k<nbcontinuum; k++){
            double dl = wavelength[k+1] - wavelength[k]; // [A]
            double e1 = energy[k  ];    // [eV]
            double e2 = energy[k+1];
            if(e2>=Eth_HI_eV)   N_LynHI += (photonFlux[k] + photonFlux[k+1])*dl;
            //else if(e2>=10.0)   N_LW  += (photonFlux[k] + photonFlux[k+1])*dl;

            // x_i     [photons/s/A/Msun]*[cm^2]*[A] = [photons*cm^2/s/Msun]
            // gamma_i [photons/s/A/Msun]*[eV]*[cm^2][A] = [eV*cm^2/s/Msun]
            // NLyn_i  [photons/s/A/Msun]*[A] = [photons/s/Msun]
            double s1HI  = PhotoionizationCrossSectionHI(e1);  // [cm^2]
            double s2HI  = PhotoionizationCrossSectionHI(e2);
            double de1HI = e1 - Eth_HI_eV;
            double de2HI = e2 - Eth_HI_eV;
            xi[ipHI]    += (photonFlux[k]*s1HI + photonFlux[k+1]*s2HI)*dl;
            if(e2>=Eth_HI_eV){
                N_Lyn[ipHI] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipHI] += (de1HI*photonFlux[k]*s1HI + de2HI*photonFlux[k+1]*s2HI)*dl;
            }

            double s1HeI  = PhotoionizationCrossSectionHeI(e1);
            double s2HeI  = PhotoionizationCrossSectionHeI(e2);
            double de1HeI = e1 - Eth_HeI_eV;
            double de2HeI = e2 - Eth_HeI_eV;
            xi[ipHeI]    += (photonFlux[k]*s1HeI + photonFlux[k+1]*s2HeI)*dl;
            if(e2>=Eth_HeI_eV){
                N_Lyn[ipHeI] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipHeI] += (de1HeI*photonFlux[k]*s1HeI + de2HeI*photonFlux[k+1]*s2HeI)*dl;
            }

            double s1HeII  = PhotoionizationCrossSectionHeII(e1);
            double s2HeII  = PhotoionizationCrossSectionHeII(e2);
            double de1HeII = e1 - Eth_HeII_eV;
            double de2HeII = e2 - Eth_HeII_eV;
            xi[ipHeII]    += (photonFlux[k]*s1HeII + photonFlux[k+1]*s2HeII)*dl;
            if(e2>=Eth_HeII_eV){
                N_Lyn[ipHeII] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipHeII] += (de1HeII*photonFlux[k]*s1HeII + de2HeII*photonFlux[k+1]*s2HeII)*dl;
            }

            double s1CI   = PhotoionizationCrossSectionCI(e1);
            double s2CI   = PhotoionizationCrossSectionCI(e2);
            double de1CI  = e1 - Eth_CI_eV;
            double de2CI  = e2 - Eth_CI_eV;
            xi[ipCI]    += (photonFlux[k]*s1CI + photonFlux[k+1]*s2CI)*dl;
            if(e2>=Eth_CI_eV){
                N_Lyn[ipCI] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipCI] += (de1CI*photonFlux[k]*s1CI + de2CI*photonFlux[k+1]*s2CI)*dl;
            }

            double s1CII  = PhotoionizationCrossSectionCII(e1);
            double s2CII  = PhotoionizationCrossSectionCII(e2);
            double de1CII = e1 - Eth_CII_eV;
            double de2CII = e2 - Eth_CII_eV;
            xi[ipCII]    += (photonFlux[k]*s1CII + photonFlux[k+1]*s2CII)*dl;
            if(e2>=Eth_CII_eV){
                N_Lyn[ipCII] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipCII] += (de1CII*photonFlux[k]*s1CII + de2CII*photonFlux[k+1]*s2CII)*dl;
            }

            double s1OI   = PhotoionizationCrossSectionOI(e1);
            double s2OI   = PhotoionizationCrossSectionOI(e2);
            double de1OI  = e1 - Eth_OI_eV;
            double de2OI  = e2 - Eth_OI_eV;
            xi[ipOI]    += (photonFlux[k]*s1OI + photonFlux[k+1]*s2OI)*dl;
            if(e2>=Eth_OI_eV){
                N_Lyn[ipOI] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipOI] += (de1OI*photonFlux[k]*s1OI + de2OI*photonFlux[k+1]*s2OI)*dl;
            }

            double s1OII  = PhotoionizationCrossSectionOII(e1);
            double s2OII  = PhotoionizationCrossSectionOII(e2);
            double de1OII = e1 - Eth_OII_eV;
            double de2OII = e2 - Eth_OII_eV;
            xi[ipOII]    += (photonFlux[k]*s1OII + photonFlux[k+1]*s2OII)*dl;
            if(e2>=Eth_OII_eV){
                N_Lyn[ipOII] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipOII] += (de1OII*photonFlux[k]*s1OII + de2OII*photonFlux[k+1]*s2OII)*dl;
            }

            double s1SiI  = PhotoionizationCrossSectionSiI(e1);
            double s2SiI  = PhotoionizationCrossSectionSiI(e2);
            double de1SiI = e1 - Eth_SiI_eV;
            double de2SiI = e2 - Eth_SiI_eV;
            xi[ipSiI]    += (photonFlux[k]*s1SiI + photonFlux[k+1]*s2SiI)*dl;
            if(e2>=Eth_SiI_eV){
                N_Lyn[ipSiI] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipSiI] += (de1SiI*photonFlux[k]*s1SiI + de2SiI*photonFlux[k+1]*s2SiI)*dl;
            }

            double s1SiII  = PhotoionizationCrossSectionSiII(e1);
            double s2SiII  = PhotoionizationCrossSectionSiII(e2);
            double de1SiII = e1 - Eth_SiII_eV;
            double de2SiII = e2 - Eth_SiII_eV;
            xi[ipSiII]    += (photonFlux[k]*s1SiII + photonFlux[k+1]*s2SiII)*dl;
            if(e2>=Eth_SiII_eV){
                N_Lyn[ipSiII] += (photonFlux[k] + photonFlux[k+1])*dl;
                gamma[ipSiII] += (de1SiII*photonFlux[k]*s1SiII + de2SiII*photonFlux[k+1]*s2SiII)*dl;
            }
        }

        N_LynHI       *= 0.5;
        for(int ip=0; ip<ipNpi; ip++){
            N_Lyn[ip] *= 0.5;           // [photons/s/Msun]
            xi[ip]    *= 0.5;           // [photons*cm^2/s/Msun]
            gamma[ip] *= 0.5*EV_TO_ERG; // [erg*cm^2/s/Msun]
        }

        if(rates != NULL){
            fprintf(rates,"%f ",age);
            for(int ip = 0; ip < ipNpi; ip++) fprintf(rates,"%e ",xi[ip]);
            for(int ip = 0; ip < ipNpi; ip++) fprintf(rates,"%e ",gamma[ip]);
            for(int ip = 0; ip < ipNpi; ip++) fprintf(rates,"%e ",N_Lyn[ip]);
            fputc('\n',rates);
        }

        if(age > 30.0) continue;
        for(int ip = 0; ip < ipNpi; ip++){
            if(N_Lyn[ip] == 0.0) continue;
            if(N_LynHI <= 0.0 || !isfinite(N_LynHI)) goto done;
            table->photonRatio[ip] += N_Lyn[ip]/N_LynHI;
            /* Every species uses the incident H-ionizing photon flux. */
            table->ionization[ip] += xi[ip]/N_LynHI;
            table->heating[ip] += gamma[ip]/N_LynHI;
        }
        youngSpectra++;
    }
    if(youngSpectra == 0) goto done;
    for(int ip = 0; ip < ipNpi; ip++){
        table->photonRatio[ip] /= youngSpectra;
        table->ionization[ip] /= youngSpectra;
        table->heating[ip] /= youngSpectra;
        if(!isfinite(table->photonRatio[ip]) || !isfinite(table->ionization[ip]) ||
           !isfinite(table->heating[ip])) goto done;
    }
    ok = true;

done:
    if(!CloseSpectrumFile(input)) ok = false;
    if(!CloseSpectrumFile(spectrum)) ok = false;
    if(!CloseSpectrumFile(rates)) ok = false;
    if(!CloseSpectrumFile(scenario)) ok = false;
    if(!ok) fprintf(stderr,"CEQLib could not load PEGASE2 spectrum %s/%s\n",directory,name);
    return ok;
}
