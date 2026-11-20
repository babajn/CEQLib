#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
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

#define nZ_PI   7
const int iZs = 4;  // for solar metallicity
static double Zsol = SOLAR_METALLICITY;
static double Z_PI[] = {0.0001,0.0004,0.004,0.008,0.02,0.05,0.1};
static double PhotoionizationRateToLyn[nZ_PI][ipNpi];
static double PhotoionizationHeatingToLyn[nZ_PI][ipNpi];
static double PhotoionizationNpiToLyn[nZ_PI][ipNpi];

static void ReadSpeactraPEAGSE(const int iZ, const char *fname);

void InitializePhotochemistory(int Check){

    if(CEQLibRunParameters.UseHIIRegion==true){
        for(int iZ=0; iZ<nZ_PI; iZ++){
            for(int ip=0; ip<ipNpi; ip++){
                PhotoionizationRateToLyn[iZ][ip]    = 
                PhotoionizationHeatingToLyn[iZ][ip] = 
                PhotoionizationNpiToLyn[iZ][ip]     = 0.e0;
            }
        }

        ReadSpeactraPEAGSE(0,"spectra_Kroupa_Z0.0001.dat");
        ReadSpeactraPEAGSE(1,"spectra_Kroupa_Z0.0004.dat");
        ReadSpeactraPEAGSE(2,"spectra_Kroupa_Z0.0040.dat");
        ReadSpeactraPEAGSE(3,"spectra_Kroupa_Z0.0080.dat");
        ReadSpeactraPEAGSE(4,"spectra_Kroupa_Z0.0200.dat");
        ReadSpeactraPEAGSE(5,"spectra_Kroupa_Z0.0500.dat");
        ReadSpeactraPEAGSE(6,"spectra_Kroupa_Z0.1000.dat");

        if(Check==1){
            sprintf(CEQFileName,"%s/PhotoionizationNormalizedNLyn.dat",CEQOutDir);
            FILE *output = fopen(CEQFileName,"w");
            for(int iZ=0; iZ<nZ_PI; iZ++){
                fprintf(output,"%f ",Z_PI[iZ]);
                for(int ip=0; ip<ipNpi; ip++){
                    fprintf(output,"%e ",PhotoionizationNpiToLyn[iZ][ip]    );
                    fprintf(output,"%e ",PhotoionizationRateToLyn[iZ][ip]   );
                    fprintf(output,"%e ",PhotoionizationHeatingToLyn[iZ][ip]);
                }
                fprintf(output,"\n");
            }
            fclose(output);
        }

        //exit(1);
    }

    return ;
}

//// Photo-ionization rate [1/s]
static double SimplifiedPhotoionizationRate(double NLyn, int iZ, int ip){
    if(NLyn==0.0)   return 0.e0;
    return PhotoionizationRateToLyn[iZ][ip]*NLyn;   // [1/s]
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
        G0 *= exp(3.0*Av);
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
    return PhotoionizationHeatingToLyn[iZ][ip]*NLyn; // [erg/s]
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

static void ReadSpeactraPEAGSE(const int iZ, const char *fname){

	char inputfile[200], outputfile[200];
    char DataDir[] = "../data/PEGASE2/";
	sprintf(inputfile,"%s/%s",DataDir,fname);
	FILE *FRV = fopen(inputfile,"r");
	if(FRV == NULL){
		fprintf(stderr,"cannot open this file\n");
		fprintf(stderr,"inputfile = %s\n",inputfile);
		fprintf(stderr,"abort\n");
		exit(1);
	}
	sprintf(outputfile,"out/%s",fname);
	FILE *synthspecflux  = fopen(outputfile,"w");
	sprintf(outputfile,"out/Photoionization_%s",fname);
	FILE *photoionflux   = fopen(outputfile,"w");
	FILE *scenariofile   = fopen("out/PEGASE_evolutionary_scenario.txt","w");
	//FILE *galaxy_mass    = fopen("out/PEGASE_mass_of_galaxy.txt","w");
	//FILE *star_mass      = fopen("out/PEGASE_mass_of_star.txt","w");
	//FILE *WD_mass        = fopen("out/PEGASE_mass_of_WD.txt","w");
	//FILE *NS_mass        = fopen("out/PEGASE_mass_of_NS.txt","w");
	//FILE *substar_mass   = fopen("out/PEGASE_mass_of_substar.txt","w");
	//FILE *gas_mass       = fopen("out/PEGASE_mass_of_gas.txt","w");
	//FILE *Z_ISM          = fopen("out/PEGASE_metal_ISM.txt","w");
	//FILE *Z_star_mass    = fopen("out/PEGASE_metal_starmass.txt","w");
	//FILE *Z_star_BL      = fopen("out/PEGASE_metal_starBL.txt","w");
	//FILE *BolometricL    = fopen("out/PEGASE_bolometric_luminosity.txt","w");
	//FILE *opt_depth      = fopen("out/PEGASE_optically_depth.txt","w");
	//FILE *dust_BL        = fopen("out/PEGASE_dust_BL_ratio.txt","w");
	//FILE *SFRate         = fopen("out/PEGASE_SFR.txt","w");
	//FILE *Ly_photon      = fopen("out/PEGASE_Ly_photon.txt","w");
	//FILE *SN2_rate       = fopen("out/PEGASE_SNII_rate.txt","w");
	//FILE *SN1a_rate      = fopen("out/PEGASE_SNIa_rate.txt","w");	
	//FILE *time_star_mass = fopen("out/PEGASE_age_of_star_mass.txt","w");	
	//FILE *time_star_BL   = fopen("out/PEGASE_age_of_star_BL.txt","w");	

	/*Evolutionary Scenarioの読み込み*/
	char scenarios[200];	
	int  nbtimesteps, nbcontinuum, nblines;
	double wave_conti[1500], ener_conti[1500];
    double flux_conti[1500], nphot_conti[1500];
    double wave_lines[100], flux_lines[100];
	double age, M_galaxy, M_star, M_WD, M_NS, M_substar, M_gas, metal_ISM, metal_star_mass, metal_star_BL;
	double BL, tau_V, dust, SFR, Ly, SN2, SN1a, age_star_mass, age_star_BL;
	for(int j=1;j<=27;j++){
		fgets(scenarios, 200, FRV);
		fputs(scenarios, scenariofile);
	}
	/*timestep、波長サンプル数の読み込み*/
	fscanf(FRV,"%d %d %d",&nbtimesteps,&nbcontinuum,&nblines);
	/*波長の読み込み*/
	for(int j=1; j<=nbcontinuum; j++)   fscanf(FRV,"%lf",&wave_conti[j]);
	for(int j=1; j<=nblines; j++)       fscanf(FRV,"%lf",&wave_lines[j]);
	/*nbtimesteps回分読み込む*/
    int Ntimes = 0;
	for(int j=1; j<=nbtimesteps; j++){
		/*propertyの読み込み*/
		fscanf(FRV,"%lf",&age);
		fscanf(FRV,"%lf",&M_galaxy);
		fscanf(FRV,"%lf",&M_star);
		fscanf(FRV,"%lf",&M_WD);
		fscanf(FRV,"%lf",&M_NS);
		fscanf(FRV,"%lf",&M_substar);
		fscanf(FRV,"%lf",&M_gas);
		fscanf(FRV,"%lf",&metal_ISM);
		fscanf(FRV,"%lf",&metal_star_mass);
		fscanf(FRV,"%lf",&metal_star_BL);
		fscanf(FRV,"%lf",&BL);
		fscanf(FRV,"%lf",&tau_V);
		fscanf(FRV,"%lf",&dust);
		fscanf(FRV,"%lf",&SFR);
		fscanf(FRV,"%lf",&Ly);
		fscanf(FRV,"%lf",&SN2);
		fscanf(FRV,"%lf",&SN1a);
		fscanf(FRV,"%lf",&age_star_mass);
		fscanf(FRV,"%lf",&age_star_BL);
		/*fluxの読み込み*/
		wave_conti[0] = wave_lines[0]= 0.0;
		for(int k=1;k<=nbcontinuum; k++)    fscanf(FRV,"%le",&flux_conti[k]);
		for(int k=1;k<=nblines; k++)        fscanf(FRV,"%le",&flux_lines[k]);
		
        for(int k=1; k<=nbcontinuum; k++){
            // [A] -> [erg]
            ener_conti[k]  = PLANK_H_CGS*LIGHT_C_CGS/(wave_conti[k]*ANGSTROM_CGS);
            // [erg/s/A] -> [photons/s/A]
            //nphot_conti[k] = 4.0*PI*flux_conti[k]/ener_conti[k];
            nphot_conti[k] = flux_conti[k]/ener_conti[k];
           
            // [erg] -> [eV]
            ener_conti[k] *= ERG_TO_EV;
        }

		/*波長vs.fluxファイル*/
		for(int k=1; k<=nbcontinuum; k++){
            fprintf(synthspecflux,"%lf ",wave_conti[k] );    // [A]
            fprintf(synthspecflux,"%lf ",ener_conti[k] );    // [eV]
            fprintf(synthspecflux,"%le ",flux_conti[k] );    // [erg/s/A/Msun]
            fprintf(synthspecflux,"%le ",nphot_conti[k]);    // [photons/s/A/Msun]
            fprintf(synthspecflux,"\n");
            //fprintf(stderr,"lamb = %f E = %f Flux=%e\n",wave_conti[k],ener_conti[k],flux_conti[k]);
        }
        fprintf(synthspecflux,"\n");
    
        // photoionization rate [cm^2/s/Msun]
        double N_LynHI = 0.e0;
        double N_Lyn[ipNpi], xi[ipNpi], gamma[ipNpi];
        for(int ip=0; ip<ipNpi; ip++)
            N_Lyn[ip] = xi[ip] = gamma[ip] = 0.e0;
		for(int k=1; k<=nbcontinuum-1; k++){
            double dl = wave_conti[k+1] - wave_conti[k]; // [A]
            double e1 = ener_conti[k  ];    // [eV]
            double e2 = ener_conti[k+1];
            if(e2>=Eth_HI_eV)   N_LynHI += (nphot_conti[k] + nphot_conti[k+1])*dl; 
            //else if(e2>=10.0)   N_LW  += (nphot_conti[k] + nphot_conti[k+1])*dl; 

            // x_i     [photons/s/A/Msun]*[cm^2]*[A] = [photons*cm^2/s/Msun]
            // gamma_i [photons/s/A/Msun]*[eV]*[cm^2][A] = [eV*cm^2/s/Msun]
            // NLyn_i  [photons/s/A/Msun]*[A] = [photons/s/Msun]
            double s1HI  = PhotoionizationCrossSectionHI(e1);  // [cm^2]
            double s2HI  = PhotoionizationCrossSectionHI(e2);
            double de1HI = e1 - Eth_HI_eV;
            double de2HI = e2 - Eth_HI_eV;
            xi[ipHI]    += (nphot_conti[k]*s1HI + nphot_conti[k+1]*s2HI)*dl;
            if(e2>=Eth_HI_eV){
                N_Lyn[ipHI] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipHI] += (de1HI*nphot_conti[k]*s1HI + de2HI*nphot_conti[k+1]*s2HI)*dl;
            }

            double s1HeI  = PhotoionizationCrossSectionHeI(e1);
            double s2HeI  = PhotoionizationCrossSectionHeI(e2);
            double de1HeI = e1 - Eth_HeI_eV;
            double de2HeI = e2 - Eth_HeI_eV;
            xi[ipHeI]    += (nphot_conti[k]*s1HeI + nphot_conti[k+1]*s2HeI)*dl;
            if(e2>=Eth_HeI_eV){
                N_Lyn[ipHeI] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipHeI] += (de1HeI*nphot_conti[k]*s1HeI + de2HeI*nphot_conti[k+1]*s2HeI)*dl;
            }

            double s1HeII  = PhotoionizationCrossSectionHeII(e1);
            double s2HeII  = PhotoionizationCrossSectionHeII(e2);
            double de1HeII = e1 - Eth_HeII_eV;
            double de2HeII = e2 - Eth_HeII_eV;
            xi[ipHeII]    += (nphot_conti[k]*s1HeII + nphot_conti[k+1]*s2HeII)*dl;
            if(e2>=Eth_HeII_eV){
                N_Lyn[ipHeII] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipHeII] += (de1HeII*nphot_conti[k]*s1HeII + de2HeII*nphot_conti[k+1]*s2HeII)*dl;
            }

            double s1CI   = PhotoionizationCrossSectionCI(e1);
            double s2CI   = PhotoionizationCrossSectionCI(e2);
            double de1CI  = e1 - Eth_CI_eV;
            double de2CI  = e2 - Eth_CI_eV;
            xi[ipCI]    += (nphot_conti[k]*s1CI + nphot_conti[k+1]*s2CI)*dl;
            if(e2>=Eth_CI_eV){
                N_Lyn[ipCI] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipCI] += (de1CI*nphot_conti[k]*s1CI + de2CI*nphot_conti[k+1]*s2CI)*dl;
            }

            double s1CII  = PhotoionizationCrossSectionCII(e1);
            double s2CII  = PhotoionizationCrossSectionCII(e2);
            double de1CII = e1 - Eth_CII_eV;
            double de2CII = e2 - Eth_CII_eV;
            xi[ipCII]    += (nphot_conti[k]*s1CII + nphot_conti[k+1]*s2CII)*dl;
            if(e2>=Eth_CII_eV){
                N_Lyn[ipCII] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipCII] += (de1CII*nphot_conti[k]*s1CII + de2CII*nphot_conti[k+1]*s2CII)*dl;
            }

            double s1OI   = PhotoionizationCrossSectionOI(e1);
            double s2OI   = PhotoionizationCrossSectionOI(e2);
            double de1OI  = e1 - Eth_OI_eV;
            double de2OI  = e2 - Eth_OI_eV;
            xi[ipOI]    += (nphot_conti[k]*s1OI + nphot_conti[k+1]*s2OI)*dl;
            if(e2>=Eth_OI_eV){
                N_Lyn[ipOI] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipOI] += (de1OI*nphot_conti[k]*s1OI + de2OI*nphot_conti[k+1]*s2OI)*dl;
            }

            double s1OII  = PhotoionizationCrossSectionOII(e1);
            double s2OII  = PhotoionizationCrossSectionOII(e2);
            double de1OII = e1 - Eth_OII_eV;
            double de2OII = e2 - Eth_OII_eV;
            xi[ipOII]    += (nphot_conti[k]*s1OII + nphot_conti[k+1]*s2OII)*dl;
            if(e2>=Eth_OII_eV){
                N_Lyn[ipOII] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipOII] += (de1OII*nphot_conti[k]*s1OII + de2OII*nphot_conti[k+1]*s2OII)*dl;
            }

            double s1SiI  = PhotoionizationCrossSectionSiI(e1);
            double s2SiI  = PhotoionizationCrossSectionSiI(e2);
            double de1SiI = e1 - Eth_SiI_eV;
            double de2SiI = e2 - Eth_SiI_eV;
            xi[ipSiI]    += (nphot_conti[k]*s1SiI + nphot_conti[k+1]*s2SiI)*dl;
            if(e2>=Eth_SiI_eV){
                N_Lyn[ipSiI] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipSiI] += (de1SiI*nphot_conti[k]*s1SiI + de2SiI*nphot_conti[k+1]*s2SiI)*dl;
            }

            double s1SiII  = PhotoionizationCrossSectionSiII(e1);
            double s2SiII  = PhotoionizationCrossSectionSiII(e2);
            double de1SiII = e1 - Eth_SiII_eV;
            double de2SiII = e2 - Eth_SiII_eV;
            xi[ipSiII]    += (nphot_conti[k]*s1SiII + nphot_conti[k+1]*s2SiII)*dl;
            if(e2>=Eth_SiII_eV){
                N_Lyn[ipSiII] += (nphot_conti[k] + nphot_conti[k+1])*dl; 
                gamma[ipSiII] += (de1SiII*nphot_conti[k]*s1SiII + de2SiII*nphot_conti[k+1]*s2SiII)*dl;
            }
        }
       
        N_LynHI       *= 0.5;
        for(int ip=0; ip<ipNpi; ip++){
            N_Lyn[ip] *= 0.5;           // [photons/s/Msun]
            xi[ip]    *= 0.5;           // [photons*cm^2/s/Msun]
            gamma[ip] *= 0.5*EV_TO_ERG; // [erg*cm^2/s/Msun]
        }
        fprintf(photoionflux,"%f ",age          ); // #01 
        fprintf(photoionflux,"%e ",xi[ipHI     ]); // #02
        fprintf(photoionflux,"%e ",xi[ipHeI    ]); // #03
        fprintf(photoionflux,"%e ",xi[ipHeII   ]); // #04
        fprintf(photoionflux,"%e ",xi[ipCI     ]); // #05
        fprintf(photoionflux,"%e ",xi[ipCII    ]); // #06
        fprintf(photoionflux,"%e ",xi[ipOI     ]); // #07
        fprintf(photoionflux,"%e ",xi[ipOII    ]); // #08
        fprintf(photoionflux,"%e ",xi[ipSiI    ]); // #09
	    fprintf(photoionflux,"%e ",xi[ipSiII   ]); // #10
        fprintf(photoionflux,"%e ",gamma[ipHI  ]); // #11
        fprintf(photoionflux,"%e ",gamma[ipHeI ]); // #12
        fprintf(photoionflux,"%e ",gamma[ipHeII]); // #13
        fprintf(photoionflux,"%e ",gamma[ipCI  ]); // #14
        fprintf(photoionflux,"%e ",gamma[ipCII ]); // #15
        fprintf(photoionflux,"%e ",gamma[ipOI  ]); // #16
        fprintf(photoionflux,"%e ",gamma[ipOII ]); // #17
        fprintf(photoionflux,"%e ",gamma[ipSiI ]); // #18
	    fprintf(photoionflux,"%e ",gamma[ipSiII]); // #19
        fprintf(photoionflux,"%e ",N_Lyn[ipHI  ]); // #20
        fprintf(photoionflux,"%e ",N_Lyn[ipHeI ]); // #21
        fprintf(photoionflux,"%e ",N_Lyn[ipHeII]); // #22
        fprintf(photoionflux,"%e ",N_Lyn[ipCI  ]); // #23
        fprintf(photoionflux,"%e ",N_Lyn[ipCII ]); // #24
        fprintf(photoionflux,"%e ",N_Lyn[ipOI  ]); // #25
        fprintf(photoionflux,"%e ",N_Lyn[ipOII ]); // #26
        fprintf(photoionflux,"%e ",N_Lyn[ipSiI ]); // #27
	    fprintf(photoionflux,"%e ",N_Lyn[ipSiII]); // #28
	    fprintf(photoionflux,"\n");

        if(age>30.0)    continue;
        for(int ip=0; ip<ipNpi; ip++){
            if(N_Lyn[ip]==0.0)  continue;
            // 
            PhotoionizationNpiToLyn[iZ][ip]     += N_Lyn[ip]/N_LynHI ; 
            // [cm^2]
            PhotoionizationRateToLyn[iZ][ip]    += xi[ip]   /N_Lyn[ip]; 
            // [erg*cm^2]
            PhotoionizationHeatingToLyn[iZ][ip] += gamma[ip]/N_Lyn[ip];
        }
        Ntimes ++;
    }
	fclose(FRV);
	fclose(scenariofile);
	fclose(synthspecflux);
	fclose(photoionflux);

    //fprintf(stderr," Z=%f: \n",Z_PI[iZ]);
    for(int ip=0; ip<ipNpi; ip++){
        PhotoionizationNpiToLyn[iZ][ip]     /= (double)Ntimes; 
        PhotoionizationRateToLyn[iZ][ip]    /= (double)Ntimes; 
        PhotoionizationHeatingToLyn[iZ][ip] /= (double)Ntimes;
        //fprintf(stderr," ip=%d: Npi/NLynHI = %e xi/NLyn = %e Gi/NLyn = %e\n",ip,
        //        PhotoionizationNpiToLyn[iZ][ip],
        //        PhotoionizationRateToLyn[iZ][ip],
        //        PhotoionizationHeatingToLyn[iZ][ip]);
    }

    return ;
}
