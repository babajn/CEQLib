#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <assert.h>
#include "Astro.h"
#include "Constants.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "GrainSurfaceChemistry.h"

#define G0_Draine   (1.7)
#define Zsol        (SOLAR_METALLICITY)

//// Wolfire et al. (1995) Table 2 
#define HeliumAbundance     (pow(10.0,-1.00))
#define CarbonAbundance     (pow(10.0,-3.52))
#define OxygenAbundance     (pow(10.0,-3.34))
#define SiliconAbundance    (pow(10.0,-5.45))
#define MagnesiumAbundance  (pow(10.0,-4.84))
#define IronAbundance       (pow(10.0,-6.15))

int main_CoolingFunctionTest(void);
int main_WolfireTest(void);
int main_WolfireMilkyWayTest(void);
int main_HIIregionTest(void);

static void main_CoolingFunctionTestBody(double nH, double G0, double Nt, double Z, double X0, char *fname);
static unsigned wolfireFailures;
static void main_WolfireTestBody(double G0, double Nt, double Z, double X0, char *fname);

int main(int argc, char *argv[]){
    
    if(mkdir("./out",0777) != 0 && errno != EEXIST){
        perror("./out");
        return EXIT_FAILURE;
    }
    
    //return main_CoolingFunctionTest();
    //return main_HIIregionTest();
    return main_WolfireTest();
    //return main_WolfireMilkyWayTest();

    return ;
}

static CEQLibStatus SetInitialHydroDataWolfireTest(double G0, double Nt, double Z, double X0, StructCEQData *PhydroPP){
    
    // tentative fractional abundances
    double fHII   = 0.90, fHI   = 0.01, fH2  = 1.0-fHI  -fHII  ;
    double fHeIII = 0.e0, fHeII = 0.99, fHeI = 1.0-fHeII-fHeIII;
    double fCIII  = 0.e0, fCII  = 0.99, fCI  = 1.0-fCII-fCIII, fCO = 0.e0;
    double fOIII  = 0.e0, fOII  = 0.90, fOI  = 1.0-fOII-fOIII;
    double fSiIII = 0.e0, fSiII = 0.99, fSiI = 1.0-fSiII-fSiIII;
    double nH = PhydroPP->n[ipHYDROGEN];

    PhydroPP->Z  = Z;
    PhydroPP->G0 = G0;
    PhydroPP->Nt = Nt;  // [H/cm2]
    PhydroPP->X0 = X0;
    PhydroPP->Zd = Z*ApproximateDustToGasRatio(PhydroPP->T);
    PhydroPP->Td = ApproximateDustTemperature(PhydroPP->Nt,PhydroPP->G0);

    // hydrogem [/cc]
    //PhydroPP->n[ipHYDROGEN]     = nH;
    PhydroPP->nion[ipHYDROGEN]  = fHII*PhydroPP->n[ipHYDROGEN];
    PhydroPP->natom[ipHYDROGEN] = fHI *PhydroPP->n[ipHYDROGEN];
    PhydroPP->nmol[ipH2]        = fH2 *PhydroPP->n[ipHYDROGEN];
    // helium
    PhydroPP->n[ipHELIUM]       = (HeliumAbundance )*nH;
    PhydroPP->nHeIII            = fHeIII*PhydroPP->n[ipHELIUM];
    PhydroPP->nion[ipHELIUM]    = fHeII *PhydroPP->n[ipHELIUM];
    PhydroPP->natom[ipHELIUM]   = fHeI  *PhydroPP->n[ipHELIUM];
    // carbon
    PhydroPP->n[ipCARBON]       = (CarbonAbundance )*nH*PhydroPP->Z/Zsol;
    PhydroPP->nCIII             = fCIII*PhydroPP->n[ipCARBON];
    PhydroPP->nion[ipCARBON]    = fCII *PhydroPP->n[ipCARBON];
    PhydroPP->natom[ipCARBON]   = fCI  *PhydroPP->n[ipCARBON];
    PhydroPP->nmol[ipCO]        = fCO  *PhydroPP->n[ipCARBON];
    // oxygen
    PhydroPP->n[ipOXYGEN]       = (OxygenAbundance )*nH*PhydroPP->Z/Zsol;
    PhydroPP->nOIII             = fOIII*PhydroPP->n[ipOXYGEN];
    PhydroPP->nion[ipOXYGEN]    = fOII *PhydroPP->n[ipOXYGEN];
    PhydroPP->natom[ipOXYGEN]   = fOI  *PhydroPP->n[ipOXYGEN];
    // silicon 
    PhydroPP->n[ipSILICON]      = (SiliconAbundance )*nH*PhydroPP->Z/Zsol;
    PhydroPP->nSiIII            = fSiIII*PhydroPP->n[ipSILICON];
    PhydroPP->nion[ipSILICON]   = fSiII *PhydroPP->n[ipSILICON];
    PhydroPP->natom[ipSILICON]  = fSiI  *PhydroPP->n[ipSILICON];
    // electron
    PhydroPP->nel               = CEQLibGetElectronNumberDensity(PhydroPP);

    CEQLibStatus status = CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(PhydroPP,1.e-3);
    if(status != CEQ_SUCCESS) return status;
    PhydroPP->Rho = CEQLibGetTotalMassDensity(PhydroPP);
    PhydroPP->U = CEQLibConvertTtoU(PhydroPP);
    return CEQ_SUCCESS;
}

static void main_WolfireTestBody(double G0, double Nt, double Z, double X0, char *fname){
    
    //CEQLibInitialization(1);
   
    int Nsample = 100;
    double dt = 1.e+3*YEAR_CGS; // [sec] 
    double distance = 10.0*PC_CGS; // [cm]
    double NLyn = 0.e0;
    //double NLyn = 1.e+49/(4.0*PI*SQ(distance)); // [1/s/cm^2]
    double lognHmin = -2 + log10(X0), lognHmax = 3 + log10(X0);
    //double lognHmin = -2 + log10(X0), lognHmax = 2 + log10(X0);
    double dlognH = (lognHmax-lognHmin)/Nsample; 
    double Tinit  = 1.e+1;
    if(NLyn!=0.e0)  Tinit = 1.e+4;
    
    FILE *out = fopen(fname,"w");
    if(out == NULL){
        perror(fname);
        wolfireFailures++;
        return;
    }
    for(int i=0; i<Nsample; i++){
        StructCEQData PhydroPP = {0};

        PhydroPP.dt     = dt;
        PhydroPP.Du     = 0.0;
        PhydroPP.n[ipHYDROGEN] = pow(10.0,lognHmin + i*dlognH);
        PhydroPP.T             = Tinit;    // [K]
        PhydroPP.NLyn          = NLyn;
        CEQLibStatus status = SetInitialHydroDataWolfireTest(G0,Nt,Z,X0,&PhydroPP);
        if(status != CEQ_SUCCESS){
            fprintf(stderr,"Wolfire initial state failed: %s index=%d nH=%g status=%d\n",
                    fname,i,PhydroPP.n[ipHYDROGEN],status);
            wolfireFailures++;
            continue;
        }

        //if(CEQLibRunParameters.UseDebugMode==true){
            fprintf(stderr,"Index=%d: \n",i);
            fprintf(stderr," T = %e" ,PhydroPP.T);
            fprintf(stderr," nH = %e",PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," Nt = %e",PhydroPP.Nt);
            fprintf(stderr," Z = %e" ,PhydroPP.Z);
            fprintf(stderr,"\n");
            fprintf(stderr," G0 = %e",PhydroPP.G0);
            fprintf(stderr," X0 = %e",PhydroPP.X0);
            fprintf(stderr," NLyn = %e",PhydroPP.NLyn);
            fprintf(stderr,"\n");
        //}

        double t = 0.0;
        fprintf(stderr,"  (U,T) = (%e,%e) (mu,gamm)=(%f,%f)\n",
                PhydroPP.U,PhydroPP.T,PhydroPP.MeanMolecularWeight,PhydroPP.Gamma);
        for(int it=0; t<500.0*MEGAYEAR_CGS; it++){
            
            status = CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&PhydroPP,1.e-4,&PhydroPP.Du);
            if(status != CEQ_SUCCESS) break;
            PhydroPP.U += PhydroPP.dt*PhydroPP.Du;
            t += PhydroPP.dt;

            // if L_net*dt << U, the gas reaches thermal equilibrium
            double L_net = fabs(PhydroPP.L_tot - PhydroPP.H_tot)/PhydroPP.Rho;
            if(L_net*PhydroPP.dt<1.e-5*PhydroPP.U)   break;
        }
        if(status != CEQ_SUCCESS){
            fprintf(stderr,"Wolfire step failed: %s index=%d nH=%g T=%g t=%g Myr status=%d\n",
                    fname,i,PhydroPP.n[ipHYDROGEN],PhydroPP.T,t/MEGAYEAR_CGS,status);
            wolfireFailures++;
            continue;
        }
        fprintf(stderr," ==> t = %f Myr: (U,T) = (%e,%e) (mu,gamm)=(%f,%f)\n",
                t/MEGAYEAR_CGS,PhydroPP.U,PhydroPP.T,PhydroPP.MeanMolecularWeight,PhydroPP.Gamma);

        if(CEQLibRunParameters.UseDebugMode==true){
            fprintf(stderr," ==>");
            fprintf(stderr," ne/nH = %e ",PhydroPP.nel             /PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," Hii/H = %e ",PhydroPP.nion[ipHYDROGEN]/PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," H2 /H = %e ",PhydroPP.nmol[ipH2      ]/PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Heiii/He = %e ",PhydroPP.nHeIII         /PhydroPP.n[ipHELIUM]);
            fprintf(stderr," Heii /He = %e ",PhydroPP.nion[ipHELIUM ]/PhydroPP.n[ipHELIUM]);
            fprintf(stderr," Hei  /He = %e ",PhydroPP.natom[ipHELIUM]/PhydroPP.n[ipHELIUM]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Ciii/C = %e ",PhydroPP.nCIII           /PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," Cii /C = %e ",PhydroPP.nion[ipCARBON  ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," Ci  /C = %e ",PhydroPP.natom[ipCARBON ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," CO /C  = %e ",PhydroPP.nmol[ipCO      ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Oiii/O = %e ",PhydroPP.nOIII           /PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr," Oii /O = %e ",PhydroPP.nion[ipOXYGEN  ]/PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr," Oi  /O = %e ",PhydroPP.natom[ipOXYGEN ]/PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," SiIII/Si = %e ",PhydroPP.nSiIII          /PhydroPP.n[ipSILICON]);
            fprintf(stderr," SiII /Si = %e ",PhydroPP.nion [ipSILICON]/PhydroPP.n[ipSILICON]);
            fprintf(stderr," SiI  /Si = %e ",PhydroPP.natom[ipSILICON]/PhydroPP.n[ipSILICON]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L_total     = %e ",PhydroPP.L_tot    );
            fprintf(stderr," H_total     = %e ",PhydroPP.H_tot    );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(Hi      ) = %e ",PhydroPP.L_HI     );
            fprintf(stderr," L(Hei     ) = %e ",PhydroPP.L_HeI    );
            fprintf(stderr," L(Heii    ) = %e ",PhydroPP.L_HeII   );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L_res       = %e ",PhydroPP.L_res    );
            fprintf(stderr," L_ff        = %e ",PhydroPP.L_ff     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L([Cii]158) = %e ",PhydroPP.L_CII158 );
            fprintf(stderr," L([Oi ]063) = %e ",PhydroPP.L_OI63   );
            fprintf(stderr," L([SiII]35) = %e ",PhydroPP.L_SiII35 );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(Dust Rec) = %e ",PhydroPP.L_DustRec);
            fprintf(stderr," H(PE      ) = %e ",PhydroPP.H_PE     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," H(PI      ) = %e ",PhydroPP.H_PI     );
            fprintf(stderr," H(CR      ) = %e ",PhydroPP.H_CR     );
            fprintf(stderr," H(XR      ) = %e ",PhydroPP.H_XR     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(H2      ) = %e ",PhydroPP.L_H2     );
            fprintf(stderr," H(H2      ) = %e ",PhydroPP.H_H2     );
            fprintf(stderr," L(CO      ) = %e ",PhydroPP.L_CO     );
            fprintf(stderr,"\n");
        }
       
        /* #01 */ fprintf(out,"%e ",PhydroPP.Pos[0]);   // kpc
        /* #02 */ fprintf(out,"%e ",PhydroPP.Pos[1]);
        /* #03 */ fprintf(out,"%e ",PhydroPP.Pos[2]);
        /* #04 */ fprintf(out,"%e ",PhydroPP.Vel[0]);   // km/s
        /* #05 */ fprintf(out,"%e ",PhydroPP.Vel[1]);
        /* #06 */ fprintf(out,"%e ",PhydroPP.Vel[2]);
        /* #07 */ fprintf(out,"%e ",PhydroPP.Mass  );   // Msun
        /* #08 */ fprintf(out,"%e ",PhydroPP.Rho   );   // g/cm^3
        /* #09 */ fprintf(out,"%e ",PhydroPP.T     );   // K
        /* #10 */ fprintf(out,"%e ",PhydroPP.G0    );
        /* #11 */ fprintf(out,"%e ",PhydroPP.Nt    );   // H/cm^2
        /* #12 */ fprintf(out,"%e ",PhydroPP.Z     );

        /* #13 */ fprintf(out,"%e ",PhydroPP.n[ipHYDROGEN]    );
        /* #14 */ fprintf(out,"%e ",PhydroPP.n[ipHELIUM  ]    );
        /* #15 */ fprintf(out,"%e ",PhydroPP.n[ipCARBON  ]    );
        /* #16 */ fprintf(out,"%e ",PhydroPP.n[ipOXYGEN  ]    );
        /* #17 */ fprintf(out,"%e ",PhydroPP.nel              );
        /* #18 */ fprintf(out,"%e ",PhydroPP.nion [ipHYDROGEN]);
        /* #19 */ fprintf(out,"%e ",PhydroPP.natom[ipHYDROGEN]);
        /* #10 */ fprintf(out,"%e ",PhydroPP.nHeIII           );
        /* #21 */ fprintf(out,"%e ",PhydroPP.nion [ipHELIUM  ]);
        /* #22 */ fprintf(out,"%e ",PhydroPP.natom[ipHELIUM  ]);
        /* #23 */ fprintf(out,"%e ",PhydroPP.nCIII            );
        /* #24 */ fprintf(out,"%e ",PhydroPP.nion [ipCARBON  ]);
        /* #25 */ fprintf(out,"%e ",PhydroPP.natom[ipCARBON  ]);
        /* #26 */ fprintf(out,"%e ",PhydroPP.nOIII            );
        /* #27 */ fprintf(out,"%e ",PhydroPP.nion [ipOXYGEN  ]);
        /* #28 */ fprintf(out,"%e ",PhydroPP.natom[ipOXYGEN  ]);

            // cooling rate [erg/s/cc], i.e. angular-integrated emission coefficient (= 4*pi*j)
        /* #29 */ fprintf(out,"%e ",PhydroPP.L_ff     );
        /* #30 */ fprintf(out,"%e ",PhydroPP.L_HI     );
        /* #31 */ fprintf(out,"%e ",PhydroPP.L_HeI    );
        /* #32 */ fprintf(out,"%e ",PhydroPP.L_HeII   );
        /* #33 */ fprintf(out,"%e ",PhydroPP.L_CII158 );
        /* #34 */ fprintf(out,"%e ",PhydroPP.L_CI610  );
        /* #35 */ fprintf(out,"%e ",PhydroPP.L_OI63   );
        /* #36 */ fprintf(out,"%e ",PhydroPP.L_DustRec);
        /* #37 */ fprintf(out,"%e ",PhydroPP.L_H2     );
        /* #38 */ fprintf(out,"%e ",PhydroPP.L_CO     );
            // heating rate [erg/s/cc]
        /* #39 */ fprintf(out,"%e ",PhydroPP.H_PE     );
        /* #40 */ fprintf(out,"%e ",PhydroPP.H_CR     );
        /* #41 */ fprintf(out,"%e ",PhydroPP.H_XR     );
        /* #42 */ fprintf(out,"%e ",PhydroPP.H_H2     );
       
            // total cooling / heating rates [erg/s/cc]
        /* #43 */ fprintf(out,"%e ",PhydroPP.L_tot    );
        /* #44 */ fprintf(out,"%e ",PhydroPP.H_tot    );

            // emissivity (4*pi*j) & aborpsion coefficiencies of HI 21cm
        /* #45 */ fprintf(out,"%e ",PhydroPP.L_HI21cm );
        /* #46 */ fprintf(out,"%e ",PhydroPP.k_HI21cm );
        
        fprintf(out,"\n");
        fflush(NULL);
    }
    fclose(out);

    return ; 

}


int main_WolfireTest(void){
        
    fprintf(stderr,"#%s():Line%d:%s\n",__FUNCTION__,__LINE__,__FILE__);
    
    CEQLibiSetRunParamDebugMode(true);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(true);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    CEQLibiSetRunParamH2Formation(false,1.e0);
    CEQLibiSetRunParamCOFormation(false,1.e0);
    
    CEQLibInitialization(1);
   
    char fname[100];
    //FILE *out   = fopen("./out/PostProcessMap.dat","w");
    double G0 = G0_Draine;
    double Nt = 1.e+19;
    double Z  = Zsol; 
    double X0 = 1.0;
   
#if 1 
    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
#endif

#if 1   // For Fig.4 of Wolfire et al. 1995
    G0 = G0_Draine;
    Nt = 1.e+18;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);

    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = G0_Draine;
    Nt = 1.e+20;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
#endif

#if 1   // For Fig.6 of Wolfire et al. 1995
    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = 0.01*Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = 0.03*Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = 0.1*Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);

    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = 0.3*Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);

    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = 3.0*Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
#endif

#if 1   // For Fig.7 of Wolfire et al. 1995 
    G0 = 0.3*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);

    G0 = 0.5*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = 3.0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = 10.0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    G0 = 100.0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);

    G0 = 1000.0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e.dat",Nt,G0/G0_Draine,Z/Zsol);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
#endif

#if 1   // For Fig.8 of Wolfire et al. 1995 
    X0 = 0.1;
    G0 = X0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e_Xo%.1f.dat",Nt,G0/G0_Draine,Z/Zsol,X0);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    X0 = 1.0;
    G0 = X0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e_Xo%.1f.dat",Nt,G0/G0_Draine,Z/Zsol,X0);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    X0 = 10.0;
    G0 = X0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e_Xo%.1f.dat",Nt,G0/G0_Draine,Z/Zsol,X0);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    X0 = 100.0;
    G0 = X0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e_Xo%.1f.dat",Nt,G0/G0_Draine,Z/Zsol,X0);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
    
    X0 = 1000.0;
    G0 = X0*G0_Draine;
    Nt = 1.e+19;
    Z  = Zsol; 
    sprintf(fname,"./out/PostProcessMap_Nt%02.1e_ISRF%.1f_Zs%02.1e_Xo%.1f.dat",Nt,G0/G0_Draine,Z/Zsol,X0);
    main_WolfireTestBody(G0,Nt,Z,X0,fname);
#endif

    return wolfireFailures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

    
int main_WolfireMilkyWayTest(void){
    
    CEQLibiSetRunParamDebugMode(true);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(true);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamH2Formation(false,1.e0);
    CEQLibiSetRunParamCOFormation(false,1.e0);
    
    CEQLibInitialization(1);

    //// Milky Way
    //G0 = 10.0*G0_Draine;
    //Nt = 1.e+20;
    //Z  = Zsol; 
    //sprintf(fname,"./out/PostProcessMap_Nt%02.1e_Go%02.1e_Zs%02.1e.dat",Nt,G0,Z/Zsol);
    //main_WolfireTestBody(G0,Nt,Z,fname);
 
    return wolfireFailures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main_HIIregionTest(void){
    
    CEQLibiSetRunParamDebugMode(true);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(false);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    CEQLibiSetRunParamH2Formation(false,1.e0);
    CEQLibiSetRunParamCOFormation(false,1.e0);
    
    CEQLibInitialization(1);
   
    char fname[100];
    double G0 = G0_Draine;
    double Nt = 1.e+19;
    double Z  = Zsol; 
    double X0 = 1.0;
    
    sprintf(fname,"./out/PostProcessMap.dat",Nt,G0/G0_Draine,Z/Zsol);
   
    //int Nsample_nH  = 5;
    //double lognHmin = -2, lognHmax = 3;
    //int Nsample_nH  = 3;
    //double lognHmin = -2, lognHmax = 2;
    int Nsample_nH  = 1;
    //double lognHmin = 1, lognHmax = 1;
    double lognHmin = 2, lognHmax = 2;
    double dlognH = (lognHmax-lognHmin)/Nsample_nH; 
    int Nsample_tmp = 100;
    double logTmin = 2.0, logTmax = 6.0;
    double dlogT  = (logTmax-logTmin)/Nsample_tmp; 
    int Nsample = Nsample_nH * Nsample_tmp;
    
    double distance = 10.0*PC_CGS; // [cm]
    //double NLyn = 0.e0;
    double NLyn = 1.e+49/(4.0*PI*SQ(distance)); // [1/s/cm^2]

    double dt = 1.e+3*YEAR_CGS; // [sec] 
    StructCEQData PhydroPP;
    FILE *out   = fopen(fname,"w");
    for(int inH =0; inH <Nsample_nH ; inH ++){
    for(int itmp=0; itmp<Nsample_tmp; itmp++){
        int i = inH + itmp*Nsample_nH;
        //int inH  = i/Nsample_tmp;
        //int itmp = i/Nsample_nH ;

        PhydroPP.t             = 0.0;
        PhydroPP.dt            = dt;
        PhydroPP.Du            = 0.0;
        //PhydroPP.n[ipHYDROGEN] = pow(10.0,lognHmin + i*dlognH);
        //PhydroPP.T             = Tinit;    // [K]
        PhydroPP.n[ipHYDROGEN] = pow(10.0,lognHmin + inH *dlognH);
        PhydroPP.T             = pow(10.0,logTmin  + itmp*dlogT );
        PhydroPP.NLyn          = NLyn;
        CEQLibStatus status = SetInitialHydroDataWolfireTest(G0,Nt,Z,X0,&PhydroPP);
        if(status != CEQ_SUCCESS){
            fprintf(stderr,"Initial state failed: %s index=%d status=%d\n",fname,i,status);
            fclose(out);
            return EXIT_FAILURE;
        }

        fprintf(stderr,"Index=%d: \n",i);
        fprintf(stderr," T = %e" ,PhydroPP.T);
        fprintf(stderr," nH = %e",PhydroPP.n[ipHYDROGEN]);
        fprintf(stderr," Nt = %e",PhydroPP.Nt);
        fprintf(stderr," Z  = %e" ,PhydroPP.Z);
        fprintf(stderr," Zd = %e" ,PhydroPP.Zd);
        fprintf(stderr,"\n");
        fprintf(stderr," G0 = %e",PhydroPP.G0);
        fprintf(stderr," X0 = %e",PhydroPP.X0);
        fprintf(stderr," NLyn = %e",PhydroPP.NLyn);
        fprintf(stderr,"\n");

        //fprintf(stderr,"  (U,T) = (%e,%e) (mu,gamm)=(%f,%f)\n",
        //        PhydroPP.U,PhydroPP.T,PhydroPP.MeanMolecularWeight,PhydroPP.Gamma);
        //for(int it=0; PhydroPP.t<500.0*MEGAYEAR_CGS; it++){
        //    
        //    PhydroPP.Du = CEQLibEnergyEquationSolverChemcalEquilibriumCoolingHeating(&PhydroPP,1.e-4);
        //    PhydroPP.U += PhydroPP.dt*PhydroPP.Du;
        //    PhydroPP.t += PhydroPP.dt;
    
        //    //CEQLibChemicalReactionEquilibrium(PhydroPP,1.e-3);
        //    //CEQLibCalculateCoolingHeatingRates(PhydroPP);

        //    // if L_net*dt << U, the gas reaches thermal equilibrium
        //    double L_net = fabs(PhydroPP.L_tot - PhydroPP.H_tot)/PhydroPP.Rho;
        //    if(L_net*PhydroPP.dt<1.e-5*PhydroPP.U)   break;
        //}
        //fprintf(stderr," ==> t = %f Myr: (U,T) = (%e,%e) (mu,gamm)=(%f,%f)\n",
        //        PhydroPP.t/MEGAYEAR_CGS,PhydroPP.U,PhydroPP.T,PhydroPP.MeanMolecularWeight,PhydroPP.Gamma);

        if(CEQLibRunParameters.UseDebugMode==true){
            fprintf(stderr," ==>");
            fprintf(stderr," ne/nH = %e ",PhydroPP.nel             /PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," Hii/H = %e ",PhydroPP.nion[ipHYDROGEN]/PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," H2 /H = %e ",PhydroPP.nmol[ipH2      ]/PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Heiii/He = %e ",PhydroPP.nHeIII         /PhydroPP.n[ipHELIUM]);
            fprintf(stderr," Heii /He = %e ",PhydroPP.nion[ipHELIUM ]/PhydroPP.n[ipHELIUM]);
            fprintf(stderr," Hei  /He = %e ",PhydroPP.natom[ipHELIUM]/PhydroPP.n[ipHELIUM]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Ciii/C = %e ",PhydroPP.nCIII           /PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," Cii /C = %e ",PhydroPP.nion[ipCARBON  ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," Ci  /C = %e ",PhydroPP.natom[ipCARBON ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," CO /C  = %e ",PhydroPP.nmol[ipCO      ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Oiii/O = %e ",PhydroPP.nOIII           /PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr," Oii /O = %e ",PhydroPP.nion[ipOXYGEN  ]/PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr," Oi  /O = %e ",PhydroPP.natom[ipOXYGEN ]/PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," SiIII/Si = %e ",PhydroPP.nSiIII          /PhydroPP.n[ipSILICON]);
            fprintf(stderr," SiII /Si = %e ",PhydroPP.nion [ipSILICON]/PhydroPP.n[ipSILICON]);
            fprintf(stderr," SiI  /Si = %e ",PhydroPP.natom[ipSILICON]/PhydroPP.n[ipSILICON]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L_total     = %e ",PhydroPP.L_tot    );
            fprintf(stderr," H_total     = %e ",PhydroPP.H_tot    );
            fprintf(stderr," L_net       = %e ",PhydroPP.L_tot-PhydroPP.H_tot);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(Hi      ) = %e ",PhydroPP.L_HI     );
            fprintf(stderr," L(Hei     ) = %e ",PhydroPP.L_HeI    );
            fprintf(stderr," L(Heii    ) = %e ",PhydroPP.L_HeII   );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L_res       = %e ",PhydroPP.L_res    );
            fprintf(stderr," L_ff        = %e ",PhydroPP.L_ff     );
            //fprintf(stderr,"\n");
            //fprintf(stderr," ==>");
            //fprintf(stderr," L([Oiii]  ) = %e ",PhydroPP.L_OIII   );
            //fprintf(stderr," L([Oii ]  ) = %e ",PhydroPP.L_OII    );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L([Cii]158) = %e ",PhydroPP.L_CII158 );
            fprintf(stderr," L([Oi ]063) = %e ",PhydroPP.L_OI63   );
            fprintf(stderr," L([SiII]35) = %e ",PhydroPP.L_SiII35 );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(Dust Rec) = %e ",PhydroPP.L_DustRec);
            fprintf(stderr," H(PE      ) = %e ",PhydroPP.H_PE     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," H(PI      ) = %e ",PhydroPP.H_PI     );
            fprintf(stderr," H(CR      ) = %e ",PhydroPP.H_CR     );
            fprintf(stderr," H(XR      ) = %e ",PhydroPP.H_XR     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(H2      ) = %e ",PhydroPP.L_H2     );
            fprintf(stderr," H(H2      ) = %e ",PhydroPP.H_H2     );
            fprintf(stderr," L(CO      ) = %e ",PhydroPP.L_CO     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(HI21cm  ) = %e ",PhydroPP.L_HI21cm );
            fprintf(stderr," L(Halpha  ) = %e ",PhydroPP.L_Halpha );
            fprintf(stderr,"\n");
        }
       
        /* #01 */ fprintf(out,"%e ",PhydroPP.Pos[0]);   // kpc
        /* #02 */ fprintf(out,"%e ",PhydroPP.Pos[1]);
        /* #03 */ fprintf(out,"%e ",PhydroPP.Pos[2]);
        /* #04 */ fprintf(out,"%e ",PhydroPP.Vel[0]);   // km/s
        /* #05 */ fprintf(out,"%e ",PhydroPP.Vel[1]);
        /* #06 */ fprintf(out,"%e ",PhydroPP.Vel[2]);
        /* #07 */ fprintf(out,"%e ",PhydroPP.Mass  );   // Msun
        /* #08 */ fprintf(out,"%e ",PhydroPP.Rho   );   // g/cm^3
        /* #09 */ fprintf(out,"%e ",PhydroPP.T     );   // K
        /* #10 */ fprintf(out,"%e ",PhydroPP.G0    );
        /* #11 */ fprintf(out,"%e ",PhydroPP.Nt    );   // H/cm^2
        /* #12 */ fprintf(out,"%e ",PhydroPP.Z     );

        /* #13 */ fprintf(out,"%e ",PhydroPP.n[ipHYDROGEN]    );
        /* #14 */ fprintf(out,"%e ",PhydroPP.n[ipHELIUM  ]    );
        /* #15 */ fprintf(out,"%e ",PhydroPP.n[ipCARBON  ]    );
        /* #16 */ fprintf(out,"%e ",PhydroPP.n[ipOXYGEN  ]    );
        /* #17 */ fprintf(out,"%e ",PhydroPP.nel              );
        /* #18 */ fprintf(out,"%e ",PhydroPP.nion [ipHYDROGEN]);
        /* #19 */ fprintf(out,"%e ",PhydroPP.natom[ipHYDROGEN]);
        /* #20 */ fprintf(out,"%e ",PhydroPP.nHeIII           );
        /* #21 */ fprintf(out,"%e ",PhydroPP.nion [ipHELIUM  ]);
        /* #22 */ fprintf(out,"%e ",PhydroPP.natom[ipHELIUM  ]);
        /* #23 */ fprintf(out,"%e ",PhydroPP.nCIII            );
        /* #24 */ fprintf(out,"%e ",PhydroPP.nion [ipCARBON  ]);
        /* #25 */ fprintf(out,"%e ",PhydroPP.natom[ipCARBON  ]);
        /* #26 */ fprintf(out,"%e ",PhydroPP.nOIII            );
        /* #27 */ fprintf(out,"%e ",PhydroPP.nion [ipOXYGEN  ]);
        /* #28 */ fprintf(out,"%e ",PhydroPP.natom[ipOXYGEN  ]);

            // cooling rate [erg/s/cc], i.e. angular-integrated emission coefficient (= 4*pi*j)
        /* #29 */ fprintf(out,"%e ",PhydroPP.L_res    );
        /* #30 */ fprintf(out,"%e ",PhydroPP.L_ff     );
        /* #31 */ fprintf(out,"%e ",PhydroPP.L_HI     );
        /* #32 */ fprintf(out,"%e ",PhydroPP.L_HeI    );
        /* #33 */ fprintf(out,"%e ",PhydroPP.L_HeII   );
        /* #34 */ fprintf(out,"%e ",PhydroPP.L_CII158 );
        /* #35 */ fprintf(out,"%e ",PhydroPP.L_OIII   );
        /* #36 */ fprintf(out,"%e ",PhydroPP.L_OII    );
        /* #37 */ fprintf(out,"%e ",PhydroPP.L_OI63   );
        /* #38 */ fprintf(out,"%e ",PhydroPP.L_DustRec);
        /* #39 */ fprintf(out,"%e ",PhydroPP.L_H2     );
        /* #40 */ fprintf(out,"%e ",PhydroPP.L_CO     );
            // heating rate [erg/s/cc]
        /* #41 */ fprintf(out,"%e ",PhydroPP.H_PI     );
        /* #42 */ fprintf(out,"%e ",PhydroPP.H_PE     );
        /* #43 */ fprintf(out,"%e ",PhydroPP.H_CR     );
        /* #44 */ fprintf(out,"%e ",PhydroPP.H_XR     );
        /* #45 */ fprintf(out,"%e ",PhydroPP.H_H2     );
       
            // total cooling / heating rates [erg/s/cc]
        /* #46 */ fprintf(out,"%e ",PhydroPP.L_tot    );
        /* #47 */ fprintf(out,"%e ",PhydroPP.H_tot    );

            // emissivity (4*pi*j) & aborpsion coefficiencies of HI 21cm
        /* #48 */ fprintf(out,"%e ",PhydroPP.L_HI21cm );
        /* #49 */ fprintf(out,"%e ",PhydroPP.k_HI21cm );
            
            // emissivity (4*pi*j) of H_alpha emission
        /* #50 */ fprintf(out,"%e ",PhydroPP.L_Halpha );
        
        fprintf(out,"\n");
        fflush(NULL);
    }
    }
    fclose(out);

    return ; 

}

int main_CoolingFunctionTest(void){
    
    char fname[100];
    
    ///////////// CIE Cooling Fuction
    CEQLibiSetRunParamDebugMode(true);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(false);
    CEQLibiSetRunParamGrainSurfaceProcesses(false);
    CEQLibiSetRunParamCosmicRayIonization(false);
    CEQLibiSetRunParamSoftXrayIonization(false);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    CEQLibiSetRunParamH2Formation(false,1.e0);
    CEQLibiSetRunParamCOFormation(false,1.e0);
    
    CEQLibInitialization(1);

    double X0 = 1.0;

    double ModelnH[] = {1.e-2,1.e-1,1.e+0,1.e+1,1.e+2};
    double ModelZ[]  = {0.01,0.03,0.1,0.3,1.0,3.0};
    double ModelNt[] = {1.e+18,1.e+19,1.e+20,1.e+21};
    double ModelG0[] = {0.01,0.1,1.0,10.0,100.0,1000.0};
    int nnH = sizeof(ModelnH)/sizeof(double);
    int nZ  = sizeof(ModelZ )/sizeof(double);
    int nNt = sizeof(ModelNt)/sizeof(double);
    int nG0 = sizeof(ModelG0)/sizeof(double);

    for(int iz=0; iz<nZ; iz++){
        double Z  = ModelZ[iz]*Zsol; 
        double nH = 10.0;
        double G0 = 1.e-20; 
        double Nt = 1.e+18;
        sprintf(fname,"./out/CIECoolingFunction_nH%02.1e_Nt%02.1e_Zs%02.1e.dat",nH,Nt,Z/Zsol);
        main_CoolingFunctionTestBody(nH,G0,Nt,Z,X0,fname);
    }
   
    ///////////// PDR Cooling Function
    CEQLibiSetRunParamDebugMode(true);
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(false);
    CEQLibiSetRunParamGrainSurfaceProcesses(false);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    CEQLibiSetRunParamH2Formation(true,1.e0);
    CEQLibiSetRunParamCOFormation(false,1.e0);
    
    CEQLibInitialization(1);

    for(int iz=0; iz<nZ; iz++){
        double Z  = ModelZ[iz]*Zsol; 
        for(int ig=0; ig<nG0; ig++){
            double G0 = ModelG0[ig]*G0_Draine;
            for(int iN=0; iN<nNt; iN++){
                double Nt = ModelNt[iN];
                for(int in=0; in<nnH; in++){
                    double nH = ModelnH[in];
                    sprintf(fname,"./out/PDRCoolingFunction_nH%02.1e_Nt%02.1e_ISRF%02.1e_Zs%02.1e.dat",nH,Nt,G0/G0_Draine,Z/Zsol);
                    main_CoolingFunctionTestBody(nH,G0,Nt,Z,X0,fname);
                }
            }
        }
    }

    return ;
} 

static void main_CoolingFunctionTestBody(double nH, double G0, double Nt, double Z, double X0, char *fname){
    
    //int Nsample_nH  = 5;
    //double lognHmin = -2, lognHmax = 3;
    //int Nsample_nH  = 3;
    //double lognHmin = -2, lognHmax = 2;
    int Nsample_nH  = 1;
    //double lognHmin = 1, lognHmax = 1;
    //double lognHmin = 2, lognHmax = 2;
    //double dlognH = (lognHmax-lognHmin)/Nsample_nH; 
    int Nsample_tmp = 100;
    double logTmin = 1.0, logTmax = 8.0;
    double dlogT  = (logTmax-logTmin)/Nsample_tmp; 
    int Nsample = Nsample_nH * Nsample_tmp;
    
    double distance = 10.0*PC_CGS; // [cm]
    double NLyn = 0.e0;
    //double NLyn = 1.e+49/(4.0*PI*SQ(distance)); // [1/s/cm^2]

    double dt = 1.e+3*YEAR_CGS; // [sec] 
    StructCEQData PhydroPP;
    FILE *out   = fopen(fname,"w");
    for(int itmp=0; itmp<Nsample_tmp; itmp++){
        int i = itmp;
        //int i = inH + itmp*Nsample_nH;
        //int inH  = i/Nsample_tmp;
        //int itmp = i/Nsample_nH ;

        PhydroPP.t             = 0.0;
        PhydroPP.dt            = dt;
        PhydroPP.Du            = 0.0;
        //PhydroPP.n[ipHYDROGEN] = pow(10.0,lognHmin + i*dlognH);
        //PhydroPP.T             = Tinit;    // [K]
        PhydroPP.n[ipHYDROGEN] = nH; //pow(10.0,lognHmin + inH *dlognH);
        PhydroPP.T             = pow(10.0,logTmin  + itmp*dlogT );
        PhydroPP.NLyn          = NLyn;
        CEQLibStatus status = SetInitialHydroDataWolfireTest(G0,Nt,Z,X0,&PhydroPP);
        if(status != CEQ_SUCCESS){
            fprintf(stderr,"Initial state failed: %s index=%d status=%d\n",fname,i,status);
            fclose(out);
            return;
        }

        fprintf(stderr,"Index=%d: \n",i);
        fprintf(stderr," T = %e" ,PhydroPP.T);
        fprintf(stderr," nH = %e",PhydroPP.n[ipHYDROGEN]);
        fprintf(stderr," Nt = %e",PhydroPP.Nt);
        fprintf(stderr," Z  = %e" ,PhydroPP.Z);
        fprintf(stderr," Zd = %e" ,PhydroPP.Zd);
        fprintf(stderr,"\n");
        fprintf(stderr," G0 = %e",PhydroPP.G0);
        fprintf(stderr," X0 = %e",PhydroPP.X0);
        fprintf(stderr," NLyn = %e",PhydroPP.NLyn);
        fprintf(stderr,"\n");

        //fprintf(stderr,"  (U,T) = (%e,%e) (mu,gamm)=(%f,%f)\n",
        //        PhydroPP.U,PhydroPP.T,PhydroPP.MeanMolecularWeight,PhydroPP.Gamma);
        //for(int it=0; PhydroPP.t<500.0*MEGAYEAR_CGS; it++){
        //    
        //    PhydroPP.Du = CEQLibEnergyEquationSolverChemcalEquilibriumCoolingHeating(&PhydroPP,1.e-4);
        //    PhydroPP.U += PhydroPP.dt*PhydroPP.Du;
        //    PhydroPP.t += PhydroPP.dt;
    
        //    //CEQLibChemicalReactionEquilibrium(PhydroPP,1.e-3);
        //    //CEQLibCalculateCoolingHeatingRates(PhydroPP);

        //    // if L_net*dt << U, the gas reaches thermal equilibrium
        //    double L_net = fabs(PhydroPP.L_tot - PhydroPP.H_tot)/PhydroPP.Rho;
        //    if(L_net*PhydroPP.dt<1.e-5*PhydroPP.U)   break;
        //}
        //fprintf(stderr," ==> t = %f Myr: (U,T) = (%e,%e) (mu,gamm)=(%f,%f)\n",
        //        PhydroPP.t/MEGAYEAR_CGS,PhydroPP.U,PhydroPP.T,PhydroPP.MeanMolecularWeight,PhydroPP.Gamma);

        if(CEQLibRunParameters.UseDebugMode==true){
            fprintf(stderr," ==>");
            fprintf(stderr," ne/nH = %e ",PhydroPP.nel             /PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," Hii/H = %e ",PhydroPP.nion[ipHYDROGEN]/PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr," H2 /H = %e ",PhydroPP.nmol[ipH2      ]/PhydroPP.n[ipHYDROGEN]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Heiii/He = %e ",PhydroPP.nHeIII         /PhydroPP.n[ipHELIUM]);
            fprintf(stderr," Heii /He = %e ",PhydroPP.nion[ipHELIUM ]/PhydroPP.n[ipHELIUM]);
            fprintf(stderr," Hei  /He = %e ",PhydroPP.natom[ipHELIUM]/PhydroPP.n[ipHELIUM]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Ciii/C = %e ",PhydroPP.nCIII           /PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," Cii /C = %e ",PhydroPP.nion[ipCARBON  ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," Ci  /C = %e ",PhydroPP.natom[ipCARBON ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr," CO /C  = %e ",PhydroPP.nmol[ipCO      ]/PhydroPP.n[ipCARBON  ]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," Oiii/O = %e ",PhydroPP.nOIII           /PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr," Oii /O = %e ",PhydroPP.nion[ipOXYGEN  ]/PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr," Oi  /O = %e ",PhydroPP.natom[ipOXYGEN ]/PhydroPP.n[ipOXYGEN  ]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," SiIII/Si = %e ",PhydroPP.nSiIII          /PhydroPP.n[ipSILICON]);
            fprintf(stderr," SiII /Si = %e ",PhydroPP.nion [ipSILICON]/PhydroPP.n[ipSILICON]);
            fprintf(stderr," SiI  /Si = %e ",PhydroPP.natom[ipSILICON]/PhydroPP.n[ipSILICON]);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L_total     = %e ",PhydroPP.L_tot    );
            fprintf(stderr," H_total     = %e ",PhydroPP.H_tot    );
            fprintf(stderr," L_net       = %e ",PhydroPP.L_tot-PhydroPP.H_tot);
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(Hi      ) = %e ",PhydroPP.L_HI     );
            fprintf(stderr," L(Hei     ) = %e ",PhydroPP.L_HeI    );
            fprintf(stderr," L(Heii    ) = %e ",PhydroPP.L_HeII   );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L_res       = %e ",PhydroPP.L_res    );
            fprintf(stderr," L_ff        = %e ",PhydroPP.L_ff     );
            //fprintf(stderr,"\n");
            //fprintf(stderr," ==>");
            //fprintf(stderr," L([Oiii]  ) = %e ",PhydroPP.L_OIII   );
            //fprintf(stderr," L([Oii ]  ) = %e ",PhydroPP.L_OII    );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L([Cii]158) = %e ",PhydroPP.L_CII158 );
            fprintf(stderr," L([Oi ]063) = %e ",PhydroPP.L_OI63   );
            fprintf(stderr," L([SiII]35) = %e ",PhydroPP.L_SiII35 );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(Dust Rec) = %e ",PhydroPP.L_DustRec);
            fprintf(stderr," H(PE      ) = %e ",PhydroPP.H_PE     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," H(PI      ) = %e ",PhydroPP.H_PI     );
            fprintf(stderr," H(CR      ) = %e ",PhydroPP.H_CR     );
            fprintf(stderr," H(XR      ) = %e ",PhydroPP.H_XR     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(H2      ) = %e ",PhydroPP.L_H2     );
            fprintf(stderr," H(H2      ) = %e ",PhydroPP.H_H2     );
            fprintf(stderr," L(CO      ) = %e ",PhydroPP.L_CO     );
            fprintf(stderr,"\n");
            fprintf(stderr," ==>");
            fprintf(stderr," L(HI21cm  ) = %e ",PhydroPP.L_HI21cm );
            fprintf(stderr," L(Halpha  ) = %e ",PhydroPP.L_Halpha );
            fprintf(stderr,"\n");
        }
       
        /* #01 */ fprintf(out,"%e ",PhydroPP.Pos[0]);   // kpc
        /* #02 */ fprintf(out,"%e ",PhydroPP.Pos[1]);
        /* #03 */ fprintf(out,"%e ",PhydroPP.Pos[2]);
        /* #04 */ fprintf(out,"%e ",PhydroPP.Vel[0]);   // km/s
        /* #05 */ fprintf(out,"%e ",PhydroPP.Vel[1]);
        /* #06 */ fprintf(out,"%e ",PhydroPP.Vel[2]);
        /* #07 */ fprintf(out,"%e ",PhydroPP.Mass  );   // Msun
        /* #08 */ fprintf(out,"%e ",PhydroPP.Rho   );   // g/cm^3
        /* #09 */ fprintf(out,"%e ",PhydroPP.T     );   // K
        /* #10 */ fprintf(out,"%e ",PhydroPP.G0    );
        /* #11 */ fprintf(out,"%e ",PhydroPP.Nt    );   // H/cm^2
        /* #12 */ fprintf(out,"%e ",PhydroPP.Z     );

        /* #13 */ fprintf(out,"%e ",PhydroPP.n[ipHYDROGEN]    );
        /* #14 */ fprintf(out,"%e ",PhydroPP.n[ipHELIUM  ]    );
        /* #15 */ fprintf(out,"%e ",PhydroPP.n[ipCARBON  ]    );
        /* #16 */ fprintf(out,"%e ",PhydroPP.n[ipOXYGEN  ]    );
        /* #17 */ fprintf(out,"%e ",PhydroPP.nel              );
        /* #18 */ fprintf(out,"%e ",PhydroPP.nion [ipHYDROGEN]);
        /* #19 */ fprintf(out,"%e ",PhydroPP.natom[ipHYDROGEN]);
        /* #20 */ fprintf(out,"%e ",PhydroPP.nHeIII           );
        /* #21 */ fprintf(out,"%e ",PhydroPP.nion [ipHELIUM  ]);
        /* #22 */ fprintf(out,"%e ",PhydroPP.natom[ipHELIUM  ]);
        /* #23 */ fprintf(out,"%e ",PhydroPP.nCIII            );
        /* #24 */ fprintf(out,"%e ",PhydroPP.nion [ipCARBON  ]);
        /* #25 */ fprintf(out,"%e ",PhydroPP.natom[ipCARBON  ]);
        /* #26 */ fprintf(out,"%e ",PhydroPP.nOIII            );
        /* #27 */ fprintf(out,"%e ",PhydroPP.nion [ipOXYGEN  ]);
        /* #28 */ fprintf(out,"%e ",PhydroPP.natom[ipOXYGEN  ]);

            // cooling rate [erg/s/cc], i.e. angular-integrated emission coefficient (= 4*pi*j)
        /* #29 */ fprintf(out,"%e ",PhydroPP.L_res    );
        /* #30 */ fprintf(out,"%e ",PhydroPP.L_ff     );
        /* #31 */ fprintf(out,"%e ",PhydroPP.L_HI     );
        /* #32 */ fprintf(out,"%e ",PhydroPP.L_HeI    );
        /* #33 */ fprintf(out,"%e ",PhydroPP.L_HeII   );
        /* #34 */ fprintf(out,"%e ",PhydroPP.L_CII158 );
        /* #35 */ fprintf(out,"%e ",PhydroPP.L_OIII   );
        /* #36 */ fprintf(out,"%e ",PhydroPP.L_OII    );
        /* #37 */ fprintf(out,"%e ",PhydroPP.L_OI63   );
        /* #38 */ fprintf(out,"%e ",PhydroPP.L_DustRec);
        /* #39 */ fprintf(out,"%e ",PhydroPP.L_H2     );
        /* #40 */ fprintf(out,"%e ",PhydroPP.L_CO     );
            // heating rate [erg/s/cc]
        /* #41 */ fprintf(out,"%e ",PhydroPP.H_PI     );
        /* #42 */ fprintf(out,"%e ",PhydroPP.H_PE     );
        /* #43 */ fprintf(out,"%e ",PhydroPP.H_CR     );
        /* #44 */ fprintf(out,"%e ",PhydroPP.H_XR     );
        /* #45 */ fprintf(out,"%e ",PhydroPP.H_H2     );
       
            // total cooling / heating rates [erg/s/cc]
        /* #46 */ fprintf(out,"%e ",PhydroPP.L_tot    );
        /* #47 */ fprintf(out,"%e ",PhydroPP.H_tot    );

            // emissivity (4*pi*j) & aborpsion coefficiencies of HI 21cm
        /* #48 */ fprintf(out,"%e ",PhydroPP.L_HI21cm );
        /* #49 */ fprintf(out,"%e ",PhydroPP.k_HI21cm );
            
            // emissivity (4*pi*j) of H_alpha emission
        /* #50 */ fprintf(out,"%e ",PhydroPP.L_Halpha );
        
        fprintf(out,"\n");
        fflush(NULL);
    }
    fclose(out);

    return ; 

}


