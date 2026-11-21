#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include "Astro.h"
#include "Constants.h"
#include "Utilities.h"
#include "NewtonRaphson.h"
#include "ChemicalEquilibriumCooling.h"
#include "ChemicalElements.h"
#include "EquilibriumAbundance.h"
#include "CollisionalGasPhaseChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#include "FineStructureLines.h"
#include "HyperFineStructureLines.h"
#include "CosmicRay.h"
#include "Xray.h"
#include "Bremsstrahlung.h"
#include "CoolingCIE.h"
#include "HIIregionLineCooling.h"

struct StructCEQLibParameters_tag CEQLibRunParameters;

//#ifdef LSODA_SUBCYCLE_COOLING_SOLVER
//#include "lsoda.h"
//enum{
//    iEner,  //0
//    ineq
//};
//static int istate;
//static int itol;
//static int itask;
//static int iopt;
//static int jt;
//static int iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
//static double rwork1, rwork5, rwork6, rwork7;
//static void InitializeLSODA(int MaximunInteration);
//static void RateEquations(double t, double *n, double *ndot, struct StructPhydroPP_tag *PhydroPPwork);
//static double EnergyEquationSolverChemcalEquilibriumCoolingHeatingLSODASubcycle(struct StructPhydroPP_tag *PhydroPP, double Accuracy);
//#endif

#define NumberOfTranslationalDegreesOfFreedom (3.0)
#define NumberOfFreedomForMolecules          (5.0) // translation + rotation


static bool EnsureDiagnosticDirectory(const char *path){
    if(mkdir(path,0755) == 0) return true;
    struct stat info;
    return errno == EEXIST && stat(path,&info) == 0 &&
           S_ISDIR(info.st_mode) && access(path,W_OK | X_OK) == 0;
}

CEQLibStatus CEQLibTryInitialization(int Check, const char *pegaseDirectory){
    if(Check == 1 && (!EnsureDiagnosticDirectory(CEQ_DIAGNOSTIC_DIRECTORY) ||
            (CEQLibRunParameters.UseHIIRegion && !EnsureDiagnosticDirectory("out")))){
        fprintf(stderr,"CEQLib could not create its diagnostic output directories\n");
        return CEQ_INITIALIZATION_FAILED;
    }
    const CEQLibStatus status = CEQInitializePhotochemistry(Check,pegaseDirectory);
    if(status != CEQ_SUCCESS) return status;

    InitializeChemicalElementTable(Check);
    InitializeCollisionalIonization(Check);
    InitializeMoleculeChemistory(Check);
    InitializeFineStructureLines(Check);
    InitializeGrainSurfaceChemistry(Check);
    InitializeCosmicRayIonization(Check);
    InitializeHyperFineStructureLines(Check);
    InitializeSoftXrayIonization(Check);
    InitializeCoolingLinesHIIregion(Check);
    InitializeCoolingTableCIE(Check);
    
    if(Check==1)
        fprintf(stderr,"... done.\n");

#ifdef LSODA_SUBCYCLE_COOLING_SOLVER
    InitializeLSODA(10000);
#endif

    return CEQ_SUCCESS;
}

void CEQLibInitialization(int Check){
    if(CEQLibTryInitialization(Check,NULL) != CEQ_SUCCESS) exit(EXIT_FAILURE);
}

void CEQLibiSetRunParamDebugMode(bool Use){
    CEQLibRunParameters.UseDebugMode = Use;
}

void CEQLibiSetRunParamRecombinationCaseB(bool Use){
    CEQLibRunParameters.UseRecombinationCaseB = Use;
}

void CEQLibiSetRunParamHIIRegion(bool Use){
    CEQLibRunParameters.UseHIIRegion = Use;
}

void CEQLibiSetRunParamDustAttenuation(bool Use){
    CEQLibRunParameters.UseDustAttenuation = Use;
}

void CEQLibiSetRunParamGrainSurfaceProcesses(bool Use){
    CEQLibRunParameters.UseGrainSurfaceProcesses = Use;
}

void CEQLibiSetRunParamCosmicRayIonization(bool Use){
    CEQLibRunParameters.UseCosmicRayIonization = Use;
}

void CEQLibiSetRunParamSoftXrayIonization(bool Use){
    CEQLibRunParameters.UseSoftXrayIonization = Use;
}

void CEQLibiSetRunParamCIEResonanceCooling(bool Use){
    CEQLibRunParameters.UseCIEResonanceCooling = Use;
}

void CEQLibiSetRunParamH2Formation(bool Use, double Cp){
    CEQLibRunParameters.UseH2Formation = Use;
    CEQLibRunParameters.CpForH2 = Cp;
}

void CEQLibiSetRunParamCOFormation(bool Use, double Cp){
    CEQLibRunParameters.UseCOFormation = Use;
    CEQLibRunParameters.CpForCO = Cp;
}

//#ifdef LSODA_SUBCYCLE_COOLING_SOLVER
//static void RateEquations(double t, double *n, double *ndot, struct StructPhydroPP_tag *PhydroPPwork){
//   
//    // energy equation [erg/g/s] 
//    ndot[iEner] = NetCoolingRateChemcalEquilibrium(PhydroPPwork);
//
//    //iteration ++;
//
//    return ;
//}
//
//static double EnergyEquationSolverChemcalEquilibriumCoolingHeatingLSODASubcycle(struct StructPhydroPP_tag *PhydroPP, double Accuracy){
//
//    istate = 1;
//
//    double dt   = PhydroPP->dt;
//    double Uini = PhydroPP->U ; 
//    /// initialize
//	double t = 0.0E0;
//	double tout = t + dt;
//	double y[ineq+1], yi[ineq+1];
//    double rtol[ineq+1], atol[ineq+1];
//    double EnvData[1];
//
//    y[1+iEner] = PhydroPP->U;
//    for(int i=0; i<ineq+1; i++){
//        yi[i]   = y[i];
//        rtol[i] = Accuracy;
//        atol[i] = 0.e0;
//    }
//    
//    StructCEQData PhydroPPwork = *PhydroPP;
//    
//    /// LSODA
//    lsoda(RateEquations,ineq,y,&t,tout,itol,rtol,atol,itask,&istate,iopt,jt,
//          iwork1,iwork2,iwork5,iwork6,iwork7,iwork8,iwork9,
//          rwork1,rwork5,rwork6,rwork7,EnvData);
//	n_lsoda_terminate();
//
//    PhydroPP->U = y[1+iEner];
//    
//    return (PhydroPP->U-Uini)/dt;
//
//    //int status = true; //ErrorCheckForLsoda(index,t,tout,yi,y,EnvData,istate,Accuracy);
//    //if(status==true){
//    //    /// update
//    //    PhydroPP->U = y[1+iEner];
//    //}
//    //else{
//    //    //exit(1);
//    //    fprintf(stderr," Repeat with Accuracy =%e\n",0.1*Accuracy);
//    //    //if(Accuracy<=1.e-10)    exit(1);
//    //    int iter_ch = EnergyEquationSolverChemcalEquilibriumCoolingHeatingLSODASubcycle(PhydroPP,0.1*Accuracy);
//    //}
//}
//#endif

CEQLibStatus CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(StructCEQData *state, double accuracy){
    if(state == NULL) return CEQ_INVALID_INPUT;

    StructCEQData trial = *state;
    int status = CEQLibChemicalReactionEquilibrium(&trial,accuracy);
    if(status < 0) return (CEQLibStatus)status;

    trial.invT = 1.0/trial.T;
    trial.logT = log10(trial.T);
    status = CEQLibCalculateCoolingHeatingRates(&trial);
    if(status != CEQ_SUCCESS) return (CEQLibStatus)status;

    *state = trial;
    return CEQ_SUCCESS;
}

void CEQLibCalulateChemcalAbundanceCoolingHeatingRates(StructCEQData *PhydroPP){
    CEQLibStatus status = CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(PhydroPP,1.e-3);
    if(status != CEQ_SUCCESS){
        fprintf(stderr,"CEQLib could not evaluate the chemical equilibrium rates: status %d\n",status);
        return;
    }

    if(CEQLibRunParameters.UseDebugMode==true){
        double Ts = SpinTemperatureHI21cm(PhydroPP->natom[ipHYDROGEN],PhydroPP->T);
        fprintf(stderr," ==>");
        //fprintf(stderr," itr = %d\n",itr);
        fprintf(stderr," ne/nH = %e ",PhydroPP->nel             /PhydroPP->n[ipHYDROGEN]);
        fprintf(stderr," Hii/H = %e ",PhydroPP->nion[ipHYDROGEN]/PhydroPP->n[ipHYDROGEN]);
        fprintf(stderr," H2 /H = %e ",PhydroPP->nmol[ipH2      ]/PhydroPP->n[ipHYDROGEN]);
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," Heiii/He = %e ",PhydroPP->nHeIII         /PhydroPP->n[ipHELIUM]);
        fprintf(stderr," Heii /He = %e ",PhydroPP->nion[ipHELIUM ]/PhydroPP->n[ipHELIUM]);
        fprintf(stderr," Hei  /He = %e ",PhydroPP->natom[ipHELIUM]/PhydroPP->n[ipHELIUM]);
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," Ciii/C = %e ",PhydroPP->nCIII           /PhydroPP->n[ipCARBON  ]);
        fprintf(stderr," Cii /C = %e ",PhydroPP->nion[ipCARBON  ]/PhydroPP->n[ipCARBON  ]);
        fprintf(stderr," Ci  /C = %e ",PhydroPP->natom[ipCARBON ]/PhydroPP->n[ipCARBON  ]);
        fprintf(stderr," CO /C  = %e ",PhydroPP->nmol[ipCO      ]/PhydroPP->n[ipCARBON  ]);
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," Oiii/O = %e ",PhydroPP->nOIII           /PhydroPP->n[ipOXYGEN  ]);
        fprintf(stderr," Oii /O = %e ",PhydroPP->nion[ipOXYGEN  ]/PhydroPP->n[ipOXYGEN  ]);
        fprintf(stderr," Oi  /O = %e ",PhydroPP->natom[ipOXYGEN ]/PhydroPP->n[ipOXYGEN  ]);
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," SiIII/Si = %e ",PhydroPP->nSiIII          /PhydroPP->n[ipSILICON]);
        fprintf(stderr," SiII /Si = %e ",PhydroPP->nion [ipSILICON]/PhydroPP->n[ipSILICON]);
        fprintf(stderr," SiI  /Si = %e ",PhydroPP->natom[ipSILICON]/PhydroPP->n[ipSILICON]);
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," L_total     = %e ",PhydroPP->L_tot    );
        fprintf(stderr," H_total     = %e ",PhydroPP->H_tot    );
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," L(Hi      ) = %e ",PhydroPP->L_HI     );
        fprintf(stderr," L(Hei     ) = %e ",PhydroPP->L_HeI    );
        fprintf(stderr," L(Heii    ) = %e ",PhydroPP->L_HeII   );
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," L([Cii]158) = %e ",PhydroPP->L_CII158 );
        fprintf(stderr," L([Oi ]063) = %e ",PhydroPP->L_OI63   );
        fprintf(stderr," L([SiII]35) = %e ",PhydroPP->L_SiII35 );
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," L(Dust Rec) = %e ",PhydroPP->L_DustRec);
        fprintf(stderr," H(PE      ) = %e ",PhydroPP->H_PE     );
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," H(PI      ) = %e ",PhydroPP->H_PI     );
        fprintf(stderr," H(CR      ) = %e ",PhydroPP->H_CR     );
        fprintf(stderr," H(XR      ) = %e ",PhydroPP->H_XR     );
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," L(H2      ) = %e ",PhydroPP->L_H2     );
        fprintf(stderr," H(H2      ) = %e ",PhydroPP->H_H2     );
        fprintf(stderr," L(CO      ) = %e ",PhydroPP->L_CO     );
        fprintf(stderr,"\n");
        fprintf(stderr," ==>");
        fprintf(stderr," L(HI21cm  ) = %e ",PhydroPP->L_HI21cm);
        fprintf(stderr," k(HI21cm  ) = %e ",PhydroPP->k_HI21cm);
        fprintf(stderr," Ts(HI21cm ) = %e ",Ts);
        fprintf(stderr,"\n");
    }
}

int CEQLibCalculateCoolingHeatingRates(StructCEQData *PhydroPP){

    if(PhydroPP == NULL || !isfinite(PhydroPP->T) || PhydroPP->T <= 0.0) return CEQ_INVALID_INPUT;

    PhydroPP->L_res = CoolingRateCIE(PhydroPP);
    if(!isfinite(PhydroPP->L_res)) return CEQ_INVALID_INPUT;

    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double T    = PhydroPP->T;

    // heating rate [erg/s/cc]
    PhydroPP->H_PI      = HeatingPhotoionization(PhydroPP); //T,nH,ne,Z,G0,Nt);
    PhydroPP->H_PE      = HeatingDustPhotoelectic(PhydroPP); //T,nH,ne,Z,G0,Nt);
    PhydroPP->H_CR      = TotalCosmicRayHeatingRate(PhydroPP);
    PhydroPP->H_XR      = TotalSoftXrayHeatingRate(PhydroPP);
    if(CEQLibRunParameters.UseH2Formation==true){
        PhydroPP->H_H2  = HeatingH2Photodissociate(PhydroPP)
                        + HeatingH2PhotomPumping(PhydroPP)
                        + HeatingH2Formation(PhydroPP);
    }else{
        PhydroPP->H_H2  = 0.e0;
    }

    PhydroPP->H_tot     = PhydroPP->H_PI
                        + PhydroPP->H_PE
                        + PhydroPP->H_CR
                        + PhydroPP->H_XR
                        + PhydroPP->H_H2;

    // cooling rate = angular-integrated emission coefficient (= 4*pi*j) [erg/s/cc]
    PhydroPP->L_ff      = CoolingBremsstrahalung(PhydroPP);
    PhydroPP->L_HI      = CoolingCollisionalExcitation(PhydroPP)
                        + CoolingCollisionalIonization(PhydroPP)
                        + CoolingRecombination(PhydroPP);

    PhydroPP->L_HeI     = CoolingHeICollisionalExcitation(PhydroPP)
                        + CoolingHeICollisionalIonization(PhydroPP)
                        + CoolingHeIIRecombination(PhydroPP)
                        + CoolingHeIIDielectronicRecombination(PhydroPP);

    PhydroPP->L_HeII    = CoolingHeIICollisionalExcitation(PhydroPP)
                        + CoolingHeIICollisionalIonization(PhydroPP)
                        + CoolingHeIIIRecombination(PhydroPP);

    PhydroPP->L_CII158  = CoolingCII158um(PhydroPP);
    PhydroPP->L_CI610   = 0.0; //CoolingCI610um(nCI,T,ne,nHII,nHI,nH2);
    PhydroPP->L_CI370   = 0.0; //CoolingCI370um(nCI,T,ne,nHII,nHI,nH2);
    PhydroPP->L_OIII    = 0.e0;
    PhydroPP->L_OII     = 0.e0;
    double coolingOI44;
    CoolingOI(PhydroPP->natom[ipOXYGEN],T,PhydroPP->nel,PhydroPP->nion[ipHYDROGEN],
              nHI,PhydroPP->nmol[ipH2],&PhydroPP->L_OI63,&coolingOI44,&PhydroPP->L_OI146);
    PhydroPP->L_SiII35  = CoolingSiII35um(PhydroPP); //nSiII,T,ne,nHI);
    PhydroPP->L_DustRec = CoolingDustRecombination(PhydroPP); //T,nH,ne,Z,G0,Nt);
    if(CEQLibRunParameters.UseH2Formation==true){
        PhydroPP->L_H2  = 
                        + CoolingH2RovibrationalLine(PhydroPP)
                        + CoolingH2CollisionalDissociation(PhydroPP)
                        ;
    }else {
        PhydroPP->L_H2  = 0.e0;
    }
    if(CEQLibRunParameters.UseCOFormation==true){
        PhydroPP->L_CO  = CoolingCORovibrationalLine(PhydroPP); //T,nCO,nH,nHI,nH2);
    }else{
        PhydroPP->L_CO  = 0.e0;
    }
    
    PhydroPP->L_tot     = PhydroPP->L_res
                        + PhydroPP->L_ff
                        + PhydroPP->L_HI
                        + PhydroPP->L_HeI
                        + PhydroPP->L_HeII
                        + PhydroPP->L_CII158
                        + PhydroPP->L_CI610
                        + PhydroPP->L_CI370
                        + PhydroPP->L_OIII
                        + PhydroPP->L_OII
                        + PhydroPP->L_OI63
                        + PhydroPP->L_OI146
                        + PhydroPP->L_SiII35
                        + PhydroPP->L_DustRec
                        + PhydroPP->L_H2
                        + PhydroPP->L_CO;      

    PhydroPP->L_Halpha   = CoolingHalphaEmission(PhydroPP);

    // 4pi * line-integrated emissivity [erg/s/cc]
    PhydroPP->L_HI21cm   = CoolingHI21cm(nHI,T);  // 4*pi*j 
    // absorption coefficient [1/cm]
    PhydroPP->k_HI21cm   = AbsorptionCoefficientHI21cm(nHI,T);
  
#if 0 
    double Hphotion  = HeatingHIPhotoionization(T,nHI,Phydro[i]->Fion);
    double Hphot     = HeatingCarbonPhotoinozation(T,nCI,G0,Av,Nc,N2);
    // CIE plasma cooling : Resonance lines (C,O,Fe ...)
    double Lres     = ne*nH*CoolingRateCIE(logT,Metal);
    // HII region line cooling [erg/cc/s]
    double Lce      = 0.e0;  //CoolingLinesHIIregionOxygen(T,ne,nOII,nOIII);
#endif

    if(!isfinite(PhydroPP->H_tot) || !isfinite(PhydroPP->L_tot)
       || !isfinite(PhydroPP->L_Halpha) || !isfinite(PhydroPP->L_HI21cm)
       || !isfinite(PhydroPP->k_HI21cm)) return CEQ_INVALID_INPUT;
    return CEQ_SUCCESS;
}

double CEQLibGetTotalNumberDensity(StructCEQData *PhydroPP){
    double Ntot = PhydroPP->nel 
                + PhydroPP->nHeIII + PhydroPP->nCIII + PhydroPP->nOIII + PhydroPP->nSiIII;
    for(int i=0; i<NumberChemicalElement; i++)
        Ntot += PhydroPP->nion[i] +PhydroPP-> natom[i];
    for(int i=0; i<NumberMoleculeSpeacies; i++)
        Ntot += PhydroPP->nmol[i];
    return Ntot;
}

double CEQLibGetTotalMassDensity(StructCEQData *PhydroPP){
    double Mass = PhydroPP->nHeIII*ChemicalElementMass[ipHELIUM]
                + PhydroPP->nCIII *ChemicalElementMass[ipCARBON]
                + PhydroPP->nOIII *ChemicalElementMass[ipOXYGEN]
                + PhydroPP->nSiIII*ChemicalElementMass[ipSILICON];
    for(int i=0; i<NumberChemicalElement; i++)
        Mass += (PhydroPP->nion[i] + PhydroPP->natom[i])*ChemicalElementMass[i];
    for(int i=0; i<NumberMoleculeSpeacies; i++)
        Mass += PhydroPP->nmol[i]*MolecularElementMass[i];
    return Mass;
}

double CEQLibGetElectronNumberDensity(StructCEQData *PhydroPP){
    double Ntot = 2*(PhydroPP->nHeIII + PhydroPP->nCIII + PhydroPP->nOIII + PhydroPP->nSiIII);
    for(int i=0; i<NumberChemicalElement; i++)
        Ntot += PhydroPP->nion[i];
    return Ntot;
}

double CEQLibGetMeanMolecularWeight(StructCEQData *PhydroPP){
    return CEQLibGetTotalMassDensity(PhydroPP)/(PROTON_MASS_CGS*CEQLibGetTotalNumberDensity(PhydroPP));
}

double CEQLibGetMeanDegreeOfFreedom(StructCEQData *PhydroPP){
    double Ntot    = CEQLibGetTotalNumberDensity(PhydroPP);
    double Freedom = PhydroPP->nHeIII*NumberOfTranslationalDegreesOfFreedom
                   + PhydroPP->nCIII *NumberOfTranslationalDegreesOfFreedom
                   + PhydroPP->nOIII *NumberOfTranslationalDegreesOfFreedom
                   + PhydroPP->nSiIII*NumberOfTranslationalDegreesOfFreedom;
    for(int i=0; i<NumberChemicalElement; i++)
        Freedom += (PhydroPP->nion[i] + PhydroPP->natom[i])*NumberOfTranslationalDegreesOfFreedom;
    for(int i=0; i<NumberMoleculeSpeacies; i++)
        Freedom += PhydroPP->nmol[i]*NumberOfFreedomForMolecules;
    Freedom += PhydroPP->nel*NumberOfTranslationalDegreesOfFreedom;
    return 0.5*Freedom/Ntot;
}

double CEQLibConvertTtoU(StructCEQData *PhydroPP){
    return PhydroPP->DegreeOfFreedom*(BOLTZMANN_CONSTANT_CGS)*PhydroPP->T/(PhydroPP->MeanMolecularWeight*PROTON_MASS_CGS);
}

double CEQLibConvertUtoT(StructCEQData *PhydroPP){
    return PhydroPP->MeanMolecularWeight*PROTON_MASS_CGS*PhydroPP->U/(PhydroPP->DegreeOfFreedom*BOLTZMANN_CONSTANT_CGS);
}

double CEQApproximateDustToGasRatio(double T){
    return ApproximateDustToGasRatio(T);
}

double CEQApproximateDustTemperature(double Nt, double G0){
    return ApproximateDustTemperature(Nt,G0);
}
#if 0
#ifdef LSODA_SUBCYCLE_COOLING_SOLVER
static void InitializeLSODA(int MaximunInteration){

    // an indicator for the type of error control.  
    itol = 2;

    // an index specifying the task to be performed
    itask = 1;
    //itask = 2;

    // an integer flag to specify whether or not any optional inputs are
    // being used on this call. 
    //  iopt=0: no optional inputs are being used. 
    //  iopt=1: one or more optional inputs are being used. 
    iopt = 1;

    // jacobian type indicator.
    jt = 2;
    
    ///// optional inputs.
    // lower and upper half-bandwidths, respectively, of the banded
    // jacobian, excluding the main diagonal. 
	iwork1 = 0;
    iwork2 = 0;

    // flag to generate extra printing at method switchies.
    //  iwork5=0: no extra printing (default).
    //  iwork5=1: print data on each switch.
    iwork5 = 0;
    
    // maximum number of (initially defined) steps allowed during one
    // call to the solver.  the default value is 500.
    iwork6 = MaximunInteration;

    // maximum number of messages printed (per problem) warning that t+h
    // = t on a step (h= step size). this must be positive to result in
    // a non-default value. the default value is 10.
    //iwork7 = 0;
    iwork7 = 10;

    // maximum order to be allowed for the nonstiff (adams) method. the
    // default value is 12.
    //iwork8 = 0;
    iwork8 = 12;

    // maximum order to be allowed for the stiff (bdf) method. the
    // default value is 5.
    //iwork9 = 0;
    iwork9 = 5;

    // 
	rwork1 = 0.e0;

    // step size to be attempted on the first step. the default value is
    // determined by the solver.
    rwork5 = 0.e0;

    // maximum absolute step size allowed. the default value is
    // infinite.
    //rwork6 = 0.e0;
    rwork6 = 1.e+10;

    // maximum absolute step size allowed. the default value is 0.
    rwork7 = 0.e0;

    return ;
}
#endif
#endif
