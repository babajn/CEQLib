#ifndef __CHEMICALEQUILIBRIUMCOOLING_H_INCLUDED__
#define __CHEMICALEQUILIBRIUMCOOLING_H_INCLUDED__

#include "ChemicalElements.h"
#include "CEQDataStructure.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

//#include "ChemicalEquilibriumCooling.h"
//#include "EquilibriumAbundance.h"
//#include "CollisionalGasPhaseChemistry.h"
//#include "PhotochemicalGasPhaseChemistry.h"
//#include "GrainSurfaceChemistry.h"
//#include "MoleculeChemistry.h"
//#include "FineStructureLines.h"
//#include "HyperFineStructureLines.h"
//#include "CosmicRay.h"
//#include "Xray.h"
//#include "Bremsstrahlung.h"
//#include "CoolingCIE.h"
//#include "HIIregionLineCooling.h"

static char CEQOutDir[] = "./CEQCooling";
static char CEQFileName[100];

extern struct StructCEQLibParameters_tag{
    bool UseDebugMode;
    bool UseRecombinationCaseB;
    bool UseDustAttenuation;
    bool UseGrainSurfaceProcesses;
    bool UseCosmicRayIonization;
    bool UseSoftXrayIonization;
    bool UseCIEResonanceCooling;
    bool UseH2Formation;
    bool UseCOFormation;
    bool UseHIIRegion;
    double CpForH2;
    double CpForCO;
} CEQLibRunParameters;

void CEQLibInitialization(int Check);
void CEQLibiSetRunParamDebugMode(bool Use);
void CEQLibiSetRunParamRecombinationCaseB(bool Use);
void CEQLibiSetRunParamHIIRegion(bool Use);
void CEQLibiSetRunParamDustAttenuation(bool Use);
void CEQLibiSetRunParamGrainSurfaceProcesses(bool Use);
void CEQLibiSetRunParamCosmicRayIonization(bool Use);
void CEQLibiSetRunParamSoftXrayIonization(bool Use);
void CEQLibiSetRunParamCIEResonanceCooling(bool Use);
void CEQLibiSetRunParamH2Formation(bool Use, double Cp);
void CEQLibiSetRunParamCOFormation(bool Use, double Cp);

int CEQLibCalculateCoolingHeatingRates(StructCEQData *PhydroPP);
int CEQLibChemicalReactionEquilibrium(StructCEQData *PhydroPP, const double Accuracy);
void CEQLibCalulateChemcalAbundanceCoolingHeatingRates(StructCEQData *PhydroPP);
double CEQLibNetCoolingRateChemcalEquilibrium(StructCEQData *PhydroPP);
double CEQLibEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructCEQData *PhydroPP, double Accuracy);

double CEQLibGetTotalNumberDensity(StructCEQData *PhydroPP);
double CEQLibGetTotalMassDensity(StructCEQData *PhydroPP);
double CEQLibGetElectronNumberDensity(StructCEQData *PhydroPP);
double CEQLibGetMeanMolecularWeight(StructCEQData *PhydroPP);
double CEQLibGetMeanDegreeOfFreedom(StructCEQData *PhydroPP);
double CEQLibConvertTtoU(StructCEQData *PhydroPP);
double CEQLibConvertUtoT(StructCEQData *PhydroPP);
double CEQApproximateDustToGasRatio(double T);
double CEQApproximateDustTemperature(double Nt, double G0);

#endif
