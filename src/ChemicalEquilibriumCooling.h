#ifndef __CHEMICALEQUILIBRIUMCOOLING_H_INCLUDED__
#define __CHEMICALEQUILIBRIUMCOOLING_H_INCLUDED__

#include <stdbool.h>
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

#define CEQ_DIAGNOSTIC_DIRECTORY "./CEQCooling"

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
    bool UseChargeTransfer;
    double CpForH2;
    double CpForCO;
    double ChemistryAccuracy;
} CEQLibRunParameters;

typedef enum {
    CEQ_SUCCESS = 0,
    CEQ_INVALID_INPUT = -1,
    CEQ_NO_EQUILIBRIUM = -2, /* No admissible root was found in the search intervals. */
    CEQ_NO_CONVERGENCE = -3,
    CEQ_INITIALIZATION_FAILED = -4
} CEQLibStatus;

/* Initialize before evaluating rates. Set H II regions before initialization
   to load spectra from pegaseDirectory, or from ../data/PEGASE2 if it is NULL.
   A check value of 1 writes diagnostic files. Change the run parameters only
   when no evaluations are running. */
CEQLibStatus CEQLibTryInitialization(int check, const char *pegaseDirectory);
void CEQLibInitialization(int Check);
void CEQLibSetRunParamDebugMode(bool Use);
void CEQLibSetRunParamRecombinationCaseB(bool Use);
/* Use NLyn to evaluate photoionization and blend the CIE metal term with
   [O II]/[O III] cooling where H is photoionized. The explicit line term
   models only oxygen and tapers from 1.5e4 to 2.5e4 K to recover CIE.
   When off, NLyn is ignored and C keeps its FUV photoionization. */
void CEQLibSetRunParamHIIRegion(bool Use);
/* With H II regions enabled, solve O0 + H+ <-> O+ + H0 and O++ + H0 -> O+ + H+
   together with the H balance (default on). Turning it off makes H II updates
   about 20-35% faster, but leaves O III too abundant where some H I remains;
   against Cloudy, the median temperature error at solar metallicity grows
   from about 1.4% to 5.3%. */
void CEQLibSetRunParamChargeTransfer(bool Use);
/* Upper bound on the relative charge residual of the chemistry inside the
   implicit energy API, which uses min(value, 0.1*accuracy). The default 1e-8
   is far tighter than the energy accuracy needs. 1e-5 typically makes energy
   updates 10-25% faster and changes T by less than about 2e-5. The net-rate
   API solves at accuracy 1e-8 and so uses at most 1e-9, whatever the value.
   Returns CEQ_INVALID_INPUT unless 0 < value < 1. */
CEQLibStatus CEQLibSetRunParamChemistryAccuracy(double value);
void CEQLibSetRunParamDustAttenuation(bool Use);
void CEQLibSetRunParamGrainSurfaceProcesses(bool Use);
void CEQLibSetRunParamCosmicRayIonization(bool Use);
void CEQLibSetRunParamSoftXrayIonization(bool Use);
/* Extrapolate the CIE metal contribution above the table (10^8.5 K) with
   sqrt(T) scaling up to 10^9 K. Higher temperatures return CEQ_INVALID_INPUT
   from the checked cooling and heating APIs when CIE is enabled. */
void CEQLibSetRunParamCIEResonanceCooling(bool Use);
void CEQLibSetRunParamH2Formation(bool Use, double Cp);
void CEQLibSetRunParamCOFormation(bool Use, double Cp);

int CEQLibCalculateCoolingHeatingRates(StructCEQData *PhydroPP);
/* Fixed-temperature equilibrium. Negative results are CEQLibStatus errors.
   Supply the previous equilibrium in nel and nmol[ipCO] to track a branch.
   A zero CO seed chooses the lowest admissible CO density on a cold start.
   A failed solve leaves the input state intact. */
int CEQLibChemicalReactionEquilibrium(StructCEQData *PhydroPP, const double Accuracy);
/* The checked calls publish outputs and state only when they succeed. */
CEQLibStatus CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(StructCEQData *state, double accuracy);
CEQLibStatus CEQLibTryNetCoolingRateChemcalEquilibrium(StructCEQData *state, double *rate);
CEQLibStatus CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructCEQData *state, double accuracy, double *Du);
void CEQLibCalulateChemcalAbundanceCoolingHeatingRates(StructCEQData *PhydroPP);
/* Net heating per unit mass (positive for heating), in erg/g/s. NaN on failure. */
double CEQLibNetCoolingRateChemcalEquilibrium(StructCEQData *PhydroPP);
/* Updates the thermochemical state for the solved energy and returns dU/dt,
   or NaN on failure. For either energy API, the caller applies U += dt * dU/dt. */
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
