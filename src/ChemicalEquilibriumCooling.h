#ifndef __CHEMICALEQUILIBRIUMCOOLING_H_INCLUDED__
#define __CHEMICALEQUILIBRIUMCOOLING_H_INCLUDED__

#include <stdbool.h>
#include "ChemicalElements.h"
#include "EQCDataStructure.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#ifdef __cplusplus
extern "C" {
#endif

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

#define EQC_DIAGNOSTIC_DIRECTORY "./EQCooling"

extern struct StructEQCLibParameters_tag{
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
    bool UseNonEquilibriumTemperatureTables;
    double CpForH2;
    double CpForCO;
    double ChemistryAccuracy;
} EQCLibRunParameters;

typedef enum {
    EQC_SUCCESS = 0,
    EQC_INVALID_INPUT = -1,
    EQC_NO_EQUILIBRIUM = -2, /* No admissible root was found in the search intervals. */
    EQC_NO_CONVERGENCE = -3,
    EQC_INITIALIZATION_FAILED = -4
} EQCLibStatus;

typedef enum {
    EQC_CHEMISTRY_EQUILIBRIUM = 0,
    EQC_CHEMISTRY_NONEQUILIBRIUM
} EQCLibChemistryMode;

/* Parse a trimmed value from a host application's configuration. */
EQCLibStatus EQCLibTryParseChemistryMode(const char *value, EQCLibChemistryMode *mode);
/* Read 'chemistry_mode = equilibrium' or 'chemistry_mode = nonequilibrium'.
   Other settings belong to the host, and # introduces a line comment.
   A NULL path or a file without the key selects equilibrium. I/O failures
   return EQC_INITIALIZATION_FAILED and invalid or duplicate entries return
   EQC_INVALID_INPUT. Both config functions leave mode intact on failure. */
EQCLibStatus EQCLibTryLoadChemistryMode(const char *path, EQCLibChemistryMode *mode);

/* Initialize before evaluating rates. Set H II regions before initialization
   to load spectra from pegaseDirectory, or from ../data/PEGASE2 if it is NULL.
   A check value of 1 writes diagnostic files. Change the run parameters only
   when no evaluations are running. */
EQCLibStatus EQCLibTryInitialization(int check, const char *pegaseDirectory);
void EQCLibInitialization(int Check);
void EQCLibSetRunParamDebugMode(bool Use);
void EQCLibSetRunParamRecombinationCaseB(bool Use);
/* Use NLyn to evaluate photoionization and blend the CIE metal term with
   [O II]/[O III] cooling where H is photoionized. The explicit line term
   models only oxygen and tapers from 1.5e4 to 2.5e4 K to recover CIE.
   When off, NLyn is ignored and C keeps its FUV photoionization. */
void EQCLibSetRunParamHIIRegion(bool Use);
/* Enable O0 + H+ <-> O+ + H0 and O++ + H0 -> O+ + H+ (default on).
   The equilibrium solver uses these reactions when H II regions are enabled.
   Turning them off makes equilibrium H II updates
   about 20-35% faster, but leaves O III too abundant where some H I remains;
   against Cloudy, the median temperature error at solar metallicity grows
   from about 1.4% to 5.3%. The thermochemistry solver uses them whenever
   oxygen is present and this parameter is enabled. */
void EQCLibSetRunParamChargeTransfer(bool Use);
/* Interpolate smooth temperature-only coefficients in the non-equilibrium
   updates (default on). False evaluates the original rate formulas directly.
   Change the setting only while no evaluations are running. */
void EQCLibSetRunParamNonEquilibriumTemperatureTables(bool use);
/* Upper bound on the relative charge residual of the chemistry inside the
   implicit energy API, which uses min(value, 0.1*accuracy). The default 1e-8
   is far tighter than the energy accuracy needs. 1e-5 typically makes energy
   updates 10-25% faster and changes T by less than about 2e-5. The net-rate
   API solves at accuracy 1e-8 and so uses at most 1e-9, whatever the value.
   Returns EQC_INVALID_INPUT unless 0 < value < 1. */
EQCLibStatus EQCLibSetRunParamChemistryAccuracy(double value);
void EQCLibSetRunParamDustAttenuation(bool Use);
void EQCLibSetRunParamGrainSurfaceProcesses(bool Use);
void EQCLibSetRunParamCosmicRayIonization(bool Use);
void EQCLibSetRunParamSoftXrayIonization(bool Use);
/* Enable CIE metal cooling (default off). At positive Z, the C, O and Si
   coefficients are weighted by n_e and each element's total number density.
   N, Ne, Mg, S, Ca and Fe, which EQCLib does not track, add the Asplund et al.
   (2009) solar-mixture coefficient times n_e n_H Z/0.0134 (no depletion); this
   term is zero at or below 1e4 K, full from 1e5 K, and smoothstepped in log T.
   Non-equilibrium cooling mixes it in as [O II] and [O III] fade between
   15000 and 25000 K and uses Steffen interpolation in log temperature.
   Equilibrium uses log-linear interpolation of the same coefficients.
   The CIE rate extrapolates above 10^8 K with sqrt(T) scaling up to 10^9 K.
   Higher temperatures are outside the supported range when using CIE. */
void EQCLibSetRunParamCIEResonanceCooling(bool Use);
/* Cp multiplies the Glover & Jappsen (2007) grain H2 formation coefficient. */
void EQCLibSetRunParamH2Formation(bool Use, double Cp);
void EQCLibSetRunParamCOFormation(bool Use, double Cp);

int EQCLibCalculateCoolingHeatingRates(StructEQCData *PhydroPP);
/* Fixed-temperature equilibrium. Negative results are EQCLibStatus errors.
   Supply the previous equilibrium in nel and nmol[ipCO] to track a branch.
   A zero CO seed chooses the lowest admissible CO density on a cold start.
   A failed solve leaves the input state intact. */
int EQCLibChemicalReactionEquilibrium(StructEQCData *PhydroPP, const double Accuracy);
/* Advance pure hydrogen chemistry by state->dt seconds at fixed T, total
   hydrogen density and radiation fields. Start from nion[ipHYDROGEN] and
   nmol[ipH2] in cm^-3. H I and electron densities are derived from nuclear
   and charge conservation. Set the other elements and CO to zero, and
   enable H2 formation for molecular gas. accuracy controls the integration
   error in the abundances per hydrogen nucleus (0 < accuracy < 1).
   Only the hydrogen species and nel are updated on success. A zero dt is
   a no-op. Initialize the library before calling this function. */
EQCLibStatus EQCLibTryAdvanceHydrogenChemistry(StructEQCData *state, double accuracy);
/* Advance hydrogen, optional helium, atomic carbon and oxygen, and specific thermal energy
   U [erg/g] over state->dt seconds at fixed Rho [g/cm^3] and radiation fields.
   Supply positive U, Rho and hydrogen nuclear density. For He, C and O, give the
   respective totals and singly ionized densities in n and nion, and doubly
   ionized densities in nHeIII, nCIII and nOIII (all in cm^-3). For optional
   silicon, hold nion[ipSILICON] and nSiIII fixed and derive Si I from
   n[ipSILICON]. Its charge and [Si II] 35 um cooling are included.
   Nonzero CO density, or enabled CO formation with a positive timestep,
   returns EQC_INVALID_INPUT.
   Reconstruct T from U and the initial species, and nel from their charge.
   On success, update U, the species, T, the equation of state, and instantaneous
   thermal rates, including [C II] 158 um, [O II] and [O III]. The oxygen
   ion lines taper from 15000 to 25000 K independently of illumination
   when CIE resonance cooling is enabled and Z is positive.
   A zero dt is a no-op. Errors leave the input state intact. */
EQCLibStatus EQCLibTryAdvanceHydrogenThermochemistry(StructEQCData *state, double accuracy);
/* Update U, temperature, composition and thermal rates over state->dt seconds
   with the selected solver. Du is the mean change in specific energy per
   second during the step. The caller takes the updated U directly. Use the
   chosen solver's initial-state requirements and carry its accepted species
   into the next step. A zero dt is a no-op. Errors leave state intact. */
EQCLibStatus EQCLibTryAdvanceCell(StructEQCData *state, double accuracy, EQCLibChemistryMode mode);
/* The checked calls publish outputs and state only when they succeed. */
EQCLibStatus EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(StructEQCData *state, double accuracy);
EQCLibStatus EQCLibTryNetCoolingRateChemcalEquilibrium(StructEQCData *state, double *rate);
EQCLibStatus EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructEQCData *state, double accuracy, double *Du);
void EQCLibCalulateChemcalAbundanceCoolingHeatingRates(StructEQCData *PhydroPP);
/* Net heating per unit mass (positive for heating), in erg/g/s. NaN on failure. */
double EQCLibNetCoolingRateChemcalEquilibrium(StructEQCData *PhydroPP);
/* Updates the thermochemical state for the solved energy and returns dU/dt,
   or NaN on failure. For either energy API, the caller applies U += dt * dU/dt. */
double EQCLibEnergyEquationSolverChemcalEquilibriumCoolingHeating(StructEQCData *PhydroPP, double Accuracy);

double EQCLibGetTotalNumberDensity(StructEQCData *PhydroPP);
double EQCLibGetTotalMassDensity(StructEQCData *PhydroPP);
double EQCLibGetElectronNumberDensity(StructEQCData *PhydroPP);
double EQCLibGetMeanMolecularWeight(StructEQCData *PhydroPP);
double EQCLibGetMeanDegreeOfFreedom(StructEQCData *PhydroPP);
double EQCLibConvertTtoU(StructEQCData *PhydroPP);
double EQCLibConvertUtoT(StructEQCData *PhydroPP);
double EQCApproximateDustToGasRatio(double T);
double EQCApproximateDustTemperature(double Nt, double G0);

#ifdef __cplusplus
}
#endif

#endif
