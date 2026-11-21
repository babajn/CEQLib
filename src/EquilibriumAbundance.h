#include "EQCDataStructure.h"
#include "MoleculeChemistry.h"
#include "TemperatureTable.h"

int CEQSolveIonLadder(double density, double ionization0, double recombination1,
                      double ionization1, double recombination2,
                      double *neutral, double *ion, double *doubleIon);
int CEQSolveHydrogenEquilibrium(StructEQCData *state, double accuracy);

/* Ladder coefficients that do not depend on n_e: photoionization rates [1/s],
   electron-impact ionization and recombination coefficients [cm^3/s]. */
typedef struct {
    double photo0, collision0, recombination1, photo1, collision1, recombination2;
} CEQLadderCoefficients;

/* H2 balance terms that do not depend on n_e or the H fractions. */
typedef struct {
    double formation;       /* n_H times the grain-surface formation rate [1/s] */
    double doppler;         /* Doppler parameter for self-shielding [km/s] */
    double dustShielding;   /* Dust shielding of the dissociating radiation */
    H2CollisionalDissociationCoefficients collision;
} CEQMolecularHydrogenCoefficients;

typedef struct {
    double formationPerH2, photoCI, photoCO;
} CEQMolecularCarbonCoefficients;

/* Coefficients fixed by T, NLyn, G0, Nt and Zd. A fixed-temperature equilibrium
   evaluates them once and reuses them for every trial electron density. */
typedef struct {
    CEQLadderCoefficients hydrogen, helium, carbon, oxygen, silicon;
    double toOxygenII, toHydrogenI, oxygenIIIToHydrogenI;   /* O-H charge transfer [cm^3/s] */
    bool hasMolecularHydrogen;                              /* Set when H2 formation is on. */
    CEQMolecularHydrogenCoefficients molecularHydrogen;
    CEQMolecularCarbonCoefficients molecularCarbon;
} CEQRateCoefficients;
void CEQPrepareRateCoefficients(const StructEQCData *state, CEQRateCoefficients *coefficients);
/* During a non-equilibrium update the radiation, columns and elemental totals
   are fixed. Initialize their coefficients once, then refresh at each stage T. */
void CEQInitializeChemistryTemperatureTables(void);
void CEQInitializeNonEquilibriumRateCoefficients(const StructEQCData *state, CEQRateCoefficients *coefficients);
void CEQRefreshNonEquilibriumRateCoefficients(const StructEQCData *state, CEQRateCoefficients *coefficients,
                                              const CEQTemperaturePoint *point);
H2CollisionalDissociationCoefficients CEQH2DissociationCoefficientsWithTemperatureTable(
        double T, const CEQTemperaturePoint *point);

/* H ionization rate per H atom and recombination rate per H+ [1/s],
   without charge transfer. They depend on the state's n_e, not on its ions. */
typedef struct {
    double ionization, recombination;
} CEQHydrogenRates;
CEQHydrogenRates CEQHydrogenIonizationRates(const StructEQCData *state);
CEQHydrogenRates CEQHydrogenIonizationRatesWithCoefficients(const StructEQCData *state,
                                                            const CEQRateCoefficients *coefficients);
/* coefficients may be NULL; the H2 terms are then evaluated here. */
int CEQSolveHydrogenEquilibriumWithRates(StructEQCData *state, const CEQHydrogenRates *rates,
                                         const CEQRateCoefficients *coefficients, double accuracy);

/* Per-species ladder rates [1/s]: ionization of the neutral, recombination of
   the ion, ionization of the ion, and recombination of the doubly charged ion. */
typedef struct {
    double ionization0, recombination1, ionization1, recombination2;
} CEQIonLadderRates;
CEQIonLadderRates CEQHeliumIonizationRates(const StructEQCData *state,
                                         const CEQRateCoefficients *coefficients);
CEQIonLadderRates CEQCarbonIonizationRates(const StructEQCData *state,
                                         const CEQRateCoefficients *coefficients);
CEQIonLadderRates CEQOxygenIonizationRates(const StructEQCData *state,
                                         const CEQRateCoefficients *coefficients);

/* Solve the O ladder with O0/O+ <-> H+/H0 and O++ + H0 charge transfer, and
   rebalance H I/H II. With H2 formation on, every trial H ionization also
   re-solves the H2 balance, so H2, H and O all balance in the result.
   hydrogen holds the H rates without charge transfer for the state's n_e.
   The unknown is s = ln(n_HII/n_HI); success means |s - ln(I/R)| <= accuracy
   for the H ionization and recombination rates including charge transfer. */
int CEQSolveHydrogenOxygenChargeTransfer(StructEQCData *state, const CEQHydrogenRates *hydrogen,
                                         const CEQIonLadderRates *oxygen,
                                         const CEQRateCoefficients *coefficients, double accuracy);
int CEQSolveCarbonEquilibrium(StructEQCData *state, double previousCO, double accuracy);

typedef struct {
    double upperFraction;
    double logColumnAtAnchor;
    double logShieldingAtAnchor;
    double slope;
    double lowerSelfShielding;
    double upperSelfShielding;
} CEQCarbonShieldingInterval;

/* Allocate room for CEQCOShieldingPointCount() + 1 intervals. */
int CEQPrepareCarbonShieldingIntervals(const StructEQCData *state,
                                       CEQCarbonShieldingInterval *intervals, int points);
/* With checkExistingCO, keep the previously solved CO if it balances at the
   current H2 density, and update the carbon outside CO as needed.
   coefficients may be NULL to compute the rates here. */
int CEQSolveCarbonWithIntervals(StructEQCData *state, double previousCO, double accuracy,
                               const CEQCarbonShieldingInterval *intervals, int count,
                               const CEQRateCoefficients *coefficients, bool checkExistingCO);

void EquilibriumHydrogenFractionalAbundance(StructEQCData *PhydroPP);
void EquilibriumHeliumFractionalAbundance(StructEQCData *PhydroPP);
void EquilibriumCarbonFractionalAbundance(StructEQCData *PhydroPP);
void EquilibriumOxygenFractionalAbundanceChargeTransfer(StructEQCData *PhydroPP);
void EquilibriumOxygenFractionalAbundance(StructEQCData *PhydroPP);
void EquilibriumSiliconFractionalAbundance(StructEQCData *PhydroPP);

//void EquilibriumHydrogenFractionalAbundance(double T, double nH, double ne, double G0, double Ntot, double Z, double X0, 
//                                            double *fHII, double *fHI, double *fH2);
//void EquilibriumHeliumFractionalAbundance(double T, double nH, double ne, double G0, double Ntot, double Z, double *fHeII, double *fHeI);
//void EquilibriumHeliumFractionalAbundanceHeIII(double T, double nH, double ne, double G0, double Ntot, double Z, double X0, 
//                                               double *fHeIII, double *fHeII, double *fHeI);
//void EquilibriumOxygenFractionalAbundance(double T, double nHI, double nHII, double *fOII, double *fOI);
//void EquilibriumOxygenFractionalAbundanceOIII(double T, double nHI, double nHII, double nH, double ne, 
//                                              double G0, double Ntot, double Z, double X0,
//                                              double *fOIII, double *fOII, double *fOI);
//void EquilibriumCarbonFractionalAbundance(double xC, double nH, double ne, double nH2, double T, double G0, double X0, double Ntot, double N2, 
//        double Z, double *fCII, double *fCI, double *fCO);
//void EquilibriumSiliconFractionalAbundance(double T, double nH, double ne, double G0, double Ntot, double Z, double X0, 
//                                               double *fSiIII, double *fSiII, double *fSiI);
