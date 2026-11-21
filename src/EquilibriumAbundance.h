#include "CEQDataStructure.h"

int CEQSolveIonLadder(double density, double ionization0, double recombination1,
                      double ionization1, double recombination2,
                      double *neutral, double *ion, double *doubleIon);
int CEQSolveHydrogenEquilibrium(StructCEQData *state, double accuracy);

/* H ionization rate per H atom and recombination rate per H+ [1/s],
   without charge transfer. They depend on the state's n_e, not on its ions. */
typedef struct {
    double ionization, recombination;
} CEQHydrogenRates;
CEQHydrogenRates CEQHydrogenIonizationRates(const StructCEQData *state);
int CEQSolveHydrogenEquilibriumWithRates(StructCEQData *state, const CEQHydrogenRates *rates,
                                         double accuracy);

/* Per-species ladder rates [1/s]: ionization of the neutral, recombination of
   the ion, ionization of the ion, and recombination of the doubly charged ion. */
typedef struct {
    double ionization0, recombination1, ionization1, recombination2;
} CEQIonLadderRates;

/* Solve the O ladder with O0/O+ <-> H+/H0 and O++ + H0 charge transfer, and
   rebalance H I/H II among the atomic H already set in state, keeping n(H2).
   hydrogen holds the H rates without charge transfer for the state's n_e.
   The unknown is s = ln(n_HII/n_HI); success means |s - ln(I/R)| <= accuracy
   for the H ionization and recombination rates including charge transfer. */
int CEQSolveHydrogenOxygenChargeTransfer(StructCEQData *state, const CEQHydrogenRates *hydrogen,
                                         const CEQIonLadderRates *oxygen, double accuracy);
int CEQSolveCarbonEquilibrium(StructCEQData *state, double previousCO, double accuracy);

typedef struct {
    double upperFraction;
    double logColumnAtAnchor;
    double logShieldingAtAnchor;
    double slope;
    double lowerSelfShielding;
    double upperSelfShielding;
} CEQCarbonShieldingInterval;

/* Allocate room for CEQCOShieldingPointCount() + 1 intervals. */
int CEQPrepareCarbonShieldingIntervals(const StructCEQData *state,
                                       CEQCarbonShieldingInterval *intervals, int points);
int CEQSolveCarbonWithIntervals(StructCEQData *state, double previousCO, double accuracy,
                               const CEQCarbonShieldingInterval *intervals, int count);

void EquilibriumHydrogenFractionalAbundance(StructCEQData *PhydroPP);
void EquilibriumHeliumFractionalAbundance(StructCEQData *PhydroPP);
void EquilibriumCarbonFractionalAbundance(StructCEQData *PhydroPP);
void EquilibriumOxygenFractionalAbundanceChargeTransfer(StructCEQData *PhydroPP);
void EquilibriumOxygenFractionalAbundance(StructCEQData *PhydroPP);
void EquilibriumSiliconFractionalAbundance(StructCEQData *PhydroPP);

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
