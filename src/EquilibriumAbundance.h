#include "CEQDataStructure.h"

int CEQSolveIonLadder(double density, double ionization0, double recombination1,
                      double ionization1, double recombination2,
                      double *neutral, double *ion, double *doubleIon);
int CEQSolveHydrogenEquilibrium(StructCEQData *state, double accuracy);
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
