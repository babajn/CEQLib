#ifndef CEQ_MOLECULE_CHEMISTRY_H
#define CEQ_MOLECULE_CHEMISTRY_H

#include "EQCDataStructure.h"
#include "TemperatureTable.h"

void InitializeHydrogenChemistory(int Check);
void InitializeMoleculeChemistory(int Check);

double DustShieldingFactor(double Ntot, double Zd, double gamma);
double VisualExtinction(double Ntot, double Zd);
double DustShieldingFactorFUV(double Ntot, double Zd);

double H2FormationRate(double T, double Zd, double Tdust);
double H2DopplerWidth(double T);
double H2SelfShieldingFactor(double N2, double b5);
double H2PhotoDissociationRate(double G0, double Nt, double N2, double b5, double Zd);
double H2PhotoDissociationRateWithDustFactor(double G0, double N2, double b5, double dustFactor);
/* Optionally return d ln(self-shielding) / d ln(N2) for the equilibrium Jacobian. */
double CEQH2PhotoRateAndShieldingSlope(double G0, double N2, double b5, double dustFactor,
                                       double *logSlope);
typedef struct {
    double electron, atomicHydrogen, molecularHydrogen;
} H2CollisionalDissociationCoefficients;
H2CollisionalDissociationCoefficients H2DissociationCoefficients(double T);
static inline double H2DissociationFromCoefficients(H2CollisionalDissociationCoefficients coefficients,
                                                  double ne, double nHI, double nH2){
    return ne*coefficients.electron + nHI*coefficients.atomicHydrogen
           + nH2*coefficients.molecularHydrogen;
}
double H2CollisionalDissociationRate(double T, double nH, double nHI, double nH2, double ne);
    
double CoolingH2RovibrationalLineLowDensityLimitGalliPalla1998(double T, double nHI, double nH2);
double CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(double T, double nHII, double nHI, double nH2,
        double nHe, double ne);
double CoolingH2RovibrationalLineHollenbachMcKee1979(double T, double nHI, double nH2);
double CoolingH2RovibrationalLineNeufeldKaufman1993(double T, double nHI, double nH2);
double CoolingH2RovibrationalLineLeBourlot1999(double T, double nHI, double nH2);
double CoolingH2RovibrationalLine(StructEQCData *PhydroPP);
double CEQCoolingH2RovibrationalLine(const StructEQCData *state, double neutralHelium);
void CEQInitializeMolecularCoolingTemperatureTables(void);
double CEQCoolingH2RovibrationalLineWithTables(const StructEQCData *state, double neutralHelium,
                                              const CEQTemperaturePoint *point);

double CoolingH2CollisionalDissociation(StructEQCData *PhydroPP);
double CEQCoolingH2CollisionalDissociationWithCoefficients(const StructEQCData *state,
                                                         H2CollisionalDissociationCoefficients coefficients);

double HeatingH2Photodissociate(StructEQCData *PhydroPP);
double HeatingH2PhotomPumping(StructEQCData *PhydroPP);
double HeatingH2Formation(StructEQCData *PhydroPP);
double CEQHeatingH2(const StructEQCData *state);
double CEQHeatingH2WithTables(const StructEQCData *state, const CEQTemperaturePoint *point);

double COShieldingFactor(double N2, double Nco, double b5);
/* The tabulated log10 column and log10 CO self-shielding factor. */
int CEQCOShieldingPointCount(void);
void CEQCOShieldingPoint(int index, double *logColumn, double *logFactor);
double COPhotoDissociationRate(double G0, double Ntot, double Zd);
#ifdef CHEMISTRY_CO_NELSON_LANGER_1997
double COFormationRate(double nH, double xO, double G0, double Ntot, double Zd);
#endif
double CoolingCORovibrationalLine(StructEQCData *PhydroPP);

#endif
