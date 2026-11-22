#include "ChemicalEquilibriumCooling.h"

void InitializePhotochemistory(int Check);
EQCLibStatus CEQInitializePhotochemistry(int check, const char *pegaseDirectory);

double PhotoionizationCrossSectionHI(double energyEV);
double PhotoionizationCrossSectionHeI(double energyEV);
double PhotoionizationCrossSectionHeII(double energyEV);
double PhotoionizationCrossSectionCI(double energyEV);
double PhotoionizationCrossSectionCII(double energyEV);
double PhotoionizationCrossSectionOI(double energyEV);
double PhotoionizationCrossSectionOII(double energyEV);
double PhotoionizationCrossSectionSiI(double energyEV);
double PhotoionizationCrossSectionSiII(double energyEV);

double HIPhotoIonizationRate(double NLyn);
double HeIPhotoIonizationRate(double NLyn);
double HeIIPhotoIonizationRate(double NLyn);
double CIPhotoionizationRate(double NLyn);
double CIIPhotoionizationRate(double NLyn);
double OIPhotoionizationRate(double NLyn);
double OIIPhotoionizationRate(double NLyn);
double SiIPhotoionizationRate(double NLyn);
double SiIIPhotoionizationRate(double NLyn);

double CIFUVPhotoionizationRate(double G0, double Av);
double SiIFUVPhotoionizationRate(double G0, double Av);

double HeatingPhotoionization(StructEQCData *PhydroPP);
/* The same heating without Si I and Si II. The non-equilibrium update holds
   the silicon ionization fixed, so these photoionizations never occur there. */
double CEQHeatingPhotoionizationWithoutSilicon(const StructEQCData *state);
double HeatingHIPhotoionization(double NLyn);
double HeatingHeIPhotoionization(double NLyn);
double HeatingHeIIPhotoionization(double NLyn);
double HeatingCIPhotoionization(double NLyn);
double HeatingCIIPhotoionization(double NLyn);
double HeatingOIPhotoionization(double NLyn);
double HeatingOIIPhotoionization(double NLyn);
double HeatingSiIPhotoionization(double NLyn);
double HeatingSiIIPhotoionization(double NLyn);
