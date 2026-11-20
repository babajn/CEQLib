#include "CEQDataStructure.h"

void InitializePhotochemistory(int Check);

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

double HeatingPhotoionization(StructCEQData *PhydroPP);
double HeatingHIPhotoionization(double NLyn);
double HeatingHeIPhotoionization(double NLyn);
double HeatingHeIIPhotoionization(double NLyn);
double HeatingCIPhotoionization(double NLyn);
double HeatingCIIPhotoionization(double NLyn);
double HeatingOIPhotoionization(double NLyn);
double HeatingOIIPhotoionization(double NLyn);
double HeatingSiIPhotoionization(double NLyn);
double HeatingSiIIPhotoionization(double NLyn);
