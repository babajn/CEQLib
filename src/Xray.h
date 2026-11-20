#include "CEQDataStructure.h"

void InitializeSoftXrayIonization(int Check);
double SoftXrayIonizationRateHI(double xe, double Ntot);
double SoftXrayIonizationRateHeI(double xe, double Ntot);
double SoftXrayIonizationRateHeII(double xe, double Ntot);
double SoftXrayIonizationRateCI(double xe, double Ntot);
double SoftXrayIonizationRateOI(double xe, double Ntot);
double SoftXrayIonizationRateSiI(double xe, double Ntot);
double SoftXrayIonizationRateSiII(double xe, double Ntot);

double TotalSoftXrayHeatingRate(StructCEQData *PhydroPP);
double TotalSoftXrayHeatingRateFixed(double nH, double X0);
