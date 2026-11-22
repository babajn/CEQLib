#include "EQCDataStructure.h"

#define EQC_XRAY_IONIZATION_MAX_XE 0.5

void InitializeSoftXrayIonization(int Check);
double SoftXrayIonizationRateHI(double xe, double Ntot);
double SoftXrayIonizationRateHeI(double xe, double Ntot);
double SoftXrayIonizationRateHeII(double xe, double Ntot);
double SoftXrayIonizationRateCI(double xe, double Ntot);
double SoftXrayIonizationRateOI(double xe, double Ntot);
double SoftXrayIonizationRateSiI(double xe, double Ntot);
double SoftXrayIonizationRateSiII(double xe, double Ntot);

double TotalSoftXrayHeatingRate(StructEQCData *PhydroPP);
double TotalSoftXrayHeatingRateFixed(double nH, double X0);
