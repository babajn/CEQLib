#ifndef EQCLIB_XRAY_H
#define EQCLIB_XRAY_H

#include <stdbool.h>
#include "EQCDataStructure.h"

#define EQC_XRAY_IONIZATION_MAX_XE 0.5

typedef struct {
    bool ready;
    double base, attenuation;
} CEQSoftXrayColumn;

/* H I and He I use the same rate. A zero-initialized cache may be reused
   while evaluating charges at the same column. Pass NULL for a single call. */
double CEQSoftXrayIonizationRate(double xe, double Ntot, CEQSoftXrayColumn *column);
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

#endif
