#include "EQCDataStructure.h"

void InitializeGrainSurfaceChemistry(int Check);

double HIIRecombinationRateOnDust(double T, double G0, double ne, double Z);
double DIIRecombinationRateOnDust(double T, double G0, double ne, double Z);
double HeIIRecombinationRateOnDust(double T, double G0, double ne, double Z);
double CIIRecombinationRateOnDust(double T, double G0, double ne, double Z);
double OIIRecombinationRateOnDust(double T, double G0, double ne, double Z);
double SiIIRecombinationRateOnDust(double T, double G0, double ne, double Z);

//double TotalDustCoolingRate(double T, double nH, double ne, double Z, double G0);
//double TotalDustCoolingRateHIIregion(double T, double nH, double ne, double Z, double G0);
double CoolingDustRecombination(StructEQCData *PhydroPP);
//double CoolingDustRecombination(double T, double nH, double ne, double Z, double G0, double Ntot);
double HeatingDustPhotoelectic(StructEQCData *PhydroPP);
//double HeatingDustPhotoelectic(double T, double nH, double ne, double Z, double G0, double Ntot);
double CoolingGasDustEnergyTransfer(double T, double Tdust, double nH, double Z);
double HeatingGasDustEnergyTransfer(double T, double Tdust, double nH, double Z);

double ApproximateDustTemperature(double Ntot, double G0);
double EffectiveDustTemperature(const StructEQCData *state);
double ApproximateDustToGasRatio(double T);
