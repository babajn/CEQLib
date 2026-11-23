#include "EQCDataStructure.h"

void InitializeCosmicRayIonization(int Check);
double CosmicRayIonizationRateHI(double xe);
double CosmicRayIonizationRateDI(double xe);
double CosmicRayIonizationRateH2(double xe);
double CEQH2CosmicRayDissociationRate(double X0);
double CosmicRayIonizationRateHeI(double xe);
double CosmicRayIonizationRateCI(double xe);
double CosmicRayIonizationRateOI(double xe);
double CosmicRayIonizationRateSiI(double xe);
double CosmicRayIonizationRateSiII(double xe);
double CosmicRayIonizationRateCO(double xe);
//double CosmicRayHeating(double nH);
double TotalCosmicRayHeatingRate(StructEQCData *PhydroPP);
double TotalCosmicRayHeatingRateFixed(double nH, double X0);
