#include "EQCDataStructure.h"
#include "TemperatureTable.h"

void InitializeCollisionalIonization(int Check);

double HICollisionalIonizationRate(double T);
double HeICollisionalIonizationRate(double T);
double HeIICollisionalIonizationRate(double T);

double CICollisionalIonizationRate(double T);
double CIICollisionalIonizationRate(double T);
double CIIICollisionalIonizationRate(double T);
double CIVCollisionalIonizationRate(double T);
double CVCollisionalIonizationRate(double T);
double CVICollisionalIonizationRate(double T);

double NICollisionalIonizationRate(double T);
double NIICollisionalIonizationRate(double T);
double NIIICollisionalIonizationRate(double T);
double NIVCollisionalIonizationRate(double T);
double NVCollisionalIonizationRate(double T);
double NVICollisionalIonizationRate(double T);
double NVIICollisionalIonizationRate(double T);

double OICollisionalIonizationRate(double T);
double OIICollisionalIonizationRate(double T);
double OIIICollisionalIonizationRate(double T);
double OIVCollisionalIonizationRate(double T);
double OVCollisionalIonizationRate(double T);
double OVICollisionalIonizationRate(double T);
double OVIICollisionalIonizationRate(double T);
double OVIIICollisionalIonizationRate(double T);

double SiICollisionalIonizationRate(double T);
double SiIICollisionalIonizationRate(double T);
double SiIIICollisionalIonizationRate(double T);
double SiIVCollisionalIonizationRate(double T);


double HIIRecombinationRate(double T);
double HIIRecombinationRateCaseA(double T);
double HIIRecombinationRateCaseB(double T);
double HeIIRecombinationRate(double T);
double HeIIRecombinationRateCaseA(double T);
double HeIIRecombinationRateCaseB(double T);
double HeIIDielectronicRecombinationRate(double T);
double HeIIIRecombinationRate(double T);
double HeIIIRecombinationRateCaseA(double T);
double HeIIIRecombinationRateCaseB(double T);
double CIIRecombinationRate(double T);
double CIIIRecombinationRate(double T);
double OIIRecombinationRate(double T);
double OIIIRecombinationRate(double T);
double SiIIRecombinationRate(double T);
double SiIIIRecombinationRate(double T);

/// ref. Tabel 1 of Glover & Jappsen (2007)
double HeIChargeTransferRateToHII(double T);
double HeIIChargeTransferRateToHI(double T);
double OIChargeTransferRateToHII(double T);
double OIIChargeTransferRateToHI(double T);
double OIIIChargeTransferRateToHI(double T);
double OIChargeTransferRateToHeII(double T);
double CIChargeTransferRateToHII(double T);
double CIIChargeTransferRateToHI(double T);
double CIChargeTransferRateToHeII(double T);
double SiIChargeTransferRateToHII(double T);
double SiIChargeTransferRateToHeII(double T);
double SiIChargeTransferRateToCII(double T);
double SiIIChargeTransferRateToHII(double T);
double SiIIIChargeTransferRateToHI(double T);


double CoolingCollisionalIonization(StructEQCData *PhydroPP);
double CoolingHeICollisionalIonization(StructEQCData *PhydroPP);
double CoolingHeIICollisionalIonization(StructEQCData *PhydroPP);
double CoolingCollisionalExcitation(StructEQCData *PhydroPP);
double CoolingHeICollisionalExcitation(StructEQCData *PhydroPP);
double CoolingHeIICollisionalExcitation(StructEQCData *PhydroPP);
double CoolingRecombination(StructEQCData *PhydroPP);
double CoolingHeIIRecombination(StructEQCData *PhydroPP);
double CoolingHeIIIRecombination(StructEQCData *PhydroPP);
double CoolingHeIIDielectronicRecombination(StructEQCData *PhydroPP);
void CEQInitializeAtomicCoolingTemperatureTables(void);
void CEQAtomicCoolingWithTables(StructEQCData *state, const CEQTemperaturePoint *point);

