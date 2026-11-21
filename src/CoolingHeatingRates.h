#ifndef CEQ_COOLING_HEATING_RATES_H
#define CEQ_COOLING_HEATING_RATES_H

#include "EQCDataStructure.h"
#include "TemperatureTable.h"

/* Rate evaluation for the network with explicit H, He, C and O ions. */
int CEQCalculateNonEquilibriumCoolingHeatingRates(StructEQCData *state);
/* Evaluate the energy equation without emissivity and opacity diagnostics. */
int CEQCalculateNonEquilibriumEnergyRates(StructEQCData *state, const CEQTemperaturePoint *point);

#endif
