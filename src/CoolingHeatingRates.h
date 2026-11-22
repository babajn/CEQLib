#ifndef CEQ_COOLING_HEATING_RATES_H
#define CEQ_COOLING_HEATING_RATES_H

#include "EQCDataStructure.h"
#include "TemperatureTable.h"
#include "UVBackground.h"

struct CEQUVBEquilibriumStats;
struct CEQSecondaryIonizationSeed;

/* Thermochemistry entry points for a UVB evaluated once by the caller. */
EQCLibStatus CEQCalculateChemicalAbundanceCoolingHeatingRatesWithUVB(StructEQCData *state, double accuracy,
                                                                    const CEQUVBRates *uvb);
EQCLibStatus CEQCalculateChemicalAbundanceCoolingHeatingRatesWithSeed(
    StructEQCData *state, double accuracy, const CEQUVBRates *uvb, struct CEQUVBEquilibriumStats *stats,
    struct CEQSecondaryIonizationSeed *seed);
EQCLibStatus CEQSolveEquilibriumEnergyWithStats(StructEQCData *state, double accuracy, double *rate,
                                              struct CEQUVBEquilibriumStats *stats);
/* Same operation as the public equilibrium update, with per-call counters. */
EQCLibStatus CEQAdvanceEquilibriumWithStats(StructEQCData *state, double accuracy,
                                           struct CEQUVBEquilibriumStats *stats);
/* Rate evaluation for the network with explicit H, He, C and O ions. */
int CEQCalculateNonEquilibriumCoolingHeatingRates(StructEQCData *state);
int CEQCalculateNonEquilibriumCoolingHeatingRatesWithUVB(StructEQCData *state, const CEQUVBRates *uvb);
/* Evaluate the energy equation without emissivity and opacity diagnostics. */
int CEQCalculateNonEquilibriumEnergyRates(StructEQCData *state, const CEQTemperaturePoint *point,
                                         double uvbHeating);

/* Fixed-temperature reference integration for the coupled H/He and H/He/O networks. */
EQCLibStatus CEQAdvanceFixedTemperatureHydrogenHeliumChemistry(StructEQCData *state, double accuracy);
EQCLibStatus CEQAdvanceFixedTemperatureHHeOChemistry(StructEQCData *state, double accuracy);

#endif
