#include <stdbool.h>
#include "EQCDataStructure.h"

bool InitializeCoolingTableCIE(int Check);
double CoolingRateCIE(StructEQCData *PhydroPP);
double CEQNonEquilibriumCoolingRateCIE(const StructEQCData *state);
