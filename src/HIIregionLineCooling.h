#ifndef CEQ_HII_REGION_LINE_COOLING_H
#define CEQ_HII_REGION_LINE_COOLING_H

#include "EQCDataStructure.h"

enum { OxygenLevelCount = 5 };

/* Kept for callers that initialized the former runtime tables. */
void InitializeCoolingLinesHIIregion(int Check);

/* Optically thin electron-impact cooling in erg/cm^3/s, or NaN for invalid
   inputs. Levels are in increasing energy order, starting with 0. If supplied,
   lines[upper][lower] receives the power in each transition in the same units.
   Outside the atomic data grid, use the nearest effective collision strength. */
double CoolingLinesOII(double T, double ne, double nOII,
                       double lines[OxygenLevelCount][OxygenLevelCount]);
double CoolingLinesOIII(double T, double ne, double nOIII,
                        double lines[OxygenLevelCount][OxygenLevelCount]);
double CoolingLinesHIIregionOxygen(double T, double ne, double nOII, double nOIII);

double CoolingHalphaEmission(StructEQCData *PhydroPP);

#endif
