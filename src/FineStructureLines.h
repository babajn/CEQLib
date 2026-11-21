#include "EQCDataStructure.h"

void InitializeFineStructureLines(int Check);
//double TotalFineStructureLinesCoolingRate(double T, double nel, double nHII, double nHI, double nH2, 
//        double nCI, double nCII, double nOI, double nSiI, double nSiII);

double CoolingCI(double nCI, double T, double nel, double nHII, double nHI, double nH2, 
        double *CoolingCI610um, double *CoolingCI230um, double *CoolingCI370um);
double CoolingCI610um(double nCI, double T, double nel, double nHII, double nHI, double nH2);
double CoolingCI230um(double nCI, double T, double nel, double nHII, double nHI, double nH2);
double CoolingCI370um(double nCI, double T, double nel, double nHII, double nHI, double nH2);
double CoolingCII158um(StructEQCData *PhydroPP);

double CoolingOI(double nOI, double T, double nel, double nHII, double nHI, double nH2, 
        double *CoolingOI063um, double *CoolingOI044um, double *CoolingOI146um);
double CoolingOI63um(StructEQCData *PhydroPP);
double CoolingOI44um(StructEQCData *PhydroPP);
double CoolingOI146um(StructEQCData *PhydroPP);
//double CoolingOI63um(double nOI, double T, double nel, double nHII, double nHI, double nH2);
//double CoolingOI44um(double nOI, double T, double nel, double nHII, double nHI, double nH2);
//double CoolingOI146um(double nOI, double T, double nel, double nHII, double nHI, double nH2);

/* With Si I disabled, all rates are zero. Enabling the unimplemented model
   returns NaN for the total and each line, including through the wrappers. */
double CoolingSiI(double nSiI, double T, double nHI, double nHII,
        double *CoolingSiI130um, double *CoolingSiI45um, double *CoolingSiI68um);
double CoolingSiI130um(double nSiI, double T, double nHI, double nHII);
double CoolingSiI45um(double nSiI, double T, double nHI, double nHII);
double CoolingSiI68um(double nSiI, double T, double nHI, double nHII);
double CoolingSiII35um(StructEQCData *PhydroPP);
//double CoolingSiII35um(double nSiII, double T, double nel, double nHI);
