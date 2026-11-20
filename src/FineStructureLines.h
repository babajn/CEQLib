#include "CEQDataStructure.h"

void InitializeFineStructureLines(int Check);
//double TotalFineStructureLinesCoolingRate(double T, double nel, double nHII, double nHI, double nH2, 
//        double nCI, double nCII, double nOI, double nSiI, double nSiII);

double CoolingCI(double nCI, double T, double nel, double nHII, double nHI, double nH2, 
        double *CoolingCI610um, double *CoolingCI230um, double *CoolingCI370um);
double CoolingCI610um(double nCI, double T, double nel, double nHII, double nHI, double nH2);
double CoolingCI230um(double nCI, double T, double nel, double nHII, double nHI, double nH2);
double CoolingCI370um(double nCI, double T, double nel, double nHII, double nHI, double nH2);
double CoolingCII158um(StructCEQData *PhydroPP);

double CoolingOI(double nOI, double T, double nel, double nHII, double nHI, double nH2, 
        double *CoolingOI063um, double *CoolingOI044um, double *CoolingOI146um);
double CoolingOI63um(StructCEQData *PhydroPP);
double CoolingOI44um(StructCEQData *PhydroPP);
double CoolingOI146um(StructCEQData *PhydroPP);
//double CoolingOI63um(double nOI, double T, double nel, double nHII, double nHI, double nH2);
//double CoolingOI44um(double nOI, double T, double nel, double nHII, double nHI, double nH2);
//double CoolingOI146um(double nOI, double T, double nel, double nHII, double nHI, double nH2);

double CoolingSiI(double nSiI, double T, double nHI, double nHII,
        double *CoolingSiI130um, double *CoolingSiI45um, double *CoolingSiI68um);
double CoolingSiI130um(double nSiI, double T, double nHI, double nHII);
double CoolingSiI45um(double nSiI, double T, double nHI, double nHII);
double CoolingSiI68um(double nSiI, double T, double nHI, double nHII);
double CoolingSiII35um(StructCEQData *PhydroPP);
//double CoolingSiII35um(double nSiII, double T, double nel, double nHI);
