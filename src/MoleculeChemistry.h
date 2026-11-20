#include "CEQDataStructure.h"

void InitializeHydrogenChemistory(int Check);
void InitializeMoleculeChemistory(int Check);

double DustShieldingFactor(double Ntot, double Zd, double gamma);
double VisualExtinction(double Ntot, double Zd);
double DustShieldingFactorFUV(double Ntot, double Zd);

double H2FormationRate(double T, double Zd, double Tdust);
double H2SelfShieldingFactor(double N2, double b5);
double H2PhotoDissociationRate(double G0, double Nt, double N2, double b5, double Zd);
double H2CollisionalDissociationRate(double T, double nH, double nHI, double nH2, double ne);
    
double CoolingH2RovibrationalLineLowDensityLimitGalliPalla1998(double T, double nHI, double nH2);
double CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(double T, double nHII, double nHI, double nH2,
        double nHe, double ne);
double CoolingH2RovibrationalLineHollenbachMcKee1979(double T, double nHI, double nH2);
double CoolingH2RovibrationalLineNeufeldKaufman1993(double T, double nHI, double nH2);
double CoolingH2RovibrationalLineLeBourlot1999(double T, double nHI, double nH2);
double CoolingH2RovibrationalLine(StructCEQData *PhydroPP);

double CoolingH2CollisionalDissociation(StructCEQData *PhydroPP);

double HeatingH2Photodissociate(StructCEQData *PhydroPP);
double HeatingH2PhotomPumping(StructCEQData *PhydroPP);
double HeatingH2Formation(StructCEQData *PhydroPP);

double COShieldingFactor(double N2, double Nco, double b5);
double COPhotoDissociationRate(double G0, double Ntot, double Zd);
#ifdef CHEMISTRY_CO_NELSON_LANGER_1997
double COFormationRate(double nH, double xO, double G0, double Ntot, double Zd);
#endif
double CoolingCORovibrationalLine(StructCEQData *PhydroPP);

