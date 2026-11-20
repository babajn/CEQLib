//#ifdef CHEMISTRY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "config.h"
#include "Astro.h"
#include "Constants.h"
#include "ChemicalElements.h"

double ChemicalElementMass[NumberChemicalElement];
double MolecularElementMass[NumberMoleculeSpeacies];

char   ChemicalElementName[NumberChemicalElement][100];
char   MolecularElementName[NumberMoleculeSpeacies][100];

void InitializeChemicalElementTable(int Check){

    ChemicalElementMass[ipHYDROGEN] = 1.00794;
    ChemicalElementMass[ipHELIUM] = 4.002602; 
//    ChemicalElementMass[ipLITHIUM] = 6.941;
//    ChemicalElementMass[ipBERYLLIUM] = 9.012182;
//    ChemicalElementMass[ipBORON] = 10.811; 
    ChemicalElementMass[ipCARBON] = 12.0107; 
//    ChemicalElementMass[ipNITROGEN] = 14.0067; 
    ChemicalElementMass[ipOXYGEN] = 15.9994; 
//    ChemicalElementMass[ipFLUORINE] = 18.9984032; 
//    ChemicalElementMass[ipNEON] = 20.1797; 
//    ChemicalElementMass[ipSODIUM] = 22.98976928; 
//    ChemicalElementMass[ipMAGNESIUM] = 24.3050; 
//    ChemicalElementMass[ipALUMINIUM] = 26.9815386; 
    ChemicalElementMass[ipSILICON] = 28.0855; 
//    ChemicalElementMass[ipPHOSPHORUS] = 30.973762; 
//    ChemicalElementMass[ipSULPHUR] = 32.065; 
//    ChemicalElementMass[ipCHLORINE] = 35.453; 
//    ChemicalElementMass[ipARGO] = 39.948; 
//    ChemicalElementMass[ipPOTASSIUM] = 39.0983; 
//    ChemicalElementMass[ipCALCIUM] = 40.078; 
//    ChemicalElementMass[ipSCANDIUM] = 44.955912; 
//    ChemicalElementMass[ipTITANIUM] = 47.867; 
//    ChemicalElementMass[ipVANADIUM] = 50.9415; 
//    ChemicalElementMass[ipCHROMIUM] = 51.9961; 
//    ChemicalElementMass[ipMANGANESE] = 54.938045; 
//    ChemicalElementMass[ipIRON] = 55.845; 
//    ChemicalElementMass[ipCOBALT] = 58.933195; 
//    ChemicalElementMass[ipNICKEL] = 58.6934; 
//    ChemicalElementMass[ipCOPPER] = 63.546; 
//    ChemicalElementMass[ipZINC] = 65.409; 
//    ChemicalElementMass[ipKRYPTON] = 83.798; 
    
    sprintf(ChemicalElementName[ipHYDROGEN],"H");
    sprintf(ChemicalElementName[ipHELIUM]  ,"He"); 
    sprintf(ChemicalElementName[ipCARBON]  ,"C"); 
    sprintf(ChemicalElementName[ipOXYGEN]  ,"O");
//    sprintf(ChemicalElementName[ipSILICON] ,"Si"); 
//    sprintf(ChemicalElementName[ipIRON]    ,"Fe"); 

    for(int i=0; i<NumberChemicalElement; i++){
        ChemicalElementMass[i] *= PROTON_MASS_CGS;  // [g] 
        if(Check==1)
            fprintf(stderr,"i=%d: Element: %s Mass=%e[g]\n",
                    i,ChemicalElementName[i],ChemicalElementMass[i]);
    }

    MolecularElementMass[ipH2] = 2.0*ChemicalElementMass[ipHYDROGEN];
    MolecularElementMass[ipCO] = ChemicalElementMass[ipCARBON]+ChemicalElementMass[ipOXYGEN];
    //MolecularElementMass[ipO2] = 2.0*ChemicalElementMass[ipOXYGEN];
    //MolecularElementMass[ipH2O]= 2.0*ChemicalElementMass[ipHYDROGEN]+ChemicalElementMass[ipOXYGEN];
    //MolecularElementMass[ipOH] = ChemicalElementMass[ipHYDROGEN]+ChemicalElementMass[ipOXYGEN];
    sprintf(MolecularElementName[ipH2] ,"H2");
    sprintf(MolecularElementName[ipCO] ,"CO"); 
    //sprintf(MolecularElementName[ipO2] ,"O2");
    //sprintf(MolecularElementName[ipH2O],"H2O"); 
    //sprintf(MolecularElementName[ipOH] ,"OH"); 

    if(Check==1){
        for(int i=0; i<NumberMoleculeSpeacies; i++){
            fprintf(stderr,"i=%d: Molecules: %s Mass=%e[g]\n",
                    i, MolecularElementName[i],MolecularElementMass[i]);
        }
    }

    return ;
}
//#endif  // CHEMISTRY

