#ifndef __CHEMICAL_ELEMENTS_H_INCLUDED__
#define __CHEMICAL_ELEMENTS_H_INCLUDED__

//#define CHEM_TINY                (1.e-50)
//#define CHEM_TINY                ((double)CHEMICAL_RATE_TINY)
/** these are indices for some elements, on the C scale */
//#define NumberChemicalElement   (6)
enum{
    ipHYDROGEN,
    ipHELIUM,
//    ipLITHIUM,
//    ipBERYLLIUM,
//    ipBORON,
    ipCARBON,
//    ipNITROGEN,
    ipOXYGEN,
//    ipFLUORINE,
//    ipNEON,
//    ipSODIUM,
//    ipMAGNESIUM,
//    ipALUMINIUM,
    ipSILICON,
//    ipPHOSPHORUS,
//    ipSULPHUR,
//    ipCHLORINE,
//    ipARGO,
//    ipPOTASSIUM,
//    ipCALCIUM,
//    ipSCANDIUM,
//    ipTITANIUM,
//    ipVANADIUM,
//    ipCHROMIUM,
//    ipMANGANESE,
//    ipIRON,
//    ipCOBALT,
//    ipNICKEL,
//    ipCOPPER,
//    ipZINC,
//    ipKRYPTON,
    NumberChemicalElement
};

//#define NumberMoleculeSpeacies  (5)
enum{
    ipH2,
    ipCO,
    //ipO2,
    //ipH2O,
    //ipOH,
    NumberMoleculeSpeacies
};

extern double ChemicalElementMass[NumberChemicalElement];
extern char   ChemicalElementName[NumberChemicalElement][100];

extern double MolecularElementMass[NumberMoleculeSpeacies];
extern char   MolecularElementName[NumberMoleculeSpeacies][100];
#endif

void InitializeChemicalElementTable(int Check);

