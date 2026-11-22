#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "PhotoelectronDeposition.h"
#include "UVBackground.h"

static const double SolarMetallicity = ASPLUND_2009_SOLAR_METALLICITY;
static const int Metals[] = {ipCARBON, ipOXYGEN, ipSILICON};
static const char *const MetalNames[] = {"C", "O", "Si"};

static void Require(bool passed, const char *label)
{
    if(!passed){
        fprintf(stderr,"Metal photoionization verification failed: %s\n",label);
        exit(EXIT_FAILURE);
    }
}

static void RequireSuccess(const char *label, EQCLibStatus status)
{
    if(status != EQC_SUCCESS){
        fprintf(stderr,"%s failed with status %d\n",label,status);
        exit(EXIT_FAILURE);
    }
}

static FILE *OpenOutput(const char *path)
{
    FILE *file = fopen(path,"w");
    if(file == NULL){
        perror(path);
        exit(EXIT_FAILURE);
    }
    return file;
}

static StructEQCData Gas(double nH, double metallicity, bool includeMetals)
{
    StructEQCData state = {.T=1.e4,.Z=metallicity*SolarMetallicity,.Rho=1.4*PROTON_MASS_CGS*nH};
    state.n[ipHYDROGEN] = nH;
    state.n[ipHELIUM] = 0.1*nH;
    if(includeMetals){
        /* Explicit gas-phase Asplund (2009) abundances shared with Cloudy. */
        state.n[ipCARBON] = nH*pow(10.0,8.43-12.0);
        state.n[ipOXYGEN] = nH*pow(10.0,8.69-12.0);
        state.n[ipSILICON] = nH*pow(10.0,7.51-12.0);
    }
    memcpy(state.natom,state.n,sizeof state.natom);
    return state;
}

#ifndef EQC_NO_UVB_METAL_PHOTOIONIZATION
static void CheckHeating(const CEQUVBRates *rates)
{
    StructEQCData gas = Gas(1.0,1.0,true);
    gas.nel = gas.nion[ipHYDROGEN] = 0.01;
    gas.natom[ipHYDROGEN] -= gas.nion[ipHYDROGEN];
    CEQUVBRates withoutStageRates = *rates;
    memset(withoutStageRates.metalIonization,0,sizeof withoutStageRates.metalIonization);
    const CEQPhotoelectronDeposition with = CEQDepositUVBPhotoelectrons(&gas,rates);
    const CEQPhotoelectronDeposition without = CEQDepositUVBPhotoelectrons(&gas,&withoutStageRates);
    Require(with.heating == without.heating && with.otherPower == without.otherPower &&
            with.secondaryHI == without.secondaryHI && with.secondaryHeI == without.secondaryHeI,
            "stage-specific rates do not add photoelectron power");
    const double kinetic = gas.natom[ipHYDROGEN]*rates->heating[CEQ_UVB_HI] +
                           gas.natom[ipHELIUM]*rates->heating[CEQ_UVB_HEI] +
                           rates->metalPowerLow+rates->metalPowerHigh;
    const double accounted = with.heating+with.otherPower +
                             EV_TO_ERG*(13.60*with.secondaryHI+24.59*with.secondaryHeI);
    Require(fabs(kinetic-accounted) <= 32.0*DBL_EPSILON*kinetic,"photoelectron energy closure");
}

static void WriteRate(FILE *output, double z, double Z, double top, double bottom)
{
    const double columns[] = {top,bottom};
    CEQUVBRates result;
    RequireSuccess("metal photoionization slab",CEQEvaluateUVBackgroundWithMetallicity(columns,Z*SolarMetallicity,&result));
    CheckHeating(&result);
    fprintf(output,"%.17g,%.17g,%.17g,%.17g",z,Z,top,bottom);
    for(size_t i=0;i<CEQ_UVB_METAL_SPECIES;i++)
        fprintf(output,",%.17g",result.metalIonization[i]);
    fputc('\n',output);
}

static void WriteRates(void)
{
    FILE *output = OpenOutput("metal-photoionization-rates.csv");
    fputs("z,Z_solar,NHI_top,NHI_bottom,CI,CII,OI,OII,SiI,SiII\n",output);
    const double metallicity[] = {0,1,2};
    for(int z=0;z<=3;z++){
        RequireSuccess("metal photoionization redshift",EQCLibTrySetUVBackgroundRedshift(z));
        WriteRate(output,z,1.0,0.0,0.0);
        for(size_t m=0;m<sizeof metallicity/sizeof metallicity[0];m++)
            for(int exponent=18;exponent<=22;exponent++){
                const double column = pow(10.0,exponent);
                WriteRate(output,z,metallicity[m],column,column);
            }
    }
    const double offGrid[][4] = {
        {.37,.17,4.37e18,1.224e19},
        {.733,.63,7.83e20,1.67e21},
        {2.21,1.54,4.37e21,6.46e21},
        {2.21,1.83,7.84e21,9.12e21},
    };
    for(size_t i=0;i<sizeof offGrid/sizeof offGrid[0];i++){
        const double *row = offGrid[i];
        RequireSuccess("off-grid metal photoionization redshift",EQCLibTrySetUVBackgroundRedshift(row[0]));
        WriteRate(output,row[0],row[1],row[2],row[3]);
    }
    Require(fclose(output) == 0,"closing metal photoionization rates");
}
#endif

static void WriteFractions(void)
{
    FILE *output = OpenOutput("metal-photoionization-eqc-fractions.csv");
    fputs("z,nH,element,neutral,singly,doubly,xe\n",output);
    /* Cloudy includes charge exchange. Exercise the UVB-only EQCLib H–O
       coupling at fixed temperature with the H II region mode off. */
    EQCLibSetRunParamHIIRegion(false);
    EQCLibSetRunParamChargeTransfer(true);
    for(int z=0;z<=2;z+=2){
        RequireSuccess("fixed-temperature redshift",EQCLibTrySetUVBackgroundRedshift(z));
        for(int density=-6;density<=0;density++){
            const double nH = pow(10.0,0.5*density);
            StructEQCData gas = Gas(nH,1.0,true);
            RequireSuccess("metal ion fractions",EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&gas,1.e-8));
            const double secondIon[] = {gas.nCIII,gas.nOIII,gas.nSiIII};
            for(size_t m=0;m<sizeof Metals/sizeof Metals[0];m++){
                const int element = Metals[m];
                fprintf(output,"%d,%.17g,%s,%.17g,%.17g,%.17g,%.17g\n",z,nH,MetalNames[m],
                        gas.natom[element]/gas.n[element],gas.nion[element]/gas.n[element],
                        secondIon[m]/gas.n[element],gas.nel/nH);
            }
        }
    }
    Require(fclose(output) == 0,"closing EQCLib ion fractions");
}

static void WriteRegression(void)
{
    FILE *output = OpenOutput("metal-photoionization-regression.csv");
    fputs("z,nH,Z_solar,mode,T,xe,xHI,xHII,xHeI,xHeII,xHeIII,nH2,nCO,nCI,nCII,nCIII,nOI,nOII,nOIII,"
          "nSiI,nSiII,nSiIII,U,Du,DQheat,MeanMolecularWeight,DegreeOfFreedom,Gamma,Gm1,GGm1,"
          "H_UVB,H_tot,H_PI,H_PE,H_CR,H_XR,H_H2,L_tot,L_res,L_ff,L_HI,L_HeI,L_HeII,L_CII158,L_CI610,L_CI370,"
          "L_OII,L_OIII,L_OI63,L_OI146,L_SiII35,L_DustRec,L_H2,L_CO,L_Halpha,L_HI21cm,k_HI21cm\n",output);
    for(int z=0;z<=2;z+=2){
        RequireSuccess("metal-free regression redshift",EQCLibTrySetUVBackgroundRedshift(z));
        for(int m=0;m<=1;m++){
            for(int exponent=-3;exponent<=1;exponent++){
                const double nH = pow(10.0,exponent);
                for(int mode=0;mode<=1;mode++){
                    StructEQCData gas = Gas(nH,m,false);
                    gas.NHIuvb[0] = nH*10.0*PC_CGS;
                    gas.NHIuvb[1] = nH*50.0*PC_CGS;
                    RequireSuccess("metal-free chemistry",EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&gas,1.e-8));
                    gas.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&gas);
                    gas.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&gas);
                    gas.U = EQCLibConvertTtoU(&gas);
                    gas.dt = 0.001*MEGAYEAR_CGS;
                    const EQCLibChemistryMode chemistry = mode ? EQC_CHEMISTRY_NONEQUILIBRIUM : EQC_CHEMISTRY_EQUILIBRIUM;
                    RequireSuccess("metal-free update",EQCLibTryAdvanceCell(&gas,1.e-5,chemistry));
                    const double nHe = gas.n[ipHELIUM];
                    const double result[] = {gas.T,gas.nel/nH,gas.natom[ipHYDROGEN]/nH,gas.nion[ipHYDROGEN]/nH,
                        gas.natom[ipHELIUM]/nHe,gas.nion[ipHELIUM]/nHe,gas.nHeIII/nHe,gas.nmol[ipH2],gas.nmol[ipCO],
                        gas.natom[ipCARBON],gas.nion[ipCARBON],gas.nCIII,gas.natom[ipOXYGEN],gas.nion[ipOXYGEN],gas.nOIII,
                        gas.natom[ipSILICON],gas.nion[ipSILICON],gas.nSiIII,gas.U,gas.Du,gas.DQheat,
                        gas.MeanMolecularWeight,gas.DegreeOfFreedom,gas.Gamma,gas.Gm1,gas.GGm1,
                        gas.H_UVB,gas.H_tot,gas.H_PI,gas.H_PE,gas.H_CR,gas.H_XR,gas.H_H2,gas.L_tot,gas.L_res,
                        gas.L_ff,gas.L_HI,gas.L_HeI,gas.L_HeII,gas.L_CII158,gas.L_CI610,gas.L_CI370,
                        gas.L_OII,gas.L_OIII,gas.L_OI63,gas.L_OI146,gas.L_SiII35,gas.L_DustRec,gas.L_H2,gas.L_CO,
                        gas.L_Halpha,gas.L_HI21cm,gas.k_HI21cm};
                    fprintf(output,"%d,%a,%d,%d",z,nH,m,mode);
                    for(size_t i=0;i<sizeof result/sizeof result[0];i++)
                        fprintf(output,",%a",result[i]);
                    fputc('\n',output);
                }
            }
        }
    }
    Require(fclose(output) == 0,"closing metal-free regression data");
}

int main(int argc, char **argv)
{
    const bool fractions = argc == 3 && strcmp(argv[2],"--fractions") == 0;
    const bool regression = argc == 3 && strcmp(argv[2],"--regression") == 0;
    if(argc != 2 && !fractions && !regression){
        fprintf(stderr,"Usage: %s SPECTRUM_DIRECTORY [--fractions|--regression]\n",argv[0]);
        return EXIT_FAILURE;
    }
    EQCLibSetRunParamRecombinationCaseB(false);
    RequireSuccess("UVB configuration",EQCLibSetRunParamUVBackground(true,argv[1]));
    RequireSuccess("UVB initialization",EQCLibTryInitialization(0,NULL));
    if(fractions) WriteFractions();
    else if(regression) WriteRegression();
#ifndef EQC_NO_UVB_METAL_PHOTOIONIZATION
    else WriteRates();
#else
    else Require(false,"pre-feature build only supports fraction and regression modes");
#endif
    return EXIT_SUCCESS;
}
