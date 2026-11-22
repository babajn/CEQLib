#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"
#include "EquilibriumAbundance.h"

typedef struct {
    double density, temperature;
    double top, bottom;
} TestCell;

static const double Accuracy = 1.0e-9;
static const char *const Stages[] = {"HI","HII","HeI","HeII","HeIII","OI","OII","OIII","xe"};
enum { StageCount = sizeof Stages/sizeof Stages[0] };
static const TestCell Cases[] = {{0.001,10000.0,0.0,0.0}, {0.0316,8000.0,1.e18,2.e18}, {0.316,5500.0,5.e19,7.e19},
                                 {1.0,150.0,3.e20,3.e20}, {10.0,50.0,3.e21,3.e21}};

static void RequireSuccess(const char *operation, EQCLibStatus status)
{
    if(status != EQC_SUCCESS){
        fprintf(stderr,"%s failed with status %d\n",operation,status);
        exit(EXIT_FAILURE);
    }
}

static StructEQCData Gas(TestCell test)
{
    const double nH = test.density;
    StructEQCData gas = {.T=test.temperature,.dt=1.0e3*MEGAYEAR_CGS,.Z=ASPLUND_2009_SOLAR_METALLICITY,
                         .Zd=0.01,.G0=0.1,.X0=1.0,.Nt=nH*100.0*PC_CGS,
                         .NHIuvb={test.top,test.bottom}};
    gas.n[ipHYDROGEN] = nH;
    gas.n[ipHELIUM] = 0.1*nH;
    gas.n[ipOXYGEN] = nH*pow(10.0,8.69-12.0);
    gas.nion[ipHYDROGEN] = 0.2*gas.n[ipHYDROGEN];
    gas.nion[ipHELIUM] = 0.1*gas.n[ipHELIUM];
    gas.nion[ipOXYGEN] = 0.1*gas.n[ipOXYGEN];
    gas.nOIII = 0.2*gas.n[ipOXYGEN];
    for(int element=0;element<NumberChemicalElement;element++)
        gas.natom[element] = gas.n[element]-gas.nion[element];
    gas.natom[ipOXYGEN] -= gas.nOIII;
    gas.nel = EQCLibGetElectronNumberDensity(&gas);
    gas.Rho = EQCLibGetTotalMassDensity(&gas);
    return gas;
}

static void Fractions(const StructEQCData *gas, double fractions[StageCount])
{
    const double nH = gas->n[ipHYDROGEN], nHe = gas->n[ipHELIUM], nO = gas->n[ipOXYGEN];
    fractions[0] = gas->natom[ipHYDROGEN]/nH;
    fractions[1] = gas->nion[ipHYDROGEN]/nH;
    fractions[2] = gas->natom[ipHELIUM]/nHe;
    fractions[3] = gas->nion[ipHELIUM]/nHe;
    fractions[4] = gas->nHeIII/nHe;
    fractions[5] = nO > 0.0 ? gas->natom[ipOXYGEN]/nO : 0.0;
    fractions[6] = nO > 0.0 ? gas->nion[ipOXYGEN]/nO : 0.0;
    fractions[7] = nO > 0.0 ? gas->nOIII/nO : 0.0;
    fractions[8] = gas->nel/nH;
}

#ifndef EQC_UVB_CT_REGRESSION_ONLY
static const double FractionFloor = 1.0e-8;

static double ScaledDifference(double actual, double reference)
{
    return fabs(actual-reference)/fmax(reference,FractionFloor);
}

/* The signed term dropped from d n(H II)/dt by the one-way approximation. */
static double RelativeHydrogenTransfer(const StructEQCData *gas)
{
    const double nHI = gas->natom[ipHYDROGEN], nHII = gas->nion[ipHYDROGEN];
    const double oxygenToHydrogen = nHI*(gas->nion[ipOXYGEN]*OIIChargeTransferRateToHI(gas->T) +
                                         gas->nOIII*OIIIChargeTransferRateToHI(gas->T));
    const double hydrogenToOxygen = nHII*gas->natom[ipOXYGEN]*OIChargeTransferRateToHII(gas->T);
    const CEQHydrogenRates h = CEQHydrogenIonizationRates(gas);
    const double budget = fmax(nHI*h.ionization,nHII*h.recombination);
    return (oxygenToHydrogen-hydrogenToOxygen)/budget;
}

static bool CheckCell(FILE *output, int redshift, TestCell test)
{
    const StructEQCData initial = Gas(test);
    StructEQCData oneWay = initial, twoWay = initial, bdf = initial;
    EQCLibSetRunParamHIIRegion(false);
    RequireSuccess("UVB-only fixed-temperature equilibrium",EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&oneWay,Accuracy));
    EQCLibSetRunParamHIIRegion(true);
    RequireSuccess("bidirectional fixed-temperature equilibrium",EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&twoWay,Accuracy));
    EQCLibSetRunParamHIIRegion(false);
    RequireSuccess("fixed-temperature BDF",CEQAdvanceFixedTemperatureHHeOChemistry(&bdf,Accuracy));
    StructEQCData bdfExtended = bdf;
    RequireSuccess("extended fixed-temperature BDF",CEQAdvanceFixedTemperatureHHeOChemistry(&bdfExtended,Accuracy));

    double one[StageCount], two[StageCount], integrated[StageCount], extended[StageCount];
    Fractions(&oneWay,one);
    Fractions(&twoWay,two);
    Fractions(&bdf,integrated);
    Fractions(&bdfExtended,extended);
    const double netTransfer = RelativeHydrogenTransfer(&oneWay);
    const double chargeError = fabs(oneWay.nel-EQCLibGetElectronNumberDensity(&oneWay))/test.density;
    bool success = bdf.T == initial.T && bdfExtended.T == initial.T && isfinite(netTransfer) && chargeError < 1.0e-8;
    double maxOneWayError = 0.0, maxTwoWayError = 0.0;
    for(size_t i=0;i<StageCount;i++){
        const double oneWayError = ScaledDifference(one[i],extended[i]);
        const double twoWayError = ScaledDifference(two[i],extended[i]);
        const double timeError = ScaledDifference(integrated[i],extended[i]);
        maxOneWayError = fmax(maxOneWayError,oneWayError);
        maxTwoWayError = fmax(maxTwoWayError,twoWayError);
        fprintf(output,"%d,%.17g,%.17g,%.17g,%.17g,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                redshift,test.density,test.temperature,test.top,test.bottom,Stages[i],one[i],two[i],integrated[i],extended[i],
                oneWayError,twoWayError,timeError,netTransfer,chargeError);
        if(!isfinite(one[i]) || !isfinite(two[i]) || !isfinite(integrated[i]) || !isfinite(extended[i]) ||
           oneWayError > 0.01 || twoWayError > 1.0e-4 || timeError > 1.0e-4){
            fprintf(stderr,"Mismatch at z=%d, nH=%g, T=%g for %s: one-way %.3g, two-way %.3g, extended BDF %.3g\n",
                    redshift,test.density,test.temperature,Stages[i],oneWayError,twoWayError,timeError);
            success = false;
        }
    }
    printf("Fixed-temperature H–O z=%d nH=%g T=%g K: one-way/BDF %.3g, two-way/BDF %.3g, "
           "omitted H flux %.3g, charge residual %.3g per H\n",redshift,test.density,test.temperature,maxOneWayError,
           maxTwoWayError,netTransfer,chargeError);
    return success;
}

static int RunBdfComparison(void)
{
    FILE *output = fopen("uvb-charge-transfer-bdf.csv","w");
    if(output == NULL){
        perror("uvb-charge-transfer-bdf.csv");
        return EXIT_FAILURE;
    }
    fputs("z,nH,T,NHI_top,NHI_bottom,stage,one_way,two_way,bdf_1Gyr,bdf_2Gyr,"
          "one_way_scaled_error,two_way_scaled_error,bdf_elapsed_scaled_error,omitted_H_flux_ratio,charge_error_per_H\n",output);
    bool success = true;
    for(int z=0;z<=2;z+=2){
        RequireSuccess("FG20 redshift",EQCLibTrySetUVBackgroundRedshift(z));
        for(size_t i=0;i<sizeof Cases/sizeof Cases[0];i++){
            if(!CheckCell(output,z,Cases[i])) success = false;
        }
    }
    if(fclose(output) != 0){
        perror("uvb-charge-transfer-bdf.csv");
        return EXIT_FAILURE;
    }
    printf("Scaled fractional errors use a %.0e floor for near-absent ion stages.\n",FractionFloor);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
#endif

typedef struct {
    const char *name;
    bool oxygen, hiiMode, chargeTransfer;
    double stellarFlux;
    EQCLibChemistryMode mode;
} RegressionCase;

static void WriteRegressionCase(FILE *output, int z, TestCell condition, RegressionCase profile)
{
    EQCLibSetRunParamHIIRegion(profile.hiiMode);
    EQCLibSetRunParamChargeTransfer(profile.chargeTransfer);
    StructEQCData gas = Gas(condition);
    if(!profile.oxygen)
        gas.n[ipOXYGEN] = gas.natom[ipOXYGEN] = gas.nion[ipOXYGEN] = gas.nOIII = 0.0;
    gas.n[ipCARBON] = gas.natom[ipCARBON] = 1.4e-4*condition.density;
    gas.n[ipSILICON] = gas.natom[ipSILICON] = 1.7e-6*condition.density;
    gas.NLyn = profile.stellarFlux;
    gas.dt = 0.001*MEGAYEAR_CGS;
    gas.nel = EQCLibGetElectronNumberDensity(&gas);
    gas.Rho = EQCLibGetTotalMassDensity(&gas);
    if(profile.mode == EQC_CHEMISTRY_EQUILIBRIUM)
        RequireSuccess(profile.name,EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&gas,Accuracy));
    gas.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&gas);
    gas.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&gas);
    gas.U = EQCLibConvertTtoU(&gas);
    RequireSuccess(profile.name,EQCLibTryAdvanceCell(&gas,1.0e-5,profile.mode));
    double fractions[StageCount];
    Fractions(&gas,fractions);
    const double outputs[] = {
        gas.T,gas.U,gas.Du,gas.DQheat,gas.MeanMolecularWeight,gas.DegreeOfFreedom,gas.Gamma,gas.Gm1,gas.GGm1,
        gas.natom[ipCARBON],gas.nion[ipCARBON],gas.nCIII,gas.natom[ipSILICON],gas.nion[ipSILICON],gas.nSiIII,
        gas.nmol[ipH2],gas.nmol[ipCO],gas.H_tot,gas.H_UVB,gas.H_PI,gas.H_PE,gas.H_CR,gas.H_XR,gas.H_H2,
        gas.L_tot,gas.L_res,gas.L_ff,gas.L_HI,gas.L_HeI,gas.L_HeII,gas.L_CII158,gas.L_CI610,gas.L_CI370,
        gas.L_OII,gas.L_OIII,gas.L_OI63,gas.L_OI146,gas.L_SiII35,gas.L_DustRec,gas.L_H2,gas.L_CO,
        gas.L_Halpha,gas.L_HI21cm,gas.k_HI21cm,
    };
    fprintf(output,"%s,%d,%a,%a",profile.name,z,condition.density,condition.temperature);
    for(size_t i=0;i<StageCount;i++) fprintf(output,",%a",fractions[i]);
    for(size_t i=0;i<sizeof outputs/sizeof outputs[0];i++) fprintf(output,",%a",outputs[i]);
    fputc('\n',output);
}

static int RunRegression(void)
{
    FILE *output = fopen("uvb-charge-transfer-regression.csv","w");
    if(output == NULL){
        perror("uvb-charge-transfer-regression.csv");
        return EXIT_FAILURE;
    }
    fputs("profile,z,nH,T_initial",output);
    for(size_t i=0;i<StageCount;i++) fprintf(output,",%s",Stages[i]);
    fputs(",T,U,Du,DQheat,MeanMolecularWeight,DegreeOfFreedom,Gamma,Gm1,GGm1,"
          "n_CI,n_CII,n_CIII,n_SiI,n_SiII,n_SiIII,n_H2,n_CO,H_tot,H_UVB,H_PI,H_PE,H_CR,H_XR,H_H2,"
          "L_tot,L_res,L_ff,L_HI,L_HeI,L_HeII,L_CII158,L_CI610,L_CI370,L_OII,L_OIII,L_OI63,L_OI146,"
          "L_SiII35,L_DustRec,L_H2,L_CO,L_Halpha,L_HI21cm,k_HI21cm\n",output);
    const RegressionCase profiles[] = {
        {"no-oxygen-solar-Z",false,false,true,0.0,EQC_CHEMISTRY_EQUILIBRIUM},
        {"charge-transfer-off",true,false,false,0.0,EQC_CHEMISTRY_EQUILIBRIUM},
        {"HII-zero-stellar-flux",true,true,true,0.0,EQC_CHEMISTRY_EQUILIBRIUM},
        {"HII-with-stellar-flux",true,true,true,1.0e5,EQC_CHEMISTRY_EQUILIBRIUM},
        {"HII-charge-transfer-off",true,true,false,1.0e5,EQC_CHEMISTRY_EQUILIBRIUM},
        {"non-equilibrium",true,false,true,0.0,EQC_CHEMISTRY_NONEQUILIBRIUM},
    };
    for(int z=0;z<=2;z+=2){
        RequireSuccess("FG20 regression redshift",EQCLibTrySetUVBackgroundRedshift(z));
        for(size_t p=0;p<sizeof profiles/sizeof profiles[0];p++)
            for(size_t i=0;i<sizeof Cases/sizeof Cases[0];i++)
                WriteRegressionCase(output,z,Cases[i],profiles[p]);
    }
    if(fclose(output) != 0){
        perror("uvb-charge-transfer-regression.csv");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    const bool regression = (argc == 3 || argc == 4) && strcmp(argv[2],"--regression") == 0;
    if(argc != 2 && !regression){
        fprintf(stderr,"Usage: %s FG20_DIRECTORY [--regression [PEGASE_DIRECTORY]]\n",argv[0]);
        return EXIT_FAILURE;
    }
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamChargeTransfer(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamNonEquilibriumTemperatureTables(false);
    /* Load PEGASE spectra for the regression profiles with stellar flux. */
    EQCLibSetRunParamHIIRegion(regression);
    RequireSuccess("UVB configuration",EQCLibSetRunParamUVBackground(true,argv[1]));
    const char *pegase = argc == 4 ? argv[3] : getenv("PEGASE_DIR");
    RequireSuccess("UVB initialization",EQCLibTryInitialization(0,pegase && *pegase ? pegase : NULL));
    if(regression) return RunRegression();
#ifndef EQC_UVB_CT_REGRESSION_ONLY
    return RunBdfComparison();
#else
    fprintf(stderr,"This executable was built for pre-change regression only.\n");
    return EXIT_FAILURE;
#endif
}
