#define _POSIX_C_SOURCE 200809L

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "CosmicRay.h"
#include "GrainSurfaceChemistry.h"
#include "MetalXray.h"
#include "MoleculeChemistry.h"
#include "PhotoelectronDeposition.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "UVBackground.h"
#include "Xray.h"

static const double IonizedHydrogen = 0.01;
static const double SolarMetallicity = ASPLUND_2009_SOLAR_METALLICITY;

static void Require(bool passed, const char *label)
{
    if(!passed){
        fprintf(stderr,"Metal X-ray verification failed: %s\n",label);
        exit(EXIT_FAILURE);
    }
}

static FILE *OpenOutput(const char *path)
{
    FILE *output = fopen(path,"w");
    if(!output){
        perror(path);
        exit(EXIT_FAILURE);
    }
    return output;
}

static double CPUTime(void)
{
    struct timespec time;
    Require(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time) == 0,"process CPU clock");
    return (double)time.tv_sec + 1.0e-9*(double)time.tv_nsec;
}

static void WriteCrossSections(void)
{
    FILE *file = OpenOutput("metal-cross-sections.csv");
    fputs("energy_keV,sigma_HI,sigma_He_weighted,sigma_metal,sigma_total\n",file);
    const double energies[] = {0.3,0.5,1.0,2.0,5.0};
    for(size_t i=0;i<sizeof energies/sizeof energies[0];i++){
        const double energy = 1000.0*energies[i];
        const double hydrogen = PhotoionizationCrossSectionHI(energy);
        const double helium = 0.1*PhotoionizationCrossSectionHeI(energy);
        double metals, low, high;
        CEQMetalXrayMixture(energy,&metals,&low,&high);
        fprintf(file,"%.17g,%.17g,%.17g,%.17g,%.17g\n",energies[i],hydrogen,helium,metals,hydrogen+helium+metals);
    }
    Require(fclose(file) == 0,"closing cross sections");
}

static void CheckLegacyAndValidation(void)
{
    const double columns[] = {1.0e19,2.0e22};
    CEQUVBRates legacy, metallic;
    Require(CEQEvaluateUVBackground(columns,&legacy) == EQC_SUCCESS,"legacy evaluator");
    Require(CEQEvaluateUVBackgroundWithMetallicity(columns,0.0,&metallic) == EQC_SUCCESS,"metallicity zero");
    Require(memcmp(legacy.ionization,metallic.ionization,sizeof legacy.ionization) == 0 &&
            memcmp(legacy.heating,metallic.heating,sizeof legacy.heating) == 0 &&
            memcmp(legacy.heatingHigh,metallic.heatingHigh,sizeof legacy.heatingHigh) == 0 &&
            metallic.metalPowerLow == 0.0 && metallic.metalPowerHigh == 0.0,"zero metallicity rate bits");
    Require(CEQEvaluateUVBackgroundWithMetallicity(columns,NAN,&metallic) == EQC_INVALID_INPUT &&
            CEQEvaluateUVBackgroundWithMetallicity(columns,-SolarMetallicity,&metallic) == EQC_INVALID_INPUT &&
            CEQEvaluateUVBackgroundWithMetallicity(columns,INFINITY,&metallic) == EQC_INVALID_INPUT,
            "invalid metallicity");
    Require(memcmp(legacy.ionization,metallic.ionization,sizeof legacy.ionization) == 0 &&
            memcmp(legacy.heating,metallic.heating,sizeof legacy.heating) == 0,
            "failed lookup leaves rates intact");

    double gamma[3], heating[3], gammaWithMetals[3], heatingWithMetals[3], metalPower[2];
    Require(EQCLibTryGetUVBackgroundHemisphereRates(columns[0],gamma,heating) == EQC_SUCCESS,
            "legacy public hemisphere");
    Require(EQCLibTryGetUVBackgroundHemisphereRatesWithMetallicity(columns[0],0.0,gammaWithMetals,
                                                                    heatingWithMetals,metalPower) == EQC_SUCCESS,
            "metallicity-zero public hemisphere");
    Require(memcmp(gamma,gammaWithMetals,sizeof gamma) == 0 &&
            memcmp(heating,heatingWithMetals,sizeof heating) == 0 && metalPower[0] == 0.0 && metalPower[1] == 0.0,
            "public hemisphere zero-metallicity rate bits");
    Require(EQCLibTryGetUVBackgroundHemisphereRatesWithMetallicity(columns[0],SolarMetallicity,gammaWithMetals,
                                                                    heatingWithMetals,metalPower) == EQC_SUCCESS,
            "solar-metallicity public hemisphere");
    Require(gammaWithMetals[0] < gamma[0] && metalPower[0]+metalPower[1] > 0.0,
            "public hemisphere metal shielding and absorption");
}

static void WriteRate(FILE *file, double z, double metallicity, double top, double bottom)
{
    const double columns[] = {top,bottom};
    CEQUVBRates rates;
    Require(CEQEvaluateUVBackgroundWithMetallicity(columns,metallicity*SolarMetallicity,&rates) == EQC_SUCCESS,
            "metallicity-dependent slab");
    StructEQCData gas = {.T=100.0,.nel=IonizedHydrogen,.n={1.0,0.1},.natom={1.0-IonizedHydrogen,0.1},
                         .nion={IonizedHydrogen},.Z=metallicity*SolarMetallicity};
    const CEQPhotoelectronDeposition result = CEQDepositUVBPhotoelectrons(&gas,&rates);
    const double expected = gas.natom[ipHYDROGEN]*rates.heating[CEQ_UVB_HI] +
                            gas.natom[ipHELIUM]*rates.heating[CEQ_UVB_HEI] +
                            rates.metalPowerLow+rates.metalPowerHigh;
    const double accounted = result.heating+result.otherPower +
                             EV_TO_ERG*(13.60*result.secondaryHI+24.59*result.secondaryHeI);
    const double closure = expected > 0.0 ? fabs(accounted/expected-1.0) : fabs(accounted);
    Require(closure <= 32.0*DBL_EPSILON && result.otherPower >= 0.0,"energy budget");

    fprintf(file,"%.17g,%.17g,%.17g,%.17g",z,metallicity,top,bottom);
    for(size_t species=0;species<CEQ_UVB_SPECIES;species++)
        fprintf(file,",%.17g,%.17g,%.17g",rates.ionization[species],
                rates.heating[species]-rates.heatingHigh[species],rates.heatingHigh[species]);
    fprintf(file,",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",rates.metalPowerLow,rates.metalPowerHigh,
            result.heating,rates.ionization[CEQ_UVB_HI]+result.ionizationHI,result.otherPower,closure);
}

static void WriteOrdinaryChemistry(FILE *file, double z, double column, double temperature)
{
    const double Nt = column/(1.0-IonizedHydrogen);
    const double nH = Nt/(100.0*PC_CGS);
    const double nHI = (1.0-IonizedHydrogen)*nH;
    const double ne = IonizedHydrogen*nH;
    const double G0 = 1.0, Zd = 0.01;
    const double Av = VisualExtinction(Nt,Zd);
    /* Ratios with the FUV and charge-transfer coefficients give a scale
       even for channels that are not active in the disk solver setup. */
    const double rates[][7] = {
        {CIFUVPhotoionizationRate(G0,Av), ne*CIChargeTransferRateToHII(temperature),
         nHI*CIIChargeTransferRateToHI(temperature), ne*CIIRecombinationRate(temperature),
         nH*CIIRecombinationRateOnDust(temperature,G0,ne,Zd),
         CosmicRayIonizationRateCI(IonizedHydrogen)+ne*CICollisionalIonizationRate(temperature),
         SoftXrayIonizationRateCI(IonizedHydrogen,Nt)},
        {0.0,ne*OIChargeTransferRateToHII(temperature),nHI*OIIChargeTransferRateToHI(temperature),
         ne*OIIRecombinationRate(temperature),nH*OIIRecombinationRateOnDust(temperature,G0,ne,Zd),
         CosmicRayIonizationRateOI(IonizedHydrogen)+ne*OICollisionalIonizationRate(temperature),0.0},
        {SiIFUVPhotoionizationRate(G0,Av),ne*SiIChargeTransferRateToHII(temperature),0.0,
         ne*SiIIRecombinationRate(temperature),nH*SiIIRecombinationRateOnDust(temperature,G0,ne,Zd),
         CosmicRayIonizationRateSiI(IonizedHydrogen)+ne*SiICollisionalIonizationRate(temperature),
         SoftXrayIonizationRateSiI(IonizedHydrogen,Nt)},
    };
    const char *const symbols[] = {"C","O","Si"};
    for(size_t i=0;i<sizeof symbols/sizeof symbols[0];i++)
        fprintf(file,"%.17g,%.17g,%.17g,%.17g,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                z,column,temperature,nH,symbols[i],rates[i][0],rates[i][1],rates[i][2],
                rates[i][3],rates[i][4],rates[i][5],rates[i][6]);
}

int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr,"Usage: %s FG20_DIRECTORY\n",argv[0]);
        return EXIT_FAILURE;
    }
    Require(EQCLibSetRunParamUVBackground(true,argv[1]) == EQC_SUCCESS,"setting FG20 directory");
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    FILE *timing = OpenOutput("metal-timing.csv");
    fputs("event,seconds,persistent_bytes\n",timing);
    double start = CPUTime();
    Require(EQCLibTryInitialization(0,NULL) == EQC_SUCCESS,"initialization");
    fprintf(timing,"initialization,%.9g,%zu\n",CPUTime()-start,CEQUVBPersistentBytes());
    start = CPUTime();
    Require(EQCLibTrySetUVBackgroundRedshift(0.005) == EQC_SUCCESS,"redshift for timing");
    fprintf(timing,"first_set,%.9g,%zu\n",CPUTime()-start,CEQUVBPersistentBytes());
    start = CPUTime();
    for(int i=0;i<10000;i++)
        Require(EQCLibTrySetUVBackgroundRedshift(i%2 ? 0.005 : 0.015) == EQC_SUCCESS,"in-bracket redshift");
    fprintf(timing,"same_bracket_mean,%.9g,%zu\n",(CPUTime()-start)/10000.0,CEQUVBPersistentBytes());
    start = CPUTime();
    Require(EQCLibTrySetUVBackgroundRedshift(2.0) == EQC_SUCCESS,"new redshift bracket");
    fprintf(timing,"changed_bracket,%.9g,%zu\n",CPUTime()-start,CEQUVBPersistentBytes());
    Require(fclose(timing) == 0,"closing timings");

    CheckLegacyAndValidation();
    WriteCrossSections();
    FILE *rates = OpenOutput("metal-rates.csv");
    fputs("z,Z_solar,NHI_top,NHI_bottom,Gamma_HI,power_low_HI,power_high_HI,"
          "Gamma_HeI,power_low_HeI,power_high_HeI,Gamma_HeII,power_low_HeII,power_high_HeII,"
          "metal_power_low,metal_power_high,H_UVB,Gamma_HI_total,Q_other,energy_closure\n",rates);
    FILE *chemistry = OpenOutput("metal-ordinary-chemistry.csv");
    fputs("z,NHI_per_side,T,nH,element,FUV,charge_transfer_ionization,charge_transfer_recombination,"
          "electron_recombination,dust_recombination,cosmic_and_collisional_ionization,wolfire_ionization\n",
          chemistry);
    const double metallicities[] = {0.0,0.1,1.0,2.0};
    for(int z=0;z<=2;z+=2){
        Require(EQCLibTrySetUVBackgroundRedshift(z) == EQC_SUCCESS,"setting grid redshift");
        for(int exponent=19;exponent<=22;exponent++){
            const double column = pow(10.0,exponent);
            for(size_t m=0;m<sizeof metallicities/sizeof metallicities[0];m++)
                WriteRate(rates,z,metallicities[m],column,column);
            WriteOrdinaryChemistry(chemistry,z,column,100.0);
            WriteOrdinaryChemistry(chemistry,z,column,8000.0);
        }
    }
    const double offGrid[][4] = {
        {0.37,0.17,4.37e19,1.224e20},
        {0.733,0.63,7.83e20,1.67e21},
        {2.21,1.54,4.37e21,6.46e21},
        {2.21,1.83,7.84e21,9.12e21},
    };
    for(size_t i=0;i<sizeof offGrid/sizeof offGrid[0];i++){
        const double *row = offGrid[i];
        Require(EQCLibTrySetUVBackgroundRedshift(row[0]) == EQC_SUCCESS,"setting off-grid redshift");
        WriteRate(rates,row[0],row[1],row[2],row[3]);
    }
    Require(fclose(rates) == 0 && fclose(chemistry) == 0,"closing rate data");
    return EXIT_SUCCESS;
}
