#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "EquilibriumAbundance.h"
#include "PhotoelectronDeposition.h"

enum { HI, HeI, HeII, SpeciesCount };
static const char *const SpeciesNames[] = {"HI", "HeI", "HeII"};
static const double Redshifts[] = {0.0, 1.0, 2.0, 3.0};
static const double Temperature = 1.0e4;

typedef CEQUVBRates RadiationRates;

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

static void Close(FILE *file)
{
    if(fclose(file) != 0){
        perror("UVB CSV");
        exit(EXIT_FAILURE);
    }
}

static RadiationRates InSlab(double upper, double lower)
{
    const double columns[] = {upper,lower};
    RadiationRates rates;
    RequireSuccess("slab radiation",CEQEvaluateUVBackground(columns,&rates));
    return rates;
}

static StructEQCData Gas(double nH)
{
    StructEQCData state = {.T = Temperature, .n = {nH,0.1*nH}, .natom = {nH,0.1*nH}};
    return state;
}

static void WriteRates(FILE *file, double redshift)
{
    const double columns[][2] = {{0,0}, {1.0e18,1.0e18}, {1.0e19,1.0e20},
                                  {0,1.0e21}, {4.7e17,8.31e19}, {1.0e21,1.0e21}, {1.0e24,1.0e24}};
    for(size_t column = 0; column < sizeof columns/sizeof columns[0]; column++){
        const double top = columns[column][0], bottom = columns[column][1];
        const RadiationRates rates = InSlab(top,bottom);
        for(int i = 0; i < SpeciesCount; i++)
            fprintf(file,"%g,%.9e,%.9e,%s,%.16e,%.16e,%.16e,%.16e\n",redshift,top,bottom,SpeciesNames[i],
                    rates.ionization[i],rates.heating[i],rates.heating[i]-rates.heatingHigh[i],rates.heatingHigh[i]);
    }
}

typedef struct {
    double electronDensity;
    double hydrogenII;
    double heliumII;
    double heliumIII;
} IonFractions;

static IonFractions AtomicLadder(double nH, double nHe, double trialElectrons, const RadiationRates *radiation,
                                 IonFractions previous)
{
    /* At fixed electron density, the analytic adjacent-stage equilibrium
       populations follow from balancing each pair of ionization/recombination rates. */
    const double nHI = nH*(1.0-previous.hydrogenII);
    const double nHeI = nHe*(1.0-previous.heliumII-previous.heliumIII);
    const double high = nHI*radiation->heatingHigh[HI] + nHeI*radiation->heatingHigh[HeI] +
                        nHe*previous.heliumII*radiation->heatingHigh[HeII];
    const double x = fmax(previous.hydrogenII,1.0e-4);
    const double extraH = nHI > 0.0 ? 0.3908*pow(1.0-pow(x,0.4092),1.7592)*high/(13.60*EV_TO_ERG*nHI) : 0.0;
    const double neutralHe = fmax(nHeI,nHe*(1.0-previous.hydrogenII));
    const double extraHe = neutralHe > 0.0 ? 0.0554*pow(1.0-pow(x,0.4614),1.666)*high/(24.59*EV_TO_ERG*neutralHe) : 0.0;
    const double ionizeH = radiation->ionization[HI]+extraH+trialElectrons*HICollisionalIonizationRate(Temperature);
    const double recombineH = trialElectrons*HIIRecombinationRateCaseA(Temperature);
    const double heliumRatio1 = (radiation->ionization[HeI]+extraHe+trialElectrons*HeICollisionalIonizationRate(Temperature))/
        (trialElectrons*(HeIIRecombinationRateCaseA(Temperature)+HeIIDielectronicRecombinationRate(Temperature)));
    const double heliumRatio2 = (radiation->ionization[HeII]+trialElectrons*HeIICollisionalIonizationRate(Temperature))/
        (trialElectrons*HeIIIRecombinationRateCaseA(Temperature));
    const double heliumNormalization = 1.0+heliumRatio1*(1.0+heliumRatio2);
    IonFractions result = {
        .hydrogenII = ionizeH/(ionizeH+recombineH),
        .heliumII = heliumRatio1/heliumNormalization,
        .heliumIII = heliumRatio1*heliumRatio2/heliumNormalization,
    };
    result.electronDensity = nH*result.hydrogenII+nHe*(result.heliumII+2.0*result.heliumIII);
    return result;
}

static IonFractions SolveAtomicLadder(double nH, double nHe, double ne, const RadiationRates *rates)
{
    IonFractions composition = {0};
    for(int i=0;i<128;i++){
        const IonFractions next = AtomicLadder(nH,nHe,ne,rates,composition);
        if(fabs(next.hydrogenII-composition.hydrogenII) < 1.0e-14 &&
           fabs(next.heliumII-composition.heliumII) < 1.0e-14 &&
           fabs(next.heliumIII-composition.heliumIII) < 1.0e-14) return next;
        composition.hydrogenII = 0.5*(composition.hydrogenII+next.hydrogenII);
        composition.heliumII = 0.5*(composition.heliumII+next.heliumII);
        composition.heliumIII = 0.5*(composition.heliumIII+next.heliumIII);
    }
    fprintf(stderr,"Reference atomic H/He ladder did not converge at nH=%g, ne=%g\n",nH,ne);
    exit(EXIT_FAILURE);
}

static IonFractions AnalyticalEquilibrium(double nH, double nHe, RadiationRates rates)
{
    double lower = 0.0, upper = nH+2.0*nHe;
    for(int i = 0; i < 64; i++){
        const double trial = 0.5*(lower+upper);
        if(SolveAtomicLadder(nH,nHe,trial,&rates).electronDensity > trial) lower = trial;
        else upper = trial;
    }
    return SolveAtomicLadder(nH,nHe,0.5*(lower+upper),&rates);
}

static void WriteEquilibrium(FILE *file, double redshift, CEQUVBEquilibriumStats *stats)
{
    const RadiationRates thin = InSlab(0.0,0.0);
    for(int exponent = -4; exponent <= 4; exponent++){
        const double nH = pow(10.0,exponent);
        StructEQCData state = Gas(nH);
        const double nHe = state.n[ipHELIUM];
        const IonFractions expected = AnalyticalEquilibrium(nH,nHe,thin);
        StructEQCData profile = state;
        if(CEQChemicalReactionEquilibriumWithUVBStats(&profile,1.0e-8,&thin,stats) < 0)
            RequireSuccess("profiled H/He equilibrium",EQC_NO_CONVERGENCE);
        RequireSuccess("H/He equilibrium",EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8));
        fprintf(file,"%g,%.9e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",redshift,nH,
                state.nel/nH,expected.electronDensity/nH,state.nion[ipHYDROGEN]/nH,expected.hydrogenII,
                state.nion[ipHELIUM]/nHe,expected.heliumII,state.nHeIII/nHe,expected.heliumIII);
    }
}

static void WriteDeepColumns(FILE *output, double z, double maximumEnergy)
{
    RequireSuccess("redshift selection",EQCLibTrySetUVBackgroundRedshift(z));
    for(int exponent=19;exponent<=21;exponent++){
        const double column = pow(10.0,exponent);
        const RadiationRates rates = InSlab(column,column);
        for(int ionized=-3;ionized<=-1;ionized++){
            const double x = pow(10.0,ionized);
            StructEQCData gas = Gas(1.0);
            gas.nion[ipHYDROGEN] = gas.nel = x;
            gas.natom[ipHYDROGEN] = 1.0-x;
            gas.NHIuvb[0] = gas.NHIuvb[1] = column;
            const double total =  (1.0-x)*rates.heating[HI]+gas.natom[ipHELIUM]*rates.heating[HeI];
            const double high = CEQHighPhotoelectronPower(&gas,&rates);
            const CEQPhotoelectronDeposition deposit = CEQDepositUVBPhotoelectrons(&gas,&rates);
            RequireSuccess("fixed composition cooling and heating",EQCLibCalculateCoolingHeatingRates(&gas));
            if(deposit.heating != gas.H_UVB){
                fprintf(stderr,"UVB heating disagrees with photoelectron budget\n");
                exit(EXIT_FAILURE);
            }
            fprintf(output,"%g,%.9e,%g,%g,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                    z,column,x,maximumEnergy,total,high,gas.H_UVB,deposit.otherPower,
                    rates.ionization[HI]+deposit.ionizationHI,rates.ionization[HeI]+deposit.ionizationHeI,
                    deposit.ionizationHI,deposit.ionizationHeI);
        }
    }
}

static void WriteXrayComparison(FILE *file, double redshift)
{
    /* Compare heating on the same disk cell, with externally specified
       linear columns, electron fraction and X0. The chemistry is held fixed. */
    EQCLibSetRunParamSoftXrayIonization(true);
    for(int exponent = 19; exponent <= 21; exponent++){
        const double neutralColumn = pow(10.0,exponent);
        StructEQCData state = Gas(1.0);
        state.nion[ipHYDROGEN] = state.nel = 0.01;
        state.natom[ipHYDROGEN] -= state.nion[ipHYDROGEN];
        state.X0 = 1.0;
        state.Nt = neutralColumn/0.99;
        state.NHIuvb[0] = state.NHIuvb[1] = neutralColumn;
        RequireSuccess("UVB/X-ray heating",EQCLibCalculateCoolingHeatingRates(&state));
        fprintf(file,"%g,%.9e,%.9e,%.9e,%.12e,%.12e,%.12e\n",redshift,neutralColumn,state.Nt,state.nel,
                state.H_UVB,state.H_XR,state.H_UVB/state.H_XR);
    }
    EQCLibSetRunParamSoftXrayIonization(false);
}

typedef struct {
    double redshift;
    double nH;
    double metallicity;
    double isrf;
    double ionizationMultiplier;
    double maximumEnergy;
    bool uvb;
    bool softXray;
    bool includeMolecularColumn;
    bool hiiRegion;
    bool chargeTransfer;
    EQCLibChemistryMode mode;
} DiskCase;

typedef struct {
    int completed;
    EQCLibStatus status;
    double temperature;
    double uvbHeatingFraction;
} DiskResult;

enum { DiskSteps = 12 };
static const double DiskChemistryAccuracy = 1.e-5;

static const DiskCase DiskDefaults = {
    .metallicity = 0.0134, .isrf = 0.1, .ionizationMultiplier = 1.0,
    .uvb = true, .softXray = true, .includeMolecularColumn = true,
    .chargeTransfer = true, .mode = EQC_CHEMISTRY_EQUILIBRIUM
};

static StructEQCData DiskInitialState(const DiskCase *test)
{
    const double nH = test->nH;
    StructEQCData state = Gas(nH);
    state.T = 8000.0;
    state.nel = 0.5*nH;
    state.nion[ipHYDROGEN] = state.nel;
    state.n[ipCARBON] = 1.4e-4*nH;
    state.n[ipOXYGEN] = 3.2e-4*nH;
    state.n[ipSILICON] = 1.7e-6*nH;
    for(int i=0;i<NumberChemicalElement;i++) state.natom[i] = state.n[i]-state.nion[i];
    state.Rho = 1.4*PROTON_MASS_CGS*nH;
    state.Z = test->metallicity;
    state.Zd = 0.01;
    state.G0 = test->isrf;
    state.X0 = test->ionizationMultiplier;
    state.Nt = nH*100.0*PC_CGS;
    state.dt = MEGAYEAR_CGS;
    return state;
}

static void WriteChargeTransferBudget(FILE *file, const DiskCase *test, StructEQCData *state, int step)
{
    CEQUVBRates uvb;
    RequireSuccess("disk charge transfer diagnostic radiation",
                   CEQEvaluateUVBackgroundWithMetallicity(state->NHIuvb,state->Z,&uvb));
    const CEQHydrogenRates hydrogen = CEQHydrogenIonizationRates(state);
    CEQRateCoefficients coefficients;
    CEQPrepareRateCoefficients(state,&coefficients,&uvb);
    const CEQIonLadderRates oxygen = CEQOxygenIonizationRates(state,&coefficients);

    const double nHI = state->natom[ipHYDROGEN], nHII = state->nion[ipHYDROGEN];
    const double nOI = state->natom[ipOXYGEN], nOII = state->nion[ipOXYGEN], nOIII = state->nOIII;
    const double ionization = nHI*hydrogen.ionization;
    const double recombination = nHII*hydrogen.recombination;
    /* Evaluate the physical flux even when the selected solver ignores it. */
    const double transferIonization = nHI*(nOII*coefficients.toHydrogenI+nOIII*coefficients.oxygenIIIToHydrogenI);
    const double transferRecombination = nHII*nOI*coefficients.toOxygenII;
    const double transfer = transferIonization-transferRecombination;
    const double oxygenIonization = nOI*oxygen.ionization0+nOII*oxygen.ionization1;
    const double oxygenRecombination = nOII*oxygen.recombination1+nOIII*oxygen.recombination2;
    const double netOxygenIonization = oxygenIonization-oxygenRecombination;
    const double hydrogenResidual = ionization-recombination;
    const double chargeFromSpecies = EQCLibGetElectronNumberDensity(state);
    const double chargeError = state->nel-chargeFromSpecies;
    if(!isfinite(chargeError) || fabs(chargeError) > DiskChemistryAccuracy*fmax(state->nel,chargeFromSpecies)){
        fprintf(stderr,"Charge balance failed at z=%g, nH=%g, G0=%g, X0=%g, step=%d\n",
                test->redshift,test->nH,test->isrf,test->ionizationMultiplier,step);
        exit(EXIT_FAILURE);
    }
    const double turnover = fmax(ionization,recombination);
    const double relativeTransfer = turnover > 0.0 ? fabs(transfer)/turnover : 0.0;
    const double nO = state->n[ipOXYGEN];
    const char *mode = test->mode == EQC_CHEMISTRY_EQUILIBRIUM ? "equilibrium" : "nonequilibrium";
    fprintf(file,"%g,%.9g,%g,%g,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                 "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                 "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            test->redshift,test->nH,test->isrf,test->ionizationMultiplier,mode,
            test->hiiRegion,test->chargeTransfer,step,state->NHIuvb[0],state->T,state->nel/test->nH,
            ionization,recombination,transferIonization,transferRecombination,
            hydrogenResidual,transfer,hydrogenResidual+transfer,relativeTransfer,
            chargeFromSpecies/test->nH,chargeError/test->nH,
            oxygenIonization,oxygenRecombination,netOxygenIonization,netOxygenIonization-transfer,
            nHI/test->nH,nHII/test->nH,state->nmol[ipH2]/test->nH,nOI/nO,nOII/nO,nOIII/nO);
}

static DiskResult RunDiskUpdates(FILE *file, DiskCase test, bool verbose, FILE *budget)
{
    EQCLibSetRunParamSoftXrayIonization(test.softXray);
    EQCLibSetRunParamHIIRegion(test.hiiRegion);
    EQCLibSetRunParamChargeTransfer(test.chargeTransfer);
    if(test.uvb) RequireSuccess("disk redshift",EQCLibTrySetUVBackgroundRedshift(test.redshift));
    StructEQCData state = DiskInitialState(&test);
    EQCLibStatus status = EQC_SUCCESS;
    /* Start the equilibrium cases with the original optically thin equilibrium.
       The non-equilibrium cases start with half the H ionized at the same T. */
    if(test.mode == EQC_CHEMISTRY_EQUILIBRIUM)
        status = EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.e-8);
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);
    int completed = 0;
    while(status == EQC_SUCCESS && completed < DiskSteps){
        /* Apply the effective neutral column to both faces using only the
           accepted composition. The H I column alone is a comparison case. */
        const double molecularNuclei = test.includeMolecularColumn ? 2.0*state.nmol[ipH2] : 0.0;
        const double column = state.Nt*((state.natom[ipHYDROGEN]+molecularNuclei)/test.nH);
        state.NHIuvb[0] = state.NHIuvb[1] = column;
        status = EQCLibTryAdvanceCell(&state,1.e-4,test.mode);
        if(status == EQC_SUCCESS){
            completed++;
            if(budget) WriteChargeTransferBudget(budget,&test,&state,completed);
        }
    }
    const double heatingFraction = state.H_tot > 0.0 ? state.H_UVB/state.H_tot : NAN;
    const char *modeName = test.mode == EQC_CHEMISTRY_EQUILIBRIUM ? "equilibrium" : "nonequilibrium";
    fprintf(file,"%g,%.9g,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%d,%s,"
                 "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%g,%g,%g,%.17g,%.17g,%.17g",
            test.redshift,test.nH,test.softXray,test.includeMolecularColumn,completed,status,state.NHIuvb[0],
            state.T,state.nel/test.nH,heatingFraction,test.uvb,modeName,
            state.H_UVB,state.H_XR,state.H_CR,state.H_PE,state.H_tot,
            test.metallicity,test.isrf,test.ionizationMultiplier,test.maximumEnergy,
            state.natom[ipHYDROGEN]/test.nH,state.nion[ipHYDROGEN]/test.nH,state.nmol[ipH2]/test.nH);
    const double nO = state.n[ipOXYGEN];
    fprintf(file,",%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            test.hiiRegion,test.chargeTransfer,state.natom[ipOXYGEN]/nO,state.nion[ipOXYGEN]/nO,state.nOIII/nO,
            state.L_CII158,state.L_OI63,state.L_OI146,state.L_OII,state.L_OIII,state.L_SiII35,state.L_res);
    if(status != EQC_SUCCESS || verbose)
        printf("Disk z=%g nH=%g UVB=%s X-rays=%s %s column=%s: %d/%d steps, status %d, T=%.2f K, UVB heating %.1f%%\n",
               test.redshift,test.nH,test.uvb ? "on" : "off",test.softXray ? "on" : "off",modeName,
               test.includeMolecularColumn ? "HI+2H2" : "HI",completed,DiskSteps,status,state.T,100.0*heatingFraction);
    if(status == EQC_SUCCESS && (!isfinite(state.T) || state.T <= 0.0 || !isfinite(heatingFraction) ||
                                (test.uvb ? state.H_UVB <= 0.0 : state.H_UVB != 0.0) ||
                                (test.softXray ? state.H_XR <= 0.0 : state.H_XR != 0.0))){
        fprintf(stderr,"Unexpected heating or temperature in the disk update for z=%g, nH=%g (%s)\n",
                test.redshift,test.nH,modeName);
        exit(EXIT_FAILURE);
    }
    return (DiskResult){completed,status,state.T,heatingFraction};
}

static DiskResult WriteDiskUpdates(FILE *file, DiskCase test, bool verbose)
{
    return RunDiskUpdates(file,test,verbose,NULL);
}

static const char DiskHeader[] =
    "z,nH,soft_xray,include_molecular_column,completed_updates,status,N_per_side,T,xe,H_UVB_over_H_tot,uvb,chemistry_mode,"
    "H_UVB,H_XR,H_CR,H_PE,H_tot,Z,G0,X0,energy_max_eV,xHI,xHII,xH2,"
    "hii_mode,charge_transfer,xOI,xOII,xOIII,L_CII158,L_OI63,L_OI146,L_OII,L_OIII,L_SiII35,L_res\n";

static void ConfigureDiskPhysics(void)
{
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    RequireSuccess("disk chemistry accuracy",EQCLibSetRunParamChemistryAccuracy(DiskChemistryAccuracy));
}

static void CheckDiskUpdates(double metallicity)
{
    ConfigureDiskPhysics();
    FILE *file = OpenOutput("disk-updates.csv");
    fputs(DiskHeader,file);
    DiskCase test = DiskDefaults;
    test.nH = 0.01;
    test.metallicity = metallicity;
    test.includeMolecularColumn = false;
    const DiskResult diffuse = WriteDiskUpdates(file,test,true);
#ifdef EQC_LEGACY_XRAY_CUTOFF
    if(diffuse.status != EQC_NO_EQUILIBRIUM || diffuse.completed != 1){
        fprintf(stderr,"The legacy soft X-ray cutoff should fail the second diffuse disk update\n");
        exit(EXIT_FAILURE);
    }
#else
    RequireSuccess("diffuse disk update with soft X-rays",diffuse.status);
#endif
    test.softXray = false;
    RequireSuccess("diffuse disk update without soft X-rays",WriteDiskUpdates(file,test,true).status);
    test.nH = 10.0;
    for(int softXray=0; softXray<=1; softXray++){
        test.softXray = softXray;
        for(int z=0; z<=2; z+=2){
            test.redshift = z;
            test.includeMolecularColumn = false;
            const DiskResult atomic = WriteDiskUpdates(file,test,true);
            RequireSuccess("molecular disk update using the atomic column",atomic.status);
            test.includeMolecularColumn = true;
            const DiskResult effective = WriteDiskUpdates(file,test,true);
            RequireSuccess("molecular disk update using the effective column",effective.status);
            if(effective.uvbHeatingFraction >= atomic.uvbHeatingFraction ||
               effective.temperature >= atomic.temperature){
                fprintf(stderr,"The molecular column proxy should reduce fractional UVB heating and temperature in this test\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    Close(file);
}

static void CheckDiskGrid(double metallicity)
{
    FILE *file = OpenOutput("disk-grid.csv");
    fputs(DiskHeader,file);
    const EQCLibChemistryMode modes[] = {EQC_CHEMISTRY_EQUILIBRIUM,EQC_CHEMISTRY_NONEQUILIBRIUM};
    DiskCase test = DiskDefaults;
    test.metallicity = metallicity;
    int successful = 0, cases = 0;
    for(int uvb=1; uvb>=0; uvb--){
        test.uvb = uvb;
        if(!test.uvb) RequireSuccess("disabling UVB",EQCLibSetRunParamUVBackground(false,NULL));
        for(int xrays=0; xrays<=1; xrays++){
            test.softXray = xrays;
            for(int z=0; z<=2; z+=2){
                test.redshift = z;
                for(size_t i=0; i<sizeof(modes)/sizeof(modes[0]); i++){
                    test.mode = modes[i];
                    for(int density=-6; density<=2; density++){
                        test.nH = pow(10.0,0.5*density);
                        const DiskResult result = WriteDiskUpdates(file,test,false);
                        cases++;
                        if(result.status == EQC_SUCCESS) successful++;
#ifdef EQC_LEGACY_XRAY_CUTOFF
                        if(result.status != EQC_NO_EQUILIBRIUM && result.status != EQC_NO_CONVERGENCE)
                            RequireSuccess("legacy disk grid update",result.status);
#endif
                    }
                }
            }
        }
    }
    Close(file);
    printf("Disk grid: %d/%d cases completed all %d updates\n",successful,cases,DiskSteps);
#ifndef EQC_LEGACY_XRAY_CUTOFF
    if(successful != cases) exit(EXIT_FAILURE);
#endif
}

static void CheckMetalSensitivity(const char *directory)
{
    ConfigureDiskPhysics();
    RequireSuccess("sensitivity UVB configuration",EQCLibSetRunParamUVBackground(true,directory));
    FILE *file = OpenOutput("disk-500ev.csv");
    fputs(DiskHeader,file);
    const EQCLibChemistryMode modes[] = {EQC_CHEMISTRY_EQUILIBRIUM,EQC_CHEMISTRY_NONEQUILIBRIUM};
    const double multipliers[] = {0.1,1.0};
    DiskCase test = DiskDefaults;
    int successful = 0;
    for(int capped=0;capped<=1;capped++){
        test.maximumEnergy = capped ? 500.0 : 0.0;
        RequireSuccess("sensitivity energy cutoff",EQCLibSetRunParamUVBackgroundMaximumEnergy(test.maximumEnergy));
        RequireSuccess("sensitivity initialization",EQCLibTryInitialization(0,NULL));
        for(int z=0;z<=2;z+=2){
            test.redshift = z;
            for(size_t g0=0;g0<sizeof(multipliers)/sizeof(multipliers[0]);g0++){
                test.isrf = multipliers[g0];
                for(size_t x0=0;x0<sizeof(multipliers)/sizeof(multipliers[0]);x0++){
                    test.ionizationMultiplier = multipliers[x0];
                    for(size_t mode=0;mode<sizeof(modes)/sizeof(modes[0]);mode++){
                        test.mode = modes[mode];
                        for(int density=-6;density<=4;density++){
                            test.nH = pow(10.0,0.5*density);
                            const DiskResult result = WriteDiskUpdates(file,test,false);
                            RequireSuccess("500 eV sensitivity disk update",result.status);
                            successful++;
                        }
                    }
                }
            }
        }
    }
    Close(file);
    printf("500 eV disk sensitivity: %d cases completed all %d updates\n",successful,DiskSteps);
}

static void CheckMetalPhotoionizationDisk(const char *directory)
{
    ConfigureDiskPhysics();
    RequireSuccess("metal photoionization disk UVB",EQCLibSetRunParamUVBackground(true,directory));
    RequireSuccess("metal photoionization disk initialization",EQCLibTryInitialization(0,NULL));
    FILE *file = OpenOutput("disk-metal-photoionization.csv");
    fputs(DiskHeader,file);
    const EQCLibChemistryMode modes[] = {EQC_CHEMISTRY_EQUILIBRIUM,EQC_CHEMISTRY_NONEQUILIBRIUM};
    DiskCase test = DiskDefaults;
    int successful = 0;
    for(int hii=0;hii<=1;hii++){
        test.hiiRegion = hii;
        for(int transfer=0;transfer<=1;transfer++){
            test.chargeTransfer = transfer;
            for(int z=0;z<=2;z+=2){
                test.redshift = z;
                for(size_t i=0;i<sizeof(modes)/sizeof(modes[0]);i++){
                    test.mode = modes[i];
                    for(int density=-6;density<=4;density++){
                        test.nH = pow(10.0,0.5*density);
                        const DiskResult result = WriteDiskUpdates(file,test,false);
                        RequireSuccess("metal photoionization disk update",result.status);
                        successful++;
                    }
                }
            }
        }
    }
    Close(file);
    printf("Metal photoionization disk: %d cases completed all %d updates\n",successful,DiskSteps);
}

static void CheckChargeTransferDisk(const char *directory)
{
    ConfigureDiskPhysics();
    RequireSuccess("charge transfer UVB",EQCLibSetRunParamUVBackground(true,directory));
    RequireSuccess("charge transfer initialization",EQCLibTryInitialization(0,NULL));
    FILE *grid = OpenOutput("disk-charge-transfer.csv");
    FILE *budget = OpenOutput("disk-charge-transfer-budget.csv");
    fputs(DiskHeader,grid);
    fputs("z,nH,G0,X0,chemistry_mode,hii_mode,charge_transfer,step,N_per_side,T,xe,"
          "H_ionizations,H_recombinations,H_CT_ionizations,H_CT_recombinations,"
          "H_non_CT_residual,H_CT_net,H_full_residual,abs_net_CT_over_H_turnover,"
          "xe_from_species,charge_error_per_H,O_ionizations,O_recombinations,O_non_CT_net,O_full_charge_residual,"
          "xHI,xHII,xH2,xOI,xOII,xOIII\n",
          budget);
    const EQCLibChemistryMode modes[] = {EQC_CHEMISTRY_EQUILIBRIUM,EQC_CHEMISTRY_NONEQUILIBRIUM};
    const double radiation[] = {0.1,1.0};
    DiskCase test = DiskDefaults;
    int successful = 0;
    for(int z=0;z<=2;z+=2){
        test.redshift = z;
        for(size_t g=0;g<sizeof radiation/sizeof radiation[0];g++){
            test.isrf = radiation[g];
            for(size_t x=0;x<sizeof radiation/sizeof radiation[0];x++){
                test.ionizationMultiplier = radiation[x];
                for(int density=-6;density<=4;density++){
                    test.nH = pow(10.0,0.5*density);
                    for(int treatment=0;treatment<3;treatment++){
                        test.hiiRegion = treatment == 2;
                        test.chargeTransfer = treatment != 0;
                        for(size_t i=0;i<sizeof modes/sizeof modes[0];i++){
                            test.mode = modes[i];
                            const DiskResult result = RunDiskUpdates(grid,test,false,budget);
                            RequireSuccess("charge transfer disk update",result.status);
                            successful++;
                        }
                    }
                }
            }
        }
    }
    Close(grid);
    Close(budget);
    printf("Charge transfer disk: %d cases completed all %d updates\n",successful,DiskSteps);
}

int main(int argc, char **argv)
{
    const bool zeroMetallicity = argc == 3 && strcmp(argv[2],"--zero-metallicity") == 0;
    const bool metalSensitivity = argc == 3 && strcmp(argv[2],"--metal-sensitivity") == 0;
    const bool metalPhotoionizationDisk = argc == 3 && strcmp(argv[2],"--metal-photoionization-disk") == 0;
    const bool chargeTransferDisk = argc == 3 && strcmp(argv[2],"--charge-transfer-disk") == 0;
    if(argc != 2 && !zeroMetallicity && !metalSensitivity && !metalPhotoionizationDisk && !chargeTransferDisk){
        fprintf(stderr,"Usage: %s FG20_DIRECTORY [--zero-metallicity|--metal-sensitivity"
                       "|--metal-photoionization-disk|--charge-transfer-disk]\n",argv[0]);
        return EXIT_FAILURE;
    }
    if(metalSensitivity){
        CheckMetalSensitivity(argv[1]);
        return EXIT_SUCCESS;
    }
    if(metalPhotoionizationDisk){
        CheckMetalPhotoionizationDisk(argv[1]);
        return EXIT_SUCCESS;
    }
    if(chargeTransferDisk){
        CheckChargeTransferDisk(argv[1]);
        return EXIT_SUCCESS;
    }
    EQCLibSetRunParamRecombinationCaseB(false);
    RequireSuccess("UVB configuration",EQCLibSetRunParamUVBackground(true,argv[1]));
    RequireSuccess("initialization",EQCLibTryInitialization(1,NULL));

    FILE *rates = OpenOutput("rates.csv");
    FILE *equilibrium = OpenOutput("equilibrium.csv");
    FILE *xrays = OpenOutput("soft-xray.csv");
    FILE *secondary = OpenOutput("secondary.csv");
    fputs("z,NHI_top,NHI_bottom,absorber,Gamma,heat_per_absorber,power_low,power_high\n",rates);
    fputs("z,nH,xe,xe_reference,xHII,xHII_reference,xHeII,xHeII_reference,xHeIII,xHeIII_reference\n",equilibrium);
    fputs("z,NHI_per_side,Nt,xe,H_UVB,H_XR,ratio\n",xrays);
    fputs("z,NHI_per_side,xHII,energy_max_eV,P_total,P_high,H_UVB,Q_other,Gamma_HI,Gamma_HeI,Gamma2_HI,Gamma2_HeI\n",secondary);
    CEQUVBEquilibriumStats stats = {0};
    for(size_t i = 0; i < sizeof Redshifts/sizeof Redshifts[0]; i++){
        const double redshift = Redshifts[i];
        RequireSuccess("redshift selection",EQCLibTrySetUVBackgroundRedshift(redshift));
        WriteRates(rates,redshift);
        WriteEquilibrium(equilibrium,redshift,&stats);
        WriteXrayComparison(xrays,redshift);
    }
    printf("SvS fixed-temperature equilibrium: %.2f cascade iterations per charge trial (%lu iterations, %lu trials)\n",
           (double)stats.iterations/stats.chargeTrials,(unsigned long)stats.iterations,(unsigned long)stats.chargeTrials);
    Close(rates);
    Close(equilibrium);
    Close(xrays);
    WriteDeepColumns(secondary,0.0,0.0);
    WriteDeepColumns(secondary,2.0,0.0);
    RequireSuccess("UVB photon cutoff",EQCLibSetRunParamUVBackgroundMaximumEnergy(300.0));
    RequireSuccess("cutoff initialization",EQCLibTryInitialization(0,NULL));
    WriteDeepColumns(secondary,0.0,300.0);
    WriteDeepColumns(secondary,2.0,300.0);
    Close(secondary);
    RequireSuccess("restore full spectrum",EQCLibSetRunParamUVBackgroundMaximumEnergy(0.0));
    RequireSuccess("full-spectrum initialization",EQCLibTryInitialization(0,NULL));
    const double metallicity = zeroMetallicity ? 0.0 : DiskDefaults.metallicity;
    CheckDiskUpdates(metallicity);
    CheckDiskGrid(metallicity);
    return EXIT_SUCCESS;
}
