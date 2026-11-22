#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <gsl/gsl_integration.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"
#include "EquilibriumAbundance.h"
#include "MoleculeChemistry.h"
#include "PhotoelectronDeposition.h"
#include "PhotochemicalGasPhaseChemistry.h"

enum { HI, HeI, HeII, SpeciesCount };

typedef struct {
    double ionization[SpeciesCount];
    double heating[SpeciesCount];
} RadiationRates;

static const double Thresholds[] = {13.60, 24.59, 54.42};
static const double Planck = 6.62607015e-27;
static const double ElectronVolt = 1.602176634e-12;
static const double UpperEnergy = 3000.0;
static const double Pi = 3.14159265358979323846;
static const char SpectrumDirectory[] = "synthetic-uvb";
static const char SpectrumPath[] = "synthetic-uvb/fg20_cloudy.ascii";

static void CheckRelative(const char *label, double actual, double expected, double tolerance)
{
    if(!isfinite(actual) || fabs(actual-expected) > tolerance*fabs(expected)){
        fprintf(stderr,"%s: got %.12g, expected %.12g (relative tolerance %g)\n",
                label,actual,expected,tolerance);
        exit(EXIT_FAILURE);
    }
}

static double Intensity(double energy)
{
    return 2.0e-23*pow(energy/Thresholds[HI],-1.5);
}

static void WriteSpectrum(void)
{
    /* Cloudy's redshift-major J_nu format, with no nodes at the ionization edges. */
    const double hcEVAngstrom = 12398.419843320027;
    const double energies[] = {3000.0, 500.0, 70.0, 40.0, 19.0, 8.0};
    enum { EnergyCount = sizeof energies/sizeof energies[0] };
    assert(mkdir(SpectrumDirectory,0700) == 0);
    FILE *file = fopen(SpectrumPath,"w");
    assert(file != NULL);
    fprintf(file,"# synthetic power-law cgs per steradian\n20060612\n1\n1\nz\n2\n%d\nlambda\n1\nF_nu\n1\n0 2\n",
            EnergyCount);
    for(int i = 0; i < EnergyCount; i++) fprintf(file,"%.17g\n",hcEVAngstrom/energies[i]);
    for(int z = 0; z <= 2; z += 2)
        for(int i = 0; i < EnergyCount; i++) fprintf(file,"%.17g\n",(1.0+1.5*z)*Intensity(energies[i]));
    assert(fclose(file) == 0);
}

typedef struct {
    int species;
    bool heating;
} Integral;

static double UnattenuatedIntegrand(double energy, void *data)
{
    const Integral *integral = data;
    double (*const crossSections[])(double) = {PhotoionizationCrossSectionHI,
        PhotoionizationCrossSectionHeI,PhotoionizationCrossSectionHeII};
    const double photons = 4.0*Pi*Intensity(energy)/(Planck*energy);
    const double rate = photons*crossSections[integral->species](energy);
    return integral->heating ? rate*(energy-Thresholds[integral->species])*ElectronVolt : rate;
}

static double IndependentThinRateFrom(int species, bool heating, double lower, double limit)
{
    if(limit <= lower) return 0.0;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(512);
    assert(workspace != NULL);
    Integral integral = {.species = species, .heating = heating};
    gsl_function function = {.function = UnattenuatedIntegrand, .params = &integral};
    double result, error;
    assert(gsl_integration_qag(&function,lower,limit,0.0,1.0e-9,512,GSL_INTEG_GAUSS61,
                               workspace,&result,&error) == 0);
    gsl_integration_workspace_free(workspace);
    return result;
}

static double IndependentThinRate(int species, bool heating, double limit)
{
    return IndependentThinRateFrom(species,heating,Thresholds[species],limit);
}

static RadiationRates Hemisphere(double column)
{
    RadiationRates result;
    assert(EQCLibTryGetUVBackgroundHemisphereRates(column,result.ionization,result.heating) == EQC_SUCCESS);
    return result;
}

static void CheckRateTable(void)
{
    assert(EQCLibTrySetUVBackgroundRedshift(0.0) == EQC_SUCCESS);
    const RadiationRates base = Hemisphere(0.0);
    for(int i = 0; i < SpeciesCount; i++){
        CheckRelative("thin ionization, two hemispheres",2.0*base.ionization[i],
                      IndependentThinRate(i,false,UpperEnergy),4.0e-4);
        CheckRelative("thin heating, two hemispheres",2.0*base.heating[i],
                      IndependentThinRate(i,true,UpperEnergy),4.0e-4);
    }

    RadiationRates previous = base;
    for(int exponent = 16; exponent <= 22; exponent++){
        const RadiationRates current = Hemisphere(pow(10.0,exponent));
        for(int i = 0; i < SpeciesCount; i++){
            assert(current.ionization[i] > 0.0 && current.ionization[i] < previous.ionization[i]);
            assert(current.heating[i] > 0.0 && current.heating[i] < previous.heating[i]);
            assert(current.heating[i]/current.ionization[i] > previous.heating[i]/previous.ionization[i]);
        }
        previous = current;
    }

    assert(EQCLibTrySetUVBackgroundRedshift(2.0) == EQC_SUCCESS);
    const RadiationRates upper = Hemisphere(0.0);
    assert(EQCLibTrySetUVBackgroundRedshift(1.0) == EQC_SUCCESS);
    const RadiationRates middle = Hemisphere(0.0);
    assert(EQCLibTrySetUVBackgroundRedshift(-0.1) != EQC_SUCCESS);
    assert(EQCLibTrySetUVBackgroundRedshift(2.1) != EQC_SUCCESS);
    assert(EQCLibTrySetUVBackgroundRedshift(NAN) != EQC_SUCCESS);
    const RadiationRates retained = Hemisphere(0.0);
    for(int i = 0; i < SpeciesCount; i++){
        CheckRelative("upper redshift ionization",upper.ionization[i],4.0*base.ionization[i],1.0e-12);
        CheckRelative("upper redshift heating",upper.heating[i],4.0*base.heating[i],1.0e-12);
        CheckRelative("midpoint redshift ionization",middle.ionization[i],2.0*base.ionization[i],1.0e-12);
        CheckRelative("midpoint redshift heating",middle.heating[i],2.0*base.heating[i],1.0e-12);
        assert(retained.ionization[i] == middle.ionization[i]);
        assert(retained.heating[i] == middle.heating[i]);
    }

    double ionization[] = {1.0,2.0,3.0}, heating[] = {4.0,5.0,6.0};
    assert(EQCLibTryGetUVBackgroundHemisphereRates(-1.0,ionization,heating) == EQC_INVALID_INPUT);
    assert(EQCLibTryGetUVBackgroundHemisphereRates(INFINITY,ionization,heating) == EQC_INVALID_INPUT);
    for(int i = 0; i < SpeciesCount; i++){
        assert(ionization[i] == 1.0+i);
        assert(heating[i] == 4.0+i);
    }
}

static StructEQCData AtomicGas(double hydrogenDensity, double heliumDensity)
{
    StructEQCData state = {.T = 1.0e4, .n = {hydrogenDensity,heliumDensity},
                          .natom = {hydrogenDensity,heliumDensity}};
    state.Rho = EQCLibGetTotalMassDensity(&state);
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);
    return state;
}

static double PureHydrogenEquilibrium(double photoionization, double nH, double temperature)
{
    const double collisional = nH*HICollisionalIonizationRate(temperature);
    const double recombination = nH*HIIRecombinationRateCaseA(temperature);
    const double linear = photoionization-collisional;
    return 2.0*photoionization/(linear + sqrt(linear*linear +
           4.0*photoionization*(collisional+recombination)));
}

static double Heating(const StructEQCData *state, RadiationRates top, RadiationRates bottom)
{
    return state->natom[ipHYDROGEN]*(top.heating[HI]+bottom.heating[HI]) +
           state->natom[ipHELIUM]*(top.heating[HeI]+bottom.heating[HeI]) +
           state->nion[ipHELIUM]*(top.heating[HeII]+bottom.heating[HeII]);
}

static void CheckChemistry(const char *pegaseDirectory)
{
    const double accuracy = 1.0e-8;
    /* Keep the independent primary-only analytic fixture below all fast-electron edges. */
    assert(EQCLibSetRunParamUVBackgroundMaximumEnergy(110.0) == EQC_SUCCESS);
    assert(EQCLibTryInitialization(0,NULL) == EQC_SUCCESS);
    assert(EQCLibTrySetUVBackgroundRedshift(0.0) == EQC_SUCCESS);
    const RadiationRates top = Hemisphere(1.0e18), bottom = Hemisphere(0.0);
    const double gamma = top.ionization[HI]+bottom.ionization[HI];
    StructEQCData neutral = AtomicGas(1.0,0.0);
    neutral.NHIuvb[0] = 1.0e18;
    const double expected = PureHydrogenEquilibrium(gamma,neutral.n[ipHYDROGEN],neutral.T);

    StructEQCData shortStep = neutral;
    shortStep.dt = 1.0e-5/gamma;
    assert(EQCLibTryAdvanceHydrogenChemistry(&shortStep,accuracy) == EQC_SUCCESS);
    CheckRelative("initial UVB-only ionization slope",shortStep.nion[ipHYDROGEN]/shortStep.dt,gamma,1.0e-4);

    StructEQCData longStep = neutral;
    longStep.dt = 100.0/gamma;
    assert(EQCLibTryAdvanceHydrogenChemistry(&longStep,accuracy) == EQC_SUCCESS);
    CheckRelative("long BDF step reaches fixed-T equilibrium",longStep.nion[ipHYDROGEN],expected,2.0e-6);

    StructEQCData equilibrium = neutral;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&equilibrium,accuracy) == EQC_SUCCESS);
    CheckRelative("UVB-only equilibrium ionization",equilibrium.nion[ipHYDROGEN],expected,2.0e-7);
    CheckRelative("UVB-only equilibrium heating",equilibrium.H_UVB,Heating(&equilibrium,top,bottom),1.0e-12);
    assert(equilibrium.H_PI == 0.0 && equilibrium.H_UVB == equilibrium.H_tot);

    StructEQCData thermo = AtomicGas(1.0,0.1);
    thermo.NHIuvb[0] = 1.0e18;
    thermo.dt = 1.0e-3/gamma;
    StructEQCData implicit = thermo;
    assert(EQCLibTryAdvanceCell(&thermo,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM) == EQC_SUCCESS);
    assert(thermo.nion[ipHYDROGEN] > 0.0 && thermo.nion[ipHELIUM] > 0.0 && thermo.nHeIII > 0.0);
    CheckRelative("BDF UVB heating",thermo.H_UVB,Heating(&thermo,top,bottom),1.0e-12);
    assert(EQCLibTryAdvanceCell(&implicit,accuracy,EQC_CHEMISTRY_EQUILIBRIUM) == EQC_SUCCESS);
    CheckRelative("implicit thermal UVB heating",implicit.H_UVB,Heating(&implicit,top,bottom),1.0e-12);
    assert(implicit.NHIuvb[0] == 1.0e18 && implicit.NHIuvb[1] == 0.0);

    StructEQCData invalid = thermo;
    invalid.NHIuvb[1] = NAN;
    const double oldEnergy = invalid.U, oldHeating = invalid.H_UVB;
    assert(EQCLibTryAdvanceCell(&invalid,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM) == EQC_INVALID_INPUT);
    assert(invalid.U == oldEnergy && invalid.H_UVB == oldHeating);
    assert(EQCLibTryAdvanceCell(&invalid,accuracy,EQC_CHEMISTRY_EQUILIBRIUM) == EQC_INVALID_INPUT);
    assert(invalid.U == oldEnergy && invalid.H_UVB == oldHeating);

    EQCLibSetRunParamHIIRegion(true);
    assert(EQCLibTryInitialization(0,pegaseDirectory) == EQC_SUCCESS);
    assert(EQCLibTrySetUVBackgroundRedshift(0.0) == EQC_SUCCESS);
    neutral = AtomicGas(1.0,0.1);
    neutral.NLyn = 1.0e-10/HIPhotoIonizationRate(1.0);
    StructEQCData ionized = neutral;
    ionized.natom[ipHYDROGEN] = ionized.natom[ipHELIUM] = 0.0;
    ionized.nion[ipHYDROGEN] = ionized.n[ipHYDROGEN];
    ionized.nHeIII = ionized.n[ipHELIUM];
    ionized.nel = EQCLibGetElectronNumberDensity(&ionized);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&neutral,accuracy) == EQC_SUCCESS);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&ionized,accuracy) == EQC_SUCCESS);
    CheckRelative("mixed stellar and UVB H II branches",neutral.nion[ipHYDROGEN],ionized.nion[ipHYDROGEN],2.0e-7);
    CheckRelative("mixed stellar and UVB He III branches",neutral.nHeIII,ionized.nHeIII,2.0e-7);
    CheckRelative("mixed heating tally",neutral.H_tot,neutral.H_PI+neutral.H_UVB,1.0e-12);
    CheckRelative("stellar H_PI is separate",neutral.H_PI,HeatingPhotoionization(&neutral),1.0e-12);
}

static void CheckEnergyBudget(CEQPhotoelectronDeposition deposit, double total)
{
    assert(deposit.heating >= 0.0 && deposit.otherPower >= 0.0);
    assert(deposit.secondaryHI >= 0.0 && deposit.secondaryHeI >= 0.0);
    CheckRelative("photoelectron energy budget",deposit.heating+deposit.otherPower+
                  EV_TO_ERG*(Thresholds[HI]*deposit.secondaryHI+Thresholds[HeI]*deposit.secondaryHeI),
                  total,8.0e-15);
}

static void CheckSecondaryChemistry(void)
{
    const double accuracy = 1.0e-9;
    EQCLibSetRunParamHIIRegion(false);
    assert(EQCLibSetRunParamUVBackgroundMaximumEnergy(0.0) == EQC_SUCCESS);
    assert(EQCLibTryInitialization(0,NULL) == EQC_SUCCESS);
    CEQUVBRates thin;
    assert(CEQEvaluateUVBackgroundHemisphere(0.0,&thin) == EQC_SUCCESS);
    for(int i=0;i<SpeciesCount;i++)
        CheckRelative("fast electron power, two unshielded hemispheres",2.0*thin.heatingHigh[i],
                      IndependentThinRateFrom(i,true,Thresholds[i]+100.0,UpperEnergy),4.0e-4);
    for(int shielded=0;shielded<=1;shielded++){
        /* Include a redshift and columns between table nodes. */
        assert(EQCLibTrySetUVBackgroundRedshift(shielded ? 0.37 : 0.0) == EQC_SUCCESS);
        StructEQCData initial = AtomicGas(1.0,0.1);
        initial.NHIuvb[0] = shielded ? 4.7e19 : 0.0;
        initial.NHIuvb[1] = shielded ? 8.31e19 : 0.0;
        CEQUVBRates uvb;
        assert(CEQEvaluateUVBackground(initial.NHIuvb,&uvb) == EQC_SUCCESS);
        const double electronPower = initial.natom[ipHYDROGEN]*uvb.heating[HI] +
                                     initial.natom[ipHELIUM]*uvb.heating[HeI];
        CheckEnergyBudget(CEQDepositUVBPhotoelectrons(&initial,&uvb),electronPower);

        StructEQCData equilibrium = initial, fixed = initial, thermo = initial;
        fixed.dt = 100.0/uvb.ionization[HI];
        thermo.dt = 1.0e-3/uvb.ionization[HI];
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&equilibrium,accuracy) == EQC_SUCCESS);
        assert(CEQAdvanceFixedTemperatureHydrogenHeliumChemistry(&fixed,accuracy) == EQC_SUCCESS);
        assert(fixed.T == initial.T);
        CheckRelative("coupled fixed-T H II",fixed.nion[ipHYDROGEN],equilibrium.nion[ipHYDROGEN],3.0e-6);
        CheckRelative("coupled fixed-T He II",fixed.nion[ipHELIUM],equilibrium.nion[ipHELIUM],3.0e-6);
        CheckRelative("coupled fixed-T He III",fixed.nHeIII,equilibrium.nHeIII,3.0e-6);
        assert(EQCLibTryAdvanceCell(&thermo,accuracy,EQC_CHEMISTRY_NONEQUILIBRIUM) == EQC_SUCCESS);
        CheckRelative("coupled BDF heat",thermo.H_UVB,CEQDepositUVBPhotoelectrons(&thermo,&uvb).heating,1.0e-14);
    }

    const double low = 1.0e-27, high = 2.0e-27;
    StructEQCData gas = AtomicGas(1.0,0.1);
    const CEQPhotoelectronDeposition neutral = CEQDepositPhotoelectrons(&gas,low,high);
    gas.nion[ipHYDROGEN] = 1.0e-4;
    gas.natom[ipHYDROGEN] = 1.0-gas.nion[ipHYDROGEN];
    const CEQPhotoelectronDeposition floor = CEQDepositPhotoelectrons(&gas,low,high);
    assert(neutral.heating == floor.heating);
    CheckRelative("SvS lower fit limit",neutral.secondaryHI,floor.secondaryHI,8.0e-15);
    CheckEnergyBudget(floor,low+high);

    gas.nion[ipHYDROGEN] = 0.01;
    gas.natom[ipHYDROGEN] = 0.99;
    const CEQPhotoelectronDeposition atomic = CEQDepositPhotoelectrons(&gas,low,high);
    gas.natom[ipHYDROGEN] *= 1.0e-8;
    gas.nmol[ipH2] = 0.5*(0.99-gas.natom[ipHYDROGEN]);
    gas.n[ipHELIUM] = gas.natom[ipHELIUM] *= 1.0e-8;
    CEQPhotoelectronDeposition spent = CEQDepositPhotoelectrons(&gas,low,high);
    CheckRelative("H I depletion",spent.secondaryHI,1.0e-8*atomic.secondaryHI,1.0e-12);
    CheckRelative("trace helium",spent.secondaryHeI,1.0e-8*atomic.secondaryHeI/0.99,1.0e-12);
    CheckEnergyBudget(spent,low+high);
    gas.natom[ipHYDROGEN] = gas.natom[ipHELIUM] = 0.0;
    spent = CEQDepositPhotoelectrons(&gas,low,high);
    assert(spent.secondaryHI == 0.0 && spent.secondaryHeI == 0.0);
    gas.nmol[ipH2] = 0.0;
    gas.nion[ipHYDROGEN] = nextafter(1.0,0.0);
    gas.natom[ipHYDROGEN] = 1.0-gas.nion[ipHYDROGEN];
    spent = CEQDepositPhotoelectrons(&gas,low,high);
    CheckRelative("SvS fully ionized approach",spent.heating,low+0.9971*high,8.0e-15);
    assert(spent.secondaryHI <= 1.0e-10*atomic.secondaryHI && spent.secondaryHeI == 0.0);
    gas.nion[ipHYDROGEN] = 1.0;
    gas.natom[ipHYDROGEN] = 0.0;
    spent = CEQDepositPhotoelectrons(&gas,low,high);
    assert(spent.secondaryHI == 0.0 && spent.secondaryHeI == 0.0);
    assert(spent.heating == low+0.9971*high);
    CheckEnergyBudget(spent,low+high);
    puts("SvS energy conservation, target limits and coupled H/He fixed-temperature BDF checks passed.");
}

static void CheckEnergyCutoff(void)
{
    EQCLibSetRunParamHIIRegion(false);
    assert(EQCLibSetRunParamUVBackgroundMaximumEnergy(30.0) == EQC_SUCCESS);
    assert(EQCLibTryInitialization(0,NULL) == EQC_SUCCESS);
    assert(EQCLibTrySetUVBackgroundRedshift(0.0) == EQC_SUCCESS);
    const RadiationRates limited = Hemisphere(0.0);
    for(int i = 0; i < SpeciesCount; i++){
        CheckRelative("limited spectrum ionization",2.0*limited.ionization[i],
                      IndependentThinRate(i,false,30.0),4.0e-4);
        CheckRelative("limited spectrum heating",2.0*limited.heating[i],
                      IndependentThinRate(i,true,30.0),4.0e-4);
    }
}

static void CheckMolecularCarbonPhotoionization(void)
{
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamHIIRegion(false);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamCOFormation(true,1.0);

    const double nH = 100.0, accuracy = 1.0e-8;
    StructEQCData state = AtomicGas(nH,0.1*nH);
    state.T = 60.0;
    state.G0 = 0.1;
    state.Nt = 1.0e21;
    state.Zd = 0.01;
    state.n[ipCARBON] = state.natom[ipCARBON] = 1.4e-4*nH;
    state.n[ipOXYGEN] = state.natom[ipOXYGEN] = 3.2e-4*nH;
    CEQUVBRates uvb;
    assert(CEQEvaluateUVBackground(state.NHIuvb,&uvb) == EQC_SUCCESS);
    assert(EQCLibChemicalReactionEquilibrium(&state,accuracy) == EQC_SUCCESS);
    assert(state.nmol[ipH2] > 0.0 && state.nmol[ipCO] > 0.0 && state.nCIII > 0.0);

    const double photoCI = 3.1e-10*state.G0*DustShieldingFactor(state.Nt,state.Zd,3.0)
                         + uvb.metalIonization[CEQ_UVB_CI];
    const double formation = 4.0e-16*pow(state.T/300.0,-0.2)*state.nmol[ipH2];
    const double formationFlux = formation*state.nion[ipCARBON];
    const double doppler = hypot(12.90*sqrt(1.0e-4*state.T),3.0);
    const double coShield = COShieldingFactor(state.Nt*state.nmol[ipH2]/nH,state.Nt*state.nmol[ipCO]/nH,doppler);
    const double dissociation = COPhotoDissociationRate(state.G0,state.Nt,state.Zd)*coShield;
    CheckRelative("UVB molecular C I balance",photoCI*state.natom[ipCARBON],formationFlux,accuracy);
    CheckRelative("UVB molecular C III balance",uvb.metalIonization[CEQ_UVB_CII]*state.nion[ipCARBON],
                  state.nel*CIIIRecombinationRate(state.T)*state.nCIII,accuracy);
    CheckRelative("UVB molecular CO balance",dissociation*state.nmol[ipCO],formationFlux,2.0*accuracy);
    const StructEQCData original = state;
    assert(CEQSolveCarbonEquilibrium(&state,state.nmol[ipCO],accuracy) == EQC_SUCCESS);
    CheckRelative("standalone UVB CO",state.nmol[ipCO],original.nmol[ipCO],2.0*accuracy);
    CheckRelative("standalone UVB molecular C III",state.nCIII,original.nCIII,2.0*accuracy);
    EQCLibRunParameters = saved;
    puts("UVB C I, C III and CO molecular fluxes balance.");
}

static void CheckMetalNonEquilibriumPhotoionization(void)
{
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamChargeTransfer(false);
    StructEQCData state = AtomicGas(1.0,0.1);
    state.T = 100.0;
    state.n[ipCARBON] = state.natom[ipCARBON] = 1.4e-4;
    state.n[ipOXYGEN] = state.natom[ipOXYGEN] = 3.2e-4;
    state.n[ipSILICON] = 0x1p-20;
    state.natom[ipSILICON] = 0.625*state.n[ipSILICON];
    state.nion[ipSILICON] = 0.25*state.n[ipSILICON];
    state.nSiIII = 0.125*state.n[ipSILICON];
    state.nel = EQCLibGetElectronNumberDensity(&state);
    state.Rho = EQCLibGetTotalMassDensity(&state);
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);

    CEQUVBRates uvb;
    assert(CEQEvaluateUVBackground(state.NHIuvb,&uvb) == EQC_SUCCESS);
    state.dt = 3.0e-4/uvb.ionization[HI];
    const StructEQCData initial = state;
    assert(EQCLibTryAdvanceCell(&state,1.0e-10,EQC_CHEMISTRY_NONEQUILIBRIUM) == EQC_SUCCESS);
    assert(state.T < 150.0 && state.n[ipSILICON] == initial.n[ipSILICON]);
    assert(state.natom[ipSILICON] == initial.natom[ipSILICON]);
    assert(state.nion[ipSILICON] == initial.nion[ipSILICON] && state.nSiIII == initial.nSiIII);
    /* With no stellar, FUV, cosmic-ray or soft X-ray field, the short step
       follows the leading terms Γ_0 dt and Γ_0 Γ_1 dt²/2. */
    const double carbon = uvb.metalIonization[CEQ_UVB_CI]*state.dt;
    const double oxygen = uvb.metalIonization[CEQ_UVB_OI]*state.dt;
    CheckRelative("BDF UVB C II fraction",state.nion[ipCARBON]/state.n[ipCARBON],carbon,5.0e-3);
    CheckRelative("BDF UVB O II fraction",state.nion[ipOXYGEN]/state.n[ipOXYGEN],oxygen,5.0e-3);
    CheckRelative("BDF UVB C III fraction",state.nCIII/state.n[ipCARBON],
                  0.5*carbon*uvb.metalIonization[CEQ_UVB_CII]*state.dt,1.0e-2);
    CheckRelative("BDF UVB O III fraction",state.nOIII/state.n[ipOXYGEN],
                  0.5*oxygen*uvb.metalIonization[CEQ_UVB_OII]*state.dt,1.0e-2);
    EQCLibRunParameters = saved;
    puts("UVB BDF ionizes C and O while preserving the Si composition.");
}

static void CheckFailedSpectrumKeepsRates(const RadiationRates *previous)
{
    assert(EQCLibTryInitialization(0,NULL) == EQC_INITIALIZATION_FAILED);
    const RadiationRates retained = Hemisphere(1.0e18);
    for(int i = 0; i < SpeciesCount; i++){
        assert(retained.ionization[i] == previous->ionization[i]);
        assert(retained.heating[i] == previous->heating[i]);
    }
}

static void CheckBadSpectrum(void)
{
    assert(EQCLibTrySetUVBackgroundRedshift(1.0) == EQC_SUCCESS);
    const RadiationRates previous = Hemisphere(1.0e18);

    /* The same valid table under the HM12 filename cannot replace FG20. */
    const char hm12Path[] = "synthetic-uvb/hm12_galaxy.ascii";
    assert(rename(SpectrumPath,hm12Path) == 0);
    fputs("[expected test diagnostics] Missing FG20 spectrum with only HM12 present.\n",stderr);
    CheckFailedSpectrumKeepsRates(&previous);
    assert(rename(hm12Path,SpectrumPath) == 0);

    FILE *file = fopen(SpectrumPath,"w");
    assert(file != NULL);
    fputs("20060612\n1\n1\nz\n1\n2\nlambda\n1\nF_nu\n1\n0\n10 100\n1e-23 nan\n",file);
    assert(fclose(file) == 0);
    fputs("[expected test diagnostics] Malformed UVB intensity.\n",stderr);
    CheckFailedSpectrumKeepsRates(&previous);
    assert(EQCLibSetRunParamUVBackground(false,NULL) == EQC_SUCCESS);
}

int main(int argc, char **argv)
{
    assert(argc == 2);
    assert(EQCLibTrySetUVBackgroundRedshift(0.0) != EQC_SUCCESS);
    WriteSpectrum();
    EQCLibSetRunParamRecombinationCaseB(false);
    assert(EQCLibSetRunParamUVBackground(true,SpectrumDirectory) == EQC_SUCCESS);
    assert(EQCLibTryInitialization(0,NULL) == EQC_SUCCESS);
    CheckRateTable();
    CheckChemistry(argv[1]);
    CheckSecondaryChemistry();
    CheckMolecularCarbonPhotoionization();
    CheckMetalNonEquilibriumPhotoionization();
    CheckEnergyCutoff();
    CheckBadSpectrum();
    puts("Synthetic UVB rates, shielding, redshift, and equilibrium/non-equilibrium tests passed.");
    return EXIT_SUCCESS;
}
