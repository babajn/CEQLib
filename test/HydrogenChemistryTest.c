#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "CoolingCIE.h"
#include "CoolingHeatingRates.h"
#include "FineStructureLines.h"
#include "HIIregionLineCooling.h"
#include "HyperFineStructureLines.h"
#include "MoleculeChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"

static const double accuracy = 1.0e-8;

static StructEQCData HydrogenState(double nH, double temperature, double ionized, double molecular)
{
    return (StructEQCData){
        .t = YEAR_CGS,
        .U = 1.0e12,
        .T = temperature,
        .Rho = ChemicalElementMass[ipHYDROGEN]*nH,
        .n = {[ipHYDROGEN] = nH},
        .natom = {[ipHYDROGEN] = (1.0-ionized-molecular)*nH},
        .nion = {[ipHYDROGEN] = ionized*nH},
        .nmol = {[ipH2] = 0.5*molecular*nH},
        .nel = ionized*nH,
    };
}

static double IonizedFraction(const StructEQCData *state)
{
    return state->nion[ipHYDROGEN]/state->n[ipHYDROGEN];
}

static double MolecularFraction(const StructEQCData *state)
{
    return 2.0*state->nmol[ipH2]/state->n[ipHYDROGEN];
}

static void CheckConservation(const StructEQCData *state, const StructEQCData *initial)
{
    const double nH = initial->n[ipHYDROGEN];
    const double nHI = state->natom[ipHYDROGEN];
    const double nHII = state->nion[ipHYDROGEN];
    const double nH2 = state->nmol[ipH2];
    assert(nHI >= 0.0 && nHII >= 0.0 && nH2 >= 0.0);
    assert(fabs(nHI+nHII+2.0*nH2-nH) <= 8.0*DBL_EPSILON*nH);
    assert(state->n[ipHYDROGEN] == nH);

    const double nHe = initial->n[ipHELIUM];
    const double nHeI = state->natom[ipHELIUM];
    const double nHeII = state->nion[ipHELIUM];
    assert(nHeI >= 0.0 && nHeII >= 0.0 && state->nHeIII >= 0.0);
    assert(fabs(nHeI+nHeII+state->nHeIII-nHe) <= 8.0*DBL_EPSILON*nHe);
    assert(state->n[ipHELIUM] == nHe);

    const double nC = initial->n[ipCARBON];
    const double nCI = state->natom[ipCARBON];
    const double nCII = state->nion[ipCARBON];
    assert(nCI >= 0.0 && nCII >= 0.0 && state->nCIII >= 0.0);
    assert(fabs(nCI+nCII+state->nCIII-nC) <= 8.0*DBL_EPSILON*nC);
    assert(state->n[ipCARBON] == nC);

    const double nO = initial->n[ipOXYGEN];
    const double nOI = state->natom[ipOXYGEN];
    const double nOII = state->nion[ipOXYGEN];
    assert(nOI >= 0.0 && nOII >= 0.0 && state->nOIII >= 0.0);
    assert(fabs(nOI+nOII+state->nOIII-nO) <= 8.0*DBL_EPSILON*nO);
    const double charge = nHII+nHeII+nCII+nOII+2.0*(state->nHeIII+state->nCIII+state->nOIII);
    assert(state->nel >= 0.0 && fabs(state->nel-charge) <= 8.0*DBL_EPSILON*charge);
    assert(state->n[ipOXYGEN] == nO);

    assert(state->dt == initial->dt && state->t == initial->t && state->Rho == initial->Rho);
    assert(state->NLyn == initial->NLyn && state->G0 == initial->G0 && state->Nt == initial->Nt);
    assert(state->Td == initial->Td && state->Zd == initial->Zd && state->X0 == initial->X0);
}

static void Advance(StructEQCData *state, double duration)
{
    state->dt = duration;
    const StructEQCData initial = *state;
    const EQCLibStatus status = EQCLibTryAdvanceHydrogenChemistry(state,accuracy);
    if(status != EQC_SUCCESS){
        fprintf(stderr,"Hydrogen chemistry failed: status=%d T=%g n_H=%g dt=%g s NLyn=%g G0=%g\n",
                status,initial.T,initial.n[ipHYDROGEN],duration,initial.NLyn,initial.G0);
    }
    assert(status == EQC_SUCCESS);
    CheckConservation(state,&initial);
    assert(state->T == initial.T && state->U == initial.U);
}

static void CheckIonizingLightSwitch(void)
{
    const double nH = 10.0;
    StructEQCData state = HydrogenState(nH,1000.0,0.0,0.0);
    const double recombination = HIIRecombinationRate(state.T)*nH;
    const double time = 1.0/recombination;
    const double crossSection = HIPhotoIonizationRate(1.0);
    assert(crossSection > 0.0);
    state.NLyn = recombination/crossSection;
    StructEQCData equilibrium = state;
    assert(EQCLibChemicalReactionEquilibrium(&equilibrium,1.0e-10) >= 0);
    const double target = IonizedFraction(&equilibrium);

    /* H I and electrons are reconstructed from the integrated species. */
    state.natom[ipHYDROGEN] = 0.0;
    state.nel = nH;
    Advance(&state,0.25*time);
    const double lit = IonizedFraction(&state);
    assert(lit > 0.0 && lit < 0.9*target);
    assert(state.nmol[ipH2] == 0.0);

    Advance(&state,30.0*time);
    const double steady = IonizedFraction(&state);
    assert(fabs(steady-target) < 20.0*accuracy);

    state.NLyn = 0.0;
    StructEQCData split = state;
    const double duration = 0.5*time;
    /* Collisional ionization is negligible at this temperature. */
    const double expected = steady/(1.0+recombination*duration*steady);
    Advance(&state,duration);
    for(int step = 0; step < 8; step++) Advance(&split,duration/8.0);
    const double dark = IonizedFraction(&state);
    const double darkSplit = IonizedFraction(&split);
    assert(dark > 0.0 && dark < steady);
    assert(fabs(dark-expected) < 20.0*accuracy);
    assert(fabs(darkSplit-expected) < 20.0*accuracy);
    assert(fabs(dark-darkSplit) < 20.0*accuracy);
    printf("Ionizing light: x_HII=%.8f after turn-on, %.8f at equilibrium, %.8f after turn-off (reference %.8f)\n",
           lit,steady,dark,expected);
}

static void CheckMolecularHydrogen(void)
{
    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamDustAttenuation(true);

    const double nH = 1.0e3;
    StructEQCData state = HydrogenState(nH,100.0,0.03,0.05);
    state.Td = 20.0;
    state.Zd = SOLAR_METALLICITY;
    state.Nt = 1.0e18;
    state.G0 = 0.5;
    state.NLyn = 0.25*nH*HIIRecombinationRate(state.T)/HIPhotoIonizationRate(1.0);
    const double time = 1.0/(2.0*nH*H2FormationRate(state.T,state.Zd,state.Td));
    StructEQCData equilibrium = state;
    assert(EQCLibChemicalReactionEquilibrium(&equilibrium,1.0e-10) >= 0);
    const double weakTarget = MolecularFraction(&equilibrium);

    Advance(&state,0.25*time);
    const double forming = MolecularFraction(&state);
    assert(forming > 0.05 && forming < 0.9*weakTarget);
    Advance(&state,40.0*time);
    const double weakFuv = MolecularFraction(&state);
    assert(fabs(weakFuv-weakTarget) < 20.0*accuracy);
    const double weakIonized = IonizedFraction(&state);
    assert(fabs(weakIonized-IonizedFraction(&equilibrium)) < 20.0*accuracy);

    state.G0 = 10.0;
    equilibrium = state;
    assert(EQCLibChemicalReactionEquilibrium(&equilibrium,1.0e-10) >= 0);
    const double strongTarget = MolecularFraction(&equilibrium);

    Advance(&state,0.1*time);
    const double dissociating = MolecularFraction(&state);
    assert(dissociating < weakFuv && dissociating > strongTarget);
    Advance(&state,40.0*time);
    const double strongFuv = MolecularFraction(&state);
    assert(fabs(strongFuv-strongTarget) < 20.0*accuracy);
    const double strongIonized = IonizedFraction(&state);
    assert(fabs(strongIonized-IonizedFraction(&equilibrium)) < 20.0*accuracy);
    assert(strongIonized > weakIonized);
    printf("FUV: 2*n_H2/n_H=%.8f while forming, %.8f at weak FUV, %.8f after FUV increase, %.8f at strong FUV\n",
           forming,weakFuv,dissociating,strongFuv);
    printf("At fixed ionizing flux, x_HII=%.8f at weak FUV and %.8f at strong FUV\n",
           weakIonized,strongIonized);
}

static void CheckFullyMolecularStart(void)
{
    StructEQCData state = HydrogenState(10.0,100.0,0.0,1.0);
    state.G0 = 1.0e3;
    const double dissociation = H2PhotoDissociationRate(state.G0,state.Nt,0.0,3.0,state.Zd);
    /* At 100 K with zero dust abundance and column, photodissociation dominates. */
    Advance(&state,1.0/dissociation);
    assert(fabs(MolecularFraction(&state)-exp(-1.0)) < 20.0*accuracy);
    Advance(&state,1.0e5/dissociation);
    assert(MolecularFraction(&state) < accuracy);
}

static void CheckUnchanged(StructEQCData *state, EQCLibStatus expected)
{
    unsigned char before[sizeof(*state)];
    memcpy(before,state,sizeof(*state));
    assert(EQCLibTryAdvanceHydrogenChemistry(state,accuracy) == expected);
    assert(memcmp(before,state,sizeof(*state)) == 0);
}

static void CheckInvalidInputs(void)
{
    StructEQCData state = HydrogenState(10.0,100.0,0.25,0.5);
    state.nel = 0.0;
    CheckUnchanged(&state,EQC_SUCCESS);

    state.dt = -1.0;
    CheckUnchanged(&state,EQC_INVALID_INPUT);

    state.dt = YEAR_CGS;
    state.n[ipHELIUM] = 1.0;
    CheckUnchanged(&state,EQC_INVALID_INPUT);

    state.n[ipHELIUM] = 0.0;
    state.nion[ipHYDROGEN] = 7.5;
    CheckUnchanged(&state,EQC_INVALID_INPUT);

    state.nion[ipHYDROGEN] = 2.5;
    EQCLibSetRunParamH2Formation(false,1.0);
    CheckUnchanged(&state,EQC_INVALID_INPUT);
}

static void SetInternalEnergy(StructEQCData *state)
{
    state->MeanMolecularWeight = EQCLibGetMeanMolecularWeight(state);
    state->DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(state);
    state->U = EQCLibConvertTtoU(state);
}

static void CheckRelative(double actual, double expected, double tolerance)
{
    assert(isfinite(actual) && isfinite(expected));
    assert(fabs(actual-expected) <= tolerance*fabs(expected));
}

static void AdvanceThermochemistry(StructEQCData *state, double duration)
{
    state->dt = duration;
    const StructEQCData initial = *state;
    const EQCLibStatus status = EQCLibTryAdvanceHydrogenThermochemistry(state,accuracy);
    if(status != EQC_SUCCESS){
        fprintf(stderr,"Hydrogen thermochemistry failed: status=%d U=%g n_H=%g dt=%g s NLyn=%g G0=%g\n",
                status,initial.U,initial.n[ipHYDROGEN],duration,initial.NLyn,initial.G0);
    }
    assert(status == EQC_SUCCESS);
    CheckConservation(state,&initial);
    assert(state->Du == initial.Du && state->DQheat == initial.DQheat && state->Z == initial.Z);
    assert(state->U > 0.0 && state->T > 0.0);
    CheckRelative(state->U,EQCLibConvertTtoU(state),1.0e-12);
    CheckRelative(state->MeanMolecularWeight,EQCLibGetMeanMolecularWeight(state),1.0e-12);
    CheckRelative(state->DegreeOfFreedom,EQCLibGetMeanDegreeOfFreedom(state),1.0e-12);
    CheckRelative(state->Gamma,1.0+1.0/state->DegreeOfFreedom,1.0e-12);
    CheckRelative(state->Gm1,state->Gamma-1.0,1.0e-12);
    CheckRelative(state->GGm1,state->Gamma*state->Gm1,1.0e-12);
    CheckRelative(state->invT,1.0/state->T,1.0e-12);
    CheckRelative(state->logT,log10(state->T),1.0e-12);

    StructEQCData rates = *state;
    assert(EQCLibCalculateCoolingHeatingRates(&rates) == EQC_SUCCESS);
    const double rateTolerance = EQCLibRunParameters.UseNonEquilibriumTemperatureTables ? 1.0e-8 : 1.0e-12;
    CheckRelative(state->H_tot,rates.H_tot,rateTolerance);
    /* The non-equilibrium update uses instantaneous lines without aggregate CIE. */
    const double oxygenII = CoolingLinesOII(state->T,state->nel,state->nion[ipOXYGEN],NULL);
    const double oxygenIII = CoolingLinesOIII(state->T,state->nel,state->nOIII,NULL);
    const double carbonII = CEQCoolingCII158um(state);
    double molecularCooling = 0.0;
    if(state->nmol[ipH2] > 0.0){
        molecularCooling = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(
            state->T,state->nion[ipHYDROGEN],state->natom[ipHYDROGEN],state->nmol[ipH2],
            state->natom[ipHELIUM],state->nel) + CoolingH2CollisionalDissociation(&rates);
    }
    const double cooling = rates.L_tot-rates.L_res-rates.L_OII-rates.L_OIII-rates.L_OI63-rates.L_OI146
                         - rates.L_CII158-rates.L_H2+oxygenII+oxygenIII+state->L_OI63+state->L_OI146
                         + carbonII+molecularCooling;
    assert(state->L_res == 0.0);
    CheckRelative(state->L_OII,oxygenII,rateTolerance);
    CheckRelative(state->L_OIII,oxygenIII,rateTolerance);
    CheckRelative(state->L_tot,cooling,rateTolerance);
    CheckRelative(state->H_PI,rates.H_PI,1.0e-12);
    CheckRelative(state->L_HI,rates.L_HI,rateTolerance);
    CheckRelative(state->L_HeI,rates.L_HeI,rateTolerance);
    CheckRelative(state->L_HeII,rates.L_HeII,rateTolerance);
    CheckRelative(state->L_CII158,carbonII,rateTolerance);
    CheckRelative(state->H_H2,rates.H_H2,rateTolerance);
    CheckRelative(state->L_H2,molecularCooling,rateTolerance);
}

static void CheckThermochemicalLightSwitch(void)
{
    const double nH = 10.0;
    const double photoionization = 1.0e-9;
    StructEQCData initial = HydrogenState(nH,8000.0,0.0,0.0);
    initial.NLyn = photoionization/HIPhotoIonizationRate(1.0);
    initial.Du = 123.0;
    initial.DQheat = 456.0;
    initial.Td = 20.0;
    SetInternalEnergy(&initial);
    assert(EQCLibCalculateCoolingHeatingRates(&initial) == EQC_SUCCESS);

    StructEQCData state = initial;
    state.T = 100.0;
    state.MeanMolecularWeight = state.DegreeOfFreedom = 0.0;
    state.H_tot = state.L_tot = 0.0;
    StructEQCData shortStep = state;
    const double shortDuration = 1.0e-6/photoionization;
    AdvanceThermochemistry(&shortStep,shortDuration);
    assert(shortStep.T > initial.T && shortStep.T < 1.0001*initial.T);
    CheckRelative((shortStep.U-initial.U)/shortDuration,(initial.H_tot-initial.L_tot)/initial.Rho,1.0e-4);
    CheckRelative(IonizedFraction(&shortStep)/shortDuration,photoionization,1.0e-4);

    StructEQCData split = state;
    const double duration = 0.8/photoionization;
    AdvanceThermochemistry(&state,duration);
    for(int step = 0; step < 8; step++) AdvanceThermochemistry(&split,duration/8.0);
    CheckRelative(state.U,split.U,100.0*accuracy);
    CheckRelative(state.T,split.T,100.0*accuracy);
    CheckRelative(IonizedFraction(&state),IonizedFraction(&split),100.0*accuracy);
    const double litTemperature = state.T;
    const double litIonized = IonizedFraction(&state);
    const double litEnergy = state.U;
    assert(litTemperature > initial.T && litEnergy > initial.U && litIonized > 0.1);
    assert(state.H_PI > 0.0 && litIonized < 0.95);

    state.NLyn = 0.0;
    AdvanceThermochemistry(&state,0.25/(nH*HIIRecombinationRate(litTemperature)));
    const double darkIonized = IonizedFraction(&state);
    assert(state.T < litTemperature && state.U < litEnergy);
    assert(darkIonized > 0.0 && darkIonized < litIonized);
    assert(state.H_PI == 0.0 && state.L_tot > state.H_tot);
    printf("Thermochemistry: T=%.2f K, x_HII=%.8f after turn-on; T=%.2f K, x_HII=%.8f after turn-off\n",
           litTemperature,litIonized,state.T,darkIonized);
}

/* Without H2 formation, each output must be accepted as the next input. In hot,
   mostly ionized gas, the implicit solver used to leave roundoff H2 in some outputs. */
static void CheckRepeatedAtomicThermochemistry(void)
{
    EQCLibSetRunParamH2Formation(false,1.0);
    for(int it = 0; it <= 26; it++){
        for(int ix = 0; ix < 4; ix++){
            StructEQCData state = HydrogenState(1.0,1.0e4*pow(1.3,it),0.1+0.3*ix,0.0);
            SetInternalEnergy(&state);
            for(int step = 0; step < 2; step++){
                AdvanceThermochemistry(&state,1.0e4*YEAR_CGS);
                assert(state.nmol[ipH2] == 0.0);
            }
            Advance(&state,1.0e4*YEAR_CGS);
            assert(state.nmol[ipH2] == 0.0);
        }
    }
}

/* The update holds Si ionization fixed, so Si photoionization must not heat
   the gas. Before the fix, neutral Si heated this cell to about 8e4 K. */
static void CheckFixedSiliconPhotoionizationHeating(void)
{
    EQCLibSetRunParamH2Formation(false,1.0);
    const double nH = 10.0;
    const double photoionization = 1.0e-8;
    StructEQCData free = HydrogenState(nH,8000.0,0.0,0.0);
    free.NLyn = photoionization/HIPhotoIonizationRate(1.0);
    free.dt = 1.0e3*YEAR_CGS;
    SetInternalEnergy(&free);

    StructEQCData silicon = free;
    const double nSi = 3.548134e-6*nH;
    silicon.n[ipSILICON] = silicon.natom[ipSILICON] = nSi;
    silicon.Rho = EQCLibGetTotalMassDensity(&silicon);
    SetInternalEnergy(&silicon);

    for(int step = 0; step < 100; step++){
        assert(EQCLibTryAdvanceHydrogenThermochemistry(&free,accuracy) == EQC_SUCCESS);
        assert(EQCLibTryAdvanceHydrogenThermochemistry(&silicon,accuracy) == EQC_SUCCESS);
        assert(silicon.natom[ipSILICON] == nSi && silicon.nion[ipSILICON] == 0.0 && silicon.nSiIII == 0.0);
        const double hydrogenHeating = silicon.natom[ipHYDROGEN]*HeatingHIPhotoionization(silicon.NLyn);
        CheckRelative(silicon.H_PI,hydrogenHeating,1.0e-12);
        /* The equilibrium heating formula still counts the Si I photoionization. */
        assert(HeatingPhotoionization(&silicon) > silicon.H_PI);
    }
    CheckRelative(silicon.T,free.T,1.0e-2);
}

static void CheckMolecularThermochemistry(void)
{
    EQCLibSetRunParamH2Formation(true,1.0);
    StructEQCData state = HydrogenState(10.0,500.0,0.0,0.8);
    state.G0 = 1.0e3;
    state.Zd = SOLAR_METALLICITY;
    state.Td = 20.0;
    SetInternalEnergy(&state);
    const double initialU = state.U;
    const double dissociation = H2PhotoDissociationRate(state.G0,state.Nt,0.0,3.0,state.Zd);
    AdvanceThermochemistry(&state,1.0/dissociation);
    /* Photodissociation dominates formation under this unshielded field. */
    CheckRelative(MolecularFraction(&state),0.8*exp(-1.0),100.0*accuracy);
    assert(state.U > initialU && state.T > 500.0 && state.H_H2 > 0.0 && state.L_H2 > 0.0);
    const double litTemperature = state.T;
    const double litMolecular = MolecularFraction(&state);

    state.G0 = 0.0;
    AdvanceThermochemistry(&state,1.0e13);
    assert(state.T < litTemperature && state.L_tot > state.H_tot);
    assert(MolecularFraction(&state) > litMolecular);
    printf("Molecular thermochemistry: T=%.2f K, 2*n_H2/n_H=%.8f after FUV turn-on; "
           "T=%.2f K, 2*n_H2/n_H=%.8f after FUV turn-off\n",
           litTemperature,litMolecular,state.T,MolecularFraction(&state));
    EQCLibSetRunParamH2Formation(false,1.0);
}

static void SetOxygenAbundances(StructEQCData *state, double nO, double ionized, double doublyIonized)
{
    state->n[ipOXYGEN] = nO;
    state->nion[ipOXYGEN] = nO*ionized;
    state->nOIII = nO*doublyIonized;
    state->natom[ipOXYGEN] = nO-state->nion[ipOXYGEN]-state->nOIII;
    const double oxygenMass = nO*ChemicalElementMass[ipOXYGEN];
    state->Rho += oxygenMass;
    state->Z = oxygenMass/state->Rho;
    state->nel = EQCLibGetElectronNumberDensity(state);
}

static void SetHeliumAbundances(StructEQCData *state, double nHe, double ionized, double doublyIonized)
{
    state->n[ipHELIUM] = nHe;
    state->nion[ipHELIUM] = nHe*ionized;
    state->nHeIII = nHe*doublyIonized;
    state->natom[ipHELIUM] = nHe-state->nion[ipHELIUM]-state->nHeIII;
    state->Rho += nHe*ChemicalElementMass[ipHELIUM];
    state->nel = EQCLibGetElectronNumberDensity(state);
}

static void CheckNeutralMolecularThermochemistry(void)
{
    const bool grainProcesses = EQCLibRunParameters.UseGrainSurfaceProcesses;
    const bool dustAttenuation = EQCLibRunParameters.UseDustAttenuation;
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    const double molecular[] = {0.1, 1.0, 0.8};
    const double temperatures[] = {20.0, 20.0, 100.0};
    for(size_t i = 0; i < sizeof molecular/sizeof molecular[0]; i++){
        StructEQCData state = HydrogenState(100.0,temperatures[i],0.0,molecular[i]);
        SetHeliumAbundances(&state,10.0,0.0,0.0);
        state.G0 = 1.0e3;
        state.Nt = 2.0e20;
        state.Zd = SOLAR_METALLICITY;
        state.Td = 20.0;
        state.dt = 1.0e5*YEAR_CGS;
        state.Rho = EQCLibGetTotalMassDensity(&state);
        SetInternalEnergy(&state);
        const StructEQCData initial = state;
        assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,1.0e-6) == EQC_SUCCESS);
        CheckConservation(&state,&initial);
        assert(state.nel == 0.0 && state.nion[ipHYDROGEN] == 0.0 && state.nion[ipHELIUM] == 0.0 &&
               state.nHeIII == 0.0);
        assert(state.T > initial.T && MolecularFraction(&state) < molecular[i]);
    }
    EQCLibSetRunParamH2Formation(false,1.0);
    EQCLibSetRunParamGrainSurfaceProcesses(grainProcesses);
    EQCLibSetRunParamDustAttenuation(dustAttenuation);
}

static void SetCarbonAbundances(StructEQCData *state, double nC, double ionized, double doublyIonized)
{
    state->n[ipCARBON] = nC;
    state->nion[ipCARBON] = nC*ionized;
    state->nCIII = nC*doublyIonized;
    state->natom[ipCARBON] = nC-state->nion[ipCARBON]-state->nCIII;
    const double carbonMass = nC*ChemicalElementMass[ipCARBON];
    const double metalMass = state->Z*state->Rho+carbonMass;
    state->Rho += carbonMass;
    state->Z = metalMass/state->Rho;
    state->nel = EQCLibGetElectronNumberDensity(state);
}

static void CheckHeliumThermochemistry(void)
{
    const double nH = 10.0;
    const double nHe = 0.1*nH;
    const double photoionization = 1.0e-7;
    StructEQCData state = HydrogenState(nH,8000.0,0.0,0.0);
    SetHeliumAbundances(&state,nHe,0.0,0.0);
    SetInternalEnergy(&state);
    state.NLyn = photoionization/HIPhotoIonizationRate(1.0);
    const double initialU = state.U;
    const double heating = nH*HeatingHIPhotoionization(state.NLyn)+nHe*HeatingHeIPhotoionization(state.NLyn);
    state.natom[ipHELIUM] = 0.0;
    state.nel = 1.0;

    StructEQCData newlyLit = state;
    const double shortIllumination = 1.0e-6/photoionization;
    AdvanceThermochemistry(&newlyLit,shortIllumination);
    CheckRelative(newlyLit.nion[ipHELIUM]/shortIllumination,nHe*HeIPhotoIonizationRate(state.NLyn),1.0e-4);
    CheckRelative((newlyLit.U-initialU)/shortIllumination,heating/state.Rho,1.0e-4);

    AdvanceThermochemistry(&state,50.0/photoionization);
    assert(state.nion[ipHELIUM] > 0.5*nHe && state.nHeIII > 0.01*nHe);
    assert(state.T > 8000.0 && state.U > initialU && state.H_PI > 0.0);
    const StructEQCData lit = state;
    const double litHeating = lit.natom[ipHYDROGEN]*HeatingHIPhotoionization(lit.NLyn)
                            + lit.natom[ipHELIUM]*HeatingHeIPhotoionization(lit.NLyn)
                            + lit.nion[ipHELIUM]*HeatingHeIIPhotoionization(lit.NLyn);
    CheckRelative(lit.H_PI,litHeating,1.0e-12);
    CheckRelative(lit.L_Halpha,CoolingHalphaEmission(&state),1.0e-12);
    CheckRelative(lit.L_HI21cm,CoolingHI21cm(lit.natom[ipHYDROGEN],lit.T),1.0e-12);
    CheckRelative(lit.k_HI21cm,AbsorptionCoefficientHI21cm(lit.natom[ipHYDROGEN],lit.T),1.0e-12);
    const double particlesPerVolume = nH+nHe+lit.nel;
    CheckRelative(lit.MeanMolecularWeight,lit.Rho/(PROTON_MASS_CGS*particlesPerVolume),1.0e-12);

    state.NLyn = 0.0;
    StructEQCData initialDark = state;
    assert(EQCLibCalculateCoolingHeatingRates(&initialDark) == EQC_SUCCESS);
    assert(initialDark.H_tot == 0.0 && initialDark.L_HeI > 0.0 && initialDark.L_HeII > 0.0);
    const double recombination = lit.nel*HeIIIRecombinationRate(lit.T);
    const double heIIIChange = lit.nel*HeIICollisionalIonizationRate(lit.T)*lit.nion[ipHELIUM] - recombination*lit.nHeIII;
    const double heIIRecombination = HeIIRecombinationRate(lit.T)+HeIIDielectronicRecombinationRate(lit.T);
    const double heIIChange = lit.nel*(HeICollisionalIonizationRate(lit.T)*lit.natom[ipHELIUM]
                                     - heIIRecombination*lit.nion[ipHELIUM]) - heIIIChange;
    const double shortDark = 1.0e-5/recombination;
    AdvanceThermochemistry(&state,shortDark);
    assert(state.H_PI == 0.0);
    CheckRelative((state.nHeIII-lit.nHeIII)/shortDark,heIIIChange,2.0e-3);
    CheckRelative((state.nion[ipHELIUM]-lit.nion[ipHELIUM])/shortDark,heIIChange,2.0e-3);
    CheckRelative((lit.U-state.U)/shortDark,initialDark.L_tot/lit.Rho,2.0e-3);

    AdvanceThermochemistry(&state,1.0/(nH*HIIRecombinationRate(lit.T)));
    const double litCharge = lit.nion[ipHELIUM]+2.0*lit.nHeIII;
    const double darkCharge = state.nion[ipHELIUM]+2.0*state.nHeIII;
    assert(state.natom[ipHELIUM] > lit.natom[ipHELIUM] && state.nHeIII < lit.nHeIII);
    assert(darkCharge < litCharge && state.T < lit.T && state.U < lit.U);
    printf("Helium thermochemistry: T=%.2f K, x_HeII=%.8f, x_HeIII=%.8f after turn-on; "
           "T=%.2f K, x_HeII=%.8f, x_HeIII=%.8f after turn-off\n",
           lit.T,lit.nion[ipHELIUM]/nHe,lit.nHeIII/nHe,state.T,state.nion[ipHELIUM]/nHe,state.nHeIII/nHe);
}

static void CheckCarbonThermochemistry(void)
{
    const double nH = 10.0;
    const double nC = 3.0e-4*nH;
    const double photoionization = 1.0e-7;
    StructEQCData state = HydrogenState(nH,8000.0,0.0,0.0);
    SetCarbonAbundances(&state,nC,0.0,0.0);
    SetInternalEnergy(&state);
    state.NLyn = photoionization/HIPhotoIonizationRate(1.0);
    state.natom[ipCARBON] = 0.0;

    StructEQCData newlyLit = state;
    const double carbonPhotoionization = CIPhotoionizationRate(state.NLyn);
    const double shortIllumination = 1.0e-6/carbonPhotoionization;
    AdvanceThermochemistry(&newlyLit,shortIllumination);
    CheckRelative(newlyLit.nion[ipCARBON]/shortIllumination,nC*carbonPhotoionization,1.0e-4);

    AdvanceThermochemistry(&state,5.0/photoionization);
    assert(state.nion[ipCARBON] > 0.1*nC && state.nCIII > 0.01*nC);
    assert(state.T > 8000.0 && state.H_PI > 0.0 && state.L_CII158 > 0.0);
    const StructEQCData lit = state;
    const double litHeating = lit.natom[ipHYDROGEN]*HeatingHIPhotoionization(lit.NLyn)
                            + lit.natom[ipCARBON]*HeatingCIPhotoionization(lit.NLyn)
                            + lit.nion[ipCARBON]*HeatingCIIPhotoionization(lit.NLyn);
    CheckRelative(lit.H_PI,litHeating,1.0e-12);

    state.NLyn = 0.0;
    const double recombination = lit.nel*CIIIRecombinationRate(lit.T);
    const double cIIIChange = lit.nel*CIICollisionalIonizationRate(lit.T)*lit.nion[ipCARBON] - recombination*lit.nCIII;
    const double cIIChange = lit.nel*(CICollisionalIonizationRate(lit.T)*lit.natom[ipCARBON]
                                    - CIIRecombinationRate(lit.T)*lit.nion[ipCARBON]) - cIIIChange;
    const double shortDark = 1.0e-5/recombination;
    AdvanceThermochemistry(&state,shortDark);
    assert(state.H_PI == 0.0 && state.L_CII158 > 0.0);
    CheckRelative((state.nCIII-lit.nCIII)/shortDark,cIIIChange,2.0e-3);
    CheckRelative((state.nion[ipCARBON]-lit.nion[ipCARBON])/shortDark,cIIChange,2.0e-3);

    AdvanceThermochemistry(&state,1.0/(nH*HIIRecombinationRate(lit.T)));
    assert(state.natom[ipCARBON] > lit.natom[ipCARBON] && state.nCIII < lit.nCIII);
    assert(state.nion[ipCARBON]+2.0*state.nCIII < lit.nion[ipCARBON]+2.0*lit.nCIII);
    assert(state.T < lit.T && state.L_CII158 < lit.L_CII158);
    printf("Carbon thermochemistry: T=%.2f K, x_CII=%.8f, x_CIII=%.8f, L_CII=%.6e after turn-on; "
           "T=%.2f K, x_CII=%.8f, x_CIII=%.8f, L_CII=%.6e after turn-off\n",
           lit.T,lit.nion[ipCARBON]/nC,lit.nCIII/nC,lit.L_CII158,
           state.T,state.nion[ipCARBON]/nC,state.nCIII/nC,state.L_CII158);
}

static void CheckCarbonFineStructureCooling(void)
{
    const double nH = 100.0;
    const double nC = 3.0e-4*nH;
    StructEQCData state = HydrogenState(nH,200.0,0.0,0.0);
    SetCarbonAbundances(&state,nC,0.0,0.0);
    SetInternalEnergy(&state);
    state.G0 = 1.0e3;
    AdvanceThermochemistry(&state,1.0/CIFUVPhotoionizationRate(state.G0,0.0));
    CheckRelative(state.nion[ipCARBON]/nC,1.0-exp(-1.0),2.0e-6);
    CheckRelative(state.nel,state.nion[ipCARBON],1.0e-12);
    assert(state.nCIII == 0.0 && state.L_CII158 > 0.0);

    state.G0 = 0.0;
    const double initialU = state.U;
    const double initialTemperature = state.T;
    const double initialCooling = CEQCoolingCII158um(&state);
    const double shortDark = 1.0e6;
    AdvanceThermochemistry(&state,shortDark);
    CheckRelative((initialU-state.U)/shortDark,initialCooling/state.Rho,1.0e-4);
    AdvanceThermochemistry(&state,1.0e12);
    assert(state.T < initialTemperature && state.L_CII158 < initialCooling);
}

static void CheckCarbonCMBExchange(void)
{
    StructEQCData state = HydrogenState(1.0e4,2.5,0.0,0.0);
    SetCarbonAbundances(&state,3.0,1.0,0.0);
    const double colder = CEQCoolingCII158um(&state);
    state.T = 3.0;
    const double warmer = CEQCoolingCII158um(&state);
    assert(colder < 0.0 && warmer > 0.0);

    state.T = CMB_TEMPERATURE;
    const double atCMB = CEQCoolingCII158um(&state);
    assert(fabs(atCMB) <= 1.0e-10*fmin(-colder,warmer));
}

static void CheckMixedThermochemistry(void)
{
    EQCLibSetRunParamH2Formation(true,1.0);
    StructEQCData state = HydrogenState(100.0,8000.0,0.5,0.05);
    SetHeliumAbundances(&state,10.0,0.5,0.3);
    SetOxygenAbundances(&state,0.05,0.2,0.3);
    SetCarbonAbundances(&state,0.03,0.35,0.25);
    SetInternalEnergy(&state);
    const StructEQCData initial = state;
    AdvanceThermochemistry(&state,1.0e8);
    assert(state.nHeIII < initial.nHeIII && state.nOIII < initial.nOIII);
    assert(state.nCIII < initial.nCIII && state.L_CII158 > 0.0);
    assert(state.nmol[ipH2] > 0.0 && state.L_H2 > 0.0);
    assert(state.L_HeI > 0.0 && state.L_HeII > 0.0 && state.L_OII > 0.0);
    EQCLibSetRunParamH2Formation(false,1.0);
}

static void CheckCIEInterpolation(void)
{
    EQCLibSetRunParamCIEResonanceCooling(true);
    StructEQCData state = HydrogenState(10.0,1.0e4,1.0,0.0);
    const double metallicities[] = {0.001,0.03,1.0,4.0};
    for(size_t i=0;i<sizeof metallicities/sizeof metallicities[0];i++){
        state.Z = metallicities[i]*SOLAR_METALLICITY;
        for(int j=0;j<=90;j++){
            state.T = pow(10.0,4.0+0.05*j);
            /* Converting the grid through T and back to log T can round by one ULP. */
            CheckRelative(CEQNonEquilibriumCoolingRateCIE(&state),CoolingRateCIE(&state),2.0e-10);
        }
        state.T = 1.0e9;
        CheckRelative(CEQNonEquilibriumCoolingRateCIE(&state),CoolingRateCIE(&state),5.0e-12);
    }

    /* Compare the one-sided slopes at a knot within the blending interval. */
    state.Z = SOLAR_METALLICITY;
    const double step = 1.0e-6;
    double linear[3], smooth[3];
    for(int i=0;i<3;i++){
        state.T = pow(10.0,4.35+(i-1)*step);
        linear[i] = log10(CoolingRateCIE(&state));
        smooth[i] = log10(CEQNonEquilibriumCoolingRateCIE(&state));
    }
    assert(fabs(linear[0]-2.0*linear[1]+linear[2]) > 0.01*step);
    assert(fabs(smooth[0]-2.0*smooth[1]+smooth[2]) < 0.001*step);
    EQCLibSetRunParamCIEResonanceCooling(false);
}

static void CheckOxygenCoolingTransition(void)
{
    const double temperatures[] = {14999.0,15000.0,15001.0,20000.0,24999.0,25000.0,25001.0,1.0e7};
    const double weights[] = {1.0,1.0,0.999999970002,0.5,0.000000029998,0.0,0.0,0.0};
    EQCLibSetRunParamCIEResonanceCooling(true);
    for(size_t i = 0; i < sizeof temperatures/sizeof temperatures[0]; i++){
        StructEQCData state = HydrogenState(10.0,temperatures[i],0.98,0.0);
        SetOxygenAbundances(&state,5.0e-3,0.2,0.7);
        const double oxygenII = CoolingLinesOII(state.T,state.nel,state.nion[ipOXYGEN],NULL);
        const double oxygenIII = CoolingLinesOIII(state.T,state.nel,state.nOIII,NULL);
        const double cie = CEQNonEquilibriumCoolingRateCIE(&state);
        assert(oxygenII > 0.0 && oxygenIII > 0.0 && cie > 0.0);
        assert(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) == EQC_SUCCESS);
        CheckRelative(state.L_OII,weights[i]*oxygenII,1.0e-8);
        CheckRelative(state.L_OIII,weights[i]*oxygenIII,1.0e-8);
        CheckRelative(state.L_res,(1.0-weights[i])*cie,2.0e-8);
        const double totalWithCIE = state.L_tot;
        const double oxygenAndCIE = state.L_OII+state.L_OIII+state.L_res;

        EQCLibSetRunParamCIEResonanceCooling(false);
        assert(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) == EQC_SUCCESS);
        assert(state.L_res == 0.0);
        CheckRelative(state.L_OII,oxygenII,1.0e-8);
        CheckRelative(state.L_OIII,oxygenIII,1.0e-8);
        CheckRelative(totalWithCIE,state.L_tot-state.L_OII-state.L_OIII+oxygenAndCIE,1.0e-12);
        EQCLibSetRunParamCIEResonanceCooling(true);

        state.Z = 0.0;
        assert(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) == EQC_SUCCESS);
        assert(state.L_res == 0.0);
        CheckRelative(state.L_OII,oxygenII,1.0e-8);
        CheckRelative(state.L_OIII,oxygenIII,1.0e-8);
    }
    EQCLibSetRunParamCIEResonanceCooling(false);
}

static void CheckIlluminatedOxygenCoolingTransition(void)
{
    EQCLibSetRunParamCIEResonanceCooling(true);
    StructEQCData state = HydrogenState(10.0,17600.0,0.98,0.0);
    SetOxygenAbundances(&state,5.0e-3,0.2,0.7);
    state.invT = 1.0/state.T;
    state.logT = log10(state.T);
    SetInternalEnergy(&state);
    assert(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) == EQC_SUCCESS);
    assert(state.H_PI == 0.0 && state.L_res > 0.0 && state.L_OII > 0.0 && state.L_OIII > 0.0);
    const StructEQCData dark = state;

    /* At fixed composition, changing the flux only switches photoionization heating. */
    state.NLyn = 1.0e-8/HIPhotoIonizationRate(1.0);
    assert(CEQCalculateNonEquilibriumCoolingHeatingRates(&state) == EQC_SUCCESS);
    assert(state.H_PI > 0.0);
    assert(state.L_res == dark.L_res && state.L_OII == dark.L_OII && state.L_OIII == dark.L_OIII);
    printf("Illuminated transition at %.0f K: [O II]+[O III]=%.3e, CIE=%.3e erg cm^-3 s^-1\n",
           state.T,state.L_OII+state.L_OIII,state.L_res);
    state.dt = YEAR_CGS;
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_SUCCESS);
    assert(state.T > 15000.0 && state.T < 25000.0);
    assert(state.L_res > 0.0 && state.L_OII > 0.0 && state.L_OIII > 0.0);
    EQCLibSetRunParamCIEResonanceCooling(false);
}

static void CheckOxygenCoolingAfterTurnOff(void)
{
    EQCLibSetRunParamCIEResonanceCooling(true);
    const double nH = 10.0;
    const double nO = 5.0e-4*nH;
    const double photoionization = 1.0e-7;
    StructEQCData state = HydrogenState(nH,8000.0,0.0,0.0);
    SetOxygenAbundances(&state,nO,0.0,0.0);
    SetInternalEnergy(&state);
    state.NLyn = photoionization/HIPhotoIonizationRate(1.0);
    state.natom[ipOXYGEN] = 0.0; /* The update reconstructs neutral oxygen. */

    StructEQCData newlyLit = state;
    const double shortIllumination = 1.0e-6/photoionization;
    AdvanceThermochemistry(&newlyLit,shortIllumination);
    CheckRelative(newlyLit.nion[ipOXYGEN]/shortIllumination,nO*OIPhotoionizationRate(state.NLyn),1.0e-4);
    assert(newlyLit.L_OI63 > 0.0 && newlyLit.L_OI146 > 0.0);

    AdvanceThermochemistry(&state,3.0/photoionization);
    assert(state.nion[ipOXYGEN] > 0.1*nO && state.nOIII > 0.01*nO);
    assert(state.L_OII > 0.0 && state.L_OIII > 0.0);

    state.NLyn = 0.0;
    const double duration = 1.0e5;
    StructEQCData initialDark = state;
    assert(EQCLibCalculateCoolingHeatingRates(&initialDark) == EQC_SUCCESS);
    assert(initialDark.L_res > 0.0);
    const double otherCooling = initialDark.L_tot-initialDark.L_res-initialDark.L_OII-initialDark.L_OIII
                              - initialDark.L_OI63-initialDark.L_OI146+state.L_OI63+state.L_OI146;
    const double litOII = state.L_OII;
    const double litOIII = state.L_OIII;
    const double netCooling = (otherCooling+state.L_OII+state.L_OIII-initialDark.H_tot)/state.Rho;
    AdvanceThermochemistry(&state,duration);
    assert(state.H_PI == 0.0 && state.L_OII > 0.0 && state.L_OIII > 0.0);
    CheckRelative(state.L_OII,litOII,1.0e-2);
    CheckRelative(state.L_OIII,litOIII,1.0e-2);
    CheckRelative((initialDark.U-state.U)/duration,netCooling,2.0e-3);
    const double earlyOII = state.L_OII;
    const double earlyOIII = state.L_OIII;
    const double earlyHI = state.natom[ipHYDROGEN];
    const double earlyOI = state.natom[ipOXYGEN];
    const double earlyDoublyIonized = state.nOIII;

    AdvanceThermochemistry(&state,1.0e13);
    assert(state.natom[ipHYDROGEN] > earlyHI && state.natom[ipOXYGEN] > earlyOI);
    assert(state.natom[ipOXYGEN] > 0.75*nO);
    assert(state.nOIII < earlyDoublyIonized);
    assert(state.L_OII < earlyOII && state.L_OIII < earlyOIII);
    printf("Oxygen cooling after turn-off [erg cm^-3 s^-1]: O II %.6e -> %.6e; O III %.6e -> %.6e\n",
           earlyOII,state.L_OII,earlyOIII,state.L_OIII);
    EQCLibSetRunParamCIEResonanceCooling(false);
}

static void CheckHydrogenOxygenChargeTransfer(void)
{
    StructEQCData initial = HydrogenState(100.0,8000.0,0.4,0.0);
    SetOxygenAbundances(&initial,2.0,0.2,0.3);
    SetInternalEnergy(&initial);
    const double toOII = OIChargeTransferRateToHII(initial.T)*initial.natom[ipOXYGEN]*initial.nion[ipHYDROGEN];
    const double toOI = OIIChargeTransferRateToHI(initial.T)*initial.nion[ipOXYGEN]*initial.natom[ipHYDROGEN];
    const double fromOIII = OIIIChargeTransferRateToHI(initial.T)*initial.nOIII*initial.natom[ipHYDROGEN];
    const double duration = 1.0e4;

    EQCLibSetRunParamHIIRegion(false);
    StructEQCData off = initial;
    EQCLibSetRunParamChargeTransfer(false);
    AdvanceThermochemistry(&off,duration);
    StructEQCData on = initial;
    EQCLibSetRunParamChargeTransfer(true);
    AdvanceThermochemistry(&on,duration);
    EQCLibSetRunParamHIIRegion(true);

    /* The difference between the short-step derivatives isolates charge transfer. */
    const double tolerance = 4.0e-3;
    CheckRelative((on.nion[ipHYDROGEN]-off.nion[ipHYDROGEN])/duration,-toOII+toOI+fromOIII,tolerance);
    CheckRelative((on.nion[ipOXYGEN]-off.nion[ipOXYGEN])/duration,toOII-toOI+fromOIII,tolerance);
    CheckRelative((on.nOIII-off.nOIII)/duration,-fromOIII,tolerance);
}

static void CheckThermochemistryInputs(void)
{
    StructEQCData state = HydrogenState(10.0,100.0,0.25,0.0);
    state.T = state.nel = 0.0;
    unsigned char before[sizeof(state)];
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_SUCCESS);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state.dt = YEAR_CGS;
    state.Rho = 0.0;
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state.Rho = 10.0*ChemicalElementMass[ipHYDROGEN];
    state.U = 0.0;
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state = HydrogenState(10.0,100.0,0.25,0.0);
    SetOxygenAbundances(&state,0.01,0.25,0.5);
    state.natom[ipOXYGEN] = 0.0;
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_SUCCESS);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state.dt = YEAR_CGS;
    CheckUnchanged(&state,EQC_INVALID_INPUT);
    state.nion[ipOXYGEN] = 0.75*state.n[ipOXYGEN];
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state = HydrogenState(10.0,100.0,0.25,0.0);
    SetHeliumAbundances(&state,1.0,0.25,0.5);
    state.natom[ipHELIUM] = 0.0;
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_SUCCESS);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state.dt = YEAR_CGS;
    CheckUnchanged(&state,EQC_INVALID_INPUT);
    state.nHeIII = 0.8*state.n[ipHELIUM];
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state = HydrogenState(10.0,100.0,0.25,0.0);
    SetCarbonAbundances(&state,0.01,0.25,0.5);
    state.natom[ipCARBON] = 0.0;
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_SUCCESS);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state.dt = YEAR_CGS;
    CheckUnchanged(&state,EQC_INVALID_INPUT);
    state.nCIII = 0.8*state.n[ipCARBON];
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(state)) == 0);

    state.nCIII = 0.5*state.n[ipCARBON];
    state.nmol[ipCO] = 0.001;
    memcpy(before,&state,sizeof(state));
    assert(EQCLibTryAdvanceHydrogenThermochemistry(&state,accuracy) == EQC_INVALID_INPUT);
    assert(memcmp(before,&state,sizeof(state)) == 0);
}

int main(int argc, char **argv)
{
    assert(argc == 2);
    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamCIEResonanceCooling(false);
    assert(EQCLibTryInitialization(0,argv[1]) == EQC_SUCCESS);

    CheckIonizingLightSwitch();
    CheckMolecularHydrogen();
    CheckFullyMolecularStart();
    CheckInvalidInputs();
    CheckThermochemicalLightSwitch();
    CheckRepeatedAtomicThermochemistry();
    CheckFixedSiliconPhotoionizationHeating();
    CheckMolecularThermochemistry();
    CheckNeutralMolecularThermochemistry();
    CheckHeliumThermochemistry();
    CheckCarbonThermochemistry();
    CheckCarbonFineStructureCooling();
    CheckCarbonCMBExchange();
    CheckMixedThermochemistry();
    CheckCIEInterpolation();
    CheckOxygenCoolingTransition();
    CheckIlluminatedOxygenCoolingTransition();
    CheckOxygenCoolingAfterTurnOff();
    CheckHydrogenOxygenChargeTransfer();
    CheckThermochemistryInputs();
    puts("Hydrogen chemistry and thermochemistry checks passed.");
    return 0;
}
