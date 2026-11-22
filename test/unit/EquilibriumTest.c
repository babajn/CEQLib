#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "ChemicalEquilibriumCooling.h"
#include "Astro.h"
#include "Constants.h"
#include "CollisionalGasPhaseChemistry.h"
#include "CosmicRay.h"
#include "EquilibriumAbundance.h"
#include "FineStructureLines.h"
#include "HIIregionLineCooling.h"
#include "MoleculeChemistry.h"
#include "GrainSurfaceChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"
#include "Xray.h"

static StructEQCData GasState(double nH, double temperature, double column)
{
    StructEQCData state = {
        .dt = 10.0 * YEAR_CGS,
        .T = temperature,
        .Nt = column,
        .Z = ASPLUND_2009_SOLAR_METALLICITY,
        .Zd = SOLAR_METALLICITY,
        .G0 = 1.7,
        .X0 = 1.0,
        .n = {nH, 0.1*nH, 3.019952e-4*nH, 4.570882e-4*nH, 3.548134e-6*nH},
    };
    state.Td = EQCApproximateDustTemperature(state.Nt,state.G0);
    for(int i = 0; i < NumberChemicalElement; i++){
        state.natom[i] = state.n[i];
        state.Rho += state.n[i] * ChemicalElementMass[i];
    }
    state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
    state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    state.U = EQCLibConvertTtoU(&state);
    return state;
}

static void CheckAtomicEquilibrium(void)
{
    const double nH = 0.01584893;
    StructEQCData state = GasState(nH,10083.56,1.0e18);
    state.nel = 0.9*nH;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-9) == EQC_SUCCESS);
    /* This equilibrium lies inside the soft X-ray taper. */
#ifdef EQC_LEGACY_XRAY_CUTOFF
    const double expectedElectronFraction = 0.27455813;
#else
    const double expectedElectronFraction = 0.24120160; /* Solve with 1e-11 charge tolerance. */
#endif
    assert(fabs(state.nel/nH - expectedElectronFraction) < 1.0e-6);
    assert(fabs(EQCLibGetElectronNumberDensity(&state) - state.nel) <= 1.0e-9 * state.nel);
    assert(state.nmol[ipH2] == 0.0 && state.nmol[ipCO] == 0.0);
    assert(isfinite(state.H_tot) && isfinite(state.L_tot));
}

static void CheckDustThermalRateLimits(void)
{
    StructEQCData state = GasState(10.0,100.0,0.0);
    state.nel = 1.0;
    /* Reference rates from the original fits, with half weight at the log midpoint. */
    const double temperatures[] = {100.0,1.0e4,sqrt(1.0e9)};
    const double cooling[] = {1.2360981601193132e-26,1.5673988492951167e-24,2.4602999223636057e-24};
    const double heating[] = {8.2712456944609747e-25,1.3077417109213918e-24,9.8383559248390100e-25};
    for(size_t i = 0; i < sizeof temperatures/sizeof temperatures[0]; i++){
        state.T = temperatures[i];
        assert(fabs(CoolingDustRecombination(&state)/cooling[i]-1.0) < 1.0e-12);
        assert(fabs(HeatingDustPhotoelectic(&state)/heating[i]-1.0) < 1.0e-12);
    }
    state.T = 1.0e5;
    assert(CoolingDustRecombination(&state) == 0.0 && HeatingDustPhotoelectic(&state) == 0.0);
    state.T = 1.0e7;
    assert(CoolingDustRecombination(&state) == 0.0 && HeatingDustPhotoelectic(&state) == 0.0);
}

static void CheckH2RateReferences(void)
{
    /* Glover & Jappsen (2007), Table 3, reaction 60 at Cp=1 and solar dust. */
    const double temperatures[] = {100.0,1000.0,300.0};
    const double dustTemperatures[] = {10.0,10.0,125.0};
    const double formation[] = {1.7652006169282126e-17,7.730961606263105e-18,1.9837262541543977e-19};
    for(size_t i=0; i<sizeof formation/sizeof formation[0]; i++)
        assert(fabs(H2FormationRate(temperatures[i],SOLAR_METALLICITY,dustTemperatures[i])/formation[i]-1.0) < 1.0e-13);

    /* Draine & Bertoldi (1996), eq. 37 divided by its value at zero column. */
    const double columns[] = {5.0e14,1.0e16,1.0e19,1.0e16};
    const double widths[] = {3.0,3.0,3.0,7.1};
    const double shielding[] = {0.5675483826058047,0.024026417882459576,0.00021947801119929764,0.07384785830882051};
    const double thin = 3.3e-11*1.7;
    assert(H2SelfShieldingFactor(0.0,3.0) == 1.0);
    double slope;
    assert(fabs(CEQH2PhotoRateAndShieldingSlope(1.7,0.0,3.0,1.0,&slope)/thin-1.0) < 1.0e-14);
    assert(slope == 0.0);
    const double step = 1.0e-5, dust = 0.42;
    for(size_t i=0; i<sizeof columns/sizeof columns[0]; i++){
        const double column = columns[i], width = widths[i], factor = shielding[i];
        assert(fabs(H2SelfShieldingFactor(column,width)/factor-1.0) < 1.0e-13);
        const double rate = CEQH2PhotoRateAndShieldingSlope(1.7,column,width,dust,&slope);
        assert(fabs(rate/(thin*factor*dust)-1.0) < 1.0e-13);
        const double above = H2SelfShieldingFactor(column*exp(step),width);
        const double below = H2SelfShieldingFactor(column*exp(-step),width);
        const double numericalSlope = log(above/below)/(2.0*step);
        assert(fabs(slope-numericalSlope) < 1.0e-9);
    }

    /* Removing the zero-temperature width leaves the thermal width of H2. */
    const double width = H2DopplerWidth(1.0e4);
    const double turbulentWidth = H2DopplerWidth(0.0);
    const double thermal = sqrt(width*width-turbulentWidth*turbulentWidth);
    assert(fabs(thermal-9.083) < 0.005);
}

static void CheckHeIIIRecombinationCooling(void)
{
    /* Hydrogenic scaling: Lambda_HeIII(T) = 8 Lambda_HII(T/4) for case A and case B. */
    const bool caseB[] = {true, false};
    for(size_t c = 0; c < sizeof(caseB)/sizeof(caseB[0]); c++){
        EQCLibSetRunParamRecombinationCaseB(caseB[c]);
        const double temperatures[] = {1.0e4, 3.0e4, 1.0e5, 1.0e6};
        for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
            StructEQCData helium = {.T = temperatures[i], .nel = 1.0, .nHeIII = 1.0};
            StructEQCData hydrogen = {.T = 0.25*temperatures[i], .nel = 1.0, .nion = {[ipHYDROGEN] = 1.0}};
            const double expected = 8.0*CoolingRecombination(&hydrogen);
            /* The He II threshold is 4.0018 times the H threshold, so the scaling holds to ~1e-3. */
            assert(fabs(CoolingHeIIIRecombination(&helium)-expected) <= 2.0e-3*expected);
        }
    }
    EQCLibSetRunParamRecombinationCaseB(true);
}

static void CheckIonizedGasEquationOfState(void)
{
    const double nH = 1.0;
    StructEQCData state = {
        .T = 1.0e6,
        .Rho = nH * ChemicalElementMass[ipHYDROGEN],
        .n = {[ipHYDROGEN] = nH},
        .nel = nH,
    };
    assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-9) == EQC_SUCCESS);
    assert(state.nion[ipHYDROGEN] > 0.999 * nH);
    assert(fabs(state.DegreeOfFreedom - 1.5) < 1.0e-12);
    assert(fabs(state.Gamma - 5.0/3.0) < 1.0e-12);

    const double expectedU = 1.5 * (nH + state.nel) * BOLTZMANN_CONSTANT_CGS * state.T / state.Rho;
    assert(fabs(EQCLibConvertTtoU(&state) - expectedU) < 1.0e-12 * expectedU);
    state.U = expectedU;
    assert(fabs(EQCLibConvertUtoT(&state) - state.T) < 1.0e-12 * state.T);
}

static void CheckRejectedRateUpdates(StructEQCData *state, EQCLibStatus expected)
{
    unsigned char original[sizeof(*state)];
    memcpy(original,state,sizeof(original));
    double rate = 42.0;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(state,1.0e-8) == expected);
    assert(memcmp(original,state,sizeof(original)) == 0);
    assert(EQCLibTryNetCoolingRateChemcalEquilibrium(state,&rate) == expected);
    assert(rate == 42.0 && memcmp(original,state,sizeof(original)) == 0);
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(state,1.0e-8,&rate) == expected);
    assert(rate == 42.0 && memcmp(original,state,sizeof(original)) == 0);
    assert(isnan(EQCLibNetCoolingRateChemcalEquilibrium(state)));
    assert(memcmp(original,state,sizeof(original)) == 0);
}

static void CheckFailureStatus(void)
{
    StructEQCData state = GasState(0.001,13242.50,1.0e18);
    state.X0 = -1.0;
    CheckRejectedRateUpdates(&state,EQC_INVALID_INPUT);
}

static void CheckSoftXrayRates(void)
{
    const double column = 1.0e20;
    const double fractions[] = {1.0e-4,0.01,0.1,0.3,0.5,1.0};
    /* Reference values from the Wolfire fits before changing the xe cutoff. */
    const double ionization[] = {1.8686041205017047e-18,1.6249656997064237e-18,7.0502336080631525e-19,
                                 3.9547076063063309e-19,2.9052818252857640e-19,0.0};
    const double heating[] = {2.4604560942189459e-29,8.4246736965932879e-29,1.2521300789597018e-28,
                              1.4368194436389935e-28,1.5144896889158104e-28,1.6080426300601547e-28};
    for(size_t i=0; i<sizeof(fractions)/sizeof(fractions[0]); i++){
        const double xe = fractions[i];
        StructEQCData state = {.Nt=column, .X0=1.0, .nel=xe, .n={1.0}};
        assert(fabs(TotalSoftXrayHeatingRate(&state)/heating[i]-1.0) < 1.0e-12);
        double expected = ionization[i];
#ifndef EQC_LEGACY_XRAY_CUTOFF
        if(xe == 0.3) expected *= 0.5;
        if(xe >= 0.5) expected = 0.0;
#endif
        const double rate = SoftXrayIonizationRateHI(xe,column);
        assert(fabs(rate-expected) <= 1.0e-12*expected);
        assert(SoftXrayIonizationRateHeI(xe,column) == rate);
    }
    assert(SoftXrayIonizationRateHI(nextafter(0.5,INFINITY),column) == 0.0);
#ifndef EQC_LEGACY_XRAY_CUTOFF
    const double step = 1.0e-6;
    const double atStart = SoftXrayIonizationRateHI(0.1,column);
    const double belowStart = SoftXrayIonizationRateHI(0.1-step,column);
    const double aboveStart = SoftXrayIonizationRateHI(0.1+step,column);
    const double leftSlope = (atStart-belowStart)/step;
    const double rightSlope = (aboveStart-atStart)/step;
    assert(fabs(rightSlope-leftSlope) <= 1.0e-3*fabs(leftSlope));

    double previous = atStart;
    for(int i=1; i<=32; i++){
        const double xe = 0.1 + (0.5-0.1)*i/32;
        const double rate = SoftXrayIonizationRateHI(xe,column);
        assert(rate < previous && rate >= 0.0);
        previous = rate;
    }

    /* The last two samples approach zero quadratically at the upper join. */
    const double nearEnd = SoftXrayIonizationRateHI(0.5-step,column);
    const double nearerEnd = SoftXrayIonizationRateHI(0.5-0.5*step,column);
    assert(nearerEnd > 0.0 && fabs(nearEnd/nearerEnd-4.0) < 1.0e-3);
#endif
}

static void CheckSoftXrayChargeTransition(void)
{
    const StructEQCData initial = GasState(0.001,13242.50,1.0e18);
    StructEQCData state = initial;
#ifdef EQC_LEGACY_XRAY_CUTOFF
    CheckRejectedRateUpdates(&state,EQC_NO_EQUILIBRIUM);
#else
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == EQC_SUCCESS);
    assert(state.nel > 0.1*state.n[ipHYDROGEN] && state.nel < 0.5*state.n[ipHYDROGEN]);
    assert(fabs(EQCLibGetElectronNumberDensity(&state)-state.nel) <= 1.0e-8*state.nel);
    state = initial;
    double rate;
    assert(EQCLibTryNetCoolingRateChemcalEquilibrium(&state,&rate) == EQC_SUCCESS);
    assert(isfinite(rate));
    state = initial;
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&rate) == EQC_SUCCESS);
    assert(isfinite(rate));
#endif
}

static void CheckZeroElectronSeed(void)
{
    StructEQCData state = GasState(1.0,100.0,1.0e19);
    state.G0 = state.Z = state.Zd = 0.0;
    for(int i = 1; i < NumberChemicalElement; i++) state.n[i] = state.natom[i] = 0.0;
    EQCLibSetRunParamCosmicRayIonization(false);
    assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-8) >= 0);
    assert(fabs(state.nel - 0.00240514) < 1.0e-7);
    EQCLibSetRunParamCosmicRayIonization(true);
}

static void CheckLowIonizationXraySearch(void)
{
    StructEQCData state = {.T = 100.0, .Nt = 1.0e24, .X0 = 1.0, .nel = 1.0, .n = {1.0e3}};
    EQCLibSetRunParamCosmicRayIonization(false);
    EQCLibSetRunParamGrainSurfaceProcesses(false);
    /* A seed inside the fitted X-ray range keeps that search domain. */
    assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-8) == EQC_NO_EQUILIBRIUM);
    state.nel = 0.05;
    assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-8) == EQC_SUCCESS);
    assert(state.nel > 0.0 && state.nel < 1.0e-4*state.n[ipHYDROGEN]);
    assert(fabs(state.nel-state.nion[ipHYDROGEN]) <= 1.0e-8*state.nel);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
}

static void CheckAtomicThermalNumerics(void)
{
    const double accuracy = 1.0e-8;
    const double samples[][2] = {{870963.6,128.5617}, {2.290868,10004.58}};
    for(size_t i = 0; i < sizeof(samples)/sizeof(samples[0]); i++){
        StructEQCData state = GasState(samples[i][0],samples[i][1],1.0e19);
        state.G0 = 1700.0;
        state.X0 = 1000.0;
        state.nel = 0.9*state.n[ipHYDROGEN];
        state.dt = 1000.0*YEAR_CGS;
        assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,accuracy,&state.Du) == EQC_SUCCESS);
        const double change = state.dt*state.Du;
        const double source = state.dt*(state.H_tot-state.L_tot)/state.Rho;
        assert(fabs(change-source) <= accuracy*fmax(fabs(change),fabs(source)));
        state.U += change;
        assert(fabs(EQCLibConvertTtoU(&state)-state.U) < 1.0e-12*state.U);
    }
}

/* Isochoric cooling over a fixed interval, split into equal implicit steps. */
static double CoolIsochorically(double nH, double temperature, double duration, int steps)
{
    StructEQCData state = GasState(nH,temperature,1.0e20*sqrt(nH));
    state.nel = 1.0e-3*nH;
    assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-9) == EQC_SUCCESS);
    state.U = EQCLibConvertTtoU(&state);
    state.dt = duration/steps;
    for(int step = 0; step < steps; step++){
        assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-4,&state.Du) == EQC_SUCCESS);
        state.U += state.dt*state.Du;
        assert(isfinite(state.U) && state.U > 0.0);
    }
    return state.T;
}

static void CheckHighTemperatureCIEExtrapolation(void)
{
    const double tableMaximum = 1.0e8;
    const double maximumTemperature = 1.0e9;
    const double nH = 0.1;
    const double justAboveTable = tableMaximum*(1.0 + 1.0e-6);
    StructEQCData endpoint = GasState(nH,tableMaximum,1.0e19);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&endpoint,1.0e-8) == EQC_SUCCESS);
    const double endpointCoefficient = endpoint.L_res/(endpoint.nel*nH);
    assert(endpointCoefficient > 0.0);

    const double temperatures[] = {justAboveTable, maximumTemperature};
    for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
        StructEQCData state = endpoint;
        state.T = temperatures[i];
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == EQC_SUCCESS);
        const double expected = endpointCoefficient*sqrt(state.T/tableMaximum);
        const double coefficient = state.L_res/(state.nel*nH);
        assert(fabs(coefficient-expected) <= 1.0e-12*expected);
    }

    const double duration = 1.0e5*YEAR_CGS;
    const double cooled = CoolIsochorically(nH,maximumTemperature,duration,1);
    assert(cooled > tableMaximum && cooled < maximumTemperature);
    assert(CoolIsochorically(nH,justAboveTable,duration,1) < tableMaximum);

    StructEQCData outside = endpoint;
    outside.T = nextafter(maximumTemperature,INFINITY);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&outside,1.0e-8) == EQC_INVALID_INPUT);
}

static void CheckTraceAbundanceCIECooling(void)
{
    const double abundanceRatio = 1.0e-3;
    const double temperatures[] = {9999.99,1.0e4,1.0e5};
    double coefficients[3];
    for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
        StructEQCData state = {
            .T = temperatures[i], .Z = ASPLUND_2009_SOLAR_METALLICITY, .nel = 1.0,
            .n = {1.0,0.1,3.0e-4,5.0e-4,4.0e-6}
        };
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == EQC_SUCCESS);
        const double coefficient = state.L_res/(state.nel*state.n[ipHYDROGEN]);
        assert(coefficient > 0.0);
        coefficients[i] = coefficient;

        /* N, Ne, Mg, S, Ca and Fe enter only above 10^4 K, in proportion to Z. */
        StructEQCData untracked = state;
        untracked.n[ipCARBON] = untracked.n[ipOXYGEN] = untracked.n[ipSILICON] = 0.0;
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&untracked,1.0e-8) == EQC_SUCCESS);
        const double untrackedCoefficient = untracked.L_res/(untracked.nel*untracked.n[ipHYDROGEN]);
        if(temperatures[i] <= 1.0e4){
            assert(untrackedCoefficient == 0.0);
        }else{
            assert(untrackedCoefficient > 0.0 && untrackedCoefficient < coefficient);
            untracked.Z *= 0.5;
            assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&untracked,1.0e-8) == EQC_SUCCESS);
            const double half = untracked.L_res/(untracked.nel*untracked.n[ipHYDROGEN]);
            assert(fabs(half-0.5*untrackedCoefficient) <= 1.0e-12*untrackedCoefficient);
        }

        /* Scaling Z and the tracked metals together scales the CIE term. */
        state.Z *= abundanceRatio;
        state.n[ipCARBON] *= abundanceRatio;
        state.n[ipOXYGEN] *= abundanceRatio;
        state.n[ipSILICON] *= abundanceRatio;
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == EQC_SUCCESS);
        const double traceCoefficient = state.L_res/(state.nel*state.n[ipHYDROGEN]);
        const double expected = abundanceRatio*coefficient;
        assert(fabs(traceCoefficient-expected) <= 1.0e-12*expected);

        state.Z = 0.0;
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == EQC_SUCCESS);
        assert(state.L_res == 0.0);
    }
    assert(fabs(coefficients[0]-coefficients[1]) <= 1.0e-4*coefficients[1]);
}

static void CheckHotGasTimeStepConvergence(void)
{
    const double duration = 2.0e5*YEAR_CGS;

    /* When the step is shorter than the cooling time, one step matches a fine integration. */
    const double diffuse[][2] = {{0.1,1.0e6}, {0.1,1.0e7}};
    for(size_t i = 0; i < sizeof(diffuse)/sizeof(diffuse[0]); i++){
        const double coarse = CoolIsochorically(diffuse[i][0],diffuse[i][1],duration,1);
        const double fine = CoolIsochorically(diffuse[i][0],diffuse[i][1],duration,1000);
        assert(fabs(coarse-fine) <= 0.01*fine);
    }

    /* Exercise the transition through 10^4 K with steps spanning the hot cooling time.
       Scale the duration inversely with density so both runs cross the transition. */
    const double dense[] = {1.0, 10.0};
    for(size_t i = 0; i < sizeof(dense)/sizeof(dense[0]); i++){
        const double coolingDuration = 1.6e6*YEAR_CGS/dense[i];
        const double fine = CoolIsochorically(dense[i],1.0e6,coolingDuration,1000);
        const double error10 = fabs(CoolIsochorically(dense[i],1.0e6,coolingDuration,10)-fine);
        const double error100 = fabs(CoolIsochorically(dense[i],1.0e6,coolingDuration,100)-fine);
        assert(fine < 1.0e4);
        assert(error100 < error10);
        assert(error100 <= 0.05*fine);
    }
}

/* The shielded CHIMES comparison exposed two stalls at nH = 1 cm^-3. */
static void CheckShieldedCoolingHistories(void)
{
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(false);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamCOFormation(false,1.0);
    EQCLibSetRunParamChargeTransfer(true);
    assert(EQCLibSetRunParamChemistryAccuracy(1.0e-5) == EQC_SUCCESS);

    const double initialTemperatures[] = {1.0e4,1.0e7};
    const double secondsPerYear = 3.15576e7;
    for(size_t i = 0; i < sizeof initialTemperatures/sizeof initialTemperatures[0]; i++){
        StructEQCData state = GasState(1.0,initialTemperatures[i],1.0e20);
        state.Td = 10.0;
        state.nel = 1.0e-3;
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-9) == EQC_SUCCESS);
        state.U = EQCLibConvertTtoU(&state);
        double previousYear = 0.0;
        for(int step = 0; step <= 180; step++){
            const double year = pow(10.0,(double)step/20.0);
            state.dt = (year-previousYear)*secondsPerYear;
            const EQCLibStatus status = EQCLibTryAdvanceCell(&state,1.0e-4,EQC_CHEMISTRY_EQUILIBRIUM);
            if(status != EQC_SUCCESS)
                fprintf(stderr,"Shielded cooling from T=%g K failed at year %g (status %d)\n",
                        initialTemperatures[i],year,status);
            assert(status == EQC_SUCCESS);
            assert(fabs(EQCLibConvertTtoU(&state)-state.U) <= 1.0e-12*state.U);
            previousYear = year;
        }
        assert(previousYear == 1.0e9 && state.T < 3.0e3);
    }
    EQCLibRunParameters = saved;
}

static void CheckMolecularBranches(void)
{
    /* Keep the test near N(H2)=10^22 cm^-2 where CO has three roots.
       A high molecular fraction keeps the H2 root stable with the DB96 fit. */
    const double savedClumping = EQCLibRunParameters.CpForH2;
    EQCLibSetRunParamH2Formation(true,343.0);
    StructEQCData initial = GasState(25.0,50.0,2.5e22);
    initial.G0 = 588.5900374411234;
    initial.Zd = 0.020 * SOLAR_METALLICITY;
    initial.n[ipCARBON] = 0.30;
    initial.n[ipOXYGEN] = 0.46;
    initial.Td = EQCApproximateDustTemperature(initial.Nt,initial.G0);
    initial.nmol[ipH2] = 10.0;
    const double seeds[] = {0.0,0.05,0.20};
    const double roots[] = {0.00218869,0.05179728,0.29849496};
    for(size_t i = 0; i < sizeof(seeds)/sizeof(seeds[0]); i++){
        StructEQCData state = initial;
        state.nmol[ipCO] = seeds[i];
        assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-9) >= 0);
        assert(fabs(state.nmol[ipCO] - roots[i]) < 1.0e-6);
        double hydrogen = state.natom[ipHYDROGEN] + state.nion[ipHYDROGEN] + 2.0*state.nmol[ipH2];
        double carbon = state.natom[ipCARBON] + state.nion[ipCARBON] + state.nCIII + state.nmol[ipCO];
        double oxygen = state.natom[ipOXYGEN] + state.nion[ipOXYGEN] + state.nOIII + state.nmol[ipCO];
        assert(fabs(hydrogen - state.n[ipHYDROGEN]) < 1.0e-9 * state.n[ipHYDROGEN]);
        assert(fabs(carbon - state.n[ipCARBON]) < 1.0e-9 * state.n[ipCARBON]);
        assert(fabs(oxygen - state.n[ipOXYGEN]) < 1.0e-9 * state.n[ipOXYGEN]);
    }
    EQCLibSetRunParamH2Formation(true,savedClumping);
}

static void CheckOxygenLimit(void)
{
    StructEQCData state = GasState(1000.0,50.0,1.0e22);
    state.n[ipCARBON] = 0.30;
    state.n[ipOXYGEN] = 0.10;
    assert(EQCLibChemicalReactionEquilibrium(&state,1.0e-8) == EQC_NO_EQUILIBRIUM);
    assert(state.nmol[ipCO] == 0.0);
}

static void CheckMolecularXrayLocalSearch(void)
{
    StructEQCData state = GasState(1.0e3,200.0,1.0e21);
    state.nel = 1.0;
    const double accuracy = 1.0e-9;
    assert(EQCLibChemicalReactionEquilibrium(&state,accuracy) == EQC_SUCCESS);

    /* Reference densities from a more tightly converged solve with the same seed. */
    const double ne = 0.29759611107942;
    const double nH2 = 498.45211736023;
    const double nCO = 0.0038340845233887;
    assert(fabs(state.nel - ne) <= 1.0e-8 * ne);
    assert(fabs(state.nmol[ipH2] - nH2) <= 1.0e-8 * nH2);
    assert(fabs(state.nmol[ipCO] - nCO) <= 1.0e-8 * nCO);
    assert(fabs(EQCLibGetElectronNumberDensity(&state) - state.nel) <= accuracy * state.nel);
}

static void CheckEnergyAndEquationOfState(void)
{
    StructEQCData state = GasState(1000.0,6000.0,1.0e22);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-10) == EQC_SUCCESS);
    const double oldH2 = state.nmol[ipH2];
    const double oldU = EQCLibConvertTtoU(&state);
    state.U = oldU;
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&state.Du) == EQC_SUCCESS);
    assert(state.U == oldU && state.Du < 0.0);
    assert(state.nmol[ipH2] > oldH2);
    state.U += state.dt * state.Du;
    const double source = state.dt * (state.H_tot - state.L_tot) / state.Rho;
    assert(fabs(state.U - oldU - source) <= 2.0e-8 * fabs(source));
    assert(fabs(EQCLibConvertTtoU(&state) - state.U) <= 1.0e-12 * state.U);

    const double solvedTemperature = state.T;
    double rate;
    assert(EQCLibTryNetCoolingRateChemcalEquilibrium(&state,&rate) == EQC_SUCCESS);
    assert(fabs(state.T - solvedTemperature) < 1.0e-7 * solvedTemperature);
    assert(fabs(EQCLibConvertTtoU(&state) - state.U) <= 2.0e-8 * state.U);
    assert(isfinite(rate));
}

static void CheckThermalBranchAtCOBoundary(void)
{
    const double nH = pow(10.0,0.5);
    StructEQCData state = GasState(nH,pow(10.0,4.5),1.0e22*sqrt(nH));
    state.nel = 1.0e-3*nH;
    state.dt = 1.25e6*YEAR_CGS;
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-4,&state.Du) == EQC_SUCCESS);
    /* A coarser search can skip the hot root next to the CO regime change. */
    assert(state.T > 1.0e4 && state.T < 1.05e4);
}

static void CheckCIIHydrogenCollisionCooling(void)
{
    StructEQCData line = {.T = 250.0};
    line.nion[ipCARBON] = 0.01;
    line.nmol[ipH2] = 10.0;
    const double atJoin = CoolingCII158um(&line);
    line.T = nextafter(line.T,INFINITY);
    assert(atJoin > 0.0 && fabs(CoolingCII158um(&line)-atJoin) <= 1.0e-12*atJoin);

    const double nH = 1.0e3;
    const double accuracy = 1.0e-4;
    const double initialTemperatures[] = {pow(10.0,5.75),1.0e6,pow(10.0,6.25)};
    for(size_t i = 0; i < sizeof(initialTemperatures)/sizeof(initialTemperatures[0]); i++){
        StructEQCData state = GasState(nH,initialTemperatures[i],1.0e20*pow(nH,0.25));
        state.nel = 1.0e-3*nH;
        state.dt = 1.0e4*YEAR_CGS;
        for(int step = 0; step < 3; step++){
            assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,accuracy,&state.Du) == EQC_SUCCESS);
            const double change = state.dt*state.Du;
            const double source = state.dt*(state.H_tot-state.L_tot)/state.Rho;
            assert(fabs(change-source) <= accuracy*fmax(fabs(change),fabs(source)));
            state.U += change;
        }
    }
}

static void CheckH2CoolingFitLimits(void)
{
    const double self = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6000.0,0,0,1,0,0);
    const double neutral = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6000.0,0,1,1,1,0);
    const double charged = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(10000.0,1,0,1,0,1);
    const double aboveSelf = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6001.0,0,0,1,0,0);
    const double aboveNeutral = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(6001.0,0,1,1,1,0);
    const double aboveCharged = CoolingH2RovibrationalLineLowDensityLimitGloverAbel2008(10001.0,1,0,1,0,1);
    assert(self > 0.0 && neutral > self && charged > self);
    assert(fabs(aboveSelf-self) <= 1.0e-12*self);
    assert(fabs(aboveNeutral-neutral) <= 1.0e-12*neutral);
    assert(fabs(aboveCharged-charged) <= 1.0e-12*charged);
}

static void CheckDisabledFineStructureLines(void)
{
    double line1 = 1.0, line2 = 2.0, line3 = 3.0;
    assert(CoolingCI(1.0,100.0,0.0,0.0,1.0,0.0,&line1,&line2,&line3) == 0.0);
    assert(line1 == 0.0 && line2 == 0.0 && line3 == 0.0);
    assert(CoolingCI610um(1.0,100.0,0.0,0.0,1.0,0.0) == 0.0);
    assert(CoolingCI230um(1.0,100.0,0.0,0.0,1.0,0.0) == 0.0);
    assert(CoolingCI370um(1.0,100.0,0.0,0.0,1.0,0.0) == 0.0);

    line1 = 1.0;
    line2 = 2.0;
    line3 = 3.0;
    assert(CoolingSiI(1.0,100.0,1.0,0.0,&line1,&line2,&line3) == 0.0);
    assert(line1 == 0.0 && line2 == 0.0 && line3 == 0.0);
    assert(CoolingSiI130um(1.0,100.0,1.0,0.0) == 0.0);
    assert(CoolingSiI45um(1.0,100.0,1.0,0.0) == 0.0);
    assert(CoolingSiI68um(1.0,100.0,1.0,0.0) == 0.0);
}

static void CheckHydrogenFlux(const StructEQCData *state, double accuracy)
{
    const double nH = state->n[ipHYDROGEN];
    const double nHI = state->natom[ipHYDROGEN];
    const double nH2 = state->nmol[ipH2];
    const double dustTemperature = state->Td > 0.0 ?
                                   state->Td : EQCApproximateDustTemperature(state->Nt,state->G0);
    const double formed = nH*H2FormationRate(state->T,state->Zd,dustTemperature)*nHI;
    const double doppler = H2DopplerWidth(state->T);
    const double photo = H2PhotoDissociationRate(state->G0,state->Nt,state->Nt*nH2/nH,doppler,state->Zd);
    const double destroyed = nH2*(photo+H2CollisionalDissociationRate(state->T,nH,nHI,nH2,state->nel));
    assert(formed > 0.0 && destroyed > 0.0);
    assert(fabs(formed-destroyed) <= accuracy*fmax(formed,destroyed));
}

static void CheckMolecularFlux(double dustTemperature)
{
    StructEQCData state = GasState(1000.0,75.0,1.0e20);
    state.Td = dustTemperature;
    const double accuracy = 1.0e-9;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == EQC_SUCCESS);
    CheckHydrogenFlux(&state,accuracy);
    assert(state.H_H2 > 0.0 && isfinite(state.H_H2));

    state.nmol[ipCO] = 1.0e-11;
    assert(CoolingCORovibrationalLine(&state) > 0.0);
}

static void CheckMolecularHydrogenBranches(void)
{
    /* The DB96 shielding curve allows two stable H2 equilibria in this pure
       hydrogen fixture. The initial abundance selects the nearby balance. */
    const double nH = 3.0;
    const double seeds[] = {0.0,5.0e-4};
    const double roots[] = {3.64371643e-6,5.83906810e-4};
    const double accuracies[] = {1.0e-9,16.0*DBL_EPSILON};
    for(size_t i=0; i<sizeof seeds/sizeof seeds[0]; i++){
        const StructEQCData initial = {
            .T = 300.0, .Td = 10.0, .Nt = 1.0e20*sqrt(nH), .Zd = SOLAR_METALLICITY, .G0 = 1.7,
            .n = {[ipHYDROGEN] = nH}, .natom = {[ipHYDROGEN] = nH*(1.0-seeds[i])},
            .nmol = {[ipH2] = 0.5*nH*seeds[i]},
        };
        /* Both the Newton solve and the tight fallback start from the same seed. */
        for(size_t j=0; j<sizeof accuracies/sizeof accuracies[0]; j++){
            StructEQCData state = initial;
            assert(CEQSolveHydrogenEquilibrium(&state,accuracies[j]) == EQC_SUCCESS);
            assert(fabs(2.0*state.nmol[ipH2]/(nH*roots[i])-1.0) < 1.0e-6);
            CheckHydrogenFlux(&state,accuracies[j]);
        }
    }
}

static void CheckTraceHydrogen(void)
{
    const double accuracy = 1.0e-12;
    StructEQCData initial = GasState(1000.0,50.0,1.0e20);
    initial.X0 = 0.0;
    StructEQCData atomic = initial;
    atomic.G0 = 1.0e100;
    assert(CEQSolveHydrogenEquilibrium(&atomic,accuracy) == EQC_SUCCESS);
    assert(atomic.nmol[ipH2] < 1.0e-90);
    CheckHydrogenFlux(&atomic,accuracy);

    StructEQCData molecular = initial;
    molecular.G0 = 1.0e-100;
    assert(CEQSolveHydrogenEquilibrium(&molecular,accuracy) == EQC_SUCCESS);
    assert(molecular.natom[ipHYDROGEN] < 1.0e-90);
    CheckHydrogenFlux(&molecular,accuracy);
}

static void CheckOxygenLineEmissivities(void)
{
    const double tolerance = 2.0e-4;
    /* Independent five-level solve using PyNeb's Z82-WFD96/Kal09 and FFT04-SZ00/SSB14 data. */
    const struct {
        double ne;
        double oxygenII;
        double oxygenIII;
        double opticalDoublet;
        double infraredDoublet;
    } reference[] = {
        {1.0,   3.5638632548e-21, 7.0532399529e-21, 1.4991843, 0.5530553},
        {100.0, 3.4963065702e-21, 6.7333964675e-21, 1.3493476, 0.7989521},
        {1.0e4, 1.7171530434e-21, 4.9690089747e-21, 0.4172276, 7.6530238},
    };
    double oII[OxygenLevelCount][OxygenLevelCount], oIII[OxygenLevelCount][OxygenLevelCount];
    for(size_t i = 0; i < sizeof(reference)/sizeof(reference[0]); i++){
        const double ne = reference[i].ne;
        const double coolingII = CoolingLinesOII(1.0e4,ne,1.0,oII)/ne;
        const double coolingIII = CoolingLinesOIII(1.0e4,ne,1.0,oIII)/ne;
        assert(fabs(coolingII/reference[i].oxygenII-1.0) < tolerance);
        assert(fabs(coolingIII/reference[i].oxygenIII-1.0) < tolerance);
        assert(fabs((oII[1][0]/oII[2][0])/reference[i].opticalDoublet-1.0) < tolerance);    /* 3729/3726 Å */
        assert(fabs((oIII[2][1]/oIII[1][0])/reference[i].infraredDoublet-1.0) < tolerance); /* 52/88 µm */
        assert(fabs(oIII[3][2]/oIII[3][1]-2.983969) < tolerance);                         /* 5007/4959 Å */
    }

    const double lowT = pow(10.0,3.7);
    const double highT = pow(10.0,4.3);
    assert(fabs(CoolingLinesOII(lowT,100.0,1.0,NULL)/9.9799349321e-21-1.0) < tolerance);
    assert(fabs(CoolingLinesOIII(highT,100.0,1.0,NULL)/1.8028688374e-18-1.0) < tolerance);
    assert(CoolingLinesOII(1.0e4,0.0,1.0,oII) == 0.0 && oII[1][0] == 0.0);
    assert(CoolingLinesOIII(1.0e4,100.0,0.0,oIII) == 0.0 && oIII[3][2] == 0.0);
    assert(isnan(CoolingLinesOII(1.0e4,-1.0,1.0,NULL)));
}

static void CheckPhotoionizedOxygenCooling(void)
{
    StructEQCData state = GasState(100.0,1.1e4,1.0e19);
    state.NLyn = 3.0e9;
    state.X0 = state.G0 = state.Zd = 0.0;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,1.0e-8) == EQC_SUCCESS);
    assert(state.L_OII > 0.0 && state.L_OIII > 0.0);

    StructEQCData cie = state;
    EQCLibSetRunParamHIIRegion(false);
    assert(EQCLibCalculateCoolingHeatingRates(&cie) == EQC_SUCCESS);
    EQCLibSetRunParamHIIRegion(true);
    assert(cie.L_OII == 0.0 && cie.L_OIII == 0.0 && cie.L_res > 0.0);

    double lines = CoolingLinesHIIregionOxygen(state.T,state.nel,state.nion[ipOXYGEN],state.nOIII);
    const double weight = (state.L_OII+state.L_OIII)/lines;
    assert(weight > 0.99 && weight <= 1.0);
    assert(fabs(state.L_res-(1.0-weight)*cie.L_res) < 1.0e-12*cie.L_res);
    const double expected = cie.L_tot - cie.L_res + state.L_res + state.L_OII + state.L_OIII;
    assert(fabs(state.L_tot-expected) < 1.0e-12*expected);

    StructEQCData transition = state;
    transition.T = 2.0e4;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&transition,1.0e-8) == EQC_SUCCESS);
    lines = CoolingLinesHIIregionOxygen(transition.T,transition.nel,transition.nion[ipOXYGEN],transition.nOIII);
    const double transitionWeight = (transition.L_OII+transition.L_OIII)/lines;
    assert(transitionWeight > 0.0 && transitionWeight < weight);

    StructEQCData hot = state;
    hot.T = 1.0e6;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&hot,1.0e-8) == EQC_SUCCESS);
    assert(hot.L_OII == 0.0 && hot.L_OIII == 0.0 && hot.L_res > 0.0);

    StructEQCData unilluminated = state;
    unilluminated.NLyn = 0.0;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&unilluminated,1.0e-8) == EQC_SUCCESS);
    assert(unilluminated.L_OII == 0.0 && unilluminated.L_OIII == 0.0 && unilluminated.L_res > 0.0);

    state.U = EQCLibConvertTtoU(&state);
    state.dt = 1000.0*YEAR_CGS;
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,1.0e-8,&state.Du) == EQC_SUCCESS);
    assert(state.Du < 0.0 && state.T < 1.0e4 && state.L_OII > 0.0 && state.L_OIII > 0.0);
    assert(state.nCIII > 0.0); /* Free carbon stays photoionized below the CO temperature threshold. */
}

static void CheckPhotoionizationTurnOn(void)
{
    EQCLibSetRunParamH2Formation(false,1.0);
    EQCLibSetRunParamCOFormation(false,1.0);
    const double nH = 100.0;
    StructEQCData initial = GasState(nH,1.0e4,1.0e19);
    initial.G0 = initial.X0 = initial.Zd = 0.0;
    initial.nel = nH;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&initial,1.0e-10) == EQC_SUCCESS);
    assert(initial.nel < 0.1*nH);
    initial.U = EQCLibConvertTtoU(&initial);
    initial.NLyn = 1.0e-3*nH*LIGHT_C_CGS;

    StructEQCData illuminated = initial, reference = initial;
    reference.nel = nH;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&reference,1.0e-12) == EQC_SUCCESS);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&illuminated,1.0e-8) == EQC_SUCCESS);
    assert(illuminated.nion[ipHYDROGEN] > 0.99*nH);
    assert(fabs(EQCLibGetElectronNumberDensity(&illuminated)-illuminated.nel) <= 1.0e-8*illuminated.nel);
    assert(fabs(illuminated.nel-reference.nel) <= 1.0e-8*reference.nel);
    assert(fabs(illuminated.nOIII-reference.nOIII) <= 1.0e-8*reference.n[ipOXYGEN]);

    initial.dt = 1.0e4*YEAR_CGS;
    StructEQCData updated = initial, tightUpdate = initial;
    tightUpdate.nel = reference.nel;
    const double accuracy = 1.0e-8;
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&updated,accuracy,&updated.Du) == EQC_SUCCESS);
    assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&tightUpdate,1.0e-9,&tightUpdate.Du) == EQC_SUCCESS);
    const double change = updated.dt*updated.Du;
    const double source = updated.dt*(updated.H_tot-updated.L_tot)/updated.Rho;
    assert(fabs(change-source) <= accuracy*fmax(fabs(change),fabs(source)));
    assert(fabs(updated.T-tightUpdate.T) <= 1.0e-8*tightUpdate.T);
}

static double RelativeDifference(double a, double b)
{
    const double scale = fmax(fabs(a),fabs(b));
    return scale > 0.0 ? fabs(a-b)/scale : 0.0;
}

static void CheckHydrogenOxygenChargeTransfer(void)
{
    EQCLibSetRunParamH2Formation(false,1.0);
    EQCLibSetRunParamCOFormation(false,1.0);
    const double nH = 100.0, accuracy = 1.0e-10;
    StructEQCData state = GasState(nH,8000.0,1.0e19);
    state.G0 = state.X0 = state.Zd = 0.0;
    state.nel = nH;
    state.NLyn = 1.0e-4*nH*LIGHT_C_CGS;   /* log U = -4 leaves a few per cent H I. */
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == EQC_SUCCESS);
    assert(fabs(EQCLibGetElectronNumberDensity(&state)-state.nel) <= 1.0e-9*state.nel);

    const double T = state.T, ne = state.nel;
    const double nHI = state.natom[ipHYDROGEN], nHII = state.nion[ipHYDROGEN];
    const double nOI = state.natom[ipOXYGEN], nOII = state.nion[ipOXYGEN], nOIII = state.nOIII;
    assert(nHI > 0.01*nH && nHI < 0.1*nH);
    const double toOxygenII = nOI*nHII*OIChargeTransferRateToHII(T);
    const double toOxygenI = nOII*nHI*OIIChargeTransferRateToHI(T);
    const double fromOxygenIII = nOIII*nHI*OIIIChargeTransferRateToHI(T);

    /* Each ladder balances with the charge-transfer fluxes included. */
    const double hydrogenIonization = nHI*(HIPhotoIonizationRate(state.NLyn) + ne*HICollisionalIonizationRate(T))
                                    + toOxygenI + fromOxygenIII;
    const double hydrogenRecombination = nHII*ne*HIIRecombinationRate(T) + toOxygenII;
    assert(RelativeDifference(hydrogenIonization,hydrogenRecombination) <= 1.0e-8);
    const double oxygenIIIonization = nOI*(OIPhotoionizationRate(state.NLyn) + ne*OICollisionalIonizationRate(T))
                                    + toOxygenII;
    const double oxygenIIRecombination = nOII*ne*OIIRecombinationRate(T) + toOxygenI;
    assert(RelativeDifference(oxygenIIIonization,oxygenIIRecombination) <= 1.0e-8);
    const double oxygenIIIIonization = nOII*(OIIPhotoionizationRate(state.NLyn) + ne*OIICollisionalIonizationRate(T));
    const double oxygenIIIRecombination = nOIII*ne*OIIIRecombinationRate(T) + fromOxygenIII;
    assert(RelativeDifference(oxygenIIIIonization,oxygenIIIRecombination) <= 1.0e-8);

    /* At this H I fraction, O++ + H dominates the O++ recombination. */
    const double withoutTransfer = (OIIPhotoionizationRate(state.NLyn) + ne*OIICollisionalIonizationRate(T))
                                 / (ne*OIIIRecombinationRate(T));
    assert(nOIII/nOII < 0.25*withoutTransfer);

    /* Dust, cosmic rays and X-rays use the same coupled solve. */
    StructEQCData dusty = GasState(nH,8000.0,1.0e19);
    dusty.nel = nH;
    dusty.NLyn = state.NLyn;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&dusty,accuracy) == EQC_SUCCESS);
    assert(fabs(EQCLibGetElectronNumberDensity(&dusty)-dusty.nel) <= 1.0e-9*dusty.nel);
}

extern double b_turb;

/* H2 destruction against formation, from the rate functions used by the solver. */
static double MolecularHydrogenImbalance(const StructEQCData *s)
{
    const double nH = s->n[ipHYDROGEN], nHI = s->natom[ipHYDROGEN], nH2 = s->nmol[ipH2];
    const double formation = nH*H2FormationRate(s->T,s->Zd,EffectiveDustTemperature(s))*nHI;
    const double doppler = H2DopplerWidth(s->T);
    const double photo = H2PhotoDissociationRateWithDustFactor(s->G0,nH2*s->Nt/nH,doppler,
                                                               DustShieldingFactorFUV(s->Nt,s->Zd));
    const double collision = H2DissociationFromCoefficients(H2DissociationCoefficients(s->T),s->nel,nHI,nH2);
    return RelativeDifference((photo+collision)*nH2,formation);
}

static void CheckMolecularChargeTransferBalance(void)
{
    /* Partly ionized molecular gas, where coupling O and H used to leave H2
       about 9% out of balance because H I moved after the H2 solve. */
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamCOFormation(true,1.0);
    const double nH = 100.0, accuracy = 1.0e-8;
    StructEQCData state = GasState(nH,60.0,1.0e21);
    state.nel = 1.0e-3*nH;
    state.NLyn = 1.0e-5*nH*LIGHT_C_CGS;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == EQC_SUCCESS);
    assert(state.nmol[ipH2] > 0.4*nH && state.nion[ipHYDROGEN] > 1.0e-3*nH);
    assert(MolecularHydrogenImbalance(&state) <= 2.0*accuracy);

    const CEQHydrogenRates rates = CEQHydrogenIonizationRates(&state);
    const double T = state.T, nHI = state.natom[ipHYDROGEN], nHII = state.nion[ipHYDROGEN];
    const double ionization = nHI*rates.ionization
                            + state.nion[ipOXYGEN]*nHI*OIIChargeTransferRateToHI(T)
                            + state.nOIII*nHI*OIIIChargeTransferRateToHI(T);
    const double recombination = nHII*rates.recombination + state.natom[ipOXYGEN]*nHII*OIChargeTransferRateToHII(T);
    assert(RelativeDifference(ionization,recombination) <= 2.0*accuracy);

    /* CO is re-solved after charge transfer moves H2: solving CO again from the
       returned H2 must leave it unchanged. The old order left a 1e-3 offset. */
    const double nCO = 316.22776601683793;
    StructEQCData molecular = GasState(nCO,1000.0,1.0e20*sqrt(nCO));
    molecular.nel = 1.0e-3*nCO;
    molecular.NLyn = 1.0e-5*nCO*LIGHT_C_CGS;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&molecular,accuracy) == EQC_SUCCESS);
    assert(molecular.nmol[ipCO] > 0.01*molecular.n[ipCARBON]);
    assert(MolecularHydrogenImbalance(&molecular) <= 2.0*accuracy);
    StructEQCData carbon = molecular;
    assert(CEQSolveCarbonEquilibrium(&carbon,molecular.nmol[ipCO],accuracy) == EQC_SUCCESS);
    assert(RelativeDifference(carbon.nmol[ipCO],molecular.nmol[ipCO]) <= 1.0e-6);
    EQCLibRunParameters = saved;
}

static void CheckH2ShieldingUpperTail(void)
{
    const double edge = 1.016e22, delta = 1.e-6;
    const double atEdge = COShieldingFactor(edge,0.0,3.0);
    const double below = COShieldingFactor(edge/(1.0+delta),0.0,3.0);
    const double above = COShieldingFactor(edge*(1.0+delta),0.0,3.0);
    assert(RelativeDifference(atEdge,4.666e-3) <= 1.e-12);
    assert(below > atEdge && atEdge > above);
    assert(RelativeDifference(below-atEdge,atEdge-above) < 1.e-5);

    /* Match the attenuation per decade in the N_CO=0 column of Visser et al. Table 5. */
    const double nextDecade = COShieldingFactor(10.0*edge,0.0,3.0);
    assert(RelativeDifference(nextDecade/atEdge,3.938e-7/1.108e-2) <= 1.e-12);
    assert(COShieldingFactor(100.0*edge,0.0,3.0) < nextDecade);
}

static void CheckCarbonBalanceAtH2ShieldingBoundary(double nH, double metallicity, double column)
{
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamChargeTransfer(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamCOFormation(true,1.0);
    /* Both fixtures put the uncoupled H2 above the edge and the coupled H2 below. */
    const double upperH2Column = 1.016e22;
    const double accuracies[] = {1.e-8, 1.e-5};
    for(size_t i=0; i<sizeof(accuracies)/sizeof(accuracies[0]); i++){
        const double accuracy = accuracies[i];
        StructEQCData state = GasState(nH,30.0,column);
        state.nel = 1.e-3*nH;
        state.NLyn = 1.e-5*nH*LIGHT_C_CGS;
        state.Z *= metallicity;
        state.Zd *= metallicity;
        for(int element=ipCARBON; element<NumberChemicalElement; element++){
            state.n[element] *= metallicity;
            state.natom[element] = state.n[element];
        }
        assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,accuracy) == EQC_SUCCESS);

        StructEQCData uncoupled = state;
        assert(CEQSolveHydrogenEquilibrium(&uncoupled,accuracy) == EQC_SUCCESS);
        const double h2Column = state.Nt*(state.nmol[ipH2]/nH);
        const double uncoupledColumn = state.Nt*(uncoupled.nmol[ipH2]/nH);
        assert(h2Column < upperH2Column && uncoupledColumn > upperH2Column);

        const double nCO = state.nmol[ipCO];
        const double freeCarbon = state.natom[ipCARBON] + state.nion[ipCARBON] + state.nCIII;
        const double association = 4.e-16*pow(state.T/300.0,-0.2)*state.nmol[ipH2];
        const double photoCI = 3.1e-10*state.G0*DustShieldingFactor(state.Nt,state.Zd,3.0);
        const double formation = (association*photoCI/(association+photoCI))*freeCarbon;
        const double destruction = COPhotoDissociationRate(state.G0,state.Nt,state.Zd)*
                                   COShieldingFactor(h2Column,state.Nt*(nCO/nH),b_turb)*nCO;
        assert(RelativeDifference(formation,destruction) <= accuracy);
        assert(MolecularHydrogenImbalance(&state) <= accuracy);
        assert(RelativeDifference(state.nel,EQCLibGetElectronNumberDensity(&state)) <= accuracy);
    }
    EQCLibRunParameters = saved;
}

static void CheckMolecularHydrogenRoundoffAccuracy(void)
{
    /* At the smallest tolerance the thermal solver uses (16 eps), the local H2
       search started near the root must still meet it when it reports success. */
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamH2Formation(true,1.0);
    const double accuracy = 16.0*DBL_EPSILON;
    const double temperatures[] = {30.0, 100.0, 1000.0};
    const double offsets[] = {1.0e-9, 1.0e-6, 1.0e-3, -1.0e-3, 0.1};
    for(size_t i = 0; i < sizeof(temperatures)/sizeof(temperatures[0]); i++){
        StructEQCData reference = GasState(100.0,temperatures[i],1.0e21);
        reference.nel = 1.0e-4*reference.n[ipHYDROGEN];
        assert(CEQSolveHydrogenEquilibrium(&reference,1.0e-8) == EQC_SUCCESS);
        for(size_t j = 0; j < sizeof(offsets)/sizeof(offsets[0]); j++){
            StructEQCData state = reference;
            state.nmol[ipH2] = reference.nmol[ipH2]*(1.0+offsets[j]);
            if(CEQSolveHydrogenEquilibrium(&state,accuracy) == EQC_SUCCESS)
                assert(MolecularHydrogenImbalance(&state) <= 1.0e-12);
        }
    }
    EQCLibRunParameters = saved;
}

static void CheckChargeTransferSwitch(void)
{
    const double nH = 100.0, accuracy = 1.0e-10;
    StructEQCData coupled = GasState(nH,8000.0,1.0e19);
    coupled.G0 = coupled.X0 = coupled.Zd = 0.0;
    coupled.nel = nH;
    coupled.NLyn = 1.0e-4*nH*LIGHT_C_CGS;
    StructEQCData uncoupled = coupled;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&coupled,accuracy) == EQC_SUCCESS);

    EQCLibSetRunParamChargeTransfer(false);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&uncoupled,accuracy) == EQC_SUCCESS);
    EQCLibSetRunParamChargeTransfer(true);

    /* Without charge transfer, the O ladder balances on photo- and electron processes alone. */
    const double T = uncoupled.T, ne = uncoupled.nel;
    const double ionization = uncoupled.nion[ipOXYGEN]*(OIIPhotoionizationRate(uncoupled.NLyn)
                                                     + ne*OIICollisionalIonizationRate(T));
    const double recombination = uncoupled.nOIII*ne*OIIIRecombinationRate(T);
    assert(RelativeDifference(ionization,recombination) <= 1.0e-8);
    assert(uncoupled.nOIII/uncoupled.nion[ipOXYGEN] > 4.0*coupled.nOIII/coupled.nion[ipOXYGEN]);
    assert(fabs(EQCLibGetElectronNumberDensity(&uncoupled)-uncoupled.nel) <= 1.0e-9*uncoupled.nel);
}

static void CheckChemistryAccuracySetting(void)
{
    const double invalid[] = {0.0, -1.0e-5, 1.0, NAN, INFINITY};
    for(size_t i = 0; i < sizeof(invalid)/sizeof(invalid[0]); i++)
        assert(EQCLibSetRunParamChemistryAccuracy(invalid[i]) == EQC_INVALID_INPUT);
    assert(EQCLibRunParameters.ChemistryAccuracy == 1.0e-8);

    /* A looser chemistry tolerance keeps energy updates within the energy accuracy. */
    const double nH = 100.0, accuracy = 1.0e-4;
    const double fluxes[] = {0.0, 1.0e-3*nH*LIGHT_C_CGS};
    for(size_t i = 0; i < sizeof(fluxes)/sizeof(fluxes[0]); i++){
        StructEQCData state = GasState(nH,1.0e4,1.0e19);
        state.nel = 0.5*nH;
        state.NLyn = fluxes[i];
        state.dt = 1.0e4*YEAR_CGS;
        StructEQCData loose = state;
        assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&state,accuracy,&state.Du) == EQC_SUCCESS);
        assert(EQCLibSetRunParamChemistryAccuracy(1.0e-5) == EQC_SUCCESS);
        assert(EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&loose,accuracy,&loose.Du) == EQC_SUCCESS);
        assert(EQCLibSetRunParamChemistryAccuracy(1.0e-8) == EQC_SUCCESS);
        assert(fabs(loose.T-state.T) <= accuracy*state.T);
    }
}

static void CheckIonizingFluxNeedsHIIRegion(void)
{
    /* With H II region mode off, NLyn must not reach the chemistry or heating,
       even after the PEGASE2 tables are loaded. Otherwise C would lose its FUV
       photoionization and drop from C II to C I in diffuse gas. */
    const struct StructEQCLibParameters_tag saved = EQCLibRunParameters;
    EQCLibSetRunParamHIIRegion(false);
    EQCLibSetRunParamH2Formation(false,1.0);
    EQCLibSetRunParamCOFormation(false,1.0);
    StructEQCData dark = GasState(1.0,100.0,1.0e19);
    dark.nel = 1.0e-3;
    StructEQCData lit = dark;
    lit.NLyn = 1.0e-3*LIGHT_C_CGS;
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&dark,1.0e-10) == EQC_SUCCESS);
    assert(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&lit,1.0e-10) == EQC_SUCCESS);
    assert(lit.nel == dark.nel && lit.H_tot == dark.H_tot && lit.L_tot == dark.L_tot);
    for(int i = 0; i < NumberChemicalElement; i++)
        assert(lit.natom[i] == dark.natom[i] && lit.nion[i] == dark.nion[i]);
    assert(dark.nion[ipCARBON] > 0.99*dark.n[ipCARBON]);
    EQCLibRunParameters = saved;
}

static void CheckPhotoionization(const char *pegaseDirectory)
{
    const double unattenuated = CIFUVPhotoionizationRate(1.7,0.0);
    const double attenuated = CIFUVPhotoionizationRate(1.7,1.0);
    assert(fabs(attenuated/unattenuated-exp(-3.0)) < 1.0e-15);
    EQCLibSetRunParamDustAttenuation(false);
    assert(CIFUVPhotoionizationRate(1.7,1.0) == unattenuated);
    EQCLibSetRunParamDustAttenuation(true);

    EQCLibSetRunParamHIIRegion(true);
    assert(EQCLibTryInitialization(0,pegaseDirectory) == EQC_SUCCESS);
    const double rate = HIPhotoIonizationRate(1.0);
    const double referenceRate = 0x1.36bcd33264bbp-58;  /* Solar PEGASE2 spectrum. */
    assert(fabs(rate-referenceRate) <= 1.0e-12*referenceRate);
    /* Solar PEGASE2 coefficients normalized by the H-ionizing flux before age averaging. */
    const double referenceHeIIRate = 5.928733248e-21;
    const double referenceHeIIHeating = 1.043341218e-31;
    assert(fabs(HeIIPhotoIonizationRate(1.0)-referenceHeIIRate) <= 1.0e-8*referenceHeIIRate);
    assert(fabs(HeatingHeIIPhotoionization(1.0)-referenceHeIIHeating) <= 1.0e-8*referenceHeIIHeating);

    struct stat info;
    assert(stat("out",&info) != 0 && stat("EQCooling",&info) != 0);
    fputs("[expected test diagnostics] Missing and malformed PEGASE.2 input files.\n",stderr);
    assert(EQCLibTryInitialization(0,"missing-pegase") == EQC_INITIALIZATION_FAILED);
    assert(mkdir("malformed",0700) == 0);
    FILE *input = fopen("malformed/spectra_Kroupa_Z0.0001.dat","w");
    assert(input != NULL);
    for(int line = 0; line < 27; line++) fputc('\n',input);
    fputs("1 1500 0\n",input);
    assert(fclose(input) == 0);
    assert(EQCLibTryInitialization(0,"malformed") == EQC_INITIALIZATION_FAILED);
    assert(HIPhotoIonizationRate(1.0) == rate);
    CheckPhotoionizedOxygenCooling();
    CheckPhotoionizationTurnOn();
    CheckHydrogenOxygenChargeTransfer();
    CheckMolecularChargeTransferBalance();
    CheckH2ShieldingUpperTail();
    CheckCarbonBalanceAtH2ShieldingBoundary(100.0,0.3,2.0320007053e22);
    CheckCarbonBalanceAtH2ShieldingBoundary(1.0,0.07,2.58616e22);
    CheckMolecularHydrogenRoundoffAccuracy();
    CheckChargeTransferSwitch();
    CheckChemistryAccuracySetting();
    CheckIonizingFluxNeedsHIIRegion();
}

int main(int argc, char **argv)
{
    assert(argc == 2);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(false);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamH2Formation(false,1.0);
    EQCLibSetRunParamCOFormation(false,1.0);
    EQCLibInitialization(0);
    CheckDisabledFineStructureLines();

    StructEQCData probe = GasState(1.0,100.0,1.e20);
    probe.nel = 0.01*probe.n[ipHYDROGEN];
    assert(TotalCosmicRayHeatingRate(&probe) == 0.0);
    EQCLibSetRunParamCosmicRayIonization(true);
    const double cosmicRayHeating = TotalCosmicRayHeatingRate(&probe);
    assert(isfinite(cosmicRayHeating) && cosmicRayHeating > 0.0);
    CheckAtomicEquilibrium();
    CheckDustThermalRateLimits();
    CheckH2RateReferences();
    CheckHeIIIRecombinationCooling();
    CheckIonizedGasEquationOfState();
    CheckFailureStatus();
    CheckSoftXrayRates();
    CheckSoftXrayChargeTransition();
    CheckZeroElectronSeed();
    CheckLowIonizationXraySearch();
    CheckAtomicThermalNumerics();
    CheckHighTemperatureCIEExtrapolation();
    CheckTraceAbundanceCIECooling();
    CheckHotGasTimeStepConvergence();

    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamCOFormation(true,1.0);
    EQCLibSetRunParamGrainSurfaceProcesses(false);
    EQCLibSetRunParamSoftXrayIonization(false);
    CheckMolecularBranches();
    CheckOxygenLimit();

    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    CheckMolecularXrayLocalSearch();
    CheckEnergyAndEquationOfState();
    CheckThermalBranchAtCOBoundary();
    CheckCIIHydrogenCollisionCooling();
    CheckH2CoolingFitLimits();
    CheckMolecularFlux(0.0);
    CheckMolecularFlux(100.0);
    CheckMolecularHydrogenBranches();
    CheckTraceHydrogen();
    CheckOxygenLineEmissivities();
    CheckPhotoionization(argv[1]);
    CheckShieldedCoolingHistories();
    puts("Chemical equilibrium and thermal update checks passed.");
    return 0;
}
