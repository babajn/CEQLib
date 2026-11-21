#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CosmicRay.h"
#include "Xray.h"

static const double Column = 1.e19;
static const double ColdElectronFraction = 2.36e-4;
static const double WarmElectronFraction = 2.18e-2;

typedef struct {
    double xrayIonization;
    double xrayHeating;
    double cosmicRayHeating;
} RadiationRates;

static RadiationRates coldScale, warmScale;
static bool useWolfireRadiation;

/* These overrides are linked only into this benchmark. Normalize the native
   rate shapes at the two published points in Section 4.3 and Figure 10b.
   For EUV/X-rays use the same ionization normalization for H and He. */
double __real_SoftXrayIonizationRateHI(double xe, double column);
double __real_SoftXrayIonizationRateHeI(double xe, double column);
double __real_CosmicRayIonizationRateHI(double xe);
double __real_TotalSoftXrayHeatingRate(StructEQCData *state);
double __real_TotalCosmicRayHeatingRate(StructEQCData *state);

static RadiationRates ComputeRadiationScale(double xe, RadiationRates target)
{
    StructEQCData state = {.n = {1.0}, .nel = xe, .X0 = 1.0, .Nt = Column};
    return (RadiationRates){
        target.xrayIonization / __real_SoftXrayIonizationRateHI(xe, Column),
        target.xrayHeating / __real_TotalSoftXrayHeatingRate(&state),
        target.cosmicRayHeating / __real_TotalCosmicRayHeatingRate(&state)
    };
}

static double InterpolateScale(double xe, double cold, double warm)
{
    if(!useWolfireRadiation) return 1.0;
    if(xe <= ColdElectronFraction) return cold;
    if(xe >= WarmElectronFraction) return warm;
    double fraction = log(xe / ColdElectronFraction) / log(WarmElectronFraction / ColdElectronFraction);
    return cold * pow(warm / cold, fraction);
}

double __wrap_SoftXrayIonizationRateHI(double xe, double column)
{
    return __real_SoftXrayIonizationRateHI(xe, column)
           * InterpolateScale(xe, coldScale.xrayIonization, warmScale.xrayIonization);
}

double __wrap_SoftXrayIonizationRateHeI(double xe, double column)
{
    return __real_SoftXrayIonizationRateHeI(xe, column)
           * InterpolateScale(xe, coldScale.xrayIonization, warmScale.xrayIonization);
}

double __wrap_CosmicRayIonizationRateHI(double xe)
{
    double primary = __real_CosmicRayIonizationRateHI(xe);
    if(!useWolfireRadiation) return primary;
    /* H secondary yield from Wolfire et al. (1995), Eq. A1, for 35 eV. */
    xe = fmin(fmax(xe, 0.0), 1.0);
    double secondary = 0.3908 * (35.0 / 13.6 - 1.0)
                       * pow(1.0 - pow(xe, 0.4092), 1.7592) * pow(0.035, 2.313 * xe);
    return primary * (1.0 + secondary);
}

double __wrap_TotalSoftXrayHeatingRate(StructEQCData *state)
{
    double xe = state->nel / state->n[ipHYDROGEN];
    return __real_TotalSoftXrayHeatingRate(state)
           * InterpolateScale(xe, coldScale.xrayHeating, warmScale.xrayHeating);
}

double __wrap_TotalCosmicRayHeatingRate(StructEQCData *state)
{
    double xe = state->nel / state->n[ipHYDROGEN];
    return __real_TotalCosmicRayHeatingRate(state)
           * InterpolateScale(xe, coldScale.cosmicRayHeating, warmScale.cosmicRayHeating);
}

static double EvaluateBalance(const StructEQCData *cell, double temperature, StructEQCData *trial)
{
    *trial = *cell;
    trial->T = temperature;
    EQCLibStatus status = EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(trial, 1.e-8);
    if(status != EQC_SUCCESS || !isfinite(trial->H_tot) || !isfinite(trial->L_tot)
       || fmax(trial->H_tot, trial->L_tot) <= 0.0){
        fprintf(stderr, "Thermochemical evaluation failed at nH=%.8g, T=%.8g (status %d)\n",
                cell->n[ipHYDROGEN], temperature, status);
        exit(EXIT_FAILURE);
    }
    return (trial->H_tot - trial->L_tot) / fmax(trial->H_tot, trial->L_tot);
}

static double FindBalance(const StructEQCData *cell, StructEQCData *equilibrium)
{
    /* At fixed density first verify that the sampled temperature range has
       one heating-to-cooling crossing. The isobaric branches are classified
       later from the slope of the pressure-density curve. */
    const int intervals = 96;
    const double minimum = 5.0, maximum = 1.2e4;
    double previousT = minimum, lower = 0.0, upper = 0.0;
    double previous = EvaluateBalance(cell, previousT, equilibrium);
    bool heatedAtMinimum = previous > 0.0;
    int crossings = 0;
    for(int i = 1; i <= intervals; i++){
        double temperature = minimum * pow(maximum / minimum, (double)i / intervals);
        double residual = EvaluateBalance(cell, temperature, equilibrium);
        if((residual > 0.0) != (previous > 0.0)){
            lower = previousT;
            upper = temperature;
            crossings++;
        }
        previous = residual;
        previousT = temperature;
    }
    if(heatedAtMinimum && previous < 0.0 && crossings == 1){
        for(int i = 0; i < 80; i++){
            double temperature = sqrt(lower * upper);
            double residual = EvaluateBalance(cell, temperature, equilibrium);
            if(fabs(residual) < 1.e-9) return residual;
            if(residual > 0.0) lower = temperature;
            else upper = temperature;
        }
    }
    fprintf(stderr, "Expected one converged thermal root at nH=%.8g (found %d brackets)\n",
            cell->n[ipHYDROGEN], crossings);
    exit(EXIT_FAILURE);
}

static void WriteCurve(const char *model)
{
    for(int i = 0; i <= 1000; i++){
        double logNH = -2.0 + 0.005 * i;
        double nH = pow(10.0, logNH);
        /* Wolfire et al. (2003), Section 3.3, following the depletion law
           of Wolfire et al. (1995), Table 2. */
        double silicon = pow(10.0, fmin(-4.45, -5.45 - 0.45 * (logNH + 0.5)));
        StructEQCData cell = {
            .n = {[ipHYDROGEN] = nH, [ipHELIUM] = 0.10 * nH, [ipCARBON] = 1.4e-4 * nH,
                  [ipOXYGEN] = 3.2e-4 * nH, [ipSILICON] = silicon * nH},
            .nel = 1.e-3 * nH,
            .G0 = 1.7, .Nt = Column, .X0 = 1.0,
            .Zd = SOLAR_METALLICITY, .Z = ASPLUND_2009_SOLAR_METALLICITY
        };
        cell.Td = EQCApproximateDustTemperature(cell.Nt, cell.G0);
        StructEQCData equilibrium;
        double residual = FindBalance(&cell, &equilibrium);
        double pressure = EQCLibGetTotalNumberDensity(&equilibrium) * equilibrium.T;
        double xe = equilibrium.nel / nH;
        printf("%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,"
               "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
               model, nH, equilibrium.T, pressure, xe,
               equilibrium.nion[ipHYDROGEN] / nH, 2.0 * equilibrium.nmol[ipH2] / nH,
               equilibrium.H_PE, equilibrium.H_XR, equilibrium.H_CR,
               equilibrium.L_CII158, equilibrium.L_OI63 + equilibrium.L_OI146,
               equilibrium.L_HI, equilibrium.L_DustRec, equilibrium.H_tot, equilibrium.L_tot, residual,
               __wrap_SoftXrayIonizationRateHI(xe, Column), __wrap_CosmicRayIonizationRateHI(xe));
    }
}

int main(void)
{
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamHIIRegion(false);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamH2Formation(true, 1.0);
    EQCLibSetRunParamCOFormation(false, 1.0);
    if(EQCLibTryInitialization(0, NULL) != EQC_SUCCESS){
        fputs("Could not initialize EQCLib\n", stderr);
        return EXIT_FAILURE;
    }

    /* Total EUV/X-ray H ionization includes secondary electrons. Heating
       anchors from Figure 10b are per H nucleus, in erg/s. */
    coldScale = ComputeRadiationScale(ColdElectronFraction,
                                     (RadiationRates){7.5e-17, 7.67e-28, 1.07e-28});
    warmScale = ComputeRadiationScale(WarmElectronFraction,
                                     (RadiationRates){5.3e-17, 2.39e-27, 1.84e-28});

    puts("model,nH_cm3,T_K,P_over_k_K_cm3,x_e,x_HII,x_H2,H_PE,H_X,H_CR,"
         "L_CII,L_OI,L_HI,L_DustRec,heating,cooling,residual,zeta_H_X,zeta_H_CR");
    useWolfireRadiation = true;
    WriteCurve("eqclib");
    useWolfireRadiation = false;
    WriteCurve("eqclib_native_ionization");
    useWolfireRadiation = true;
    EQCLibSetRunParamH2Formation(false, 1.0);
    WriteCurve("eqclib_atomic");
    return fflush(stdout) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
