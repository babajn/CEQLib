#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "Astro.h"
#include "Constants.h"
#include "ChemicalElements.h"
#include "HIIregionLineCooling.h"
#include "HIIregionOxygenData.h"

void InitializeCoolingLinesHIIregion(int Check)
{
    (void)Check;
}

static bool SolveOxygenLevels(double balance[OxygenLevelCount][OxygenLevelCount],
                              double population[OxygenLevelCount])
{
    /* Fix the ground population to one and solve for the excited-to-ground ratios.
       The excited block is a nonsingular M-matrix, so elimination needs no pivoting. */
    for(int level = 1; level < OxygenLevelCount; level++){
        const double diagonal = balance[level][level];
        if(diagonal <= 0.0 || !isfinite(diagonal)) return false;

        for(int row = level + 1; row < OxygenLevelCount; row++){
            const double factor = balance[row][level]/diagonal;
            balance[row][0] -= factor*balance[level][0];
            for(int column = level + 1; column < OxygenLevelCount; column++)
                balance[row][column] -= factor*balance[level][column];
        }
    }

    double totalPopulation = 1.0;
    population[0] = 1.0;
    for(int level = OxygenLevelCount - 1; level > 0; level--){
        double gain = -balance[level][0];
        for(int upper = level + 1; upper < OxygenLevelCount; upper++)
            gain -= balance[level][upper]*population[upper];
        population[level] = gain/balance[level][level];
        if(population[level] < 0.0 || !isfinite(population[level])) return false;
        totalPopulation += population[level];
    }

    for(int level = 0; level < OxygenLevelCount; level++)
        population[level] /= totalPopulation;
    return true;
}

static double CoolingOxygenIon(const OxygenIonData *ion, double T, double ne, double density,
                               double lines[OxygenLevelCount][OxygenLevelCount])
{
    if(!isfinite(T) || T <= 0.0 || !isfinite(ne) || ne < 0.0 ||
       !isfinite(density) || density < 0.0)
        return NAN;
    if(lines != NULL) memset(lines,0,OxygenLevelCount*sizeof(*lines));
    if(ne == 0.0 || density == 0.0) return 0.0;

    const double logT = log10(T);
    const double *grid = ion->logTemperature;
    size_t index = 0;
    while(index + 2 < ion->temperatureCount && logT > grid[index + 1])
        index++;
    const double fraction = fmax(0.0,fmin(1.0,(logT-grid[index])/(grid[index + 1]-grid[index])));

    /* Scale all transition rates together to keep large electron densities finite. */
    const double scale = fmax(1.0,ne);
    const double electronRate = (ne/scale) * 8.629e-6/sqrt(T);
    double balance[OxygenLevelCount][OxygenLevelCount] = {{0.0}};
    double energies[OxygenTransitionCount];
    for(size_t i = 0; i < OxygenTransitionCount; i++){
        const OxygenTransition *transition = &ion->transitions[i];
        const unsigned int upper = transition->upper;
        const unsigned int lower = transition->lower;
        const OxygenLevel *levels = ion->levels;
        const double energy = PLANK_H_CGS*LIGHT_C_CGS*(levels[upper].wavenumber-levels[lower].wavenumber);
        const double *strength = transition->collisionStrength;
        const double collision = electronRate*(strength[index] + fraction*(strength[index + 1]-strength[index]));
        const double upward = collision/levels[lower].statisticalWeight * exp(-(energy/BOLTZMANN_CONSTANT_CGS)/T);
        const double downward = collision/levels[upper].statisticalWeight + transition->spontaneousRate/scale;

        balance[lower][lower] += upward;
        balance[upper][upper] += downward;
        balance[lower][upper] -= downward;
        balance[upper][lower] -= upward;
        energies[i] = energy;
    }

    double population[OxygenLevelCount];
    if(!SolveOxygenLevels(balance,population)) return NAN;

    double total = 0.0;
    for(size_t i = 0; i < OxygenTransitionCount; i++){
        const OxygenTransition *transition = &ion->transitions[i];
        const double power = density * population[transition->upper] * transition->spontaneousRate * energies[i];
        if(lines != NULL) lines[transition->upper][transition->lower] = power;
        total += power;
    }
    return total;
}

double CoolingLinesOII(double T, double ne, double nOII,
                       double lines[OxygenLevelCount][OxygenLevelCount])
{
    return CoolingOxygenIon(&OxygenII,T,ne,nOII,lines);
}

double CoolingLinesOIII(double T, double ne, double nOIII,
                        double lines[OxygenLevelCount][OxygenLevelCount])
{
    return CoolingOxygenIon(&OxygenIII,T,ne,nOIII,lines);
}

double CoolingLinesHIIregionOxygen(double T, double ne, double nOII, double nOIII)
{
    return CoolingLinesOII(T,ne,nOII,NULL) + CoolingLinesOIII(T,ne,nOIII,NULL);
}

// 4*pi* line-integrated emissivity [erg/cc/s]
// ref. eq(4)(5) of Peters et al. MNRAS, 466, 3293, (2017)
static double G13(double T){
    //if(T<4.3e+3)    return 0.0;
    if(T<=2.5e+4)
        return 0.350 + T*(-2.62e-7 + T*(-8.15e-11 + 6.19e-15*T));
    else if(T<5.e+5)
        return 0.276 + T*(4.99e-6 - T*(8.85e-12+ 7.18e-18*T)); 
    else
        return 0.e0;
}
double CoolingHalphaEmission(StructEQCData *PhydroPP){
    double T    = PhydroPP->T;
    double T4   = 1.e-4*T;
    double lnT4 = log(T4);
    double ne   = PhydroPP->nel;
    double nHI  = PhydroPP->natom[ipHYDROGEN];
    double nHII = PhydroPP->nion[ipHYDROGEN];
    // Two processes contribute to the emis- sion of Hα radiation. 
    // The first one is the recombination of ionised hydrogen with a free electron. 
    double j_a_R = 4.0*PI*2.82e-26*ne*nHII*pow(T4,-0.942-0.031*lnT4); 

    // The second contribution comes from collisional excitation of neutral hydrogen by free electrons. 
    double kT    = ERG_TO_EV*BOLTZMANN_CONSTANT_CGS*T;
    double j_a_C = 1.30e-17*ne*nHI*exp(-12.1/kT)*G13(T)/sqrt(T); 

    return j_a_R + j_a_C;   // [erg/cc/s]
}


