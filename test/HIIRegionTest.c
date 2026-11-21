#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CollisionalGasPhaseChemistry.h"
#include "HIIRegionFixture.h"
#include "PhotochemicalGasPhaseChemistry.h"

static const double logTMin = 3.5, logTMax = 4.5;

static bool Evaluate(HIIRegion region, double temperature, StructCEQData *state)
{
    *state = HIIRegionGasState(region, temperature);
    const CEQLibStatus status = CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(state, 1.0e-9);
    if(status == CEQ_SUCCESS) return true;
    fprintf(stderr, "H II chemistry failed at n_H=%g, Z/Z_sun=%g, log U=%g, T=%g K (status %d)\n",
            region.nH, region.metallicity, region.logIonization, temperature, status);
    return false;
}

static bool FindThermalEquilibrium(HIIRegion region, StructCEQData *state)
{
    /* Select a crossing from net heating to net cooling at fixed density. */
    const int intervals = 40;
    double lower = pow(10.0, logTMin);
    if(!Evaluate(region, lower, state)) return false;
    double previous = state->H_tot - state->L_tot;
    for(int i = 1; i <= intervals; i++){
        double upper = pow(10.0, logTMin + (logTMax-logTMin)*i/intervals);
        if(!Evaluate(region, upper, state)) return false;
        const double net = state->H_tot - state->L_tot;
        if(previous > 0.0 && net <= 0.0){
            for(int j = 0; j < 50; j++){
                if(!Evaluate(region, 0.5*(lower+upper), state)) return false;
                const double residual = state->H_tot - state->L_tot;
                if(fabs(residual) <= 1.0e-8 * fmax(state->H_tot,state->L_tot)) return true;
                if(residual > 0.0) lower = state->T;
                else upper = state->T;
            }
            break;
        }
        lower = upper;
        previous = net;
    }
    fprintf(stderr, "No resolved thermal equilibrium at n_H=%g, Z/Z_sun=%g, log U=%g\n",
            region.nH, region.metallicity, region.logIonization);
    return false;
}

static double Fraction(double density, double total)
{
    return total > 0.0 ? density/total : NAN;
}

static void WriteSample(FILE *file, HIIRegion region, const StructCEQData *s)
{
    const double heatingH = s->natom[ipHYDROGEN] * HeatingHIPhotoionization(s->NLyn);
    const double heatingHe = s->natom[ipHELIUM] * HeatingHeIPhotoionization(s->NLyn)
                           + s->nion[ipHELIUM] * HeatingHeIIPhotoionization(s->NLyn);
    const double heatingMetals = s->natom[ipCARBON] * HeatingCIPhotoionization(s->NLyn)
                                + s->nion[ipCARBON] * HeatingCIIPhotoionization(s->NLyn)
                               + s->natom[ipOXYGEN] * HeatingOIPhotoionization(s->NLyn)
                                + s->nion[ipOXYGEN] * HeatingOIIPhotoionization(s->NLyn)
                               + s->natom[ipSILICON] * HeatingSiIPhotoionization(s->NLyn)
                                + s->nion[ipSILICON] * HeatingSiIIPhotoionization(s->NLyn);
    const double doubleIons[NumberChemicalElement] = {0.0, s->nHeIII, s->nCIII, s->nOIII, s->nSiIII};
    fprintf(file, "%.9g,%.9g,%.9g,%.12e,%.12e,%.12e", region.nH, region.metallicity,
            region.logIonization, s->NLyn, s->T, s->nel);
    for(int i = 1; i < NumberChemicalElement; i++) fprintf(file, ",%.12e", s->n[i]/region.nH);
    for(int i = 0; i < NumberChemicalElement; i++){
        fprintf(file, ",%.12e,%.12e", Fraction(s->natom[i],s->n[i]), Fraction(s->nion[i],s->n[i]));
        if(i != ipHYDROGEN) fprintf(file, ",%.12e", Fraction(doubleIons[i],s->n[i]));
    }
    const double rates[] = {
        s->H_tot, s->H_PI, heatingH, heatingHe, heatingMetals, s->H_PE, s->H_CR, s->H_XR, s->H_H2,
        s->L_tot, s->L_res, s->L_ff, s->L_HI, s->L_HeI, s->L_HeII, s->L_OII, s->L_OIII,
        s->L_OI63, s->L_OI146, s->L_CII158, s->L_CI610, s->L_CI370, s->L_SiII35,
        s->L_DustRec, s->L_H2, s->L_CO, s->H_tot - s->L_tot,
    };
    for(size_t i = 0; i < sizeof(rates)/sizeof(rates[0]); i++) fprintf(file, ",%.12e", rates[i]);
    fputc('\n', file);
}

static void WriteSampleHeader(FILE *file)
{
    fputs("nH,Z_solar,log_U,NLyn,T,ne,He_H,C_H,O_H,Si_H,HI,HII,HeI,HeII,HeIII,CI,CII,CIII,OI,OII,OIII,SiI,SiII,SiIII,"
          "H_tot,H_PI,H_PI_H,H_PI_He,H_PI_metals,H_PE,H_CR,H_XR,H_H2,"
          "L_tot,L_res,L_ff,L_HI,L_HeI,L_HeII,L_OII,L_OIII,L_OI63,L_OI146,L_CII158,L_CI610,L_CI370,"
          "L_SiII35,L_DustRec,L_H2,L_CO,net_heating\n", file);
}

static double HeliumIIRecombination(double temperature)
{
    return HeIIRecombinationRate(temperature) + HeIIDielectronicRecombinationRate(temperature);
}

static bool WriteReactions(FILE *file, HIIRegion region, const StructCEQData *s, const char *sample)
{
    /* Evaluate n_i Gamma_i, n_i n_e C_i, n_{i+1} n_e alpha_{i+1} and the charge-transfer
       fluxes using the same rate coefficients as the solver. */
    const struct {
        const char *ion;
        double lower, upper;
        double (*photo)(double), (*collision)(double), (*recombine)(double);
    } transitions[] = {
        {"HI", s->natom[ipHYDROGEN], s->nion[ipHYDROGEN], HIPhotoIonizationRate, HICollisionalIonizationRate, HIIRecombinationRate},
        {"HeI", s->natom[ipHELIUM], s->nion[ipHELIUM], HeIPhotoIonizationRate, HeICollisionalIonizationRate, HeliumIIRecombination},
        {"HeII", s->nion[ipHELIUM], s->nHeIII, HeIIPhotoIonizationRate, HeIICollisionalIonizationRate, HeIIIRecombinationRate},
        {"CI", s->natom[ipCARBON], s->nion[ipCARBON], CIPhotoionizationRate, CICollisionalIonizationRate, CIIRecombinationRate},
        {"CII", s->nion[ipCARBON], s->nCIII, CIIPhotoionizationRate, CIICollisionalIonizationRate, CIIIRecombinationRate},
        {"OI", s->natom[ipOXYGEN], s->nion[ipOXYGEN], OIPhotoionizationRate, OICollisionalIonizationRate, OIIRecombinationRate},
        {"OII", s->nion[ipOXYGEN], s->nOIII, OIIPhotoionizationRate, OIICollisionalIonizationRate, OIIIRecombinationRate},
        {"SiI", s->natom[ipSILICON], s->nion[ipSILICON], SiIPhotoionizationRate, SiICollisionalIonizationRate, SiIIRecombinationRate},
        {"SiII", s->nion[ipSILICON], s->nSiIII, SiIIPhotoionizationRate, SiIICollisionalIonizationRate, SiIIIRecombinationRate},
    };
    /* H-O charge transfer: O + H+ -> O+ + H, O+ + H -> O + H+, O++ + H -> O+ + H+. */
    const double toOxygenII = s->natom[ipOXYGEN]*s->nion[ipHYDROGEN]*OIChargeTransferRateToHII(s->T);
    const double toOxygenI = s->nion[ipOXYGEN]*s->natom[ipHYDROGEN]*OIIChargeTransferRateToHI(s->T);
    const double toOxygenIIfromIII = s->nOIII*s->natom[ipHYDROGEN]*OIIIChargeTransferRateToHI(s->T);
    const double transferIonization[] = {toOxygenI + toOxygenIIfromIII, 0, 0, 0, 0, toOxygenII, 0, 0, 0};
    const double transferRecombination[] = {toOxygenII, 0, 0, 0, 0, toOxygenI, toOxygenIIfromIII, 0, 0};
    for(size_t i = 0; i < sizeof(transitions)/sizeof(transitions[0]); i++){
        const double photo = transitions[i].lower * transitions[i].photo(s->NLyn);
        const double collision = transitions[i].lower * s->nel * transitions[i].collision(s->T);
        const double recombination = transitions[i].upper * s->nel * transitions[i].recombine(s->T);
        const double ionization = photo + collision + transferIonization[i];
        const double loss = recombination + transferRecombination[i];
        if(!(fabs(ionization - loss) <= 1.0e-6 * fmax(ionization,loss))){
            fprintf(stderr, "%s reaction balance failed at n_H=%g, Z/Z_sun=%g, log U=%g, T=%g K\n",
                    transitions[i].ion, region.nH, region.metallicity, region.logIonization, s->T);
            return false;
        }
        fprintf(file, "%s,%.9g,%.9g,%.9g,%.12e,%s,%.12e,%.12e,%.12e,%.12e,%.12e\n", sample,
                region.nH, region.metallicity, region.logIonization, s->T, transitions[i].ion,
                photo, collision, recombination, transferIonization[i], transferRecombination[i]);
    }
    return true;
}

static bool CheckNeutralSeed(HIIRegion region, const StructCEQData *equilibrium)
{
    /* A nearly neutral seed takes the photoionized electron-density path;
       the ionized seed of the equilibrium solve uses the general search. */
    StructCEQData seeded = HIIRegionGasState(region, equilibrium->T);
    seeded.nel = 1.0e-3*region.nH;
    if(CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(&seeded, 1.0e-9) == CEQ_SUCCESS &&
       fabs(seeded.nel-equilibrium->nel) <= 1.0e-8*equilibrium->nel)
        return true;
    fprintf(stderr, "Neutral seed reached a different electron density at n_H=%g, Z/Z_sun=%g, log U=%g\n",
            region.nH, region.metallicity, region.logIonization);
    return false;
}

static bool RunModels(FILE *equilibrium, FILE *thermal, FILE *reactions)
{
    WriteSampleHeader(equilibrium);
    WriteSampleHeader(thermal);
    fputs("sample,nH,Z_solar,log_U,T,ion,photoionization,collisional_ionization,recombination,"
          "charge_transfer_ionization,charge_transfer_recombination\n", reactions);
    StructCEQData state;
    const double densities[] = {1.0, 10.0, 100.0, 1000.0};
    const double metallicities[] = {0.0, 0.1, 1.0};
    for(size_t i = 0; i < sizeof(densities)/sizeof(densities[0]); i++){
        for(size_t z = 0; z < sizeof(metallicities)/sizeof(metallicities[0]); z++){
            for(int u = 0; u <= 8; u++){
                const HIIRegion region = {densities[i], metallicities[z], -4.0 + 0.25*u};
                if(!FindThermalEquilibrium(region, &state) || !CheckNeutralSeed(region, &state) ||
                   !WriteReactions(reactions, region, &state, "equilibrium"))
                    return false;
                WriteSample(equilibrium, region, &state);
            }
        }
    }

    const HIIRegion region = {100.0, 1.0, -3.0};
    const int intervals = 120;
    for(int i = 0; i <= intervals; i++){
        const double temperature = pow(10.0, logTMin + (logTMax-logTMin)*i/intervals);
        if(!Evaluate(region, temperature, &state) || !WriteReactions(reactions, region, &state, "thermal_curve"))
            return false;
        WriteSample(thermal, region, &state);
    }
    return true;
}

int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr, "Usage: %s PEGASE2_DIRECTORY\n", argv[0]);
        return EXIT_FAILURE;
    }
    if(InitializeHIIRegionFixture(argv[1]) != CEQ_SUCCESS) return EXIT_FAILURE;

    const char *paths[] = {"equilibrium.csv", "thermal_curve.csv", "reactions.csv"};
    FILE *files[sizeof(paths)/sizeof(paths[0])] = {NULL};
    bool success = true;
    for(size_t i = 0; i < sizeof(files)/sizeof(files[0]); i++){
        files[i] = fopen(paths[i], "w");
        if(files[i] == NULL){
            perror(paths[i]);
            success = false;
            break;
        }
    }
    if(success) success = RunModels(files[0], files[1], files[2]);
    for(size_t i = 0; i < sizeof(files)/sizeof(files[0]); i++){
        if(files[i] != NULL && fclose(files[i]) != 0){
            perror(paths[i]);
            success = false;
        }
    }
    if(success) puts("H II region equilibrium, seed and reaction balance checks passed.");
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
