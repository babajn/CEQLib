#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"
#include "PhotochemicalGasPhaseChemistry.h"

enum { PhaseSamples = 121 };
static const double Accuracy = 1.0e-8;
enum { DirectRates, TableRates, RateModeCount };

typedef struct {
    double temperature, fraction, totalCooling;
} TrajectoryError;

static void InitializeEnergy(StructEQCData *state)
{
    state->nel = EQCLibGetElectronNumberDensity(state);
    state->Rho = EQCLibGetTotalMassDensity(state);
    state->MeanMolecularWeight = EQCLibGetMeanMolecularWeight(state);
    state->DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(state);
    state->U = EQCLibConvertTtoU(state);
    state->invT = 1.0/state->T;
    state->logT = log10(state->T);
}

static StructEQCData Cell(double nH, double temperature, double ionized, double molecular)
{
    StructEQCData state = {.T = temperature, .Td = 20.0};
    state.n[ipHYDROGEN] = nH;
    state.n[ipHELIUM] = 0.10*nH;
    state.n[ipCARBON] = 3.0e-4*nH;
    state.n[ipOXYGEN] = 5.0e-4*nH;
    for(int element = 0; element < NumberChemicalElement; element++)
        state.natom[element] = state.n[element];
    state.natom[ipHYDROGEN] = nH*(1.0-ionized-molecular);
    state.nion[ipHYDROGEN] = nH*ionized;
    state.nmol[ipH2] = 0.5*nH*molecular;
    return state;
}

/* Verify the inventory from the returned densities, independently of the EOS. */
static bool CheckInventory(const StructEQCData *state, const StructEQCData *initial)
{
    const double extra[NumberChemicalElement] = {
        [ipHYDROGEN] = 2.0*state->nmol[ipH2], [ipHELIUM] = state->nHeIII,
        [ipCARBON] = state->nCIII, [ipOXYGEN] = state->nOIII,
    };
    double charge = 0.0;
    for(int i = 0; i < NumberChemicalElement; i++){
        const double atom = state->natom[i], ion = state->nion[i], total = state->n[i];
        if(!(total == initial->n[i] && atom >= 0.0 && ion >= 0.0 && extra[i] >= 0.0 &&
             fabs(atom+ion+extra[i]-total) <= 1.0e-12*total)) return false;
        charge += ion + (i == ipHYDROGEN ? 0.0 : 2.0*extra[i]);
    }
    return state->Rho == initial->Rho && fabs(charge-state->nel) <= 1.0e-12*state->n[ipHYDROGEN] &&
           (EQCLibRunParameters.UseH2Formation || state->nmol[ipH2] == 0.0) &&
           isfinite(state->T) && state->T > 0.0 && isfinite(state->U) && state->U > 0.0;
}

static bool Advance(StructEQCData *state, double seconds, const char *test)
{
    const StructEQCData initial = *state;
    state->dt = seconds;
    const EQCLibStatus status = EQCLibTryAdvanceHydrogenThermochemistry(state,Accuracy);
    if(status == EQC_SUCCESS && CheckInventory(state,&initial)) return true;
    fprintf(stderr,"%s: update or inventory check failed (status %d, n_H=%g, T=%g K, dt=%g yr)\n",
            test,status,state->n[ipHYDROGEN],state->T,seconds/YEAR_CGS);
    return false;
}

static bool WriteState(FILE *file, const char *test, const char *phase, double start, double elapsed,
                       const StructEQCData *state)
{
    /* Compare at the same T and radiation using the same rate evaluator for
       both compositions. Also record the total from the native equilibrium
       rate prescriptions, which treat ion and collider abundances differently. */
    StructEQCData equilibrium = *state;
    const int chemistryStatus = EQCLibChemicalReactionEquilibrium(&equilibrium,Accuracy);
    if(chemistryStatus != EQC_SUCCESS){
        fprintf(stderr,"%s %s, t=%g yr, T=%g K: equilibrium chemistry failed (status %d)\n",
                test,phase,elapsed/YEAR_CGS,state->T,chemistryStatus);
        return false;
    }
    StructEQCData native = equilibrium;
    const int sharedStatus = CEQCalculateNonEquilibriumCoolingHeatingRates(&equilibrium);
    const int nativeStatus = EQCLibCalculateCoolingHeatingRates(&native);
    if(sharedStatus != EQC_SUCCESS || nativeStatus != EQC_SUCCESS){
        fprintf(stderr,"%s %s, t=%g yr, T=%g K: equilibrium cooling failed (shared status %d, native status %d)\n",
                test,phase,elapsed/YEAR_CGS,state->T,sharedStatus,nativeStatus);
        return false;
    }

    const double nH = state->n[ipHYDROGEN];
    fprintf(file,"%s,%s,%.12e,%.12e",test,phase,(start+elapsed)/YEAR_CGS,elapsed/YEAR_CGS);
    const double values[] = {
        nH,state->Nt,state->NLyn,state->G0,state->Zd,state->T,state->U,state->nel/nH,
        state->natom[ipHYDROGEN]/nH,state->nion[ipHYDROGEN]/nH,2.0*state->nmol[ipH2]/nH,
        state->natom[ipHELIUM]/state->n[ipHELIUM],state->nion[ipHELIUM]/state->n[ipHELIUM],state->nHeIII/state->n[ipHELIUM],
        state->natom[ipCARBON]/state->n[ipCARBON],state->nion[ipCARBON]/state->n[ipCARBON],state->nCIII/state->n[ipCARBON],
        state->natom[ipOXYGEN]/state->n[ipOXYGEN],state->nion[ipOXYGEN]/state->n[ipOXYGEN],state->nOIII/state->n[ipOXYGEN],
        state->H_tot,state->H_PI,state->H_PE,state->H_H2,state->L_tot,state->L_HI,state->L_HeI,state->L_HeII,
        state->L_CII158,state->L_OI63,state->L_OII,state->L_OIII,state->L_H2,
        equilibrium.nel/nH,2.0*equilibrium.nmol[ipH2]/nH,
        equilibrium.L_tot,equilibrium.L_CII158,equilibrium.L_OII,equilibrium.L_OIII,equilibrium.L_H2,native.L_tot,
    };
    for(size_t i = 0; i < sizeof values/sizeof values[0]; i++) fprintf(file,",%.12e",values[i]);
    return fputc('\n',file) != EOF;
}

static double MaximumFractionDifference(const StructEQCData *a, const StructEQCData *b)
{
    const double extraA[] = {2.0*a->nmol[ipH2],a->nHeIII,a->nCIII,a->nOIII};
    const double extraB[] = {2.0*b->nmol[ipH2],b->nHeIII,b->nCIII,b->nOIII};
    double difference = 0.0;
    for(int i = ipHYDROGEN; i <= ipOXYGEN; i++){
        const double ions = fmax(fabs(a->nion[i]-b->nion[i]),fabs(extraA[i]-extraB[i]));
        difference = fmax(difference,fmax(ions,fabs(a->natom[i]-b->natom[i]))/a->n[i]);
    }
    return difference;
}

static bool CheckTrajectory(const StructEQCData *direct, const StructEQCData *lookup, TrajectoryError *maximum,
                            const char *scenario, const char *phase, double elapsed)
{
    const double temperature = fabs(lookup->T-direct->T)/direct->T;
    const double fraction = fmax(MaximumFractionDifference(direct,lookup),
                                 fabs(lookup->nel-direct->nel)/direct->n[ipHYDROGEN]);
    const double cooling = lookup->L_tot == direct->L_tot ? 0.0 : fabs((lookup->L_tot-direct->L_tot)/direct->L_tot);
    if(!(temperature < 5.0e-5 && fraction < 5.0e-6 && cooling < 2.0e-4)){
        fprintf(stderr,"Temperature-table comparison failed at %s %s, t=%g yr (T rel=%g, fraction abs=%g, L_tot rel=%g)\n",
                scenario,phase,elapsed/YEAR_CGS,temperature,fraction,cooling);
        return false;
    }
    maximum->temperature = fmax(maximum->temperature,temperature);
    maximum->fraction = fmax(maximum->fraction,fraction);
    maximum->totalCooling = fmax(maximum->totalCooling,cooling);
    return true;
}

static bool RunPhase(FILE *const files[RateModeCount], StructEQCData states[RateModeCount],
                     const char *test, const char *phase, double start, double first, double end,
                     TrajectoryError *maximum)
{
    double previous = 0.0;
    for(int i = -1; i < PhaseSamples; i++){
        double elapsed = 0.0;
        if(i >= 0) elapsed = i == PhaseSamples-1 ? end : first*pow(end/first,(double)i/(PhaseSamples-1));
        for(int mode = 0; mode < RateModeCount; mode++){
            StructEQCData *state = &states[mode];
            EQCLibSetRunParamNonEquilibriumTemperatureTables(mode == TableRates);
            /* At the switch, record the jump in rates with composition and U held fixed. */
            if(i < 0){
                if(CEQCalculateNonEquilibriumCoolingHeatingRates(state) != EQC_SUCCESS){
                    fprintf(stderr,"%s %s: cannot evaluate rates at the illumination switch\n",test,phase);
                    return false;
                }
            }else if(!Advance(state,elapsed-previous,test)) return false;
            if(!WriteState(files[mode],test,phase,start,elapsed,state)) return false;
        }
        if(!CheckTrajectory(&states[DirectRates],&states[TableRates],maximum,test,phase,elapsed)) return false;
        previous = elapsed;
    }
    return true;
}

static bool RunPulses(FILE *const files[RateModeCount])
{
    for(int mode = 0; mode < RateModeCount; mode++)
        fputs("scenario,phase,time_yr,phase_time_yr,nH,Nt,NLyn,G0,Zd,T,U,xe,HI,HII,H2,HeI,HeII,HeIII,CI,CII,CIII,OI,OII,OIII,"
              "H_tot,H_PI,H_PE,H_H2,L_tot,L_HI,L_HeI,L_HeII,L_CII158,L_OI63,L_OII,L_OIII,L_H2,"
              "eq_xe,eq_H2,eq_L_tot,eq_L_CII158,eq_L_OII,eq_L_OIII,eq_L_H2,eq_native_L_tot\n",files[mode]);
    TrajectoryError maximum = {0};
    EQCLibSetRunParamH2Formation(false,1.0);
    StructEQCData ionizing = Cell(10.0,8000.0,0.0,0.0);
    InitializeEnergy(&ionizing);
    ionizing.NLyn = 1.0e-7/HIPhotoIonizationRate(1.0);
    const double ionizingDuration = 5.0e8;
    StructEQCData ionizingStates[RateModeCount] = {ionizing,ionizing};
    if(!RunPhase(files,ionizingStates,"ionizing","on",0.0,1.0e-5*YEAR_CGS,ionizingDuration,&maximum)) return false;
    StructEQCData ionized[RateModeCount];
    for(int mode = 0; mode < RateModeCount; mode++){
        ionized[mode] = ionizingStates[mode];
        const StructEQCData *state = &ionized[mode];
        if(!(state->T > 8000.0 && state->nHeIII > 0.0 && state->nCIII > 0.0 && state->nOIII > 0.0)){
            fprintf(stderr,"Ionizing pulse did not heat and ionize the cell\n");
            return false;
        }
        ionizingStates[mode].NLyn = 0.0;
    }
    if(!RunPhase(files,ionizingStates,"ionizing","off",ionizingDuration,1.0e-3*YEAR_CGS,1.0e6*YEAR_CGS,&maximum))
        return false;
    for(int mode = 0; mode < RateModeCount; mode++){
        const StructEQCData *state = &ionizingStates[mode], *before = &ionized[mode];
        if(!(state->T < before->T && state->nel < before->nel && state->L_OIII < before->L_OIII &&
             state->L_OII < before->L_OII && state->L_CII158 < before->L_CII158)){
            fprintf(stderr,"Ionizing pulse did not cool and recombine in the dark\n");
            return false;
        }
    }

    EQCLibSetRunParamH2Formation(true,1.0);
    StructEQCData fuv = Cell(100.0,100.0,0.0,0.8);
    fuv.Nt = 2.0e20;
    fuv.Zd = SOLAR_METALLICITY;
    fuv.G0 = 1.0e3;
    InitializeEnergy(&fuv);
    const double fuvDuration = 3.0e5*YEAR_CGS;
    StructEQCData fuvStates[RateModeCount] = {fuv,fuv};
    if(!RunPhase(files,fuvStates,"fuv","on",0.0,1.0e-4*YEAR_CGS,fuvDuration,&maximum)) return false;
    double irradiatedH2[RateModeCount];
    for(int mode = 0; mode < RateModeCount; mode++){
        irradiatedH2[mode] = fuvStates[mode].nmol[ipH2];
        if(!(irradiatedH2[mode] < 40.0 && fuvStates[mode].nion[ipCARBON] > 0.0)){
            fprintf(stderr,"FUV pulse did not dissociate H2 and ionize carbon\n");
            return false;
        }
        fuvStates[mode].G0 = 0.0;
    }
    if(!RunPhase(files,fuvStates,"fuv","off",fuvDuration,1.0e-2*YEAR_CGS,3.0e7*YEAR_CGS,&maximum)) return false;
    for(int mode = 0; mode < RateModeCount; mode++){
        if(!(fuvStates[mode].nmol[ipH2] > irradiatedH2[mode])){
            fprintf(stderr,"H2 did not reform after the FUV pulse\n");
            return false;
        }
    }
    printf("Temperature-table pulse differences: T relative %.6g, fractions absolute %.6g, L_tot relative %.6g\n",
           maximum.temperature,maximum.fraction,maximum.totalCooling);
    return true;
}

static bool RunRecombination(FILE *file)
{
    EQCLibSetRunParamH2Formation(false,1.0);
    StructEQCData state = {.T = 1000.0, .n = {10.0}, .nion = {10.0}, .nel = 10.0};
    const double nH = state.n[ipHYDROGEN];
    const double recombinationTime = 1.0/(nH*HIIRecombinationRate(state.T));
    double previous = 0.0;
    fputs("time_yr,time_over_trec,HII,analytic_HII,abs_error\n0,0,1,1,0\n",file);
    for(int i = 0; i < PhaseSamples; i++){
        const double normalizedTime = pow(10.0,-5.0+8.0*i/(PhaseSamples-1));
        const double time = normalizedTime*recombinationTime;
        state.dt = time-previous;
        const EQCLibStatus status = EQCLibTryAdvanceHydrogenChemistry(&state,Accuracy);
        const double exact = 1.0/(1.0+normalizedTime);
        const double fraction = state.nion[ipHYDROGEN]/nH;
        const double error = fabs(fraction-exact);
        if(status != EQC_SUCCESS || !(error < 2.0e-6) || state.nmol[ipH2] != 0.0){
            fprintf(stderr,"Analytic recombination comparison failed at t/t_rec=%g (status %d, error %g)\n",
                    normalizedTime,status,error);
            return false;
        }
        fprintf(file,"%.12e,%.12e,%.12e,%.12e,%.12e\n",time/YEAR_CGS,normalizedTime,fraction,exact,error);
        previous = time;
    }
    return true;
}

static bool RunRecombinationBoth(FILE *const files[RateModeCount])
{
    for(int mode = 0; mode < RateModeCount; mode++){
        EQCLibSetRunParamNonEquilibriumTemperatureTables(mode == TableRates);
        if(!RunRecombination(files[mode])) return false;
    }
    return true;
}

static bool CheckCaseAAndTemperatureFallback(void)
{
    const bool caseB = EQCLibRunParameters.UseRecombinationCaseB;
    const bool h2 = EQCLibRunParameters.UseH2Formation;
    const double clumping = EQCLibRunParameters.CpForH2;
    const bool tables = EQCLibRunParameters.UseNonEquilibriumTemperatureTables;
    EQCLibSetRunParamRecombinationCaseB(false);
    EQCLibSetRunParamH2Formation(false,1.0);

    const double nH = 10.0;
    const double temperature = 1200.0;
    const double duration = 0.5/(nH*HIIRecombinationRateCaseA(temperature));
    const StructEQCData hot = {.T = 2.0e8, .dt = YEAR_CGS, .n = {nH}, .natom = {0.5*nH},
                               .nion = {0.5*nH}, .nel = 0.5*nH};
    StructEQCData updated[RateModeCount] = {hot,hot};
    bool success = true;
    for(int mode = 0; mode < RateModeCount; mode++){
        StructEQCData recombination = {.T = temperature, .dt = duration, .n = {nH}, .nion = {nH}, .nel = nH};
        EQCLibSetRunParamNonEquilibriumTemperatureTables(mode == TableRates);
        const EQCLibStatus recombinationStatus = EQCLibTryAdvanceHydrogenChemistry(&recombination,Accuracy);
        const double error = fabs(recombination.nion[ipHYDROGEN]/nH-2.0/3.0);
        if(recombinationStatus != EQC_SUCCESS || !(error < 2.0e-6)){
            fprintf(stderr,"Case A recombination failed (mode %d, status %d, error %g)\n",
                    mode,recombinationStatus,error);
            success = false;
        }
        const EQCLibStatus hotStatus = EQCLibTryAdvanceHydrogenChemistry(&updated[mode],Accuracy);
        if(hotStatus != EQC_SUCCESS){
            fprintf(stderr,"Out-of-range temperature update failed (mode %d, status %d)\n",mode,hotStatus);
            success = false;
        }
    }
    const StructEQCData *direct = &updated[DirectRates], *lookup = &updated[TableRates];
    if(direct->natom[ipHYDROGEN] != lookup->natom[ipHYDROGEN] || direct->nion[ipHYDROGEN] != lookup->nion[ipHYDROGEN] ||
       direct->nmol[ipH2] != lookup->nmol[ipH2] || direct->nel != lookup->nel){
        fputs("Out-of-range temperature fallback disagrees with direct evaluation\n",stderr);
        success = false;
    }
    EQCLibSetRunParamNonEquilibriumTemperatureTables(tables);
    EQCLibSetRunParamH2Formation(h2,clumping);
    EQCLibSetRunParamRecombinationCaseB(caseB);
    return success;
}

static bool RunConvergence(FILE *const files[RateModeCount])
{
    const double densities[] = {0.03,1.0,1.0e3};
    const double temperatures[] = {30.0,100.0,8000.0,1.0e5,1.0e6};
    const double steps[] = {1.0e2,2.0e4,1.0e6};
    const char *const radiation[] = {"dark","fuv","ionizing"};
    const double ionized[] = {0.005,0.005,0.20,0.95,0.95};
    const double molecular[] = {0.30,0.30,0.05,0.02,0.0};
    for(int mode = 0; mode < RateModeCount; mode++)
        fputs("nH,T_initial,radiation,dt_yr,H2,T_single,T_split,T_rel_error,U_rel_error,fraction_abs_error\n",files[mode]);
    double maximumTemperature = 0.0, maximumFraction = 0.0;
    EQCLibSetRunParamCosmicRayIonization(true);
    for(int h2 = 0; h2 <= 1; h2++){
        EQCLibSetRunParamH2Formation(h2 != 0,1.0);
        for(size_t n = 0; n < sizeof densities/sizeof densities[0]; n++){
            const double nH = densities[n];
            for(size_t t = 0; t < sizeof temperatures/sizeof temperatures[0]; t++){
                const double temperature = temperatures[t];
                StructEQCData initial = Cell(nH,temperature,ionized[t],h2 ? molecular[t] : 0.0);
                initial.Nt = 1.0e20*sqrt(nH);
                initial.Zd = SOLAR_METALLICITY;
                initial.X0 = 1.0;
                const double singly = t < 3 ? 0.15 : 0.10;
                const double doubly = t < 3 ? 0.01 : 0.80;
                for(int i = ipHELIUM; i <= ipOXYGEN; i++){
                    initial.natom[i] = initial.n[i]*(1.0-singly-doubly);
                    initial.nion[i] = initial.n[i]*singly;
                }
                initial.nHeIII = initial.n[ipHELIUM]*doubly;
                initial.nCIII = initial.n[ipCARBON]*doubly;
                initial.nOIII = initial.n[ipOXYGEN]*doubly;
                InitializeEnergy(&initial);
                for(size_t light = 0; light < sizeof radiation/sizeof radiation[0]; light++){
                    initial.G0 = light > 0 ? 1.7 : 0.0;
                    initial.NLyn = light == 2 ? 1.0e-3*nH*LIGHT_C_CGS : 0.0;
                    for(size_t d = 0; d < sizeof steps/sizeof steps[0]; d++){
                        const double dt = steps[d]*YEAR_CGS;
                        StructEQCData singles[RateModeCount] = {initial,initial};
                        StructEQCData splits[RateModeCount] = {initial,initial};
                        for(int mode = 0; mode < RateModeCount; mode++){
                            StructEQCData *single = &singles[mode], *split = &splits[mode];
                            EQCLibSetRunParamNonEquilibriumTemperatureTables(mode == TableRates);
                            if(!Advance(single,dt,radiation[light])) return false;
                            for(int substep = 0; substep < 8; substep++)
                                if(!Advance(split,dt/8.0,radiation[light])) return false;
                            const double errorT = fabs(single->T-split->T)/fmax(single->T,split->T);
                            const double errorU = fabs(single->U-split->U)/fmax(single->U,split->U);
                            const double errorSpecies = MaximumFractionDifference(single,split);
                            if(!(errorT < 5.0e-5 && errorU < 5.0e-5 && errorSpecies < 5.0e-6)){
                                fprintf(stderr,"Step splitting failed at n_H=%g, T=%g K, %s, dt=%g yr, H2=%d (%g, %g, %g)\n",
                                        nH,temperature,radiation[light],steps[d],h2,errorT,errorU,errorSpecies);
                                return false;
                            }
                            fprintf(files[mode],"%g,%g,%s,%g,%d,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                                    nH,temperature,radiation[light],steps[d],h2,single->T,split->T,errorT,errorU,errorSpecies);
                        }
                        const double singleT = fabs(singles[TableRates].T-singles[DirectRates].T)/singles[DirectRates].T;
                        const double splitT = fabs(splits[TableRates].T-splits[DirectRates].T)/splits[DirectRates].T;
                        const double errorT = fmax(singleT,splitT);
                        const double errorSpecies = fmax(MaximumFractionDifference(&singles[TableRates],&singles[DirectRates]),
                                                         MaximumFractionDifference(&splits[TableRates],&splits[DirectRates]));
                        if(!(errorT < 5.0e-5 && errorSpecies < 5.0e-6)){
                            fprintf(stderr,"Temperature-table convergence comparison failed at n_H=%g, T=%g K, %s, dt=%g yr, H2=%d (T rel=%g, fractions abs=%g)\n",
                                    nH,temperature,radiation[light],steps[d],h2,errorT,errorSpecies);
                            return false;
                        }
                        maximumTemperature = fmax(maximumTemperature,errorT);
                        maximumFraction = fmax(maximumFraction,errorSpecies);
                    }
                }
            }
        }
    }
    printf("Temperature-table step-splitting differences: T relative %.6g, fractions absolute %.6g\n",
           maximumTemperature,maximumFraction);
    return true;
}

static bool RunOutputFiles(const char *const names[RateModeCount], bool (*run)(FILE *const files[RateModeCount]))
{
    FILE *files[RateModeCount] = {NULL};
    bool success = true;
    for(int mode = 0; mode < RateModeCount; mode++){
        files[mode] = fopen(names[mode],"w");
        if(files[mode] == NULL){
            perror(names[mode]);
            success = false;
            break;
        }
    }
    if(success) success = run(files);
    for(int mode = 0; mode < RateModeCount; mode++){
        if(files[mode] && fclose(files[mode]) != 0){
            perror(names[mode]);
            success = false;
        }
    }
    return success;
}

int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr,"Usage: %s PEGASE2_DIRECTORY\n",argv[0]);
        return EXIT_FAILURE;
    }
    EQCLibSetRunParamHIIRegion(true);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCIEResonanceCooling(false);
    if(EQCLibTryInitialization(0,argv[1]) != EQC_SUCCESS) return EXIT_FAILURE;
    if(!CheckCaseAAndTemperatureFallback()) return EXIT_FAILURE;

    const struct { const char *names[RateModeCount]; bool (*run)(FILE *const files[RateModeCount]); } tests[] = {
        {{"evolution_direct.csv", "evolution.csv"}, RunPulses},
        {{"recombination_direct.csv", "recombination.csv"}, RunRecombinationBoth},
        {{"convergence_direct.csv", "convergence.csv"}, RunConvergence},
    };
    for(size_t i = 0; i < sizeof tests/sizeof tests[0]; i++)
        if(!RunOutputFiles(tests[i].names,tests[i].run)) return EXIT_FAILURE;
    puts("Non-equilibrium pulse, analytic recombination and step-splitting checks passed for direct evaluations and tables.");
    return EXIT_SUCCESS;
}
