#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"

typedef struct {
    EQCLibChemistryMode chemistryMode;
    double nH, temperature, timestepYears, steps;
    double ionizedHydrogen, molecularHydrogen;
    double G0, NLyn, Nt, Zd;
} CellSettings;

static const double CellAccuracy = 1.0e-8;

static char *Trim(char *text)
{
    while(isspace((unsigned char)*text)) text++;
    char *end = text+strlen(text);
    while(end > text && isspace((unsigned char)end[-1])) end--;
    *end = '\0';
    return text;
}

static bool ReadSettings(const char *path, CellSettings *settings)
{
    struct { const char *name; double *value; bool seen; } numbers[] = {
        {"nH", &settings->nH, false}, {"T", &settings->temperature, false},
        {"dt_years", &settings->timestepYears, false}, {"steps", &settings->steps, false},
        {"x_HII", &settings->ionizedHydrogen, false}, {"x_H2", &settings->molecularHydrogen, false},
        {"G0", &settings->G0, false}, {"NLyn", &settings->NLyn, false},
        {"Nt", &settings->Nt, false}, {"Zd", &settings->Zd, false},
    };
    FILE *file = fopen(path,"r");
    if(!file){
        perror(path);
        return false;
    }

    char *line = NULL;
    size_t capacity = 0, lineNumber = 0;
    bool seenMode = false, valid = true;
    while(getline(&line,&capacity,file) >= 0){
        lineNumber++;
        char *comment = strchr(line,'#');
        if(comment) *comment = '\0';
        char *key = Trim(line);
        if(!*key) continue;
        char *equals = strchr(key,'=');
        if(!equals){
            fprintf(stderr,"%s:%zu: expected name = value\n",path,lineNumber);
            valid = false;
            break;
        }
        *equals = '\0';
        key = Trim(key);
        char *value = Trim(equals+1);
        if(strcmp(key,"chemistry_mode") == 0){
            valid = !seenMode && EQCLibTryParseChemistryMode(value,&settings->chemistryMode) == EQC_SUCCESS;
            if(!valid){
                fprintf(stderr,"%s:%zu: set chemistry_mode once to equilibrium or nonequilibrium\n",path,lineNumber);
                break;
            }
            seenMode = true;
            continue;
        }

        size_t i = 0, count = sizeof numbers/sizeof numbers[0];
        while(i < count && strcmp(key,numbers[i].name) != 0) i++;
        if(i == count){
            fprintf(stderr,"%s:%zu: unknown setting '%s'\n",path,lineNumber,key);
            valid = false;
            break;
        }
        char *end;
        errno = 0;
        double number = strtod(value,&end);
        if(numbers[i].seen || end == value || *end || errno == ERANGE || !isfinite(number)){
            fprintf(stderr,"%s:%zu: set %s once to a finite number\n",path,lineNumber,key);
            valid = false;
            break;
        }
        numbers[i].seen = true;
        *numbers[i].value = number;
    }
    if(valid && !feof(file)){
        perror(path);
        valid = false;
    }
    free(line);
    if(fclose(file) != 0){
        perror(path);
        valid = false;
    }
    return valid;
}

static bool ValidateSettings(const CellSettings *settings)
{
    if(settings->nH <= 0.0 || settings->temperature <= 0.0 || settings->timestepYears <= 0.0){
        fputs("nH, T, and dt_years must be positive\n",stderr);
        return false;
    }
    if(settings->steps < 1.0 || settings->steps > INT_MAX || floor(settings->steps) != settings->steps ||
       settings->timestepYears > DBL_MAX/YEAR_CGS/settings->steps){
        fputs("steps must be a positive integer and the total time must be finite\n",stderr);
        return false;
    }
    if(settings->ionizedHydrogen < 0.0 || settings->molecularHydrogen < 0.0 ||
       settings->ionizedHydrogen+settings->molecularHydrogen > 1.0){
        fputs("x_HII and x_H2 must be nonnegative and sum to at most one\n",stderr);
        return false;
    }
    if(settings->G0 < 0.0 || settings->NLyn < 0.0 || settings->Nt < 0.0 || settings->Zd < 0.0){
        fputs("G0, NLyn, Nt, and Zd must be nonnegative\n",stderr);
        return false;
    }
    return true;
}

static EQCLibStatus InitializeCell(const CellSettings *settings, StructEQCData *initial)
{
    const double nH = settings->nH;
    const bool nonequilibrium = settings->chemistryMode == EQC_CHEMISTRY_NONEQUILIBRIUM;
    StructEQCData state = {
        .T = settings->temperature, .dt = settings->timestepYears*YEAR_CGS,
        .G0 = settings->G0, .NLyn = settings->NLyn, .Nt = settings->Nt, .Zd = settings->Zd,
        .Td = EQCApproximateDustTemperature(settings->Nt,settings->G0),
    };
    state.n[ipHYDROGEN] = nH;
    state.n[ipHELIUM] = 0.10*nH;
    state.n[ipCARBON] = 3.0e-4*nH;
    state.n[ipOXYGEN] = 5.0e-4*nH;
    for(int i=0;i<NumberChemicalElement;i++) state.natom[i] = state.n[i];
    if(nonequilibrium){
        state.natom[ipHYDROGEN] = nH*(1.0-(settings->ionizedHydrogen+settings->molecularHydrogen));
        state.nion[ipHYDROGEN] = nH*settings->ionizedHydrogen;
        state.nmol[ipH2] = 0.5*nH*settings->molecularHydrogen;
        state.nel = EQCLibGetElectronNumberDensity(&state);
    }else{
        /* Start the charge search away from the zero-electron root. */
        state.nel = nH;
    }
    state.Rho = EQCLibGetTotalMassDensity(&state);
    if(!isfinite(state.Rho) || state.Rho <= 0.0) return EQC_INVALID_INPUT;
    state.Z = (state.n[ipCARBON]*ChemicalElementMass[ipCARBON]+state.n[ipOXYGEN]*ChemicalElementMass[ipOXYGEN])/state.Rho;
    state.invT = 1.0/state.T;
    state.logT = log10(state.T);
    if(nonequilibrium){
        const int status = CEQCalculateNonEquilibriumCoolingHeatingRates(&state);
        if(status != EQC_SUCCESS) return (EQCLibStatus)status;
        state.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&state);
        state.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&state);
    }else{
        const EQCLibStatus status = EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&state,0.1*CellAccuracy);
        if(status != EQC_SUCCESS) return status;
    }
    /* Use the composition at the requested temperature to set the energy. */
    state.U = EQCLibConvertTtoU(&state);
    if(!isfinite(state.U) || state.U <= 0.0) return EQC_INVALID_INPUT;
    *initial = state;
    return EQC_SUCCESS;
}

static bool WriteCell(const StructEQCData *state, const char *mode, double years)
{
    const double nH = state->n[ipHYDROGEN], nHe = state->n[ipHELIUM];
    const double nC = state->n[ipCARBON], nO = state->n[ipOXYGEN];
    const double values[] = {
        years, state->T, state->U, state->nel/nH, state->natom[ipHYDROGEN]/nH,
        state->nion[ipHYDROGEN]/nH, 2.0*state->nmol[ipH2]/nH,
        state->natom[ipHELIUM]/nHe, state->nion[ipHELIUM]/nHe, state->nHeIII/nHe,
        state->natom[ipCARBON]/nC, state->nion[ipCARBON]/nC, state->nCIII/nC,
        state->natom[ipOXYGEN]/nO, state->nion[ipOXYGEN]/nO, state->nOIII/nO,
        state->H_tot, state->L_tot, state->L_CII158, state->L_OII, state->L_OIII,
    };
    if(printf("%s",mode) < 0) return false;
    for(size_t i=0;i<sizeof values/sizeof values[0];i++)
        if(printf(",%.12e",values[i]) < 0) return false;
    return putchar('\n') != EOF;
}

int main(int argc, char *argv[])
{
    if(argc > 3 || (argc == 2 && strcmp(argv[1],"--help") == 0)){
        fputs("Usage: eqclib-cell [CONFIG_FILE [PEGASE2_DIRECTORY]]\n",argc > 3 ? stderr : stdout);
        return argc > 3 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    CellSettings settings = {
        .chemistryMode = EQC_CHEMISTRY_EQUILIBRIUM,
        .nH = 10.0, .temperature = 8000.0, .timestepYears = 1000.0, .steps = 10.0,
        .ionizedHydrogen = 0.5,
    };
    if((argc > 1 && !ReadSettings(argv[1],&settings)) || !ValidateSettings(&settings)) return EXIT_FAILURE;

    EQCLibSetRunParamHIIRegion(settings.NLyn > 0.0);
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    const bool nonequilibrium = settings.chemistryMode == EQC_CHEMISTRY_NONEQUILIBRIUM;
    const bool useH2 = settings.Zd > 0.0 || (nonequilibrium && settings.molecularHydrogen > 0.0);
    EQCLibSetRunParamH2Formation(useH2,1.0);
    if(EQCLibTryInitialization(0,argc == 3 ? argv[2] : NULL) != EQC_SUCCESS) return EXIT_FAILURE;

    StructEQCData state;
    const char *mode = nonequilibrium ? "nonequilibrium" : "equilibrium";
    const EQCLibStatus initialStatus = InitializeCell(&settings,&state);
    if(initialStatus != EQC_SUCCESS){
        fprintf(stderr,"%s: could not initialize the cell at T=%g K (status %d)\n",
                mode,settings.temperature,initialStatus);
        return EXIT_FAILURE;
    }
    fputs("chemistry_mode,time_yr,T_K,U_erg_g,xe,HI,HII,H2,HeI,HeII,HeIII,CI,CII,CIII,OI,OII,OIII,"
          "H_tot,L_tot,L_CII158,L_OII,L_OIII\n",stdout);
    if(!WriteCell(&state,mode,0.0)) return EXIT_FAILURE;
    for(int step=0;step<(int)settings.steps;step++){
        const EQCLibStatus status = EQCLibTryAdvanceCell(&state,CellAccuracy,settings.chemistryMode);
        if(status != EQC_SUCCESS){
            fprintf(stderr,"%s: step %d failed (status %d)\n",mode,step+1,status);
            return EXIT_FAILURE;
        }
        if(!WriteCell(&state,mode,(step+1)*settings.timestepYears)) return EXIT_FAILURE;
    }
    return fflush(stdout) == EOF ? EXIT_FAILURE : EXIT_SUCCESS;
}
