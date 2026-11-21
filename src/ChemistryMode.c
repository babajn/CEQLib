#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"

EQCLibStatus EQCLibTryParseChemistryMode(const char *value, EQCLibChemistryMode *mode)
{
    if(!value || !mode) return EQC_INVALID_INPUT;
    if(strcmp(value,"equilibrium") == 0) *mode = EQC_CHEMISTRY_EQUILIBRIUM;
    else if(strcmp(value,"nonequilibrium") == 0) *mode = EQC_CHEMISTRY_NONEQUILIBRIUM;
    else return EQC_INVALID_INPUT;
    return EQC_SUCCESS;
}

static char *SkipSpace(char *text)
{
    while(isspace((unsigned char)*text)) text++;
    return text;
}

EQCLibStatus EQCLibTryLoadChemistryMode(const char *path, EQCLibChemistryMode *mode)
{
    if(!mode) return EQC_INVALID_INPUT;
    if(!path){
        *mode = EQC_CHEMISTRY_EQUILIBRIUM;
        return EQC_SUCCESS;
    }
    FILE *file = fopen(path,"r");
    if(!file){
        perror(path);
        return EQC_INITIALIZATION_FAILED;
    }

    const char key[] = "chemistry_mode";
    const size_t keyLength = sizeof key-1;
    char *line = NULL;
    size_t capacity = 0, lineNumber = 0;
    bool seen = false;
    EQCLibStatus status = EQC_SUCCESS;
    EQCLibChemistryMode selected = EQC_CHEMISTRY_EQUILIBRIUM;
    while(getline(&line,&capacity,file) >= 0){
        lineNumber++;
        char *comment = strchr(line,'#');
        if(comment) *comment = '\0';
        char *entry = SkipSpace(line);
        if(strncmp(entry,key,keyLength) != 0) continue;
        const unsigned char separator = entry[keyLength];
        if(separator && separator != '=' && !isspace(separator)) continue;

        char *equals = SkipSpace(entry+keyLength);
        if(*equals != '=' || seen){
            status = EQC_INVALID_INPUT;
            break;
        }
        char *value = SkipSpace(equals+1);
        char *end = value+strlen(value);
        while(end > value && isspace((unsigned char)end[-1])) end--;
        *end = '\0';
        status = EQCLibTryParseChemistryMode(value,&selected);
        if(status != EQC_SUCCESS) break;
        seen = true;
    }
    if(status == EQC_INVALID_INPUT){
        fprintf(stderr,"%s:%zu: set chemistry_mode once to equilibrium or nonequilibrium\n",path,lineNumber);
    }else if(!feof(file)){
        perror(path);
        status = EQC_INITIALIZATION_FAILED;
    }
    free(line);
    if(fclose(file) != 0){
        perror(path);
        status = EQC_INITIALIZATION_FAILED;
    }
    if(status == EQC_SUCCESS) *mode = selected;
    return status;
}

EQCLibStatus EQCLibTryAdvanceCell(StructEQCData *state, double accuracy, EQCLibChemistryMode mode)
{
    if(!state || !isfinite(state->dt) || state->dt < 0.0 ||
       !isfinite(accuracy) || accuracy <= 0.0 || accuracy >= 1.0 ||
       (mode != EQC_CHEMISTRY_EQUILIBRIUM && mode != EQC_CHEMISTRY_NONEQUILIBRIUM)) return EQC_INVALID_INPUT;
    if(state->dt == 0.0) return EQC_SUCCESS;

    StructEQCData updated = *state;
    if(mode == EQC_CHEMISTRY_EQUILIBRIUM){
        const EQCLibStatus status =
            EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating(&updated,accuracy,&updated.Du);
        if(status != EQC_SUCCESS) return status;
        updated.U += state->dt*updated.Du;
    }else{
        const EQCLibStatus status = EQCLibTryAdvanceHydrogenThermochemistry(&updated,accuracy);
        if(status != EQC_SUCCESS) return status;
        updated.Du = (updated.U-state->U)/state->dt;
    }
    if(!isfinite(updated.U) || updated.U <= 0.0 || !isfinite(updated.Du)) return EQC_NO_CONVERGENCE;
    *state = updated;
    return EQC_SUCCESS;
}
