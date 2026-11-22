#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "UVBackground.h"

enum { TokenSize = 128 };
static const double PhotonEnergyEVAngstrom = 12398.419843320026;

/* Read one whitespace-delimited token, skipping Cloudy comment lines. */
static int ReadToken(FILE *stream, char token[TokenSize])
{
    int c;
    do {
        c = fgetc(stream);
        if(c == '#'){
            do c = fgetc(stream); while(c != EOF && c != '\n');
        }
    } while(c != EOF && isspace(c));
    if(c == EOF) return ferror(stream) ? -1 : 0;

    size_t length = 0;
    while(c != EOF && !isspace(c) && c != '#'){
        if(c == '\0' || length + 1 == TokenSize) return -1;
        token[length++] = (char)c;
        c = fgetc(stream);
    }
    if(c == '#' && ungetc(c,stream) == EOF) return -1;
    token[length] = '\0';
    return 1;
}

static bool ReadLabel(FILE *stream, const char *label)
{
    char token[TokenSize];
    return ReadToken(stream,token) == 1 && strcmp(token,label) == 0;
}

static bool ReadSize(FILE *stream, size_t *value)
{
    char token[TokenSize], *end;
    if(ReadToken(stream,token) != 1 || !isdigit((unsigned char)token[0])) return false;
    errno = 0;
    const unsigned long long parsed = strtoull(token,&end,10);
    if(*end != '\0' || errno != 0 || parsed > SIZE_MAX) return false;
    *value = (size_t)parsed;
    return true;
}

static bool ReadNumber(FILE *stream, double *value)
{
    char token[TokenSize], *end;
    if(ReadToken(stream,token) != 1) return false;
    const double parsed = strtod(token,&end);
    if(end == token || *end != '\0' || !isfinite(parsed) || parsed < 0.0) return false;
    *value = parsed;
    return true;
}

static bool ReadCloudySpectrum(FILE *stream, CEQUVBSpectrum *spectrum)
{
    size_t date, one1, one2, nz, ne;
    double wavelengthScale, intensityScale;
    if(!ReadSize(stream,&date) || date < 10000101 || date > 99991231 ||
       !ReadSize(stream,&one1) || one1 != 1 || !ReadSize(stream,&one2) || one2 != 1 ||
       !ReadLabel(stream,"z") || !ReadSize(stream,&nz) || !ReadSize(stream,&ne) ||
       nz == 0 || ne < 2 || ne > SIZE_MAX / sizeof(double) / nz ||
       !ReadLabel(stream,"lambda") || !ReadNumber(stream,&wavelengthScale) || wavelengthScale == 0.0 ||
       !ReadLabel(stream,"F_nu") || !ReadNumber(stream,&intensityScale) || intensityScale == 0.0)
        return false;

    spectrum->redshifts = malloc(nz*sizeof(double));
    spectrum->energyEV = malloc(ne*sizeof(double));
    spectrum->intensity = malloc(nz*ne*sizeof(double));
    if(!spectrum->redshifts || !spectrum->energyEV || !spectrum->intensity) return false;
    spectrum->redshiftCount = nz;
    spectrum->energyCount = ne;

    for(size_t i=0;i<nz;i++){
        if(!ReadNumber(stream,&spectrum->redshifts[i]) ||
           (i > 0 && spectrum->redshifts[i] <= spectrum->redshifts[i-1])) return false;
    }
    for(size_t i=0;i<ne;i++){
        double wavelength;
        if(!ReadNumber(stream,&wavelength) || wavelength == 0.0) return false;
        const double energy = PhotonEnergyEVAngstrom/(wavelength*wavelengthScale);
        if(!isfinite(energy) || energy == 0.0 ||
           (i > 0 && energy >= spectrum->energyEV[ne-i])) return false;
        spectrum->energyEV[ne-1-i] = energy;
    }
    for(size_t iz=0;iz<nz;iz++){
        for(size_t i=0;i<ne;i++){
            double intensity;
            if(!ReadNumber(stream,&intensity)) return false;
            intensity *= intensityScale;
            if(!isfinite(intensity)) return false;
            spectrum->intensity[iz*ne+ne-1-i] = intensity;
        }
    }
    char extra[TokenSize];
    return ReadToken(stream,extra) == 0;
}

EQCLibStatus CEQReadUVBackgroundSpectrum(const char *directory, CEQUVBSpectrum *spectrum)
{
    static const char file[] = "fg20_cloudy.ascii";
    const size_t length = strlen(directory);
    if(length > SIZE_MAX - sizeof(file) - 1) return EQC_INITIALIZATION_FAILED;
    const size_t size = length + sizeof(file) + 1;
    char *path = malloc(size);
    if(!path) return EQC_INITIALIZATION_FAILED;
    snprintf(path,size,"%s/%s",directory,file);
    FILE *stream = fopen(path,"r");
    if(!stream){
        perror(path);
        free(path);
        return EQC_INITIALIZATION_FAILED;
    }

    CEQUVBSpectrum parsed = {0};
    const bool read = ReadCloudySpectrum(stream,&parsed);
    const bool closed = fclose(stream) == 0;
    if(!read || !closed){
        fprintf(stderr,"%s: invalid UV background spectrum or read failure\n",path);
        free(parsed.redshifts);
        free(parsed.energyEV);
        free(parsed.intensity);
    }else{
        *spectrum = parsed;
    }
    free(path);
    return read && closed ? EQC_SUCCESS : EQC_INITIALIZATION_FAILED;
}
