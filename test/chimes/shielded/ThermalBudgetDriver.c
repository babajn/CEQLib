#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ThermalDiagnostics.h"

static void WriteBudget(FILE *out, const char *sample, double nH, double T, int status,
                        const ThermalRates *r)
{
    const ThermalComposition *s = &r->composition;
    const double values[] = {
        s->fH2,s->xe,r->xHI,s->xHII,s->xCII,s->xOI,s->xSiII,r->pe_g0,r->pe_psi,
        r->heat_pe,r->heat_cr,r->heat_h2,r->heat_photo,r->heat_other,r->heat,
        r->cool_h2,r->cool_h2_diss,r->cool_cii,r->cool_oi,r->cool_ci,r->cool_si2,r->cool_co,
        r->cool_dust_rec,r->cool_gas_dust,r->cool_other,r->cool,r->cool_cii_e_only,
    };
    fprintf(out,"%s,%s,%.9g,%.9g,%d",BackendName(),sample,nH,T,status);
    for(size_t i=0;i<sizeof values/sizeof values[0];i++)
        fprintf(out,",%.12g",values[i]);
    fputc('\n',out);
}

static void SampleEquilibrium(FILE *out)
{
    static const double densities[] = {3.0,10.0,30.0,100.0};
    static const double additionalTemperatures[] = {30.0,89.0,129.0,211.0,300.0,3000.0,5000.0};
    const int additionalCount = sizeof additionalTemperatures/sizeof additionalTemperatures[0];
    Cell *cell = BackendNewCell();
    for(size_t i=0;i<sizeof densities/sizeof densities[0];i++){
        const double nH = densities[i];
        for(int j=-additionalCount;j<=100;j++){
            const double T = j < 0 ? additionalTemperatures[j+additionalCount] : pow(10.0,1.5+0.025*j);
            int iterations;
            ThermalRates rates = {0};
            int status = BackendEquilibrium(cell,nH,T,&iterations);
            if(status == 0) status = BackendThermalRates(cell,nH,T,NULL,&rates);
            WriteBudget(out,"native",nH,T,status,&rates);
        }
    }
}

static int SampleFixedCompositions(FILE *out, const char *requestFile)
{
    FILE *requests = fopen(requestFile,"r");
    if(requests == NULL){
        perror(requestFile);
        return 1;
    }
    char sample[32];
    double nH, T;
    ThermalComposition state;
    int fields;
    while((fields = fscanf(requests,"%31s %lf %lf %lf %lf %lf %lf %lf %lf",
                            sample,&nH,&T,&state.fH2,&state.xe,&state.xHII,
                            &state.xCII,&state.xOI,&state.xSiII)) == 9){
        ThermalRates rates = {0};
        const int status = BackendThermalRates(NULL,nH,T,&state,&rates);
        WriteBudget(out,sample,nH,T,status,&rates);
    }
    const int malformed = fields != EOF || ferror(requests);
    if(malformed) fprintf(stderr,"Could not read the temperature and composition from %s\n",requestFile);
    fclose(requests);
    return malformed;
}

int main(int argc, char **argv)
{
    if(argc != 4 && argc != 5){
        fprintf(stderr,"Usage: %s VARIANT DATA_DIRECTORY OUTPUT.csv [FIXED_COMPOSITIONS]\n",argv[0]);
        return 2;
    }
    if(BackendInit(argv[1],argv[2]) != 0) return 1;
    FILE *out = fopen(argv[3],"w");
    if(out == NULL){
        perror(argv[3]);
        return 1;
    }
    fputs("code,sample,nH,T,status,fH2,xe,xHI,xHII,xCII,xOI,xSiII,pe_g0,pe_psi,"
          "heat_pe,heat_cr,heat_h2,heat_photo,heat_other,heat,"
          "cool_h2,cool_h2_diss,cool_cii,cool_oi,cool_ci,cool_si2,cool_co,cool_dust_rec,cool_gas_dust,"
          "cool_other,cool,cool_cii_e_only\n",out);
    int status = 0;
    if(argc == 5){
        status = SampleFixedCompositions(out,argv[4]);
    }else{
        SampleEquilibrium(out);
    }
    if(fclose(out) != 0){
        perror(argv[3]);
        return 1;
    }
    return status;
}
