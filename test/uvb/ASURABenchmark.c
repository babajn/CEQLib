#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CoolingHeatingRates.h"
#include "EquilibriumAbundance.h"

enum { WarmupSteps = 6, TimingBatches = 5, MaximumCalls = 1 << 18 };
static const double EnergyAccuracy = 1.e-4;
static const double BatchSeconds = 0.01;

static const struct {
    const char *name;
    bool enabled;
    double maximumEnergy;
} RadiationCases[] = {
    {"off",false,0.0},
    {"cut110",true,110.0},
    {"full",true,0.0}
};

static void RequireSuccess(const char *operation, EQCLibStatus status)
{
    if(status != EQC_SUCCESS){
        fprintf(stderr,"%s failed with status %d\n",operation,status);
        exit(EXIT_FAILURE);
    }
}

static StructEQCData InitialCell(double nH, double metallicity)
{
    StructEQCData cell = {
        .T = 8000.0, .dt = 1.e5*YEAR_CGS,
        .Nt = nH*100.0*PC_CGS, .Z = metallicity, .Zd = 0.01,
        .G0 = 1.0, .X0 = 1.0
    };
    cell.n[ipHYDROGEN] = nH;
    cell.n[ipHELIUM] = 0.1*nH;
    cell.n[ipCARBON] = 1.4e-4*nH;
    cell.n[ipOXYGEN] = 3.2e-4*nH;
    cell.nion[ipHYDROGEN] = 0.5*nH;
    for(int i=0;i<NumberChemicalElement;i++) cell.natom[i] = cell.n[i]-cell.nion[i];
    cell.nel = EQCLibGetElectronNumberDensity(&cell);
    cell.Rho = EQCLibGetTotalMassDensity(&cell);
    cell.MeanMolecularWeight = EQCLibGetMeanMolecularWeight(&cell);
    cell.DegreeOfFreedom = EQCLibGetMeanDegreeOfFreedom(&cell);
    cell.U = EQCLibConvertTtoU(&cell);
    return cell;
}

static void UpdateColumns(StructEQCData *cell)
{
    const double neutralFraction = (cell->natom[ipHYDROGEN]+2.0*cell->nmol[ipH2])/cell->n[ipHYDROGEN];
    cell->NHIuvb[0] = cell->NHIuvb[1] = cell->Nt*neutralFraction;
}

static double ProcessSeconds(void)
{
    struct timespec now;
    if(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&now) != 0){
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return (double)now.tv_sec + 1.e-9*(double)now.tv_nsec;
}

static double Measure(const StructEQCData *initial, EQCLibChemistryMode mode, unsigned int calls,
                      StructEQCData *result)
{
    const double start = ProcessSeconds();
    for(unsigned int i=0;i<calls;i++){
        *result = *initial;
        const EQCLibStatus status = EQCLibTryAdvanceCell(result,EnergyAccuracy,mode);
        if(status != EQC_SUCCESS) RequireSuccess("timed ASURA cell update",status);
    }
    return ProcessSeconds()-start;
}

static int CompareTimings(const void *a, const void *b)
{
    const double left = *(const double *)a, right = *(const double *)b;
    return (left > right)-(left < right);
}

static void Benchmark(FILE *output, const char *radiation, int z, const StructEQCData *initial,
                      EQCLibChemistryMode mode)
{
    const bool equilibrium = mode == EQC_CHEMISTRY_EQUILIBRIUM;
    StructEQCData result;
    unsigned int calls = 1;
    while(Measure(initial,mode,calls,&result) < BatchSeconds && calls < MaximumCalls)
        calls *= 2;
    double timing[TimingBatches];
    for(int i=0;i<TimingBatches;i++) timing[i] = Measure(initial,mode,calls,&result);
    qsort(timing,TimingBatches,sizeof timing[0],CompareTimings);

    CEQUVBEquilibriumStats stats = {0};
    if(equilibrium){
        StructEQCData profiled = *initial;
        RequireSuccess("profiled ASURA equilibrium update",
                       CEQAdvanceEquilibriumWithStats(&profiled,EnergyAccuracy,&stats));
        if(profiled.T != result.T || profiled.nel != result.nel || profiled.U != result.U ||
           profiled.H_UVB != result.H_UVB || stats.temperatureTrials == 0){
            fprintf(stderr,"Mismatch in profiled ASURA equilibrium at z=%d, nH=%.9g, UVB=%s\n",
                    z,initial->n[ipHYDROGEN],radiation);
            exit(EXIT_FAILURE);
        }
    }

    const double nH = result.n[ipHYDROGEN];
    fprintf(output,"%s,%d,%.9g,%s,%u,%.9g,%lu,%lu,%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            radiation,z,nH,equilibrium ? "equilibrium" : "nonequilibrium",
            calls,1.e6*timing[TimingBatches/2]/calls,
            (unsigned long)stats.temperatureTrials,(unsigned long)stats.chargeTrials,(unsigned long)stats.iterations,
            initial->T,initial->nel/nH,initial->NHIuvb[0],
            result.T,result.nel/nH,result.natom[ipHYDROGEN]/nH,result.nion[ipHYDROGEN]/nH,
            result.nmol[ipH2]/nH,result.H_UVB,result.Z,result.H_tot > 0.0 ? result.H_UVB/result.H_tot : NAN);
}

int main(int argc, char **argv)
{
    const bool zeroMetallicity = argc == 4 && strcmp(argv[3],"--zero-metallicity") == 0;
    const bool chargeTransfer = argc == 4 && strcmp(argv[3],"--charge-transfer") == 0;
    const bool bidirectional = argc == 4 && strcmp(argv[3],"--charge-transfer-bidirectional") == 0;
    const bool compareTransfer = chargeTransfer || bidirectional;
    if(argc != 3 && !zeroMetallicity && !compareTransfer){
        fprintf(stderr,"Usage: %s PEGASE2_DIRECTORY FG20_DIRECTORY [--zero-metallicity"
                       "|--charge-transfer|--charge-transfer-bidirectional]\n",argv[0]);
        return EXIT_FAILURE;
    }

    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(true);
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamHIIRegion(false);
    EQCLibSetRunParamChargeTransfer(true);
    EQCLibSetRunParamH2Formation(true,1.0);
    EQCLibSetRunParamCOFormation(false,1.0);
    EQCLibSetRunParamNonEquilibriumTemperatureTables(true);
    RequireSuccess("ASURA chemistry accuracy",EQCLibSetRunParamChemistryAccuracy(1.e-5));

    const char *path = bidirectional ? "asura-charge-transfer-bidirectional.csv" :
                       chargeTransfer ? "asura-charge-transfer.csv" : "asura-benchmark.csv";
    FILE *output = fopen(path,"w");
    if(output == NULL){
        perror(path);
        return EXIT_FAILURE;
    }
    fputs("uvb,z,nH,chemistry_mode,calls,cpu_us,temperature_trials,charge_trials,fixed_point_iterations,"
          "T_initial,xe_initial,N_per_side,T,xe,xHI,xHII,xH2,H_UVB,Z,H_UVB_over_H_tot\n",output);
    for(size_t radiation=0;radiation<sizeof RadiationCases/sizeof RadiationCases[0];radiation++){
        if(compareTransfer && strcmp(RadiationCases[radiation].name,"full") != 0) continue;
        const bool enabled = RadiationCases[radiation].enabled;
        RequireSuccess("ASURA radiation setup",EQCLibSetRunParamUVBackground(enabled,argv[2]));
        RequireSuccess("ASURA radiation energy cutoff",
                       EQCLibSetRunParamUVBackgroundMaximumEnergy(RadiationCases[radiation].maximumEnergy));
        RequireSuccess("ASURA initialization",EQCLibTryInitialization(0,argv[1]));
        for(int z=0;z<=2;z+=2){
            if(enabled) RequireSuccess("ASURA redshift",EQCLibTrySetUVBackgroundRedshift(z));
            for(int exponent=-6;exponent<=4;exponent++){
                if(compareTransfer){
                    /* Give each transfer treatment identical warmed input. */
                    EQCLibSetRunParamHIIRegion(false);
                    EQCLibSetRunParamChargeTransfer(false);
                }
                StructEQCData cell = InitialCell(pow(10.0,0.5*exponent),zeroMetallicity ? 0.0 : 0.0134);
                /* Relax one common starting state for both measured modes. */
                for(int step=0;step<WarmupSteps;step++){
                    UpdateColumns(&cell);
                    const EQCLibStatus status = EQCLibTryAdvanceCell(&cell,EnergyAccuracy,EQC_CHEMISTRY_EQUILIBRIUM);
                    if(status != EQC_SUCCESS){
                        fprintf(stderr,"ASURA warm-up at z=%d, nH=%.9g, UVB=%s, step=%d\n",
                                z,cell.n[ipHYDROGEN],RadiationCases[radiation].name,step+1);
                        RequireSuccess("ASURA warm-up",status);
                    }
                }
                UpdateColumns(&cell);
                if(compareTransfer){
                    EQCLibSetRunParamHIIRegion(bidirectional);
                    EQCLibSetRunParamChargeTransfer(true);
                }
                Benchmark(output,RadiationCases[radiation].name,z,&cell,EQC_CHEMISTRY_EQUILIBRIUM);
                Benchmark(output,RadiationCases[radiation].name,z,&cell,EQC_CHEMISTRY_NONEQUILIBRIUM);
            }
        }
    }
    if(fclose(output) != 0){
        perror(path);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
