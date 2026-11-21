#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ChemicalEquilibriumCooling.h"
#include "Constants.h"
#include "CompareBackend.h"

#ifdef EQC_H2_DIAGNOSTICS
#include "H2Diagnostics.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#endif
#ifdef EQC_THERMAL_DIAGNOSTICS
#include "FineStructureLines.h"
#include "ThermalDiagnostics.h"
#endif

struct Cell {
    StructEQCData gas;
    int lastMethod;   /* -1: equilibrium solve */
};

typedef struct {
    const char *name;
    EQCLibChemistryMode mode;
    double accuracy;
} Method;

/* Recommended ASURA setting (H2 on, CO off, H-O charge transfer on) with the
 * energy tolerance 1e-4, plus the non-equilibrium solver at two tolerances. */
static const Method DefaultMethods[] = {
    {"eqclib_eq", EQC_CHEMISTRY_EQUILIBRIUM, 1.0e-4},
    {"eqclib_noneq_1e-4", EQC_CHEMISTRY_NONEQUILIBRIUM, 1.0e-4},
    {"eqclib_noneq_1e-6", EQC_CHEMISTRY_NONEQUILIBRIUM, 1.0e-6},
};
static const Method COMethods[] = {
    {"eqclib_co_eq", EQC_CHEMISTRY_EQUILIBRIUM, 1.0e-4},
};
static const Method *Methods = DefaultMethods;
static int MethodCount = 3;
static const char *Name = "eqclib";

int BackendInit(const char *variant, const char *dataDirectory)
{
    bool co = false;
    double clumping = 1.0;
    if(variant == NULL || strcmp(variant, "default") == 0){
        co = false;
    }else if(strcmp(variant, "co") == 0){
        co = true;
        Methods = COMethods;
        MethodCount = 1;
        Name = "eqclib_co";
#ifdef EQC_H2_DIAGNOSTICS
    }else if(strcmp(variant, "form025") == 0){
        clumping = 0.25;
        Name = "eqclib_form025";
#endif

    }else{
        fprintf(stderr, "unknown EQCLib variant %s\n", variant);
        return 1;
    }
    EQCLibSetRunParamRecombinationCaseB(true);
    EQCLibSetRunParamDustAttenuation(true);
    EQCLibSetRunParamGrainSurfaceProcesses(true);
    EQCLibSetRunParamCosmicRayIonization(true);
    EQCLibSetRunParamSoftXrayIonization(false);   /* CHIMES has no X-ray background here */
    EQCLibSetRunParamCIEResonanceCooling(true);
    EQCLibSetRunParamHIIRegion(true);              /* NLyn = 0; enables H-O charge transfer */
    EQCLibSetRunParamH2Formation(true, clumping);
    EQCLibSetRunParamCOFormation(co, 1.0);
    EQCLibSetRunParamChargeTransfer(true);
    if(EQCLibSetRunParamChemistryAccuracy(1.0e-5) != EQC_SUCCESS) return 1;
    return EQCLibTryInitialization(0, dataDirectory) == EQC_SUCCESS ? 0 : 1;
}

const char *BackendName(void) { return Name; }
int BackendMethodCount(void) { return MethodCount; }
const char *BackendMethodName(int method) { return Methods[method].name; }

Cell *BackendNewCell(void)
{
    Cell *cell = calloc(1, sizeof *cell);
    if(cell == NULL){ perror("calloc"); exit(EXIT_FAILURE); }
    return cell;
}

void BackendCopyCell(Cell *dst, const Cell *src) { *dst = *src; }

static StructEQCData InitialGas(double nH, double T)
{
    /* Use the benchmark's gas metallicity for CIE and the original dust ratio. */
    StructEQCData s = {
        .T = T, .Nt = CmpColumn(nH), .Z = CMP_GAS_Z, .Zd = SOLAR_METALLICITY, .G0 = CMP_G0, .X0 = 1.0,
        .Td = CMP_DUST_TEMP, .nel = 1.0e-3*nH,
        .n = {nH, CMP_HE*nH, CMP_C*nH, CMP_O*nH, CMP_SI*nH},
    };
    for(int i = 0; i < NumberChemicalElement; i++){
        s.natom[i] = s.n[i];
        s.Rho += s.n[i]*ChemicalElementMass[i];
    }
    return s;
}

int BackendEquilibrium(Cell *cell, double nH, double T, int *iterations)
{
    cell->gas = InitialGas(nH, T);
    cell->lastMethod = -1;
    *iterations = 1;
    if(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&cell->gas, 1.0e-9) != EQC_SUCCESS)
        return 1;
    cell->gas.U = EQCLibConvertTtoU(&cell->gas);
    return 0;
}

int BackendAdvance(Cell *cell, double dt, int method)
{
    cell->gas.dt = dt;
    cell->lastMethod = method;
    return EQCLibTryAdvanceCell(&cell->gas, Methods[method].accuracy, Methods[method].mode) == EQC_SUCCESS ? 0 : 1;
}

static double Ratio(double a, double b) { return b > 0.0 ? a/b : NAN; }

void BackendSnapshot(Cell *cell, Snapshot *out, int computeRates)
{
    const StructEQCData *s = &cell->gas;
    const double nH = s->n[ipHYDROGEN];
    *out = (Snapshot){
        .T = s->T, .xe = s->nel/nH, .xHII = s->nion[ipHYDROGEN]/nH,
        .fH2 = 2.0*s->nmol[ipH2]/nH,
        .xHeII = Ratio(s->nion[ipHELIUM], s->n[ipHELIUM]), .xHeIII = Ratio(s->nHeIII, s->n[ipHELIUM]),
        .xCII = Ratio(s->nion[ipCARBON], s->n[ipCARBON]), .xCIII = Ratio(s->nCIII, s->n[ipCARBON]),
        .xCO = Ratio(s->nmol[ipCO], s->n[ipCARBON]),
        .xOII = Ratio(s->nion[ipOXYGEN], s->n[ipOXYGEN]), .xOIII = Ratio(s->nOIII, s->n[ipOXYGEN]),
        .xSiII = Ratio(s->nion[ipSILICON], s->n[ipSILICON]),
        .heat = NAN, .cool = NAN,
    };
    if(!computeRates) return;
    if(cell->lastMethod < 0 || Methods[cell->lastMethod].mode == EQC_CHEMISTRY_NONEQUILIBRIUM){
        out->heat = s->H_tot;
        out->cool = s->L_tot;
    }else{
        /* The equilibrium mode re-equilibrates at the final temperature. */
        StructEQCData copy = *s;
        if(EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates(&copy, 1.0e-9) == EQC_SUCCESS){
            out->heat = copy.H_tot;
            out->cool = copy.L_tot;
        }
    }
}

#ifdef EQC_H2_DIAGNOSTICS
int BackendH2Rates(Cell *cell, double nH, double T, double fH2, H2Rates *rates)
{
    const bool shared = fH2 >= 0.0;
    StructEQCData gas;
    if(shared){
        gas = InitialGas(nH,T);
        gas.nmol[ipH2] = 0.5*fH2*nH;
        gas.natom[ipHYDROGEN] = (1.0-fH2)*nH;
        gas.natom[ipCARBON] = 0.0;
        gas.nion[ipCARBON] = gas.n[ipCARBON];
        gas.nel = gas.n[ipCARBON];
        gas.invT = 1.0/T;
        gas.logT = log10(T);
        const int status = EQCLibCalculateCoolingHeatingRates(&gas);
        if(status != EQC_SUCCESS) return status;
    }else{
        gas = cell->gas;
    }

    const double nHI = gas.natom[ipHYDROGEN], nH2 = gas.nmol[ipH2];
    const double kForm = H2FormationRate(gas.T,gas.Zd,EffectiveDustTemperature(&gas));
    const double dust = DustShieldingFactorFUV(gas.Nt,gas.Zd);
    const double nH2Column = gas.Nt*nH2/nH;
    const double b5 = H2DopplerWidth(gas.T);
    const double self = H2SelfShieldingFactor(nH2Column,b5);
    const double kPhoto0 = H2PhotoDissociationRateWithDustFactor(gas.G0,0.0,b5,1.0);
    const double coll = nH2*H2CollisionalDissociationRate(gas.T,nH,nHI,nH2,gas.nel);
    const double dissCooling = CoolingH2CollisionalDissociation(&gas);
    *rates = (H2Rates){
        .fH2 = 2.0*nH2/nH, .xe = gas.nel/nH, .xHII = gas.nion[ipHYDROGEN]/nH,
        .form_grain = kForm*nH*nHI, .form_other = 0.0,
        .destroy_photo = nH2*kPhoto0*dust*self, .destroy_coll = coll, .destroy_other = 0.0,
        .k_form = kForm, .k_photo0 = kPhoto0,
        .dust_shield = dust, .self_shield = self, .shield_other = 1.0,
        .cool_h2 = fmax(0.0,gas.L_H2-dissCooling), .cool_h2_diss = dissCooling,
        .heat_h2 = HeatingH2Photodissociate(&gas)+HeatingH2PhotomPumping(&gas)+HeatingH2Formation(&gas),
        .heat = gas.H_tot, .cool = gas.L_tot, .energy = shared ? NAN : gas.U*gas.Rho,
    };
    return 0;
}
#endif

#ifdef EQC_THERMAL_DIAGNOSTICS
int BackendThermalRates(const Cell *cell, double nH, double T,
                        const ThermalComposition *fixed, ThermalRates *rates)
{
    StructEQCData gas = fixed ? InitialGas(nH,T) : cell->gas;
    if(fixed){
        gas.nmol[ipH2] = 0.5*fixed->fH2*nH;
        gas.nion[ipHYDROGEN] = fixed->xHII*nH;
        gas.natom[ipHYDROGEN] = fmax(0.0,1.0-fixed->fH2-fixed->xHII)*nH;
        gas.nel = fixed->xe*nH;
        gas.nion[ipCARBON] = fixed->xCII*gas.n[ipCARBON];
        gas.natom[ipCARBON] = (1.0-fixed->xCII)*gas.n[ipCARBON];
        gas.nion[ipOXYGEN] = (1.0-fixed->xOI)*gas.n[ipOXYGEN];
        gas.natom[ipOXYGEN] = fixed->xOI*gas.n[ipOXYGEN];
        gas.nion[ipSILICON] = fixed->xSiII*gas.n[ipSILICON];
        gas.natom[ipSILICON] = (1.0-fixed->xSiII)*gas.n[ipSILICON];
        gas.invT = 1.0/T;
        gas.logT = log10(T);
        const int status = EQCLibCalculateCoolingHeatingRates(&gas);
        if(status != EQC_SUCCESS) return status;
    }
    const double g0 = gas.G0*DustShieldingFactor(gas.Nt,gas.Zd,2.5);
    const double h2Dissociation = CoolingH2CollisionalDissociation(&gas);
    *rates = (ThermalRates){
        .composition = {
            .fH2 = 2.0*gas.nmol[ipH2]/nH, .xe = gas.nel/nH, .xHII = gas.nion[ipHYDROGEN]/nH,
            .xCII = gas.nion[ipCARBON]/gas.n[ipCARBON],
            .xOI = gas.natom[ipOXYGEN]/gas.n[ipOXYGEN],
            .xSiII = gas.nion[ipSILICON]/gas.n[ipSILICON],
        },
        .xHI = gas.natom[ipHYDROGEN]/nH,
        .pe_g0 = g0, .pe_psi = gas.nel > 0.0 ? g0*sqrt(gas.T)/gas.nel : NAN,
        .heat_pe = gas.H_PE, .heat_cr = gas.H_CR, .heat_h2 = gas.H_H2, .heat_photo = gas.H_PI,
        .heat_other = gas.H_XR, .heat = gas.H_tot,
        .cool_h2 = gas.L_H2-h2Dissociation, .cool_h2_diss = h2Dissociation,
        .cool_cii = gas.L_CII158, .cool_oi = gas.L_OI63+gas.L_OI146, .cool_ci = gas.L_CI610+gas.L_CI370,
        .cool_si2 = gas.L_SiII35, .cool_co = gas.L_CO, .cool_dust_rec = gas.L_DustRec,
        .cool_gas_dust = 0.0, .cool = gas.L_tot,
    };
    rates->cool_other = gas.L_tot - (gas.L_H2+rates->cool_cii+rates->cool_oi+rates->cool_ci
                                     +rates->cool_si2+rates->cool_co+rates->cool_dust_rec);
    gas.natom[ipHYDROGEN] = gas.nmol[ipH2] = 0.0;
    rates->cool_cii_e_only = CoolingCII158um(&gas);
    return 0;
}
#endif
