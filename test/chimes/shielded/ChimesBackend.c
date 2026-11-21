#ifdef EQC_H2_DIAGNOSTICS
#define _DEFAULT_SOURCE    /* CHIMES uses M_LN10 in its interpolation helpers. */
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chimes_proto.h"
#include "CompareBackend.h"

#ifdef EQC_H2_DIAGNOSTICS
#include "chimes_vars.h"
#include "H2Diagnostics.h"
#ifdef EQC_THERMAL_DIAGNOSTICS
#include "ThermalDiagnostics.h"
#endif
#endif

#define CHIMES_LIGHT_SPEED 2.99792458e10
#ifndef CMP_CHIMES_DOPPLER_B
#define CMP_CHIMES_DOPPLER_B 7.1
#endif

struct Cell {
    struct gasVariables gas;
};

static struct globalVariables Global;
static double ElementAbundances[10];
static double PhotonDensityScale = 1.0;
static const char *Name = "chimes";

int BackendInit(const char *variant, const char *dataDirectory)
{
    /* He, C, N, O, Ne, Mg, Si, S, Ca, Fe relative to H.  He, C, O and Si use
     * the EQCLib gas-phase values in every variant; the "full" variant adds
     * the remaining CHIMES metals at the Wiersma et al. (2009) solar values. */
    const double massFraction[10] = {0.2806, 2.07e-3, 8.36e-4, 5.49e-3, 1.41e-3,
                                     5.91e-4, 6.83e-4, 4.09e-4, 6.44e-5, 1.1e-3};
    const double atomicMass[10] = {4., 12., 14., 16., 20., 24., 28., 32., 40., 56.};
    const double XH = 1.0 - 0.0129 - 0.2806;
    int include[9] = {1, 0, 1, 0, 0, 1, 0, 0, 0};   /* C N O Ne Mg Si S Ca Fe */
    /* EQCLib runs with NLyn = 0, so by default the > 13.6 eV part of the B87
     * spectrum is removed; the "-euv" variants keep it as a control. */
    const size_t length = strlen(variant);
    const int keepEUV = length > 4 && strcmp(variant+length-4, "-euv") == 0;
    if(strncmp(variant, "matched", 7) == 0 && length == (size_t)(7+4*keepEUV)){
        Name = keepEUV ? "chimes_matched_euv" : "chimes_matched";
    }else if(strncmp(variant, "full", 4) == 0 && length == (size_t)(4+4*keepEUV)){
        Name = keepEUV ? "chimes_full_euv" : "chimes_full";
        for(int i = 0; i < 9; i++) include[i] = 1;
    }else{
        fprintf(stderr, "unknown CHIMES variant %s\n", variant);
        return 1;
    }
    for(int i = 0; i < 10; i++) ElementAbundances[i] = massFraction[i]/(atomicMass[i]*XH);
    ElementAbundances[0] = CMP_HE;
    ElementAbundances[1] = CMP_C;
    ElementAbundances[3] = CMP_O;
    ElementAbundances[6] = CMP_SI;
    for(int i = 0; i < 9; i++)
        if(!include[i]) ElementAbundances[i+1] = 0.0;

    memset(&Global, 0, sizeof Global);
    snprintf(Global.MainDataTablePath, sizeof Global.MainDataTablePath, "%s/chimes_main_data.hdf5", dataDirectory);
    snprintf(Global.PhotoIonTablePath[0], sizeof Global.PhotoIonTablePath[0], "%s/cross_sections_B87.hdf5", dataDirectory);
    snprintf(Global.EqAbundanceTablePath, sizeof Global.EqAbundanceTablePath, "None");
    Global.cellSelfShieldingOn = 1;
    Global.N_spectra = 1;
    Global.rt_update_flux = 0;
    Global.rt_use_on_the_spot_approx = 0;
    Global.redshift_dependent_UVB_index = -1;
    Global.use_redshift_dependent_eqm_tables = 0;
    Global.redshift = 0.0;
    Global.reionisation_redshift = 0.0;
    Global.StaticMolCooling = 1;
    Global.Tmol_K = 1.0e5;
    Global.grain_temperature = CMP_DUST_TEMP;
    Global.cmb_temperature = 2.725;
    Global.relativeTolerance = 1.0e-4;
    Global.absoluteTolerance = 1.0e-10;
    Global.explicitTolerance = 0.01;
    Global.rt_density_absoluteTolerance = 1.0e-10;
    Global.rt_flux_absoluteTolerance = 1.0e-10;
    for(int i = 0; i < 9; i++) Global.element_included[i] = include[i];
    Global.scale_metal_tolerances = 1;
    Global.chimes_debug = 0;
    Global.hybrid_cooling_mode = 0;
    init_chimes(&Global);

    /* Scale the B87 spectrum to G0 = CMP_G0 Habing in the 6-13.6 eV band. */
    const double G0table = chimes_table_spectra.isotropic_photon_density[0]*CHIMES_LIGHT_SPEED*
                           chimes_table_spectra.G0_parameter[0];
    PhotonDensityScale = CMP_G0/G0table;
    double rateHI = 0.0;
    for(int i = 0; i < chimes_table_photoion_euv.N_reactions[1]; i++)
        if(chimes_table_photoion_euv.reactants[i] == Global.speciesIndices[sp_HI])
            rateHI = chimes_table_spectra.isotropic_photon_density[0]*PhotonDensityScale*CHIMES_LIGHT_SPEED*
                     chimes_table_photoion_euv.sigmaPhot[i];
    fprintf(stderr, "%s: %d species, B87 table G0 = %.4f Habing, scaled by %.4f, "
            "unshielded HI photoionization %.3e s^-1 (%s)\n", Name, Global.totalNumberOfSpecies,
            G0table, PhotonDensityScale, rateHI, keepEUV ? "kept" : "removed");
    if(!keepEUV){
        for(int i = 0; i < Global.N_spectra*chimes_table_photoion_euv.N_reactions[1]; i++)
            chimes_table_photoion_euv.sigmaPhot[i] = 0.0;
        for(int i = 0; i < Global.N_spectra*chimes_table_photoion_auger_euv.N_reactions[1]; i++)
            chimes_table_photoion_auger_euv.sigmaPhot[i] = 0.0;
    }
    return 0;
}

const char *BackendName(void) { return Name; }
int BackendMethodCount(void) { return 1; }
const char *BackendMethodName(int method) { (void)method; return Name; }

Cell *BackendNewCell(void)
{
    Cell *cell = calloc(1, sizeof *cell);
    if(cell == NULL){ perror("calloc"); exit(EXIT_FAILURE); }
    allocate_gas_abundances_memory(&cell->gas, &Global);
    return cell;
}

void BackendCopyCell(Cell *dst, const Cell *src)
{
    ChimesFloat *abundances = dst->gas.abundances;
    ChimesFloat *photon = dst->gas.isotropic_photon_density;
    ChimesFloat *G0 = dst->gas.G0_parameter;
    ChimesFloat *H2 = dst->gas.H2_dissocJ;
    dst->gas = src->gas;
    dst->gas.abundances = abundances;
    dst->gas.isotropic_photon_density = photon;
    dst->gas.G0_parameter = G0;
    dst->gas.H2_dissocJ = H2;
    memcpy(abundances, src->gas.abundances, Global.totalNumberOfSpecies*sizeof *abundances);
    photon[0] = src->gas.isotropic_photon_density[0];
    G0[0] = src->gas.G0_parameter[0];
    H2[0] = src->gas.H2_dissocJ[0];
}

static void SetGas(struct gasVariables *g, double nH, double T)
{
    for(int i = 0; i < 10; i++) g->element_abundances[i] = ElementAbundances[i];
    g->nH_tot = nH;
    g->temperature = T;
    g->TempFloor = 10.1;
    g->divVel = 1.0;
    g->doppler_broad = CMP_CHIMES_DOPPLER_B;
    g->isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0]*PhotonDensityScale;
    g->G0_parameter[0] = chimes_table_spectra.G0_parameter[0];
    g->H2_dissocJ[0] = chimes_table_spectra.H2_dissocJ[0];
    g->cr_rate = CMP_CR_RATE;
    g->metallicity = 1.0;
    g->dust_ratio = 1.0;
    g->dust_boost_factor = 1.0;
    g->cell_size = CmpColumn(nH)/nH;
    g->ForceEqOn = 0;
    g->ThermEvolOn = 0;
    g->temp_floor_mode = 0;
    g->InitIonState = T >= 1.5e4 ? 1 : 0;
    g->constant_heating_rate = 0.0;
    g->light_speed_reduction_factor = 1.0;
}

int BackendEquilibrium(Cell *cell, double nH, double T, int *iterations)
{
    struct gasVariables *g = &cell->gas;
    const int n = Global.totalNumberOfSpecies;
    ChimesFloat previous[CHIMES_TOTSIZE];
    SetGas(g, nH, T);
    initialise_gas_abundances(g, &Global);
    g->hydro_timestep = 1.0e9*CMP_YEAR;
    /* Integrate at fixed T in 1 Gyr chunks until the composition is stationary. */
    for(int it = 1; it <= 40; it++){
        memcpy(previous, g->abundances, n*sizeof *previous);
        if(chimes_network(g, &Global) != 0) return 1;
        double change = 0.0;
        for(int i = 0; i < n; i++){
            if(fabs(g->abundances[i]) < 1.0e-12 && fabs(previous[i]) < 1.0e-12) continue;
            change = fmax(change, fabs(g->abundances[i]-previous[i])/fmax(fabs(previous[i]), 1.0e-12));
        }
        *iterations = it;
        if(it >= 2 && change < 1.0e-3) return 0;
    }
    return 2;   /* not stationary after 40 Gyr */
}

int BackendAdvance(Cell *cell, double dt, int method)
{
    (void)method;
    cell->gas.ThermEvolOn = 1;
    cell->gas.hydro_timestep = dt;
    const int status = chimes_network(&cell->gas, &Global);
    cell->gas.ThermEvolOn = 0;
    return status;
}

static double Abundance(const struct gasVariables *g, int species)
{
    const int index = Global.speciesIndices[species];
    return index >= 0 ? g->abundances[index] : NAN;
}

void BackendSnapshot(Cell *cell, Snapshot *out, int computeRates)
{
    const struct gasVariables *g = &cell->gas;
    const double He = g->element_abundances[0], C = g->element_abundances[1];
    const double O = g->element_abundances[3], Si = g->element_abundances[6];
    *out = (Snapshot){
        .T = g->temperature, .xe = Abundance(g, sp_elec), .xHII = Abundance(g, sp_HII),
        .fH2 = 2.0*Abundance(g, sp_H2),
        .xHeII = Abundance(g, sp_HeII)/He, .xHeIII = Abundance(g, sp_HeIII)/He,
        .xCII = Abundance(g, sp_CII)/C, .xCIII = Abundance(g, sp_CIII)/C,
        .xCO = Abundance(g, sp_CO)/C,
        .xOII = Abundance(g, sp_OII)/O, .xOIII = Abundance(g, sp_OIII)/O,
        .xSiII = Abundance(g, sp_SiII)/Si,
        .heat = NAN, .cool = NAN,
    };
    if(computeRates){
        /* CHIMES fills the photoheating coefficients only when ThermEvolOn is set. */
        const int thermal = cell->gas.ThermEvolOn;
        cell->gas.ThermEvolOn = 1;
        out->cool = chimes_cooling_rate_external(&cell->gas, &Global, 1);
        out->heat = chimes_cooling_rate_external(&cell->gas, &Global, 2);
        cell->gas.ThermEvolOn = thermal;
    }
}

#ifdef EQC_H2_DIAGNOSTICS
static double MolecularColumn(const struct gasVariables *gas, int species, double column)
{
    const int index = Global.speciesIndices[species];
    return index < 0 ? 0.0 : fmax(gas->abundances[index],0.0)*column;
}

static struct UserData PrepareDiagnosticRates(struct gasVariables *gas, struct Species_Structure *species,
                                              struct chimes_current_rates_struct *native)
{
    gas->ThermEvolOn = 1;  /* Initialize the native photoheating coefficients. */
    memset(species,0,Global.totalNumberOfSpecies*sizeof *species);
    int networkSize, atomicNetworkSize;
    set_species_structures(species,gas,&networkSize,&atomicNetworkSize,&Global);
    const double column = Global.cellSelfShieldingOn ? gas->cell_size*gas->nH_tot : 0.0;
    allocate_current_rates_memory(native,&Global);
    struct UserData data = {
        .myGasVars = gas, .myGlobalVars = &Global, .species = species,
        .chimes_current_rates = native, .network_size = networkSize, .mol_flag_index = 1,
        .HI_column = Abundance(gas,sp_HI)*column, .H2_column = Abundance(gas,sp_H2)*column,
        .HeI_column = Abundance(gas,sp_HeI)*column, .HeII_column = Abundance(gas,sp_HeII)*column,
        .CO_column = MolecularColumn(gas,sp_CO,column), .H2O_column = MolecularColumn(gas,sp_H2O,column),
        .OH_column = MolecularColumn(gas,sp_OH,column), .extinction = (double)DUSTEFFSIZE*column*gas->dust_ratio,
    };
    data.case_AB_index[0] = CASE_AB_CROSS_SECTION_HYDROGEN*data.HI_column >= 1.0;
    data.case_AB_index[1] = CASE_AB_CROSS_SECTION_HELIUM*data.HeI_column >= 1.0;
    set_initial_rate_coefficients(gas,&Global,data);
    update_rates(gas,&Global,data);
    return data;
}

/* Remove H2 that appears unchanged on both sides of a reaction. */
static double H2SpectatorFlow(const int *reactants, int reactantCount,
                              const int *products, int count, const ChimesFloat *flow)
{
    const int molecule = Global.speciesIndices[sp_H2];
    double spectator = 0.0;
    for(int reaction=0;reaction<count;reaction++){
        int consumed = 0, produced = 0;
        for(int j=0;j<reactantCount;j++)
            consumed += reactants[reactantCount*reaction+j] == molecule;
        for(int j=0;j<3;j++)
            produced += products[3*reaction+j] == molecule;
        spectator += (consumed < produced ? consumed : produced)*flow[reaction];
    }
    return spectator;
}

static int PairMatches(const int *pair, int first, int second)
{
    return (pair[0] == first && pair[1] == second) ||
           (pair[0] == second && pair[1] == first);
}

/* Positive for X + H2 -> XH + HI, negative for the reverse. */
static double H2ExchangeNet(const int *reactants, int reactantCount, const int *products,
                            int count, const ChimesFloat *flow, int atom, int hydride)
{
    int molecule = Global.speciesIndices[sp_H2], neutral = Global.speciesIndices[sp_HI];
    atom = Global.speciesIndices[atom];
    hydride = Global.speciesIndices[hydride];
    if(atom < 0 || hydride < 0) return 0.0;
    double net = 0.0;
    for(int reaction=0;reaction<count;reaction++){
        const int *input = reactants+reactantCount*reaction;
        const int *output = products+3*reaction;
        if((reactantCount == 3 && input[2] >= 0) || output[2] >= 0) continue;
        if(PairMatches(input,atom,molecule) && PairMatches(output,hydride,neutral))
            net += flow[reaction];
        if(PairMatches(input,hydride,neutral) && PairMatches(output,atom,molecule))
            net -= flow[reaction];
    }
    return net;
}

/* Read the rovibrational component from the same CHIMES tables and
 * low-density/LTE interpolation used by chimes_cooling.c. */
static double H2RovibrationalCooling(const struct gasVariables *gas)
{
    int t;
    ChimesFloat fraction;
    chimes_get_table_index(chimes_table_bins.mol_cool_Temperatures,
                           chimes_table_bins.N_mol_cool_Temperatures,
                           log10(gas->temperature),&t,&fraction);
    const int colliders[] = {sp_H2,sp_HI,sp_HII,sp_HeI,sp_elec};
    const ChimesFloat *tables[] = {
        chimes_table_cooling.H2_cool_lowDens_H2,
        chimes_table_cooling.H2_cool_lowDens_HI,
        chimes_table_cooling.H2_cool_lowDens_HII,
        chimes_table_cooling.H2_cool_lowDens_HeI,
        chimes_table_cooling.H2_cool_lowDens_elec,
    };
    double lowDensity = 0.0;
    for(size_t i=0;i<sizeof colliders/sizeof colliders[0];i++)
        lowDensity += chimes_exp10(chimes_interpol_1d(tables[i],t,fraction))*Abundance(gas,colliders[i]);
    lowDensity *= gas->nH_tot;
    if(lowDensity <= 0.0) return 0.0;
    const double lte = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_LTE,t,fraction));
    return lte/(1.0+lte/lowDensity)*Abundance(gas,sp_H2)*gas->nH_tot;
}

static double H2DissociationCooling(double nH, const struct chimes_current_rates_struct *native)
{
    return 7.2e-12f*nH*(native->H2_collis_dissoc_rate[chimes_table_H2_collis_dissoc.Heating_reaction_index]
                        +native->T_dependent_rate[chimes_table_T_dependent.H2_collis_dissoc_heating_reaction_index]);
}

static double H2Heating(const struct gasVariables *gas, const struct chimes_current_rates_struct *native)
{
    const double nH = gas->nH_tot;
    const double xHI = Abundance(gas,sp_HI);
    const double xH2 = Abundance(gas,sp_H2);
    const double critical = xHI+xH2 == 0.0 ? 0.0 :
        (xHI+xH2)/(xHI/native->H2_collis_dissoc_crit_H+xH2/native->H2_collis_dissoc_crit_H2);
    const double photo = native->H2_photodissoc_rate[0];
    return 6.4e-13f*photo*nH
        + (2.7e-11f*photo+7.16e-12f*native->H2_dust_formation_rate
           +2.93e-12f*native->T_dependent_rate[chimes_table_T_dependent.H2_form_heating_reaction_index]
           +5.65e-12f*native->constant_rate[chimes_table_constant.H2_form_heating_reaction_index])
          *nH*nH/(nH+critical);
}

int BackendH2Rates(Cell *cell, double nH, double T, double fH2, H2Rates *out)
{
    struct gasVariables gas = cell->gas;
    ChimesFloat abundances[Global.totalNumberOfSpecies];
    ChimesFloat photons[1], G0[1], dissoc[1];   /* Benchmark uses one spectrum. */
    gas.abundances = abundances;
    if(fH2 >= 0.0){
        gas.isotropic_photon_density = photons;
        gas.G0_parameter = G0;
        gas.H2_dissocJ = dissoc;
        SetGas(&gas,nH,T);
        memset(abundances,0,sizeof abundances);
        abundances[Global.speciesIndices[sp_HI]] = 1.0-fH2;
        abundances[Global.speciesIndices[sp_H2]] = 0.5*fH2;
        abundances[Global.speciesIndices[sp_HeI]] = CMP_HE;
        abundances[Global.speciesIndices[sp_CII]] = CMP_C;
        abundances[Global.speciesIndices[sp_OI]] = CMP_O;
        abundances[Global.speciesIndices[sp_SiI]] = CMP_SI;
        abundances[Global.speciesIndices[sp_elec]] = CMP_C;
    }else{
        memcpy(abundances,cell->gas.abundances,sizeof abundances);
    }
    if(gas.temperature <= 0.0 || gas.temperature > Global.Tmol_K ||
       gas.nH_tot <= 0.0 || fH2 > 1.0) return 1;
    struct Species_Structure species[Global.totalNumberOfSpecies];
    struct chimes_current_rates_struct native;
    struct UserData data = PrepareDiagnosticRates(&gas,species,&native);
    update_rate_vector(species,&gas,&Global,data);

    /* CHIMES uses one LW reaction with separate dust and molecular factors.
     * The fitted LW shielding depends on T, N(H2) and turbulent width. */
    const double density = gas.nH_tot;
    const double xH2 = Abundance(&gas,sp_H2);
    const double thin = gas.isotropic_photon_density[0]*LIGHTSPEED*gas.light_speed_reduction_factor*
                        gas.H2_dissocJ[0]*chimes_table_H2_photodissoc.rates[0];
    const double dust = exp(-chimes_table_H2_photodissoc.gamma[0]*data.extinction);
    const double photo = native.H2_photodissoc_rate[0];
    /* These four channels leave the collider (HI, e, HeI or H2) unchanged.
     * Dissociative attachment and ion/molecule exchanges remain in other. */
    const int collH2 = chimes_table_T_dependent.H2_collis_dissoc_heating_reaction_index;
    double coll = native.T_dependent_rate[collH2];
    for(int reaction=0;reaction<chimes_table_H2_collis_dissoc.N_reactions[1];reaction++)
        coll += native.H2_collis_dissoc_rate[reaction];
    const double spectator =
        H2SpectatorFlow(chimes_table_T_dependent.reactants,3,chimes_table_T_dependent.products,
                        chimes_table_T_dependent.N_reactions[1],native.T_dependent_rate)
        + H2SpectatorFlow(chimes_table_constant.reactants,2,chimes_table_constant.products,
                           chimes_table_constant.N_reactions[1],native.constant_rate);
    const struct Species_Structure *h2 = &species[Global.speciesIndices[sp_H2]];
    *out = (H2Rates){
        .fH2 = 2.0*xH2, .xe = Abundance(&gas,sp_elec), .xHII = Abundance(&gas,sp_HII),
        .form_grain = native.H2_dust_formation_rate*density,
        .form_other = fmax(h2->creation_rate-native.H2_dust_formation_rate-spectator,0.0)*density,
        .destroy_photo = photo*density, .destroy_coll = coll*density,
        .destroy_other = fmax(h2->destruction_rate-photo-coll-spectator,0.0)*density,
        .k_form = native.H2_dust_formation_rate_coefficient*gas.dust_ratio*gas.dust_boost_factor,
        .k_photo0 = thin, .dust_shield = dust, .self_shield = native.H2_photodissoc_shield_factor[0]/dust,
        .shield_other = 1.0,
        .cool_h2 = H2RovibrationalCooling(&gas),
        .cool_h2_diss = H2DissociationCooling(density,&native),
        .heat_h2 = H2Heating(&gas,&native),
        .heat = calculate_total_cooling_rate(&gas,&Global,data,2),
        .cool = calculate_total_cooling_rate(&gas,&Global,data,1),
        .energy = 1.5*BOLTZMANNCGS*gas.temperature*calculate_total_number_density(abundances,density,&Global),
        .net_oh = density*H2ExchangeNet(chimes_table_T_dependent.reactants,3,chimes_table_T_dependent.products,
                                      chimes_table_T_dependent.N_reactions[1],native.T_dependent_rate,sp_OI,sp_OH),
        .net_chplus = density*(
            H2ExchangeNet(chimes_table_T_dependent.reactants,3,chimes_table_T_dependent.products,
                          chimes_table_T_dependent.N_reactions[1],native.T_dependent_rate,sp_CII,sp_CHp)
            + H2ExchangeNet(chimes_table_constant.reactants,2,chimes_table_constant.products,
                            chimes_table_constant.N_reactions[1],native.constant_rate,sp_CII,sp_CHp)),
    };
    free_current_rates_memory(&native,&Global);
    return 0;
}

#ifdef EQC_THERMAL_DIAGNOSTICS
/* The CHIMES total-rate function refreshes all three cooling tables.
   Sum a coolant's contributions using the same units as that function. */
static double AtomicCooling(const struct gasVariables *gas,
                            const struct chimes_current_rates_struct *native, int species)
{
    const int index = Global.speciesIndices[species];
    if(index < 0) return 0.0;
    const double xe = Abundance(gas,sp_elec);
    const double nH = gas->nH_tot;
    double coefficient = 0.0;
    for(int i=0;i<chimes_table_cooling.N_coolants;i++)
        if(chimes_table_cooling.coolants[i] == index)
            coefficient += native->cooling_rate[i]*xe;
    for(int i=0;i<chimes_table_cooling.N_coolants_2d;i++)
        if(chimes_table_cooling.coolants_2d[i] == index)
            coefficient += native->cooling_rate_2d[i]*xe;
    /* The 4-D coefficient already includes the collider densities. */
    for(int i=0;i<chimes_table_cooling.N_coolants_4d;i++)
        if(chimes_table_cooling.coolants_4d[i] == index)
            coefficient += native->cooling_rate_4d[i]/nH;
    return coefficient*gas->abundances[index]*nH*nH;
}

static double PhotoionizationHeating(const struct gasVariables *gas,
                                     const struct chimes_current_rates_struct *native, int molecular)
{
    double heat = 0.0;
    for(int i=0;i<chimes_table_photoion_fuv.N_reactions[molecular];i++)
        heat += native->photoion_fuv_heat_rate[i]*gas->abundances[chimes_table_photoion_fuv.reactants[i]];
    for(int i=0;i<chimes_table_photoion_euv.N_reactions[molecular];i++)
        heat += native->photoion_euv_heat_rate[i]*gas->abundances[chimes_table_photoion_euv.reactants[i]];
    return heat*gas->nH_tot;
}

static double CosmicRayHeating(const struct gasVariables *gas,
                               const struct chimes_current_rates_struct *native, int molecular)
{
    double rate = 0.0;
    for(int i=0;i<chimes_table_cosmic_ray.N_reactions[molecular];i++)
        rate += native->cosmic_ray_rate[i];
    const double logIonized = chimes_log10(chimes_max(Abundance(gas,sp_HII),CHIMES_FLT_MIN));
    const int bins = chimes_table_bins.N_secondary_cosmic_ray_xHII;
    int index;
    ChimesFloat fraction;
    chimes_get_table_index(chimes_table_bins.secondary_cosmic_ray_xHII,bins,logIonized,&index,&fraction);
    for(int i=0;i<2;i++){
        const double secondary = chimes_exp10(chimes_interpol_2d_fix_x(
            chimes_table_cosmic_ray.secondary_ratio,i,index,fraction,bins));
        rate -= native->cosmic_ray_rate[chimes_table_cosmic_ray.secondary_base_reaction[i]]
                *secondary/(1.0+secondary);
    }
    return 3.2e-11f*rate*gas->nH_tot;
}

/* Evaluate the dust terms at the actual CHIMES table coordinates. */
static void DustThermalRates(const struct gasVariables *gas, double extinction, ThermalRates *out)
{
    const double nH = gas->nH_tot;
    const double xe = Abundance(gas,sp_elec);
    int t;
    ChimesFloat dT;
    chimes_get_table_index(chimes_table_bins.Temperatures,chimes_table_bins.N_Temperatures,
                           chimes_log10(gas->temperature),&t,&dT);
    out->cool_gas_dust = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.gas_grain_transfer,t,dT))
                       *gas->dust_ratio*(gas->temperature-Global.grain_temperature)*nH*nH;

    double G0 = 0.0;
    const double lightSpeed = LIGHTSPEED*gas->light_speed_reduction_factor;
    for(int i=0;i<Global.N_spectra;i++)
        G0 += gas->isotropic_photon_density[i]*lightSpeed*gas->G0_parameter[i];
    G0 *= chimes_exp(-extinction*G0_GAMMA);
    out->pe_g0 = G0;
    if(G0 <= 0.0 || xe <= 0.0){
        out->pe_psi = NAN;
        return;
    }
    const double psi = G0*chimes_sqrt(gas->temperature)/chimes_max(nH*xe*0.5f,CHIMES_FLT_MIN);
    out->pe_psi = psi;
    const int bins = chimes_table_bins.N_Psi;
    int p;
    ChimesFloat dPsi;
    chimes_get_table_index(chimes_table_bins.Psi,bins,chimes_log10(chimes_max(psi,CHIMES_FLT_MIN)),&p,&dPsi);
    out->heat_pe = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.photoelectric_heating,t,p,dT,dPsi,bins))
                  *G0*gas->dust_ratio*nH;
    out->cool_dust_rec = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.grain_recombination,t,p,dT,dPsi,bins))
                       *gas->dust_ratio*xe*nH*nH;
}

int BackendThermalRates(const Cell *cell, double nH, double T,
                        const ThermalComposition *fixed, ThermalRates *out)
{
    struct gasVariables gas = fixed ? (struct gasVariables){0} : cell->gas;
    ChimesFloat abundances[Global.totalNumberOfSpecies], photons[1], G0[1], dissoc[1];
    gas.abundances = abundances;
    if(fixed){
        gas.isotropic_photon_density = photons;
        gas.G0_parameter = G0;
        gas.H2_dissocJ = dissoc;
        SetGas(&gas,nH,T);
        memset(abundances,0,sizeof abundances);
        abundances[Global.speciesIndices[sp_HI]] = fmax(0.0,1.0-fixed->fH2-fixed->xHII);
        abundances[Global.speciesIndices[sp_HII]] = fixed->xHII;
        abundances[Global.speciesIndices[sp_H2]] = 0.5*fixed->fH2;
        abundances[Global.speciesIndices[sp_elec]] = fixed->xe;
        abundances[Global.speciesIndices[sp_HeI]] = CMP_HE;
        abundances[Global.speciesIndices[sp_CII]] = CMP_C*fixed->xCII;
        abundances[Global.speciesIndices[sp_CI]] = CMP_C*(1.0-fixed->xCII);
        abundances[Global.speciesIndices[sp_OI]] = CMP_O*fixed->xOI;
        abundances[Global.speciesIndices[sp_OII]] = CMP_O*(1.0-fixed->xOI);
        abundances[Global.speciesIndices[sp_SiII]] = CMP_SI*fixed->xSiII;
        abundances[Global.speciesIndices[sp_SiI]] = CMP_SI*(1.0-fixed->xSiII);
    }else{
        memcpy(abundances,cell->gas.abundances,sizeof abundances);
    }
    if(gas.temperature <= 0.0 || gas.temperature > Global.Tmol_K || gas.nH_tot <= 0.0)
        return 1;

    struct Species_Structure species[Global.totalNumberOfSpecies];
    struct chimes_current_rates_struct native;
    struct UserData data = PrepareDiagnosticRates(&gas,species,&native);
    *out = (ThermalRates){
        .composition = {
            .fH2 = 2.0*Abundance(&gas,sp_H2), .xe = Abundance(&gas,sp_elec),
            .xHII = Abundance(&gas,sp_HII),
            .xCII = Abundance(&gas,sp_CII)/gas.element_abundances[1],
            .xOI = Abundance(&gas,sp_OI)/gas.element_abundances[3],
            .xSiII = Abundance(&gas,sp_SiII)/gas.element_abundances[6],
        },
        .xHI = Abundance(&gas,sp_HI),
        .heat = calculate_total_cooling_rate(&gas,&Global,data,2),
        .cool = calculate_total_cooling_rate(&gas,&Global,data,1),
        .heat_cr = CosmicRayHeating(&gas,&native,data.mol_flag_index),
        .heat_photo = PhotoionizationHeating(&gas,&native,data.mol_flag_index),
        .heat_h2 = H2Heating(&gas,&native),
        .cool_h2 = H2RovibrationalCooling(&gas),
        .cool_h2_diss = H2DissociationCooling(gas.nH_tot,&native),
    };
    DustThermalRates(&gas,data.extinction,out);
    out->cool_cii = AtomicCooling(&gas,&native,sp_CII);
    /* The C II table depends on T and electron density, with no neutral-collider axis. */
    out->cool_cii_e_only = out->cool_cii;
    out->cool_oi = AtomicCooling(&gas,&native,sp_OI);
    out->cool_ci = AtomicCooling(&gas,&native,sp_CI);
    out->cool_si2 = AtomicCooling(&gas,&native,sp_SiII);
    if(MolecularColumn(&gas,sp_CO,1.0) > 0.0){
        /* In the thermal kernel this carbon flag gates only the CO term.
           Keep abundances, collider densities and rate coefficients fixed. */
        struct globalVariables withoutCO = Global;
        withoutCO.element_included[0] = 0;
        out->cool_co = out->cool - calculate_total_cooling_rate(&gas,&withoutCO,data,1);
    }
    out->heat_other = out->heat - (out->heat_pe+out->heat_cr+out->heat_h2+out->heat_photo);
    out->cool_other = out->cool - (out->cool_h2+out->cool_h2_diss+out->cool_cii+out->cool_oi
                                 +out->cool_ci+out->cool_si2+out->cool_co
                                 +out->cool_dust_rec+out->cool_gas_dust);
    free_current_rates_memory(&native,&Global);
    return 0;
}
#endif
#endif
