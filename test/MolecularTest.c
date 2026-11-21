/* Reaction-balance check for H2 and CO chemistry with H-O charge transfer.
 *
 * The run settings match a galaxy simulation that enables H II regions, H2,
 * CO, dust, cosmic rays and soft X-rays together. The program writes
 *   equilibrium.csv  : a density, temperature and illumination grid,
 *   column_scan.csv  : H/H2 and C+/C/CO across total column density,
 * and fails if any reaction balance or the CO consistency exceeds its limit.
 * Every balance is rebuilt here from the public rate functions, independently
 * of the solver's own residuals. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Astro.h"
#include "ChemicalEquilibriumCooling.h"
#include "CollisionalGasPhaseChemistry.h"
#include "Constants.h"
#include "CosmicRay.h"
#include "EquilibriumAbundance.h"
#include "GrainSurfaceChemistry.h"
#include "MoleculeChemistry.h"
#include "PhotochemicalGasPhaseChemistry.h"

extern double b_turb;

static const double chemistryAccuracy = 1.0e-9;
static const double balanceLimit = 1.0e-6;      /* Each balance, relative to its larger flux. */
static const double consistencyLimit = 1.0e-6;  /* CO change when re-solved from the final H2. */

typedef struct {
    double nH, T, G0, Nt, logU;   /* logU = -INFINITY for no ionizing flux. */
} Condition;

static double Relative(double a, double b)
{
    const double scale = fmax(fabs(a), fabs(b));
    return scale > 0.0 ? fabs(a-b)/scale : 0.0;
}

static StructCEQData GasState(Condition c)
{
    const double nH = c.nH, Z = SOLAR_METALLICITY;
    StructCEQData s = {
        .T = c.T, .Nt = c.Nt, .Z = Z, .Zd = Z, .G0 = c.G0, .X0 = 1.0, .nel = 1.0e-3*nH,
        .NLyn = isfinite(c.logU) ? pow(10.0, c.logU)*nH*LIGHT_C_CGS : 0.0,
        .n = {nH, 0.1*nH, 3.019952e-4*nH, 4.570882e-4*nH, 3.548134e-6*nH},
    };
    for(int i = 0; i < NumberChemicalElement; i++){
        s.natom[i] = s.n[i];
        s.Rho += s.n[i]*ChemicalElementMass[i];
    }
    return s;
}

/* Reaction fluxes [cm^-3 s^-1] built from the same rate functions as the solver. */
typedef struct {
    double h2Formation, h2Photodissociation, h2CollisionalDissociation;
    double coFormation, coPhotodissociation;
    double hIonization, hRecombination;              /* Including charge transfer. */
    double oIonization, oRecombination;              /* O I <-> O II, including charge transfer. */
    double oIIIonization, oIIIRecombination;         /* O II <-> O III, including charge transfer. */
    double transferOtoOII, transferOIItoO, transferOIIItoOII;
} Fluxes;

static Fluxes ReactionFluxes(const StructCEQData *s)
{
    const double nH = s->n[ipHYDROGEN], T = s->T, ne = s->nel, xe = ne/nH;
    const double nHI = s->natom[ipHYDROGEN], nHII = s->nion[ipHYDROGEN], nH2 = s->nmol[ipH2];
    const double nOI = s->natom[ipOXYGEN], nOII = s->nion[ipOXYGEN], nOIII = s->nOIII;
    const double nCO = s->nmol[ipCO];
    Fluxes f = {0};

    /* H2: grain-surface formation, shielded photodissociation, collisional dissociation. */
    const double doppler = sqrt(SQ(12.90*sqrt(1.e-4*T))+SQ(b_turb));
    f.h2Formation = nH*H2FormationRate(T, s->Zd, EffectiveDustTemperature(s))*nHI;
    f.h2Photodissociation = nH2*H2PhotoDissociationRateWithDustFactor(s->G0, nH2*s->Nt/nH, doppler,
                                                                      DustShieldingFactorFUV(s->Nt, s->Zd));
    f.h2CollisionalDissociation = nH2*H2DissociationFromCoefficients(H2DissociationCoefficients(T), ne, nHI, nH2);

    /* CO (Glover & Clark 2012): C+ + H2 forms CO at k*Cp*n(H2), competing with
       C photoionization; CO is photodissociated with H2 and CO shielding. */
    const double formation = 4.e-16*pow(T/300.0,-0.2)*CEQLibRunParameters.CpForCO*nH2;
    const double photoCI = 3.1e-10*s->G0*DustShieldingFactor(s->Nt, s->Zd, 3.0);
    const double effective = formation+photoCI > 0.0 ? formation*photoCI/(formation+photoCI) : 0.0;
    const double freeCarbon = s->natom[ipCARBON] + s->nion[ipCARBON] + s->nCIII;
    f.coFormation = effective*freeCarbon;
    f.coPhotodissociation = nCO*COPhotoDissociationRate(s->G0, s->Nt, s->Zd)
                          * COShieldingFactor(s->Nt*nH2/nH, s->Nt*nCO/nH, doppler);

    /* H and O ionization with O-H charge transfer. */
    f.transferOtoOII = nOI*nHII*OIChargeTransferRateToHII(T);
    f.transferOIItoO = nOII*nHI*OIIChargeTransferRateToHI(T);
    f.transferOIIItoOII = nOIII*nHI*OIIIChargeTransferRateToHI(T);
    const CEQHydrogenRates h = CEQHydrogenIonizationRates(s);
    f.hIonization = nHI*h.ionization + f.transferOIItoO + f.transferOIIItoOII;
    f.hRecombination = nHII*h.recombination + f.transferOtoOII;
    f.oIonization = nOI*(OIPhotoionizationRate(s->NLyn) + ne*OICollisionalIonizationRate(T)
                         + s->X0*CosmicRayIonizationRateOI(xe)) + f.transferOtoOII;
    f.oRecombination = nOII*(ne*OIIRecombinationRate(T) + nH*OIIRecombinationRateOnDust(T, s->G0, ne, s->Zd))
                     + f.transferOIItoO;
    f.oIIIonization = nOII*(OIIPhotoionizationRate(s->NLyn) + ne*OIICollisionalIonizationRate(T));
    f.oIIIRecombination = nOIII*ne*OIIIRecombinationRate(T) + f.transferOIIItoOII;
    return f;
}

typedef struct {
    double h2, co, hydrogen, oxygen, oxygenII, coConsistency;
} Residuals;

static Residuals Balances(const StructCEQData *s, const Fluxes *f)
{
    Residuals r = {
        .h2 = Relative(f->h2Photodissociation + f->h2CollisionalDissociation, f->h2Formation),
        .co = Relative(f->coPhotodissociation, f->coFormation),
        .hydrogen = Relative(f->hIonization, f->hRecombination),
        .oxygen = Relative(f->oIonization, f->oRecombination),
        .oxygenII = Relative(f->oIIIonization, f->oIIIRecombination),
    };
    /* CO is limited by the oxygen supply when every O atom is in CO. */
    if(s->nmol[ipCO] >= s->n[ipOXYGEN]) r.co = 0.0;
    return r;
}

static void WriteHeader(FILE *file)
{
    fputs("nH,T,G0,Nt,Av,log_U,ne,"
          "HI,HII,H2,CI,CII,CIII,CO,OI,OII,OIII,"
          "H2_formation,H2_photodissociation,H2_collisional_dissociation,"
          "CO_formation,CO_photodissociation,"
          "H_ionization,H_recombination,O_ionization,O_recombination,OII_ionization,OIII_recombination,"
          "CT_O_to_OII,CT_OII_to_O,CT_OIII_to_OII,"
          "residual_H2,residual_CO,residual_H,residual_OI,residual_OII,CO_consistency\n", file);
}

static bool WriteRow(FILE *file, Condition c, const StructCEQData *s)
{
    const Fluxes f = ReactionFluxes(s);
    Residuals r = Balances(s, &f);
    if(s->nmol[ipCO] > 0.0){
        StructCEQData carbon = *s;
        const CEQLibStatus status = CEQSolveCarbonEquilibrium(&carbon, s->nmol[ipCO], chemistryAccuracy);
        if(status != CEQ_SUCCESS){
            fprintf(stderr, "CO re-solve failed at n_H=%g, T=%g K, G0=%g, N=%g, log U=%g (status %d)\n",
                    c.nH, c.T, c.G0, c.Nt, c.logU, status);
            return false;
        }
        r.coConsistency = Relative(carbon.nmol[ipCO], s->nmol[ipCO]);
    }
    const double nH = s->n[ipHYDROGEN], nC = s->n[ipCARBON], nO = s->n[ipOXYGEN];
    fprintf(file, "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.12e", c.nH, c.T, c.G0, c.Nt,
            VisualExtinction(c.Nt, s->Zd), c.logU, s->nel);
    const double fractions[] = {
        s->natom[ipHYDROGEN]/nH, s->nion[ipHYDROGEN]/nH, 2.0*s->nmol[ipH2]/nH,
        s->natom[ipCARBON]/nC, s->nion[ipCARBON]/nC, s->nCIII/nC, s->nmol[ipCO]/nC,
        s->natom[ipOXYGEN]/nO, s->nion[ipOXYGEN]/nO, s->nOIII/nO,
    };
    for(size_t i = 0; i < sizeof(fractions)/sizeof(fractions[0]); i++) fprintf(file, ",%.12e", fractions[i]);
    const double values[] = {
        f.h2Formation, f.h2Photodissociation, f.h2CollisionalDissociation,
        f.coFormation, f.coPhotodissociation,
        f.hIonization, f.hRecombination, f.oIonization, f.oRecombination, f.oIIIonization, f.oIIIRecombination,
        f.transferOtoOII, f.transferOIItoO, f.transferOIIItoOII,
        r.h2, r.co, r.hydrogen, r.oxygen, r.oxygenII, r.coConsistency,
    };
    for(size_t i = 0; i < sizeof(values)/sizeof(values[0]); i++) fprintf(file, ",%.12e", values[i]);
    fputc('\n', file);

    const double worst = fmax(fmax(fmax(r.h2, r.co), fmax(r.hydrogen, r.oxygen)), r.oxygenII);
    if(worst <= balanceLimit && r.coConsistency <= consistencyLimit) return true;
    fprintf(stderr, "Balance failed at n_H=%g, T=%g K, G0=%g, N=%g, log U=%g: "
            "H2 %.2e, CO %.2e, H %.2e, O I %.2e, O II %.2e, CO re-solve %.2e\n",
            c.nH, c.T, c.G0, c.Nt, c.logU, r.h2, r.co, r.hydrogen, r.oxygen, r.oxygenII, r.coConsistency);
    return false;
}

static bool Solve(Condition c, StructCEQData *state)
{
    const CEQLibStatus status = CEQLibTryCalculateChemicalAbundanceCoolingHeatingRates(state, chemistryAccuracy);
    if(status == CEQ_SUCCESS) return true;
    fprintf(stderr, "Chemistry failed at n_H=%g, T=%g K, G0=%g, N=%g, log U=%g (status %d)\n",
            c.nH, c.T, c.G0, c.Nt, c.logU, status);
    return false;
}

static bool RunGrid(FILE *file)
{
    /* Densities 1-1e5 cm^-3 with N = 2e21 (n_H/100)^(1/2) cm^-2, so A_V rises with density. */
    const double temperatures[] = {30.0, 100.0, 300.0, 1000.0};
    const double logU[] = {-INFINITY, -5.0, -3.0};
    WriteHeader(file);
    bool ok = true;
    for(size_t u = 0; u < sizeof(logU)/sizeof(logU[0]); u++)
    for(size_t t = 0; t < sizeof(temperatures)/sizeof(temperatures[0]); t++)
    for(int i = 0; i <= 20; i++){
        const double nH = pow(10.0, 0.25*i);
        const Condition c = {nH, temperatures[t], 1.7, 2.0e21*sqrt(nH/100.0), logU[u]};
        StructCEQData state = GasState(c);
        if(!Solve(c, &state)) return false;
        ok = WriteRow(file, c, &state) && ok;
    }
    return ok;
}

static bool RunColumnScan(FILE *file)
{
    /* A photodissociation-region sequence without ionizing photons: each column
       starts from the previous solution, as a cloud sampled from its surface inward. */
    const double densities[] = {100.0, 1000.0};
    const double fields[] = {1.7, 17.0};
    WriteHeader(file);
    bool ok = true;
    for(size_t d = 0; d < sizeof(densities)/sizeof(densities[0]); d++)
    for(size_t g = 0; g < sizeof(fields)/sizeof(fields[0]); g++){
        StructCEQData previous = {0};
        for(int i = 0; i <= 40; i++){
            const Condition c = {densities[d], 30.0, fields[g], pow(10.0, 19.0 + 0.1*i), -INFINITY};
            StructCEQData state = GasState(c);
            if(i > 0){
                state.nel = previous.nel;
                state.nmol[ipH2] = previous.nmol[ipH2];
                state.nmol[ipCO] = previous.nmol[ipCO];
            }
            if(!Solve(c, &state)) return false;
            ok = WriteRow(file, c, &state) && ok;
            previous = state;
        }
    }
    return ok;
}

int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr, "Usage: %s PEGASE2_DIRECTORY\n", argv[0]);
        return EXIT_FAILURE;
    }
    CEQLibiSetRunParamRecombinationCaseB(true);
    CEQLibiSetRunParamDustAttenuation(true);
    CEQLibiSetRunParamGrainSurfaceProcesses(true);
    CEQLibiSetRunParamCosmicRayIonization(true);
    CEQLibiSetRunParamSoftXrayIonization(true);
    CEQLibiSetRunParamCIEResonanceCooling(true);
    CEQLibiSetRunParamHIIRegion(true);
    CEQLibiSetRunParamH2Formation(true, 1.0);
    CEQLibiSetRunParamCOFormation(true, 1.0);
    if(CEQLibTryInitialization(0, argv[1]) != CEQ_SUCCESS) return EXIT_FAILURE;

    FILE *grid = fopen("equilibrium.csv", "w");
    FILE *scan = fopen("column_scan.csv", "w");
    if(grid == NULL || scan == NULL){
        perror("MolecularTest");
        return EXIT_FAILURE;
    }
    const bool gridOk = RunGrid(grid);
    const bool scanOk = RunColumnScan(scan);
    const bool closed = fclose(grid) == 0 && fclose(scan) == 0;
    if(gridOk && scanOk && closed){
        puts("Molecular H2/CO reaction balance and CO consistency checks passed.");
        return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}
