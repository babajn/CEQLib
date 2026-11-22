"""Independent Verner shell integrals and CHIMES FG20 comparisons for C, O and Si."""

import argparse
from pathlib import Path

import h5py
import numpy as np
from scipy.special import expn

from compare import ELECTRON_SPLIT, ELECTRON_VOLT, LIGHT_SPEED, PLANCK, THRESHOLDS, cross_sections, read_cloudy_spectrum
from metal_compare import AtomicModel, deposition, log_interpolate, records, spectrum_at_redshift, write_rows


STAGES = ("CI", "CII", "OI", "OII", "SiI", "SiII")
ATOMIC_NUMBERS = (6, 6, 8, 8, 14, 14)
STAGE_THRESHOLDS = np.array([11.26, 24.38, 13.62, 35.12, 8.152, 16.35])
CHIMES_SPECIES = (7, 8, 23, 24, 57, 58)  # chimes/src/chimes_proto.h
ABUNDANCES = 10.0 ** (np.array([8.43, 8.69, 7.51]) - 12.0)
ATOMIC_DATA_DEFAULT = Path(__file__).resolve().parents[3] / "c25.00" / "data"


def columns(path):
    return [line.split("#",1)[0].split() for line in path.read_text().splitlines()
            if line.split("#",1)[0].strip()]


class Photoionization:
    def __init__(self, directory):
        # Read independently from the Cloudy c25 tables: Verner & Yakovlev
        # (1995), Verner et al. (1996), and Kaastra & Mewe (1993).
        _, angular_momenta, inner, _, *fit_data = columns(directory / "phfit.dat")
        first_end = fit_data.index(["-1", "-1", "-1"])
        partial, outer = fit_data[:first_end], fit_data[first_end+1:-1]
        self.edges = [*STAGE_THRESHOLDS]
        self.outer = []
        self.subshells = []
        self.first_inner = []
        yields = {tuple(map(int, line[:3])): sum(map(float,line[5:]))
                  for line in columns(directory / "mewe_nelectron.dat")}
        for i, atomic_number in enumerate(ATOMIC_NUMBERS):
            charge = i % 2
            electrons = atomic_number-charge
            key = [str(electrons-1), str(atomic_number-1)]
            averaged, = (list(map(float,row[2:])) for row in outer if row[:2] == key)
            self.outer.append(averaged)
            shells = []
            for row in partial:
                if row[1:3] == key:
                    shell = int(row[0])
                    parameters = list(map(float,row[3:]))
                    # Cloudy uses newer outer thresholds for some subshells.
                    # The high-energy fits here follow the authors' PHFIT2
                    # values for the specified ions and their filled shells.
                    if atomic_number == 14 and charge == 0 and shell == 3:
                        parameters[0] = 15.17
                    if atomic_number == 14 and charge == 1 and shell == 3:
                        parameters[0] = 22.40
                    shells.append((shell, int(angular_momenta[shell]), parameters,
                                   yields[atomic_number,charge+1,shell+1]))
                    self.edges.append(parameters[0])
            self.subshells.append(shells)
            last_inner_shell = int(inner[electrons-1])-1
            self.first_inner.append(next(parameters[0] for shell, _, parameters, _ in shells
                                         if shell == last_inner_shell))

    def evaluate(self, energy):
        direct = np.zeros((len(STAGES), energy.size))
        multiple = np.zeros_like(direct)
        for i, threshold in enumerate(STAGE_THRESHOLDS):
            e0, s0, ya, p, yw, y0, y1 = self.outer[i]
            x = energy / e0 - y0
            y = np.hypot(x,y1)
            valence = 1.0e-18 * s0 * ((x-1)**2 + yw**2) * y**(0.5*p-5.5) / (1+np.sqrt(y/ya))**p
            first_edge = self.first_inner[i]
            direct[i] = np.where((energy >= threshold) & (energy < first_edge), valence, 0.0)
            for _, angular, (edge, e0, s0, ya, p, yw), probability in self.subshells[i]:
                x = energy/e0
                sigma = 1.0e-18 * s0 * ((x-1)**2 + yw**2) * x**(0.5*p-5.5-angular) / (1+np.sqrt(x/ya))**p
                sigma = np.where(energy >= max(first_edge, edge), sigma, 0.0)
                direct[i] += sigma
                multiple[i] += probability*sigma
        return direct, multiple


class Slab:
    def __init__(self, energies, j, atomic, opacity, lower_energy=13.60, upper_energy=None):
        nodes, weights = np.polynomial.legendre.leggauss(24)
        lower = max(energies[0],lower_energy)
        upper = energies[-1] if upper_energy is None else min(upper_energy,energies[-1])
        boundaries = np.unique(np.r_[energies, THRESHOLDS, THRESHOLDS+ELECTRON_SPLIT,
                                     atomic.edges, opacity.edges, lower, upper, 30.0])
        boundaries = np.log(boundaries[(boundaries >= lower) & (boundaries <= upper)])
        width = np.diff(boundaries)/2
        log_energy = ((boundaries[:-1]+width)[:,None]+width[:,None]*nodes).ravel()
        energy = np.exp(log_energy)
        source = np.log(energies)
        prior = np.clip(np.searchsorted(source,log_energy)-1,0,energies.size-2)
        fraction = (log_energy-source[prior])/np.diff(source)[prior]
        intensity = log_interpolate(j[prior],j[prior+1],fraction)
        self.measure = 2*np.pi/PLANCK*intensity*(width[:,None]*weights).ravel()
        self.direct, self.multiple = atomic.evaluate(energy)
        hydrogen_helium = cross_sections(energy)
        self.opacity_hhe = hydrogen_helium[0] + 0.1*hydrogen_helium[1]
        metal = opacity.evaluate(energy)
        self.opacity_metals = metal[0]
        self.metal_power = metal[1:3]*ELECTRON_VOLT
        kinetic = np.maximum(0.0, energy-THRESHOLDS[:,None])
        power = hydrogen_helium*kinetic*ELECTRON_VOLT
        high = np.where(kinetic >= ELECTRON_SPLIT,power,0.0)
        self.hydrogen_helium = np.stack((hydrogen_helium,power-high,high),axis=1).reshape(9,-1)
        excess = np.maximum(0.0, energy - STAGE_THRESHOLDS[:,None])
        self.q_below_30 = np.where(energy < 30.0, self.direct*excess*ELECTRON_VOLT, 0.0)

    def evaluate(self, z_solar, top, bottom):
        opacity = self.opacity_hhe + z_solar*self.opacity_metals
        measure = self.measure*(expn(2,top*opacity)+expn(2,bottom*opacity))
        direct, multiple = self.direct @ measure, self.multiple @ measure
        existing = np.r_[self.hydrogen_helium @ measure,z_solar*(self.metal_power @ measure)]
        # Each element can distribute its entire abundance among I, II and
        # higher states. Assign it to whichever tracked ion heats more.
        upper_heat = z_solar*(ABUNDANCES @ np.max((self.q_below_30 @ measure).reshape(3,2),axis=1))
        return direct, multiple, upper_heat, deposition(existing)[0]


def compare_integrals(folder, get_slab):
    rows = records(folder / "metal-photoionization-rates.csv")
    comparison = []
    fractions = []
    maxima = {"grid": (0,()), "interpolation": (0,())}
    omitted_heat_max = (0,())
    for row in rows:
        z, metallicity, top, bottom = (float(row[key]) for key in ("z", "Z_solar", "NHI_top", "NHI_bottom"))
        case = (z,metallicity,top,bottom)
        reference, multiple, omitted_heat, existing_heat = get_slab(z).evaluate(metallicity,top,bottom)
        actual = np.array([float(row[stage]) for stage in STAGES])
        if not np.isfinite(actual).all() or np.any(actual < 0.0) or np.any(reference <= 0.0):
            raise ValueError(f"Invalid metal photoionization rates at {case}")
        errors = np.abs(actual/reference-1.0)
        auger_fraction = multiple/reference
        i = np.argmax(errors)
        label = "grid" if z.is_integer() and top == bottom else "interpolation"
        if errors[i] > maxima[label][0]:
            maxima[label] = errors[i], (*case,STAGES[i])
        ratio = omitted_heat/existing_heat
        omitted_heat_max = max(omitted_heat_max,(ratio,case),key=lambda pair: pair[0])
        comparison.append([*row.values(), *reference, *errors, omitted_heat, existing_heat, ratio])
        fractions.append([*case,*auger_fraction])
    write_rows(folder / "metal-photoionization-integrals.csv",list(rows[0]) +
               [f"{s}_independent" for s in STAGES] + [f"{s}_relative_error" for s in STAGES] +
               ["upper_13p6_to_30_eV_heat_per_H", "existing_heat_per_H", "upper_omitted_heat_fraction"], comparison)
    write_rows(folder / "metal-photoionization-auger.csv",["z", "Z_solar", "NHI_top", "NHI_bottom", *STAGES],fractions)
    for name, (error, case) in maxima.items():
        print(f"Metal photoionization independent {name} maximum {error:.4%} at {case}")
    print(f"13.6–30 eV omitted metal heating upper estimate {omitted_heat_max[0]:.4%} of existing UVB heating, "
          f"at {omitted_heat_max[1]} (x_HII=0.01, y_HeI=0.1)")
    print("Maximum multiple-ejection fractions over tested columns:",
          ", ".join(f"{stage} {fraction:.3%}" for stage, fraction in zip(STAGES,np.max(fractions,axis=0)[4:])))
    return rows, all(maximum[0] <= 0.01 for maximum in maxima.values())


def compare_chimes(folder, data, rows, get_euv, get_fuv):
    with h5py.File(data / "chimes_main_data.hdf5") as reactions:
        base = reactions["photoion_euv/reactants"][:]
        auger = reactions["photoion_auger_euv/reactants"][:]
        fuv = reactions["photoion_fuv/reactants"][:]
        auger_fuv = reactions["photoion_auger_fuv/reactants"][:]
    result = []
    errors = []
    flux_errors = []
    for row in rows:
        if float(row["NHI_top"]) != 0.0 or float(row["NHI_bottom"]) != 0.0:
            continue
        z = float(row["z"])
        path = data / "FG20_cross_sections" / f"z{z:.3f}_cross_sections.hdf5"
        with h5py.File(path) as tables:
            if tables["Header/rt_coupling_flag"][()] != 0:
                raise ValueError(f"{path}: expected optically thin full-EUV normalization")
            flux = LIGHT_SPEED * float(tables["isotropic_photon_density"][0])
            flux_errors.append(abs(2*sum(get_euv(z).measure)/flux-1.0))
            below_13p6 = get_fuv(z).direct @ (2*get_fuv(z).measure)
            for i, (stage, species) in enumerate(zip(STAGES,CHIMES_SPECIES)):
                actual = float(row[stage])
                indexes = np.flatnonzero(base == species)
                if len(indexes) == 0:
                    other = np.flatnonzero(fuv == species).item()
                    # The species-based FUV group also contains the shortward
                    # contribution. There is no independently tabulated split.
                    full_band = flux*(float(tables["photoion_fuv/sigmaPhot"][other]) +
                                      sum(tables["photoion_auger_fuv/sigmaPhot"][:][auger_fuv == species]))
                    full_error = (actual+below_13p6[i])/full_band-1.0
                    print(f"CHIMES z={z:g} {stage}: no EUV-only entry; separate full-band comparison {full_error:+.3%}")
                    result.append([z,stage,actual,"","","","","",other,full_band,below_13p6[i],full_error])
                    continue
                reaction = indexes.item()
                single = flux * float(tables["photoion_euv/sigmaPhot"][reaction])
                extra = flux * sum(tables["photoion_auger_euv/sigmaPhot"][:][auger == species])
                # Assigning all events to one stage in EQCLib means one
                # must compare the total neutral/ion disappearance rate.
                single_error = actual/single-1.0
                total_error = actual/(single+extra)-1.0
                errors.append(abs(total_error))
                result.append([z,stage,actual,reaction,single,extra,single_error,total_error,"","","",""])
                print(f"CHIMES z={z:g} {stage}: EQCLib {actual:.6e}; single {single:.6e}, {single_error:+.3%}; "
                      f"single + Auger {single+extra:.6e}, {total_error:+.3%}")
    write_rows(folder / "metal-photoionization-chimes.csv", ["z", "stage", "EQCLib", "CHIMES_EUV_reaction",
               "CHIMES_single", "CHIMES_Auger", "relative_error_vs_single", "relative_error_vs_total",
               "CHIMES_FUV_reaction", "CHIMES_full_band", "independent_below_13p6_eV",
               "relative_error_full_band"],result)
    print(f"Maximum CHIMES difference for available EUV stages: {max(errors):.3%}; "
          f"EUV photon flux normalization differs by at most {max(flux_errors):.4%}")


def impact_estimate(folder):
    # Use the high-energy Lotz scaling ln(E/I)/(E I). Attribute all
    # valence electrons to each species' lowest threshold for an indicative
    # upper estimate, using two electrons for He and four/six/four for C/O/Si.
    rows = []
    thresholds = np.r_[THRESHOLDS[:2],STAGE_THRESHOLDS[::2]]
    electrons = np.array([1,2,4,6,4])
    for energy in (100,300,1000,10000):
        cross_section = electrons*np.log(energy/thresholds)/thresholds
        metal = ABUNDANCES*cross_section[2:]
        for x_HII in (0.01,0.5,0.9):
            x_HeI = 0.1*(1.0-x_HII)
            reference = (1.0-x_HII)*cross_section[0]+x_HeI*cross_section[1]
            rows.append([energy,x_HII,x_HeI,*(metal/reference),sum(metal)/reference])
    print(f"Indicative C/O/Si to H/He impact-ionization ratio at x_HII=0.01: "
          f"{max(row[-1] for row in rows if row[1] == 0.01):.3%} maximum (100 eV to 10 keV)")
    write_rows(folder / "metal-photoionization-impact-estimate.csv",
               ["electron_eV", "x_HII", "x_HeI_per_H", "C_to_HHe", "O_to_HHe", "Si_to_HHe", "total_to_HHe"], rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum",type=Path)
    parser.add_argument("chimes_data",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--atomic-data",type=Path,default=ATOMIC_DATA_DEFAULT)
    args = parser.parse_args()
    redshifts, energy, spectra = read_cloudy_spectrum(args.spectrum)
    photoionization = Photoionization(args.atomic_data)
    opacity = AtomicModel()
    slabs = {}
    fuv_slabs = {}

    def get_slab(z):
        if z not in slabs:
            slabs[z] = Slab(energy,spectrum_at_redshift(redshifts,spectra,z),photoionization,opacity)
        return slabs[z]

    def get_fuv(z):
        if z not in fuv_slabs:
            fuv_slabs[z] = Slab(energy,spectrum_at_redshift(redshifts,spectra,z),photoionization,opacity,
                                lower_energy=min(STAGE_THRESHOLDS),upper_energy=13.60)
        return fuv_slabs[z]

    rows, success = compare_integrals(args.output,get_slab)
    compare_chimes(args.output,args.chimes_data,rows,get_slab,get_fuv)
    impact_estimate(args.output)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
