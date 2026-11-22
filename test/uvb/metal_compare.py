"""Independent FG20 transport and deposition integrals for the metal UVB tests."""

import argparse
import csv
from pathlib import Path
import re

import numpy as np
from scipy.special import expn

from compare import (
    ELECTRON_SPLIT,
    ELECTRON_VOLT,
    PLANCK,
    SMALL_BAND_FRACTION,
    SPECIES,
    THRESHOLDS,
    cross_sections,
    read_cloudy_spectrum,
    relative_error,
)


IONIZED_HYDROGEN = 0.01
ABSORBERS = np.array([1.0 - IONIZED_HYDROGEN, 0.1, 0.0])
ATOMIC_METALS = ("C", "O", "Si")
BAND_EDGES = (100.0, 300.0, 2000.0)
BAND_NAMES = ("below_100_eV", "100_to_300_eV", "300_eV_to_2_keV", "above_2_keV")
RATE_COLUMNS = tuple(f"{part}_{name}" for name in SPECIES for part in ("Gamma", "power_low", "power_high"))
METAL_POWER_COLUMNS = ("metal_power_low", "metal_power_high")
ENERGY_COLUMNS = ("H_UVB", "Gamma_HI_total", "Q_other")
MAX_RELATIVE_ERROR = 0.01
BCM_AVOGADRO = 6.022045e23
ATOMIC_COEFFICIENTS = Path(__file__).resolve().parents[2] / "src" / "MetalXrayData.h"


def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_rows(path, header, rows):
    with path.open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


class AtomicModel:
    def __init__(self):
        # Published atomic constants are shared input data. Opacity, energy
        # loss and the subsequent spectral integrals are calculated in Python.
        data = ATOMIC_COEFFICIENTS.read_text()
        blocks = re.split(r"\[CEQ_METAL_XRAY_([A-Z]+)\]\s*=\s*{", data)[1:]
        self.elements = {}
        self.edges = [30.0]
        for name, block in zip(blocks[::2], blocks[1::2]):
            if name == "ELEMENTS":
                continue
            metadata, polynomials = block.split("{", 1)
            meta = np.fromstring(metadata.strip().rstrip(","), sep=",")
            assert meta.size == 9
            count = int(meta[-1])
            fits = []
            for upper, terms, coefficients in re.findall(r"{\s*([^,]+),\s*(\d+),\s*{\s*([^{}]+)}\s*}",
                                                         polynomials):
                coefficients = np.fromstring(coefficients, sep=",")
                assert coefficients.size == int(terms)
                fits.append((float(upper), coefficients))
            assert len(fits) == count
            abundance, mass, first, second, yield_k, alpha_k, edge_l, edge_k, _ = meta
            loss_k = first + (1.0 - yield_k) * second + yield_k * alpha_k
            self.elements[name] = (abundance, mass, first, first + second, loss_k, edge_l, edge_k, fits)
            self.edges.extend(upper for upper, _ in fits)
            self.edges.extend((edge_l, edge_k, first + ELECTRON_SPLIT, first + second + ELECTRON_SPLIT,
                               loss_k + ELECTRON_SPLIT))
        assert len(self.elements) == 15

    def evaluate(self, energy):
        # Below 30 eV the metal opacity is outside the BCM92 range. Above
        # 10 keV the last value of sigma E^3 continues the high-energy tail.
        log_fit_energy = np.log(np.clip(energy,30.0,10000.0))
        cross_sections = {}
        result = np.zeros((3,energy.size))
        for name, (abundance, mass, loss_outer, loss_l, loss_k, edge_l, edge_k, fits) in self.elements.items():
            boundaries = [upper for upper, _ in fits[:-1]]
            interval = np.searchsorted(boundaries,energy,side="right")
            polynomial = np.zeros_like(energy)
            for i, (_, coefficients) in enumerate(fits):
                mask = interval == i
                polynomial[mask] = np.polynomial.polynomial.polyval(log_fit_energy[mask],coefficients)
            sigma_atom = np.where(energy >= 30.0, (mass/BCM_AVOGADRO)*np.exp(polynomial)/energy**3, 0.0)
            cross_sections[name] = sigma_atom
            lower_loss = np.where((edge_l > 0.0) & (energy >= edge_l), loss_l, loss_outer)
            kinetic = np.maximum(0.0, energy - np.where(energy >= edge_k, loss_k, lower_loss))
            opacity = sigma_atom * abundance
            power = opacity * kinetic
            result[0] += opacity
            result[1] += np.where(kinetic < ELECTRON_SPLIT, power, 0.0)
            result[2] += np.where(kinetic >= ELECTRON_SPLIT, power, 0.0)
        return np.vstack((result, [cross_sections[metal.upper()] for metal in ATOMIC_METALS]))


def log_interpolate(left, right, fraction):
    value = (1.0 - fraction) * left + fraction * right
    positive = (left > 0.0) & (right > 0.0)
    t = np.broadcast_to(fraction, value.shape)[positive]
    value[positive] = np.exp((1.0 - t) * np.log(left[positive]) + t * np.log(right[positive]))
    return value


def spectrum_at_redshift(redshifts, intensities, redshift):
    upper = np.searchsorted(redshifts, redshift)
    if redshifts[upper] == redshift:
        return intensities[upper]
    lower = upper - 1
    fraction = (redshift - redshifts[lower]) / (redshifts[upper] - redshifts[lower])
    return log_interpolate(intensities[lower], intensities[upper], fraction)


class Integral:
    def __init__(self, energies, intensities, atomic):
        nodes, weights = np.polynomial.legendre.leggauss(24)
        minimum = max(energies[0], THRESHOLDS[0])
        maximum = energies[-1]
        edges = np.unique(np.r_[energies, THRESHOLDS, THRESHOLDS + ELECTRON_SPLIT, atomic.edges,
                                BAND_EDGES, minimum, maximum])
        edges = np.log(edges[(edges >= minimum) & (edges <= maximum)])
        midpoints = 0.5 * (edges[:-1] + edges[1:])
        widths = 0.5 * np.diff(edges)
        log_energy = (midpoints[:, None] + widths[:, None] * nodes).ravel()
        energy = np.exp(log_energy)
        source_log_energy = np.log(energies)
        left = np.searchsorted(source_log_energy, log_energy) - 1
        left = np.clip(left, 0, energies.size - 2)
        position = (log_energy - source_log_energy[left]) / np.diff(source_log_energy)[left]
        intensity = log_interpolate(intensities[left], intensities[left + 1], position)
        self.photon_measure = (2.0 * np.pi / PLANCK) * intensity * (widths[:, None] * weights).ravel()

        hydrogen_helium = cross_sections(energy)
        atomic_samples = atomic.evaluate(energy)
        self.base_opacity = hydrogen_helium[0] + 0.1 * hydrogen_helium[1]
        self.metal_opacity = atomic_samples[0]
        self.metal_kinetic_kernel = atomic_samples[1:3] * ELECTRON_VOLT
        self.metal_atomic_cross_sections = atomic_samples[3:]
        k_edges = np.array([atomic.elements[metal.upper()][6] for metal in ATOMIC_METALS])
        self.metal_atomic_above_k_cross_sections = self.metal_atomic_cross_sections * (energy >= k_edges[:, None])
        self.bands = np.searchsorted(BAND_EDGES, energy, side="right")
        kinetic = np.maximum(energy - THRESHOLDS[:, None], 0.0)
        power = hydrogen_helium * kinetic * ELECTRON_VOLT
        low = np.where(kinetic < ELECTRON_SPLIT, power, 0.0)
        high = power - low
        self.hydrogen_helium_kernels = np.stack((hydrogen_helium, low, high), axis=1).reshape(9, -1)

    def measure(self, metallicity, top, bottom):
        optical_depth_per_column = self.base_opacity + metallicity * self.metal_opacity
        return self.photon_measure * (expn(2, top * optical_depth_per_column) +
                                      expn(2, bottom * optical_depth_per_column))

    def integrate(self, metallicity, top, bottom, band=None):
        measure = self.measure(metallicity, top, bottom)
        if band is not None:
            measure *= self.bands == band
        elemental = self.hydrogen_helium_kernels @ measure
        metals = metallicity * (self.metal_kinetic_kernel @ measure)
        return np.r_[elemental, metals]

    def direct_metal_rates(self, metallicity, column):
        measure = self.measure(metallicity, column, column)
        return self.metal_atomic_cross_sections @ measure, self.metal_atomic_above_k_cross_sections @ measure


def deposition(values):
    primaries = values[:9].reshape(3, 3)
    low_power, high_power = ABSORBERS @ primaries[:, 1:] + values[9:]
    x = IONIZED_HYDROGEN
    heat_share = 0.9971 * (1.0 - (1.0 - x**0.2663)**1.3163)
    hi_share = 0.3908 * (1.0 - x**0.4092)**1.7592
    hei_share = 0.0554 * (1.0 - x**0.4614)**1.666
    heating = low_power + heat_share * high_power
    gamma = primaries[0, 0] + hi_share * high_power / (ABSORBERS[0] * THRESHOLDS[0] * ELECTRON_VOLT)
    other = (1.0 - heat_share - hi_share - hei_share) * high_power
    return heating, gamma, other


def compare_rates(folder, get_slab):
    rows = records(folder / "metal-rates.csv")
    fields = RATE_COLUMNS + METAL_POWER_COLUMNS + ENERGY_COLUMNS
    direct_max = (0.0, None)
    interp_max = (0.0, None)
    closure = 0.0
    compared = []
    for row in rows:
        z, metallicity, top, bottom = (float(row[key]) for key in ("z", "Z_solar", "NHI_top", "NHI_bottom"))
        case = (z, metallicity, top, bottom)
        integrated = get_slab(z).integrate(metallicity, top, bottom)
        expected = dict(zip(fields, [*integrated, *deposition(integrated)]))
        band_floors = {}
        for i, name in enumerate(SPECIES):
            floor = SMALL_BAND_FRACTION * sum(integrated[3 * i + 1:3 * i + 3])
            band_floors[f"power_low_{name}"] = band_floors[f"power_high_{name}"] = floor
        for name in METAL_POWER_COLUMNS:
            band_floors[name] = SMALL_BAND_FRACTION * sum(integrated[-2:])
        errors = {name: relative_error(float(row[name]), expected[name], band_floors.get(name, 0.0))
                  for name in fields}
        for key, error in errors.items():
            diagnostic = abs(error), (key, *case)
            if (z in (0.0, 2.0) and metallicity in (0.0, 0.1, 1.0, 2.0) and
                    top == bottom and top in (1.0e19, 1.0e20, 1.0e21, 1.0e22)):
                direct_max = max(direct_max, diagnostic, key=lambda item: item[0])
            else:
                interp_max = max(interp_max, diagnostic, key=lambda item: item[0])
        closure = max(closure,float(row["energy_closure"]))
        compared.append([*row.values(), *(expected[name] for name in fields), *(errors[name] for name in fields)])
    write_rows(folder / "metal-comparison.csv", list(rows[0]) + [name + "_integral" for name in fields] +
               [name + "_relative_error" for name in fields], compared)
    print(f"Metal UVB independent integration grid maximum {direct_max[0]:.6%} at {direct_max[1]}")
    print(f"Metal UVB independent integration interpolation maximum {interp_max[0]:.6%} at {interp_max[1]}")
    print(f"Metal UVB budget closure {closure:.3e}; negligible band denominator floor "
          f"{SMALL_BAND_FRACTION:g} of the absorber power")
    passed = max(direct_max[0], interp_max[0]) <= MAX_RELATIVE_ERROR and closure <= 32*np.finfo(float).eps
    return passed, rows


def cross_section_comparison(folder, atomic):
    # Morrison & McCammon (1983), Table 2; their historical abundances differ
    # from Asplund et al. (2009). E is in keV and the polynomial has cm²/H units.
    mm83 = {
        0.3: (78.1, 18.8, 4.3),
        0.5: (71.4, 66.8, -51.4),
        1.0: (120.6, 169.3, -47.7),
        2.0: (202.7, 104.7, -17.0),
        5.0: (433.9, -2.4, 0.75),
    }
    rows = records(folder / "metal-cross-sections.csv")
    result = []
    for row in rows:
        keV = float(row["energy_keV"])
        c0, c1, c2 = mm83[keV]
        reference = (c0 + keV * (c1 + keV * c2)) * 1.0e-24 / keV**3
        actual = float(row["sigma_total"])
        difference = actual / reference - 1.0
        independent_metal = atomic.evaluate(np.array([1000.0 * keV]))[0,0]
        atomic_difference = float(row["sigma_metal"]) / independent_metal - 1.0
        print(f"Cross sections ({keV:g} keV): EQCLib {actual:.6e}, MM83 {reference:.6e}, "
              f"difference {difference:+.2%}")
        result.append([*row.values(), reference, difference, independent_metal, atomic_difference])
    write_rows(folder / "metal-cross-section-comparison.csv", list(rows[0]) +
               ["sigma_Morrison_McCammon_1983", "relative_difference_MM83",
                "sigma_metal_Python", "relative_difference_metal_Python"], result)
    return all(abs(row[-1]) <= 1.0e-9 for row in result)


def write_depth_tables(folder, rows, get_slab):
    index = {(float(row["z"]), float(row["NHI_top"]), float(row["Z_solar"])): row for row in rows
             if row["NHI_top"] == row["NHI_bottom"]}
    totals = []
    bands = []
    for z in (0.0, 2.0):
        for column in (1.0e19, 1.0e20, 1.0e21, 1.0e22):
            without = index[z, column, 0.0]
            solar = index[z, column, 1.0]
            heat_without, heat_solar = float(without["H_UVB"]), float(solar["H_UVB"])
            gamma_without, gamma_solar = float(without["Gamma_HI_total"]), float(solar["Gamma_HI_total"])
            totals.append((z, column, heat_without, heat_solar, heat_solar / heat_without,
                           gamma_without, gamma_solar, gamma_solar / gamma_without))
            print(f"Deep slab z={z:g}, N={column:.0e}: heat {heat_without:.4e} -> {heat_solar:.4e}, "
                  f"total HI {gamma_without:.4e} -> {gamma_solar:.4e}")
            for metallicity in (0.0, 1.0):
                for band, label in enumerate(BAND_NAMES):
                    moments = get_slab(z).integrate(metallicity,column,column,band)
                    heating, gamma, _ = deposition(moments)
                    bands.append((z, column, metallicity, label, heating, gamma, sum(moments[-2:])))
    write_rows(folder / "metal-deep-comparison.csv",
               ["z", "NHI_per_side", "heating_H_Z0", "heating_H_Zsolar", "heating_ratio",
                "Gamma_HI_total_Z0", "Gamma_HI_total_Zsolar", "Gamma_HI_total_ratio"], totals)
    write_rows(folder / "metal-deep-bands.csv",
               ["z", "NHI_per_side", "Z_solar", "photon_band", "heating_per_H", "Gamma_HI_total",
                "metal_power_per_H"], bands)


def compare_metal_chemistry(folder, get_slab):
    rows = records(folder / "metal-ordinary-chemistry.csv")
    comparisons = []
    maxima = {element: np.zeros(4) for element in ATOMIC_METALS}
    for row in rows:
        z, column = float(row["z"]), float(row["NHI_per_side"])
        element = row["element"]
        all_energy, above_k = get_slab(z).direct_metal_rates(1.0,column)
        index = ATOMIC_METALS.index(element)
        xray, xray_above_k = all_energy[index], above_k[index]
        fuv = float(row["FUV"])
        ct_ionization = float(row["charge_transfer_ionization"])
        ct_recombination = float(row["charge_transfer_recombination"])
        recombination = float(row["electron_recombination"]) + float(row["dust_recombination"])
        # The ASURA benchmark enables H–O charge transfer. The C and Si charge
        # transfer coefficients and the Si FUV coefficient are diagnostic only.
        active_fuv = fuv if element == "C" else 0.0
        ordinary_ionization = (active_fuv + float(row["cosmic_and_collisional_ionization"]) +
                               float(row["wolfire_ionization"]))
        asura_ionization = ordinary_ionization + (ct_ionization if element == "O" else 0.0)
        asura_recombination = recombination + (ct_recombination if element == "O" else 0.0)
        asura_scale = max(asura_ionization, asura_recombination)
        ratios = (xray / asura_ionization, xray_above_k / asura_ionization,
                  xray / max(ordinary_ionization, recombination), xray_above_k / asura_scale)
        rates = (fuv, ct_ionization, ct_recombination, recombination)
        comparison = [xray / rate if rate > 0.0 else "" for rate in rates]
        comparisons.append([*row.values(), xray, xray_above_k, *comparison, *ratios])
        maxima[element] = np.maximum(maxima[element], ratios)
    write_rows(folder / "metal-chemistry-comparison.csv", list(rows[0]) +
               ["FG20_neutral_absorption_s-1", "FG20_absorption_above_K_edge_s-1", "ratio_FUV",
                "ratio_charge_transfer_ionization", "ratio_charge_transfer_recombination", "ratio_recombination",
                "ratio_ASURA_neutral_ionization", "ratio_above_K_edge_ASURA_neutral_ionization",
                "ratio_ladder_scale_without_charge_transfer", "ratio_above_K_edge_ASURA_ladder_scale"], comparisons)
    for element, (ionization, hard_ionization, ladder_without_transfer, hard_ladder) in maxima.items():
        print(f"Neutral {element} FG20 absorption / ordinary ASURA ionization: max {ionization:.3e}, "
              f"above K edge {hard_ionization:.3e}; all energies / ladder scale without CT "
              f"{ladder_without_transfer:.3e}, above K edge / ASURA ladder scale {hard_ladder:.3e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    redshifts, energies, intensities = read_cloudy_spectrum(args.spectrum)
    atomic = AtomicModel()
    slabs = {}

    def get_slab(z):
        if z not in slabs:
            slabs[z] = Integral(energies, spectrum_at_redshift(redshifts,intensities,z), atomic)
        return slabs[z]

    passed, rows = compare_rates(args.folder, get_slab)
    cross_sections_passed = cross_section_comparison(args.folder, atomic)
    write_depth_tables(args.folder, rows, get_slab)
    compare_metal_chemistry(args.folder, get_slab)
    return 0 if passed and cross_sections_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
