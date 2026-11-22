import argparse
import csv
import math
from pathlib import Path

import h5py
import numpy as np
from scipy.special import expn


SPECIES = ("HI", "HeI", "HeII")
THRESHOLDS = np.array((13.60, 24.59, 54.42))
PLANCK = 6.62607015e-27
ELECTRON_VOLT = 1.602176634e-12
EQCLIB_ELECTRON_VOLT = 1.602177335e-12  # src/Constants.h, for the C energy balance check
HC_EV_ANGSTROM = 12398.419843320027
LIGHT_SPEED = 2.99792458e10
ELECTRON_SPLIT = 100.0
RELATIVE_LIMIT = 0.01
SMALL_BAND_FRACTION = 1.0e-7


def read_cloudy_spectrum(path):
    with path.open() as stream:
        lines = (line.strip() for line in stream)
        lines = [line for line in lines if line and not line.startswith("#")]
    assert lines[0] == "20060612" and lines[3] == "z"
    assert lines[6] == "lambda" and lines[8] == "F_nu"
    nz, ne = map(int, lines[4:6])
    wavelengths_scale, intensity_scale = float(lines[7]), float(lines[9])
    values = np.fromstring(" ".join(lines[10:]), sep=" ")
    assert values.size == nz + ne + nz * ne
    redshifts = values[:nz]
    energies = HC_EV_ANGSTROM / (values[nz:nz + ne][::-1] * wavelengths_scale)
    intensities = intensity_scale * values[nz + ne:].reshape(nz, ne)[:, ::-1]
    return redshifts, energies, intensities


def cross_sections(energy):
    # Verner et al. (1996), Table 1. Values at and above the threshold only.
    values = np.empty((len(SPECIES), energy.size))
    for i, (e0, s0, ya, p, yw, y0, y1) in enumerate((
        (0.4298, 5.475e-14, 32.88, 2.963, 0, 0, 0),
        (13.61, 9.492e-16, 1.469, 3.188, 2.039, 0.4434, 2.136),
        (1.720, 1.369e-14, 32.88, 2.963, 0, 0, 0),
    )):
        x = energy / e0 - y0
        y = np.hypot(x, y1)
        values[i] = s0 * (np.square(x - 1) + yw**2) * y**(0.5 * p - 5.5) * (1 + np.sqrt(y / ya))**(-p)
        values[i, energy < THRESHOLDS[i]] = 0.0
    return values


class ReferenceSlab:
    def __init__(self, source_energy, source_intensity, maximum_energy=0.0):
        # Gaussian quadrature in ln(E) on every tabulated interval, with
        # breakpoints at the ionization and photoelectron edges.
        nodes, weights = np.polynomial.legendre.leggauss(16)
        lower = max(THRESHOLDS[0], source_energy[0])
        upper = min(maximum_energy, source_energy[-1]) if maximum_energy > 0.0 else source_energy[-1]
        edges = np.unique(np.r_[source_energy, THRESHOLDS, THRESHOLDS + ELECTRON_SPLIT, lower, upper])
        edges = np.log(edges[(edges >= lower) & (edges <= upper)])
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = 0.5 * np.diff(edges)
        log_energy = (centers[:, None] + widths[:, None] * nodes).ravel()
        energy = np.exp(log_energy)
        log_source_energy = np.log(source_energy)
        source_index = np.clip(np.searchsorted(log_source_energy, log_energy) - 1, 0, source_energy.size - 2)
        interp = (log_energy - log_source_energy[source_index]) / np.diff(log_source_energy)[source_index]
        left, right = source_intensity[source_index], source_intensity[source_index + 1]
        intensity = (1 - interp) * left + interp * right
        positive = (left > 0) & (right > 0)
        intensity[positive] = np.exp((1 - interp[positive]) * np.log(left[positive]) +
                                     interp[positive] * np.log(right[positive]))
        sigma = cross_sections(energy)
        self.opacity = sigma[0] + 0.1 * sigma[1]
        self.ionization = 2 * np.pi / PLANCK * intensity * sigma * (widths[:, None] * weights).ravel()
        kinetic_energy = np.maximum(0, energy - THRESHOLDS[:, None])
        power = self.ionization * kinetic_energy * ELECTRON_VOLT
        is_high = kinetic_energy >= ELECTRON_SPLIT
        self.power_low = np.where(is_high, 0.0, power)
        self.power_high = np.where(is_high, power, 0.0)

    def rates(self, top, bottom):
        transmission = expn(2, top * self.opacity) + expn(2, bottom * self.opacity)
        return self.ionization @ transmission, self.power_low @ transmission, self.power_high @ transmission

    def high_energy_fraction(self, column):
        _, power_low, power_high = self.rates(column, column)
        total = power_low + power_high
        return np.divide(power_high, total, out=np.zeros_like(total), where=total > 0.0)


def relative_error(value, reference, floor=0.0):
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"Invalid nonnegative quantity: {value}")
    scale = max(reference, floor)
    if scale == 0.0:
        return 0.0 if value == reference else math.inf
    return (value - reference) / scale


def chimes_reaction_indices(data):
    # CHIMES species identifiers: 1=H I, 4=He I, 5=He II.
    with h5py.File(data / "chimes_main_data.hdf5") as rates:
        ions = rates["photoion_euv/reactants"][:]
        return [int(np.flatnonzero(ions == number).item()) for number in (1, 4, 5)]


def chimes_photoionization(data, z, reactions):
    path = data / "FG20_cross_sections" / f"z{z:.3f}_cross_sections.hdf5"
    with h5py.File(path) as rates:
        flux = LIGHT_SPEED * rates["isotropic_photon_density"][0]
        return flux * rates["photoion_euv/sigmaPhot"][:][reactions]


def compare_rates(output, spectra, chimes_data):
    redshifts, energies, intensities = spectra
    slabs = {}
    maximum = 0.0
    split_maximum = 0.0
    split_closure = 0.0
    chimes_maximum = 0.0
    thin_comparisons = []
    reaction_indices = chimes_reaction_indices(chimes_data)
    with (output / "rates.csv").open() as actual_stream, (output / "comparison.csv").open("w") as comparison_stream:
        actual = csv.DictReader(actual_stream)
        comparison = csv.writer(comparison_stream)
        comparison.writerow(["z", "NHI_top", "NHI_bottom", "absorber", "Gamma", "Gamma_integral",
                             "Gamma_relative_error", "heat_per_absorber", "heat_integral", "heat_relative_error",
                             "Gamma_CHIMES", "CHIMES_relative_difference",
                             "power_low", "power_low_integral", "power_low_relative_error",
                             "power_high", "power_high_integral", "power_high_relative_error"])
        for row in actual:
            z, top, bottom = (float(row[key]) for key in ("z", "NHI_top", "NHI_bottom"))
            index = SPECIES.index(row["absorber"])
            if z not in slabs:
                source = np.flatnonzero(np.isclose(redshifts, z, rtol=0, atol=1e-12)).item()
                slabs[z] = (ReferenceSlab(energies, intensities[source]), chimes_photoionization(chimes_data, z, reaction_indices))
            slab, chimes = slabs[z]
            reference_ionization, reference_low, reference_high = slab.rates(top, bottom)
            low, high = reference_low[index], reference_high[index]
            power = low + high
            ionization, heating = float(row["Gamma"]), float(row["heat_per_absorber"])
            power_low, power_high = float(row["power_low"]), float(row["power_high"])
            gamma_error = relative_error(ionization, reference_ionization[index])
            heating_error = relative_error(heating, power)
            # A vanishing band can be below the precision of 1 - a stored high fraction.
            low_error = relative_error(power_low, low, SMALL_BAND_FRACTION * power)
            high_error = relative_error(power_high, high, SMALL_BAND_FRACTION * power)
            maximum = max(maximum, abs(gamma_error), abs(heating_error))
            split_maximum = max(split_maximum, abs(low_error), abs(high_error))
            split_closure = max(split_closure, abs(relative_error(power_low + power_high, heating)))
            chimes_value = chimes_error = ""
            if top == bottom == 0.0:
                chimes_value = chimes[index]
                chimes_error = ionization / chimes_value - 1.0
                chimes_maximum = max(chimes_maximum, abs(chimes_error))
                thin_comparisons.append((z, SPECIES[index], ionization, chimes_value, chimes_error))
            comparison.writerow([z, top, bottom, SPECIES[index], ionization, reference_ionization[index], gamma_error,
                                 heating, power, heating_error, chimes_value, chimes_error,
                                 power_low, low, low_error, power_high, high, high_error])
    print(f"FG20 independent integration: maximum relative difference {maximum:.3%} (limit 1%)")
    print(f"FG20 power split: maximum relative difference {split_maximum:.3%} (limit 1%, "
          f"denominator floor {SMALL_BAND_FRACTION:g} of total power), closure {split_closure:.3e}")
    for z, absorber, ionization, reference, error in thin_comparisons:
        print(f"CHIMES z={z:g} {absorber}: EQCLib {ionization:.6e}, CHIMES {reference:.6e}, difference {error:+.2%}")
    print(f"CHIMES optically thin: maximum relative difference {chimes_maximum:.3%} (limit 3%)")
    passed = (maximum <= RELATIVE_LIMIT and split_maximum <= RELATIVE_LIMIT and
              split_closure <= 5.0e-12 and chimes_maximum <= 0.03)
    return passed, {z: slab for z, (slab, _) in slabs.items()}


def compare_heat_diagnostics(output, slabs):
    expected = {}
    for z, slab in slabs.items():
        for column in (0.0, 1.0e19, 1.0e20, 1.0e21):
            for species, fraction in zip(SPECIES, slab.high_energy_fraction(column)):
                expected[z, column, species] = fraction
    maximum = 0.0
    with (output / "UVBackgroundRates.dat").open() as stream:
        next(stream)
        for row in stream:
            z, column, species, _, _, fraction = row.split()
            key = float(z), float(column), species
            if key in expected:
                difference = abs(float(fraction) - expected.pop(key))
                maximum = max(maximum, difference)
    print(f"Photoelectron power fractions at K >= 100 eV: maximum absolute difference {maximum:.3e} (limit 1e-3)")
    return not expected and maximum <= 1.0e-3


def compare_secondary(output, spectra, slabs):
    redshifts, energies, intensities = spectra
    uncapped_and_capped = {0.0, 300.0}
    expected_grid = {(z, column, x, cap)
                     for z in (0.0, 2.0)
                     for column in (1.0e19, 1.0e20, 1.0e21)
                     for x in (0.001, 0.01, 0.1)
                     for cap in uncapped_and_capped}
    quantities = ["P_total", "P_high", "H_UVB", "Q_other", "Gamma_HI", "Gamma_HeI",
                  "Gamma2_HI", "Gamma2_HeI"]
    reference_slabs = {(z, 0.0): slab for z, slab in slabs.items()}
    maximum = 0.0
    closure = 0.0
    anchors = []
    with (output / "secondary.csv").open() as actual_stream, (output / "secondary-comparison.csv").open("w") as dest:
        rows = csv.DictReader(actual_stream)
        comparison = csv.writer(dest)
        comparison.writerow(rows.fieldnames + [name + "_integral" for name in quantities] +
                             ["H_UVB_over_legacy", "Gamma_HI_over_primary", "Gamma_HeI_over_primary"])
        for row in rows:
            z, column, x, cap = (float(row[key]) for key in ("z", "NHI_per_side", "xHII", "energy_max_eV"))
            case = (z, column, x, cap)
            if case not in expected_grid:
                raise ValueError(f"Unknown or repeated deposition case: {case}")
            expected_grid.remove(case)
            if (z, cap) not in reference_slabs:
                source = np.flatnonzero(np.isclose(redshifts, z, rtol=0.0, atol=1.0e-12)).item()
                reference_slabs[z, cap] = ReferenceSlab(energies, intensities[source], cap)
            primary, low, high = reference_slabs[z, cap].rates(column, column)
            # This fixture has nH=1, fully neutral He/H=0.1, and no H2.
            absorber_densities = np.array([1.0 - x, 0.1, 0.0])
            low_power = low @ absorber_densities
            high_power = high @ absorber_densities
            heat_fraction = 0.9971 * (1.0 - (1.0 - x**0.2663)**1.3163)
            hydrogen_fraction = 0.3908 * (1.0 - x**0.4092)**1.7592
            helium_fraction = 0.0554 * (1.0 - x**0.4614)**1.666
            secondary = (high_power * np.array([hydrogen_fraction, helium_fraction]) /
                         (THRESHOLDS[:2] * ELECTRON_VOLT * absorber_densities[:2]))
            expected = [low_power + high_power, high_power, low_power + heat_fraction * high_power,
                        (1.0 - heat_fraction - hydrogen_fraction - helium_fraction) * high_power,
                        *(primary[:2] + secondary), *secondary]
            actual = [float(row[key]) for key in quantities]
            maximum = max(maximum, *(abs(relative_error(value, ref)) for value, ref in zip(actual, expected)))
            accounted = actual[2] + actual[3] + ((THRESHOLDS[:2] * EQCLIB_ELECTRON_VOLT) *
                                                 absorber_densities[:2]) @ np.asarray(actual[6:8])
            closure = max(closure, abs(relative_error(accounted, actual[0])))
            heating_ratio = actual[2] / actual[0]
            ionization_ratios = [total / (total - new) for total, new in zip(actual[4:6], actual[6:8])]
            comparison.writerow([*row.values(), *expected, heating_ratio, *ionization_ratios])
            if (z, column, x) == (0.0, 1.0e20, 0.01):
                high_fraction = actual[1] / actual[0]
                if cap == 0.0:
                    anchors.append(abs(high_fraction - 0.9927) <= 0.001 and
                                   abs(heating_ratio - 0.3705) <= 0.003 and
                                   abs(ionization_ratios[0] - 19.42) <= 0.2)
                print(f"FG20 SvS at z=0, NHI=1e20 per side, xHII=.01, Emax={cap or 'full'}: "
                      f"fast power fraction {high_fraction:.3%}, H_UVB/legacy {heating_ratio:.4f}, "
                      f"total/primary HI {ionization_ratios[0]:.3f}, HeI {ionization_ratios[1]:.3f}")
    print(f"FG20 independent SvS deposition: maximum relative difference {maximum:.6%} (limit 1%), "
          f"energy closure {closure:.3e}")
    return (not expected_grid and anchors == [True] and maximum <= RELATIVE_LIMIT and
            closure <= 5.0e-12)


def compare_equilibrium(output):
    maximum = 0.0
    with (output / "equilibrium.csv").open() as stream:
        for row in csv.DictReader(stream):
            for species in ("xe", "xHII", "xHeII", "xHeIII"):
                actual, expected = float(row[species]), float(row[species + "_reference"])
                if not math.isfinite(actual) or not math.isfinite(expected):
                    raise ValueError(f"Invalid equilibrium for {species}: {row}")
                maximum = max(maximum, abs(actual - expected))
    print(f"H/He fixed-temperature reference balance: maximum absolute fraction difference {maximum:.3e}")
    return maximum <= 2.0e-7


def main():
    parser = argparse.ArgumentParser(description="Compare EQCLib with an independent FG20 integration and CHIMES.")
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("chimes_data", type=Path)
    parser.add_argument("output", type=Path)
    options = parser.parse_args()
    spectra = read_cloudy_spectrum(options.spectrum)
    rates_ok, slabs = compare_rates(options.output, spectra, options.chimes_data)
    heat_ok = compare_heat_diagnostics(options.output, slabs)
    secondary_ok = compare_secondary(options.output, spectra, slabs)
    equilibrium_ok = compare_equilibrium(options.output)
    return 0 if rates_ok and heat_ok and secondary_ok and equilibrium_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
