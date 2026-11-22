"""Compare single-zone metal ion fractions with Cloudy c25 and identical EUV spectra."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np

from compare import read_cloudy_spectrum
from metal_compare import log_interpolate, spectrum_at_redshift


ELEMENTS = {"C": "carbon", "O": "oxygen", "Si": "silicon"}
ION_STAGES = {"C": 7, "O": 9, "Si": 15}
OTHER_ELEMENTS = ("lithium", "beryllium", "boron", "nitrogen", "fluorine", "neon", "sodium",
                  "magnesium", "aluminium", "phosphorus", "sulphur", "chlorine", "argon", "potassium",
                  "calcium", "scandium", "titanium", "vanadium", "chromium", "manganese", "iron",
                  "cobalt", "nickel", "copper", "zinc")
REDSHIFTS = (0, 2)
LOG_DENSITIES = tuple(i / 2 for i in range(-6, 1))
HYDROGEN_EDGE = 13.60
RYDBERG_EV = 13.605693122994
NORMALIZATION_RYD = 2.0


def intensity_at(energies, intensities, energy):
    index = np.searchsorted(energies, energy)
    if energies[index] == energy:
        return intensities[index]
    fraction = math.log(energy / energies[index-1]) / math.log(energies[index] / energies[index-1])
    return log_interpolate(intensities[index-1:index], intensities[index:index+1], fraction)[0]


def cloudy_sed(output, energies, intensities):
    above = energies > HYDROGEN_EDGE
    edge_intensity = intensity_at(energies, intensities, HYDROGEN_EDGE)
    energies_out = np.r_[HYDROGEN_EDGE, energies[above]] / RYDBERG_EV
    intensities_out = np.r_[edge_intensity, intensities[above]]
    # Cloudy interpolates in log space. Use a numerically negligible positive
    # value wherever the original table contains a zero.
    floor = np.min(intensities_out[intensities_out > 0.0]) * 1.0e-30
    np.savetxt(output, np.c_[energies_out, np.maximum(intensities_out, floor)], fmt="%.16e")
    return math.log10(4.0 * math.pi * intensity_at(energies, intensities, NORMALIZATION_RYD * RYDBERG_EV))


def run_cloudy(executable, working, name, redshift, log_density, energies, intensities):
    work = working / f"{name}-{redshift}-{log_density:g}"
    work.mkdir()
    log_flux = cloudy_sed(work / "euv.sed", energies, intensities)
    commands = [
        f"title {name} EUV z={redshift}, log(nH)={log_density:g}",
        'table SED "euv.sed"',
        f"f(nu) = {log_flux:.14g} at {NORMALIZATION_RYD:.8g} Ryd",
        f"hden {log_density:g}",
        "constant temperature 10000 K",
        "stop zone 1",
        "set dr 0",
        "no molecules",
        "element abundance helium -1.00",
        "element abundance carbon -3.57",
        "element abundance oxygen -3.31",
        "element abundance silicon -4.49",
    ]
    commands += [f"element {element} off" for element in OTHER_ELEMENTS]
    commands += [f'save element {element} last "{symbol}.dat"' for symbol, element in ELEMENTS.items()]
    model = subprocess.run([executable], cwd=work, input="\n".join(commands)+"\n", text=True,
                           capture_output=True, timeout=180)
    if model.returncode != 0 or "Cloudy ends: 1 zone" not in model.stdout:
        raise RuntimeError(f"Cloudy failed for {name}, z={redshift}, log(nH)={log_density}\n"
                           f"{model.stdout[-5000:]}\n{model.stderr}")
    fractions = {}
    for symbol in ELEMENTS:
        values = np.loadtxt(work / f"{symbol}.dat", ndmin=2)
        if values.shape[0] != 1 or not np.isfinite(values).all() or np.any(values[0,1:] < 0.0):
            raise ValueError(f"Unexpected Cloudy one-zone output for {work.name}, {symbol}")
        # The columns following the ion ladder contain excited levels and
        # molecular populations, which cannot be added to its ion fractions.
        stages = values[0,1:1+ION_STAGES[symbol]]
        if not math.isclose(sum(stages), 1.0, abs_tol=3.0e-3):
            raise ValueError(f"Cloudy {symbol} stage fractions sum to {sum(stages):g}, not unity")
        fractions[symbol] = stages[:3], sum(stages[3:])
    return name, redshift, log_density, fractions


def run_eqclib(executable, source, work):
    work.mkdir()
    (work / "fg20_cloudy.ascii").symlink_to(source)
    subprocess.run([executable, work, "--fractions"], cwd=work, check=True, timeout=180)
    with (work / "metal-photoionization-eqc-fractions.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    return {(int(row["z"]), round(math.log10(float(row["nH"])), 2), row["element"]):
            np.array([float(row[part]) for part in ("neutral", "singly", "doubly")]) for row in rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cloudy_directory", type=Path)
    parser.add_argument("fg20_spectrum", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("reference_executable", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    cloudy = (args.cloudy_directory / "source" / "cloudy.exe").resolve(strict=True)
    spectra = {"HM12": (args.cloudy_directory / "data" / "hm12_galaxy.ascii").resolve(strict=True),
               "FG20": args.fg20_spectrum.resolve(strict=True)}
    eqclib = args.reference_executable.resolve(strict=True)

    with TemporaryDirectory(prefix="eqclib-cloudy-metal-") as temporary:
        work = Path(temporary)
        reference = {name: run_eqclib(eqclib, source, work / f"eqclib-{name}") for name, source in spectra.items()}
        cases = []
        for name, source in spectra.items():
            redshifts, energies, intensities = read_cloudy_spectrum(source)
            for z in REDSHIFTS:
                j = spectrum_at_redshift(redshifts, intensities, z)
                cases.extend((cloudy, work, name, z, density, energies, j) for density in LOG_DENSITIES)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            outputs = list(pool.map(lambda case: run_cloudy(*case), cases))

    header = ["spectrum", "z", "nH", "element", "EQCLib_I", "EQCLib_II", "EQCLib_III",
              "Cloudy_I", "Cloudy_II", "Cloudy_III", "Cloudy_IV_and_above", "max_stage_absolute_difference"]
    with (args.output_directory / "metal-photoionization-cloudy.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for name, z, log_density, fractions in outputs:
            for symbol, (cloudy_fractions, stages_above) in fractions.items():
                eqc_fractions = reference[name][z, log_density, symbol]
                difference = max(abs(eqc_fractions - cloudy_fractions))
                writer.writerow([name, z, 10**log_density, symbol, *eqc_fractions, *cloudy_fractions,
                                 stages_above, difference])
                print(f"Cloudy {name} z={z}, nH={10**log_density:.4g}, {symbol}: "
                      f"max |Δ(I–III)|={difference:.4g}, Cloudy IV+={stages_above:.4g}")


if __name__ == "__main__":
    main()
