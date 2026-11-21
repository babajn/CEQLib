"""Run Cloudy C23.01 on the EQCLib dust-free H II equilibrium grid."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import math
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np


TEST_DIR = Path(__file__).resolve().parents[2]
SPECTRUM = TEST_DIR.parent / "data/PEGASE2/spectra_Kroupa_Z0.0200.dat"
KEY = ("nH", "Z_solar", "log_U")
ABUNDANCES = ("He_H", "C_H", "O_H", "Si_H")
ELEMENTS = ("helium", "carbon", "oxygen", "silicon")
OTHER_ELEMENTS = (
    "lithium beryllium boron nitrogen fluorine neon sodium magnesium aluminum phosphorus "
    "sulphur chlorine argon potassium calcium scandium titanium vanadium chromium manganese "
    "iron cobalt nickel copper zinc"
).split()
ION_COLUMNS = {
    "HI": "HI", "HII": "HII", "HeI": "HeI", "HeII": "HeII", "HeIII": "HeIII",
    "OI": "O1", "OII": "O2", "OIII": "O3", "OIV": "O4", "OV": "O5", "OVI": "O6",
}


def write_sed(destination):
    with SPECTRUM.open() as source:
        for _ in range(27):
            next(source)
        values = np.fromstring(source.read(), sep=" ")
    steps, count, lines = map(int, values[:3])
    wavelength = values[3:3 + count]
    epochs = values[3 + count + lines:].reshape(steps, 19 + count + lines)
    flux = epochs[epochs[:, 0] <= 30, 19:19 + count]

    # Use exactly the hydrogen-ionizing trapezoids and constants in
    # ReadPEGASESpectrum before averaging the normalized F_lambda over ages.
    energy = 6.62606876e-27 * 2.997924580e10 / (wavelength * 1e-8)
    photons = flux / energy
    ionizing = energy[1:] / 1.602177335e-12 >= 13.60
    integrals = 0.5 * (photons[:, :-1] + photons[:, 1:]) * np.diff(wavelength)
    hydrogen_photons = integrals[:, ionizing].sum(axis=1)
    shape = (flux / hydrogen_photons[:, None]).mean(axis=0) * wavelength  # nu F_nu
    if not np.isfinite(shape).all() or np.any(shape <= 0.0):
        raise ValueError("Invalid PEGASE2 continuum shape")
    with destination.open("w") as output:
        output.write("# PEGASE2 Kroupa Z=0.0200, 0–30 Myr, averaged after normalization by Q(H)\n")
        for i, (wave, intensity) in enumerate(zip(wavelength, shape)):
            units = " units Angstrom nuFnu" if i == 0 else ""
            print(f"{wave:.9e} {intensity:.12e}{units}", file=output)


def model_input(sample):
    commands = [
        'table SED "pegase2.sed"',
        f"ionization parameter {sample['log_U']:.12g}",
        f"hden {math.log10(sample['nH']):.12g}",
    ]
    for element, column in zip(ELEMENTS, ABUNDANCES):
        abundance = sample[column]
        setting = f"abundance {abundance:.12e} linear" if abundance > 0 else "off"
        commands.append(f"element {element} {setting}")
    commands.extend(f"element {element} off" for element in OTHER_ELEMENTS)
    commands.extend([
        "no molecules",
        "case B no photoionization no pdest",
        "no induced processes",
        "stop zone 1",
        "set dr 0",  # 1 cm, so the zone is optically thin to incident radiation
        "set weakheatcool 0",
        'save overview last "zone.ovr"',
        'save cooling last "zone.cool"',
    ])
    return "\n".join(commands) + "\n"


def zone_table(path):
    header, *lines = path.read_text().splitlines()
    rows = [line.split("\t") for line in lines if line.strip() and not line.startswith("#")]
    if len(rows) != 1 or any(line.startswith("#>>>>") for line in lines):
        raise ValueError(f"{path}: expected one converged Cloudy zone")
    return header.lstrip("#").split("\t"), rows[0]


def run_model(cloudy, root, index, sample):
    directory = root / f"zone_{index:03d}"
    directory.mkdir(exist_ok=True)
    shutil.copyfile(root / "pegase2.sed", directory / "pegase2.sed")
    (directory / "zone.in").write_text(model_input(sample))
    computed = subprocess.run([cloudy, "-r", "zone"], cwd=directory, timeout=120,
                              capture_output=True, text=True)
    if computed.returncode:
        raise ValueError(f"{directory}: Cloudy exited {computed.returncode}. "
                         f"See zone.out. {computed.stderr.strip()}")
    log = (directory / "zone.out").read_text(errors="replace")
    version = re.search(r"^\s*Cloudy (23\.01)\s*$", log, flags=re.MULTILINE)
    flux = re.search(r"Ion pht flx:\s*([\d.Ee+-]+)", log)
    if not version or not flux or "Cloudy ends: 1 zone," not in log or "Cloudy exited OK" not in log:
        raise ValueError(f"{directory}: expected a successful Cloudy C23.01 one-zone model")

    names, values = zone_table(directory / "zone.ovr")
    overview = dict(zip(names, map(float, values), strict=True))
    _, cooling = zone_table(directory / "zone.cool")
    heating, total_cooling = map(float, cooling[2:4])
    if (not math.isclose(overview["hden"], sample["nH"], rel_tol=1e-5)
            or not math.isclose(float(flux[1]), sample["NLyn"], rel_tol=5e-4)
            or overview["AV(extend)"] != 0 or overview["Tau912"] >= 1e-6
            or not math.isclose(heating, total_cooling, rel_tol=5e-3)):
        raise ValueError(f"{directory}: unmatched density, ionizing flux, opacity, or thermal balance")

    # Cloudy reports signed net collisional line cooling per species as a
    # fraction of total cooling. Include the summed species plus any separate
    # transitions of that ion in the same cooling ledger.
    by_ion = {2: 0.0, 3: 0.0}
    ions_found = set()
    for label, fraction in zip(cooling[5::2], cooling[6::2], strict=True):
        match = re.fullmatch(r"O\s+([23])\s+\S+", label.strip())
        if match:
            ion = int(match[1])
            by_ion[ion] += total_cooling * float(fraction)
            ions_found.add(ion)
    has_oxygen = sample["O_H"] > 0
    if has_oxygen and ions_found != set(by_ion):
        raise ValueError(f"{directory}: missing O II/O III in the cooling ledger")

    row = {field: sample[field] for field in (*KEY, *ABUNDANCES, "NLyn")}
    row.update({"cloudy_version": version[1], "incident_flux": float(flux[1]),
                "T": overview["Te"], "ne": overview["eden"],
                "H_tot": heating, "L_tot": total_cooling, "tau912": overview["Tau912"]})
    row.update({field: overview[name] if has_oxygen or not field.startswith("O") else math.nan
                for field, name in ION_COLUMNS.items()})
    row.update({"L_OII": by_ion[2], "L_OIII": by_ion[3]})
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cloudy", required=True, type=Path, help="path to the built cloudy.exe")
    parser.add_argument("--equilibrium", type=Path, default=TEST_DIR / "out/hii/equilibrium.csv")
    parser.add_argument("--work-dir", type=Path, default=TEST_DIR / "out/hii/cloudy")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    try:
        with args.equilibrium.open(newline="") as source:
            samples = [{field: float(row[field]) for field in (*KEY, *ABUNDANCES, "NLyn")}
                       for row in csv.DictReader(source)]
        if not samples or args.jobs < 1:
            raise ValueError("Expected at least one equilibrium sample and a positive job count")
        args.work_dir.mkdir(parents=True, exist_ok=True)
        write_sed(args.work_dir / "pegase2.sed")
        cloudy = args.cloudy.resolve(strict=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            models = [pool.submit(run_model, cloudy, args.work_dir, i, sample)
                      for i, sample in enumerate(samples, 1)]
            rows = []
            for i, model in enumerate(models, 1):
                rows.append(model.result())
                if i % 12 == 0 or i == len(models):
                    print(f"Cloudy: {i}/{len(models)} converged", flush=True)
        output = args.work_dir / "cloudy_equilibrium.csv"
        with output.open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=rows[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(output)
    except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
