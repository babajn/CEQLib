import csv
import io
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def run_cell(program, config, settings):
    config.write_text(settings, encoding="utf-8")
    result = subprocess.run([program, config], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    rows = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(rows) == 2, len(rows)
    return rows


def check_temperature(rows, initial):
    assert float(rows[0]["T_K"]) == initial
    assert math.isclose(float(rows[1]["T_K"]), initial, rel_tol=1.0e-7)


def main():
    program = Path(sys.argv[1]).resolve()
    step = "steps = 1\ndt_years = 1e-6\n"
    with TemporaryDirectory() as directory:
        config = Path(directory) / "cell.conf"
        atomic = step + "T = 8000\nx_HII = 0.5\n"
        equilibrium = run_cell(program, config, atomic)
        assert equilibrium[0]["chemistry_mode"] == "equilibrium"
        check_temperature(equilibrium, 8000.0)
        assert float(equilibrium[0]["HII"]) < 0.01
        assert equilibrium == run_cell(program, config, step + "chemistry_mode = equilibrium\nx_HII = 0\n")
        assert equilibrium == run_cell(program, config, step + "x_HII = 0\nx_H2 = 0.8\n")
        nonequilibrium = run_cell(program, config, "chemistry_mode = nonequilibrium\n" + atomic)
        assert float(nonequilibrium[0]["HII"]) == 0.5

        # At this grid point, the solar CIE metal term is the solar table entry
        # minus the primordial entry, multiplied by n_e * n_H = 100 cm^-6.
        hot = run_cell(program, config, step + "chemistry_mode = nonequilibrium\nnH = 10\n"
                       f"T = {10.0**4.4:.17g}\nx_HII = 1\n")
        expected_metal_cooling = 100.0 * (10.0**-21.88 - 10.0**-22.26)
        total_cooling = float(hot[0]["L_tot"])
        assert expected_metal_cooling <= total_cooling <= 1.05 * expected_metal_cooling, total_cooling

        irradiated = step + "T = 100\nnH = 100\nG0 = 1000\nNt = 2e20\nZd = 0.02\nx_HII = 0\n"
        molecular = run_cell(program, config, irradiated + "x_H2 = 0.8\n")
        check_temperature(molecular, 100.0)
        assert float(molecular[0]["H2"]) < 0.01
        assert molecular == run_cell(program, config, irradiated)
        molecular_nonequilibrium = run_cell(program, config,
                                           "chemistry_mode = nonequilibrium\n" + irradiated + "x_H2 = 0.8\n")
        assert float(molecular_nonequilibrium[0]["H2"]) == 0.8
    print("Cell runner initialization checks passed.")


if __name__ == "__main__":
    main()
