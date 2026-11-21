import csv
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys

CASES = ("hydrogen", "recombining", "molecular", "hot")
SOLVERS = ("EQCLib", "CHIMES")
STEPS = (10, 1000, 100000)


def read_csv(path):
    with path.open(newline="") as file:
        return list(csv.DictReader(file))


def summarize(directory):
    timings = defaultdict(list)
    for row in read_csv(directory / "timings.csv"):
        key = row["case"], int(row["dt_yr"]), row["solver"]
        timings[key].append(1e6 * float(row["cpu_s"]) / int(row["calls"]))
    endpoints = {
        (row["case"], int(row["time_yr"]), row["solver"], float(row["rtol"])): row
        for row in read_csv(directory / "steps.csv")
    }

    summaries = []
    for case in CASES:
        for step in STEPS:
            summary = {"case": case, "dt_yr": step}
            for solver in SOLVERS:
                prefix = solver.lower()
                samples = timings[case, step, solver]
                if len(samples) != 5 or min(samples) <= 0:
                    raise ValueError(f"Invalid timing batches for {case} {step} {solver}")
                summary[f"{prefix}_cpu_us"] = median(samples)
                summary[f"{prefix}_min_us"] = min(samples)
                summary[f"{prefix}_max_us"] = max(samples)
                nominal = endpoints[case, step, solver, 1e-6]
                tight = endpoints[case, step, solver, 1e-8]
                for field in ("T_K", "x_HII", "f_H2", "x_e", "x_OIII", "cooling"):
                    summary[f"{prefix}_{field}"] = float(tight[field])
                summary[f"{prefix}_T_tolerance_change_pct"] = 100 * abs(float(nominal["T_K"]) / float(tight["T_K"]) - 1)
                summary[f"{prefix}_x_HII_tolerance_change"] = abs(float(nominal["x_HII"]) - float(tight["x_HII"]))
            summary["chimes_over_eqclib_cpu"] = summary["chimes_cpu_us"] / summary["eqclib_cpu_us"]
            summary["eqclib_vs_chimes_T_pct"] = 100 * (summary["eqclib_T_K"] / summary["chimes_T_K"] - 1)
            summaries.append(summary)

    with (directory / "summary.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=summaries[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)


if __name__ == "__main__":
    summarize(Path(sys.argv[1]))
