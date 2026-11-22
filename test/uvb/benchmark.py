import argparse
from collections import defaultdict
import csv
from pathlib import Path
import statistics


def measure(path):
    batches = defaultdict(list)
    with path.open() as stream:
        for row in csv.DictReader(stream):
            key = row["scenario"], row["method"]
            batches[key].append(1e6 * float(row["elapsed_s"]) / int(row["calls"]))
    return {key: statistics.median(batch) for key, batch in batches.items()}


def main():
    parser = argparse.ArgumentParser(description="Summarize the existing thermochemistry benchmark's UVB CPU timings.")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    baseline = measure(args.output / "benchmark-off.csv")
    illuminated = measure(args.output / "benchmark-on.csv")
    assert baseline and baseline.keys() == illuminated.keys()
    with (args.output / "benchmark-summary.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scenario", "method", "off_us_per_cell", "on_us_per_cell", "relative_difference"])
        for (scenario, method), off in baseline.items():
            on = illuminated[scenario, method]
            change = on / off - 1.0
            writer.writerow([scenario, method, off, on, change])
            print(f"{scenario}, {method}: {off:.3f} → {on:.3f} µs per cell ({change:+.1%})")


if __name__ == "__main__":
    main()
