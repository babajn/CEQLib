"""Compare single-step thermal updates with H2, CO and H–O charge transfer switches."""

import argparse
import csv
from pathlib import Path

import matplotlib

if __name__ == "__main__":
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import FixedLocator
import numpy as np

from plot_style import BLUE, GRAY, ORANGE, PLOT_STYLE, SKY

ROOT = Path(__file__).resolve().parents[1]
MODES = ("h2_co_ct", "h2_co", "h2_ct", "h2", "atomic_ct", "atomic")
MOLECULES = (r"H$_2$ + CO", r"H$_2$ only", "No molecules")
COMPARISONS = (
    ("CO off", "h2_co_ct", "h2_ct"),
    (r"H$_2$ off (CO off)", "h2_ct", "atomic_ct"),
    (r"H$_2$ + CO off", "h2_co_ct", "atomic_ct"),
    ("H–O off (molecules)", "h2_co_ct", "h2_co"),
    (r"H–O off (H$_2$ only)", "h2_ct", "h2"),
    ("H–O off (atomic)", "atomic_ct", "atomic"),
    ("All three off", "h2_co_ct", "atomic"),
)
METRICS = {
    "T": ("Temperature", "% of reference"),
    "delta_U_over_U0": ("Energy change / initial energy", "pp"),
    "H2": (r"$2n_{\rm H_2}/n_{\rm H}$", "pp"),
    "CO": ("CO / C", "pp"),
    "HII": ("H II / H", "pp"),
    "OII": ("O II / O", "pp"),
    "OIII": ("O III / O", "pp"),
}


def by_case(data, mode, year):
    rows = data[(data["mode"] == mode) & (data["dt_years"] == year)]
    return np.sort(rows, order=("block", "case"))


def paired(data, year, ref, changed):
    base = by_case(data, ref, year)
    compared = by_case(data, changed, year)
    mask = (base["block"] == base["block"].min()) & (base["status"] == 0) & (compared["status"] == 0)
    return base[mask], compared[mask]


def differences(base, changed, key):
    delta = changed[key] - base[key]
    if key == "T":
        delta = delta / base[key]
    return 100 * delta


def trial_times(data, steps):
    timing = {}
    for year in steps:
        rows = [by_case(data, mode, year) for mode in MODES]
        ids = rows[0][["block", "case"]]
        if not ids.size or np.unique(ids).size != ids.size:
            raise ValueError(f"Missing or duplicate case identifiers for dt={year:g}")
        if any(not np.array_equal(row[["block", "case"]], ids) for row in rows):
            raise ValueError(f"Unmatched case identifiers for dt={year:g}")
        common = np.logical_and.reduce([(r["status"] == 0) & np.isfinite(r["cpu_us"]) for r in rows])
        blocks = np.unique(ids["block"])
        if not np.all(common.reshape(len(blocks), -1) == common.reshape(len(blocks), -1)[0]):
            raise ValueError(f"Successful case sets changed across repeated blocks for dt={year:g}")
        if not common.any():
            raise ValueError(f"No common converged cases for dt={year:g}")
        for mode, row in zip(MODES, rows):
            means = [row["cpu_us"][common & (ids["block"] == b)].mean() for b in blocks]
            first = ids["block"] == blocks[0]
            timing[year, mode] = (np.median(means), min(means), max(means),
                                  np.count_nonzero(first & (row["status"] == 0)),
                                  np.count_nonzero(common & first), np.count_nonzero(first))
    return timing


def plot_timings(steps, timing):
    fig, axes = plt.subplots(1, len(steps), figsize=(12, 4), layout="constrained", squeeze=False)
    for year, ax in zip(steps, axes[0]):
        for offset, modes, color, label in ((-0.19, MODES[::2], BLUE, "H–O on"),
                                              (0.19, MODES[1::2], SKY, "H–O off")):
            typical, low, high, *_ = np.asarray([timing[year, m] for m in modes]).T
            ax.bar(np.arange(3) + offset, typical, width=0.38, color=color, label=label,
                   yerr=[typical-low, high-typical], capsize=2)
        *_, common, total = timing[year, MODES[0]]
        ax.set(title=f"{year:g} yr · {common}/{total} common cases",
               xticks=np.arange(3), xticklabels=MOLECULES, ylabel="CPU time per update [µs]")
    axes[0, 0].legend()
    return fig


def summarize(data, steps):
    rows = []
    for year in steps:
        for label, ref, mode in COMPARISONS:
            base, changed = paired(data, year, ref, mode)
            for key in METRICS:
                delta = np.abs(differences(base, changed, key))
                stats = np.percentile(delta, [50, 95, 100]) if delta.size else [np.nan]*3
                rows.append((year, label, ref, mode, key, len(delta), *stats))
    return rows


def plot_effects(summary, year):
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.6), sharey=True, layout="constrained")
    fig.suptitle(f"Effect of switching off reactions · {year:g} yr")
    for ax, (key, (title, units)) in zip(axes.flat, METRICS.items()):
        stats = np.array([row[6:] for row in summary if row[0] == year and row[4] == key])
        y = np.arange(len(stats))
        ax.hlines(y, stats[:, 0], stats[:, 2], color=GRAY, alpha=0.5)
        for i, (marker, color, label) in enumerate((("o", BLUE, "Median"), ("s", ORANGE, "95th percentile"),
                                                  ("|", GRAY, "Maximum"))):
            ax.scatter(stats[:, i], y, marker=marker, color=color, s=25, label=label, zorder=3)
        ax.set_xscale("symlog", linthresh=1e-6)
        ax.xaxis.set_major_locator(FixedLocator([0, 1e-4, 1e-2, 1, 100, 1e4, 1e6]))
        ax.set(xlim=(0, None), title=title, xlabel=f"Absolute difference [{units}]",
               yticks=y, yticklabels=[c[0] for c in COMPARISONS])
    axes[0, 0].set_ylim(len(COMPARISONS)-0.5, -0.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    axes[1, 3].legend(handles, labels, loc="center")
    axes[1, 3].text(0.5, 0.28, "pp = percentage points", ha="center", va="center",
                    transform=axes[1, 3].transAxes, fontsize=9, color=GRAY)
    axes[1, 3].axis("off")
    return fig


def plot_temperature(data, year):
    indices = (0, 2, 3)
    conditions = [paired(data, year, *COMPARISONS[i][1:]) for i in indices]
    finite_values = [np.abs(differences(base, mode, "T")) for base, mode in conditions if base.size]
    bound = max(0.01, *(v.max() for v in finite_values))
    norm = SymLogNorm(linthresh=0.01, vmin=-bound, vmax=bound)
    fig, axes = plt.subplots(3, 3, figsize=(10, 7), sharex=True, sharey=True, layout="constrained")
    for row, (base, changed), i in zip(axes, conditions, indices):
        delta = differences(base, changed, "T")
        for col, (ax, flux) in enumerate(zip(row, (-np.inf, -5, -3))):
            select = base["log_U"] == flux
            dots = ax.scatter(base["nH"][select], base["T0"][select], c=delta[select],
                              cmap="RdBu_r", norm=norm, s=29, marker="s")
            if i == indices[0]:
                ax.set_title("No ionizing flux" if col == 0 else rf"$\log_{{10}} U={flux:g}$")
            ax.set(xscale="log", yscale="log")
        row[0].set_ylabel(COMPARISONS[i][0] + "\n" + r"Initial $T$ [K]")
    for ax in axes[-1]:
        ax.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    fig.colorbar(dots, ax=axes, label=r"$(T_{\rm off}/T_{\rm ref}-1)\times100$ [%]")
    fig.suptitle(f"Temperature change after one update · {year:g} yr")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "test/out/co/molecular_settings.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "plot/out/co/settings")
    parser.add_argument("--dt", type=float, default=1e4, help="Time step to visualize in the result plots [yr]")
    args = parser.parse_args()
    try:
        data = np.genfromtxt(args.input, names=True, delimiter=",", dtype=None, encoding="utf-8", ndmin=1)
        steps = np.unique(data["dt_years"])
        if args.dt not in steps:
            raise ValueError(f"dt={args.dt:g} is not present in {args.input}")
        timing = trial_times(data, steps)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize(data, steps)
    with (args.output_dir / "comparison_summary.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("dt_years", "comparison", "reference", "changed", "quantity", "paired_cases",
                         "median_abs_percent_or_pp", "p95_abs_percent_or_pp", "max_abs_percent_or_pp"))
        writer.writerows(summary)
    with (args.output_dir / "timing_summary.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("dt_years", "mode", "median_us", "min_us", "max_us", "converged", "common", "total"))
        writer.writerows((*key, *values) for key, values in timing.items())
    with plt.rc_context(PLOT_STYLE):
        figures = {"timings": plot_timings(steps, timing),
                   "effects": plot_effects(summary, args.dt),
                   "temperature_map": plot_temperature(data, args.dt)}
        for name, fig in figures.items():
            fig.savefig(args.output_dir / f"{name}.png", dpi=180, bbox_inches="tight")
            plt.close(fig)
    print(f"Molecular settings plots and summary CSVs written to {args.output_dir}")


if __name__ == "__main__":
    main()
