"""Plot and summarize the shielded EQCLib–CHIMES grid comparison.

From the EQCLib root, run python3 plot/plot_chimes_comparison.py.
Use --data-dir and --out-dir to select directories for the CSVs and figures.
"""
import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_style import BLACK, GRAY, BLUE, ORANGE, GREEN, SKY, PURPLE, PLOT_STYLE

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DIR = ROOT / "test/out/chimes/shielded"

STYLE = {
    "chimes_full": dict(color=BLACK, ls="-", label="CHIMES (full, 157 species)"),
    "chimes_matched": dict(color=GRAY, ls="--", label="CHIMES (H, He, C, O, Si)"),
    "chimes_full_euv": dict(color=PURPLE, ls=":", label="CHIMES full + B87 EUV (control)"),
    "eqclib": dict(color=BLUE, ls="-", label="EQCLib equilibrium"),
    "eqclib_eq": dict(color=BLUE, ls="-", label="EQCLib equilibrium"),
    "eqclib_co": dict(color=SKY, ls="-.", label="EQCLib equilibrium + CO"),
    "eqclib_co_eq": dict(color=SKY, ls="-.", label="EQCLib equilibrium + CO"),
    "eqclib_noneq_1e-4": dict(color=ORANGE, ls="-", label="EQCLib non-eq (tol 1e-4)"),
    "eqclib_noneq_1e-6": dict(color=GREEN, ls=":", label="EQCLib non-eq (tol 1e-6)"),
}


def read(path):
    with open(path) as f:
        return [{k: (v if k in ("code", "method") else float(v)) for k, v in row.items()}
                for row in csv.DictReader(f)]


def key(nH, T):
    return (round(math.log10(nH), 2), round(math.log10(T), 2))


def load_equilibrium(data_dir):
    out = {}
    for name in ("eqclib", "eqclib_co", "chimes_matched", "chimes_full", "chimes_full_euv"):
        path = os.path.join(data_dir, f"eq_{name}.csv")
        if os.path.exists(path):
            out[name] = {key(r["nH"], r["T0"]): r for r in read(path)}
    return out


def grid(table, column):
    keys = sorted(table)
    lognH = sorted({k[0] for k in keys})
    logT = sorted({k[1] for k in keys})
    z = np.full((len(logT), len(lognH)), np.nan)
    for i, t in enumerate(logT):
        for j, n in enumerate(lognH):
            r = table.get((n, t))
            if r is None or r["status"] != 0:
                continue
            z[i, j] = column(r) if callable(column) else r[column]
    return np.array(lognH), np.array(logT), z


def plot_profiles(eq, out_dir):
    densities = (-1.0, 0.0, 1.0, 2.0, 3.0)
    quantities = (("xe", "electron fraction $x_e$", True),
                  ("fH2", r"$2n_{\rm H_2}/n_{\rm H}$", True),
                  ("cool", r"cooling $\Lambda$ [erg cm$^{-3}$ s$^{-1}$]", True),
                  ("heat", r"heating $\Gamma$ [erg cm$^{-3}$ s$^{-1}$]", True))
    fig, axes = plt.subplots(len(quantities), len(densities), figsize=(17, 12), sharex=True, sharey="row")
    for j, ln in enumerate(densities):
        for i, (col, label, logy) in enumerate(quantities):
            ax = axes[i, j]
            for name, table in eq.items():
                pts = sorted((k[1], r[col]) for k, r in table.items() if k[0] == ln and r["status"] == 0)
                if not pts:
                    continue
                x, y = zip(*pts)
                y = np.where(np.array(y) > 0, y, np.nan)
                ax.plot(x, y, **{k: v for k, v in STYLE[name].items() if k != "label"},
                        label=STYLE[name]["label"])
            if logy:
                ax.set_yscale("log")
            if i == 0:
                ax.set_title(rf"$n_{{\rm H}} = 10^{{{ln:g}}}$ cm$^{{-3}}$")
            if j == 0:
                ax.set_ylabel(label)
            if i == len(quantities)-1:
                ax.set_xlabel(r"$\log_{10} T$ [K]")
    axes[1, 0].set_ylim(1e-8, 2)
    axes[0, 0].set_ylim(1e-7, 2)
    axes[0, 0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Fixed-temperature chemical equilibrium (G0 = 1.7, no > 13.6 eV photons, CR 1.8e-17 s$^{-1}$, "
                 r"$N_{\rm H} = 10^{20} n_{\rm H}^{1/2}$ cm$^{-2}$)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "equilibrium_profiles.png"), dpi=130)
    plt.close(fig)


def plot_ratio_maps(eq, out_dir, reference="chimes_full"):
    ref = eq[reference]
    panels = (("eqclib", "xe", r"$x_e$"), ("eqclib", "cool", r"$\Lambda$"), ("eqclib", "heat", r"$\Gamma$"),
              ("eqclib_co", "cool", r"$\Lambda$ (EQCLib + CO)"))
    fig, axes = plt.subplots(1, len(panels), figsize=(18, 4.6), sharey=True)
    for ax, (name, col, label) in zip(axes, panels):
        n, t, a = grid(eq[name], col)
        _, _, b = grid(ref, col)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.log10(a/b)
        m = ax.pcolormesh(n, t, z, cmap="RdBu_r", vmin=-1, vmax=1, shading="nearest")
        ax.set_title(f"log10 EQCLib/CHIMES: {label}")
        ax.set_xlabel(r"$\log_{10} n_{\rm H}$ [cm$^{-3}$]")
        ax.grid(False)
    axes[0].set_ylabel(r"$\log_{10} T$ [K]")
    fig.colorbar(m, ax=axes, shrink=0.9, label="dex")
    fig.savefig(os.path.join(out_dir, "equilibrium_ratio_maps.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def thermal_equilibria(table):
    """Temperatures where heat = cool along each density column (log-linear interpolation)."""
    roots = []
    by_density = defaultdict(list)
    for (n, t), r in table.items():
        if r["status"] == 0 and r["heat"] > 0 and r["cool"] > 0:
            by_density[n].append((t, math.log10(r["cool"]/r["heat"])))
    for n, pts in by_density.items():
        pts.sort()
        for (t0, f0), (t1, f1) in zip(pts, pts[1:]):
            if f0 == 0.0:
                roots.append((n, t0, f1 > f0))
            elif f0*f1 < 0:
                roots.append((n, t0 + (t1-t0)*f0/(f0-f1), f1 > f0))
    return roots


def plot_thermal_equilibrium(eq, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for name, table in eq.items():
        roots = thermal_equilibria(table)
        style = STYLE[name]
        for stable, marker in ((True, "o"), (False, "x")):
            pts = sorted((n, t) for n, t, s in roots if s == stable)
            if not pts:
                continue
            n, t = map(np.array, zip(*pts))
            lab = style["label"] if stable else None
            axes[0].plot(n, t, marker, color=style["color"], ms=4, label=lab, mfc="none" if not stable else style["color"])
            axes[1].plot(n, n + t, marker, color=style["color"], ms=4, mfc="none" if not stable else style["color"])
    axes[0].set_ylabel(r"$\log_{10} T_{\rm eq}$ [K]")
    axes[1].set_ylabel(r"$\log_{10} P/k$ [K cm$^{-3}$] (= $n_{\rm H} T$)")
    for ax in axes:
        ax.set_xlabel(r"$\log_{10} n_{\rm H}$ [cm$^{-3}$]")
    axes[0].legend(fontsize=8)
    axes[0].set_title(r"Thermal equilibrium $\Gamma = \Lambda$ (filled: stable, open: unstable)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "thermal_equilibrium.png"), dpi=130)
    plt.close(fig)
    return {name: thermal_equilibria(table) for name, table in eq.items()}


def load_timing(data_dir):
    rows = []
    for name in ("eqclib", "eqclib_co", "chimes_matched", "chimes_full"):
        path = os.path.join(data_dir, f"timing_{name}.csv")
        if os.path.exists(path):
            rows += read(path)
    return rows


def plot_timing(rows, out_dir):
    methods = [m for m in STYLE if any(r["method"] == m for r in rows)]
    years = sorted({r["dt_yr"] for r in rows})
    temps = sorted({r["T0"] for r in rows})
    fig, axes = plt.subplots(1, len(years), figsize=(16, 4.6), sharey=True)
    for ax, dt in zip(axes, years):
        for m in methods:
            med, lo, hi, xs = [], [], [], []
            for T in temps:
                v = [r["cpu_us_median"] for r in rows if r["method"] == m and r["dt_yr"] == dt and r["T0"] == T
                     and r["status"] == 0]
                if v:
                    xs.append(T)
                    med.append(np.median(v)); lo.append(np.min(v)); hi.append(np.max(v))
            if xs:
                st = STYLE[m]
                ax.plot(xs, med, marker="o", ms=3, color=st["color"], ls=st["ls"], label=st["label"])
                ax.fill_between(xs, lo, hi, color=st["color"], alpha=0.10)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"dt = {dt:g} yr")
        ax.set_xlabel(r"initial $T$ [K]")
    axes[0].set_ylabel(r"CPU time per cell update [$\mu$s]")
    axes[0].legend(fontsize=8)
    fig.suptitle("Single-core CPU time per isochoric update (line: median over 13 densities, band: min-max)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "timing.png"), dpi=130)
    plt.close(fig)


def timing_summary(rows):
    methods = [m for m in STYLE if any(r["method"] == m for r in rows)]
    years = sorted({r["dt_yr"] for r in rows})
    lines = ["| method | dt [yr] | success | median us | p90 us | mean us | Chimes-full / method (median of per-cell ratios) |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    ref = {(r["nH"], r["T0"], r["dt_yr"]): r["cpu_us_median"] for r in rows
           if r["method"] == "chimes_full" and r["status"] == 0}
    for m in methods:
        for dt in years:
            sel = [r for r in rows if r["method"] == m and r["dt_yr"] == dt]
            ok = [r for r in sel if r["status"] == 0]
            v = np.array([r["cpu_us_median"] for r in ok])
            ratios = [ref[(r["nH"], r["T0"], dt)]/r["cpu_us_median"] for r in ok if (r["nH"], r["T0"], dt) in ref]
            lines.append(f"| {m} | {dt:g} | {len(ok)}/{len(sel)} | {np.median(v):.1f} | {np.percentile(v, 90):.1f} | "
                         f"{v.mean():.1f} | {np.median(ratios):.1f} |")
    return "\n".join(lines)


def speed_ratios(rows):
    """CHIMES time / EQCLib time on cells where both succeeded, aggregated two ways."""
    t = {(r["method"], r["nH"], r["T0"], r["dt_yr"]): r["cpu_us_median"] for r in rows if r["status"] == 0}
    years = sorted({r["dt_yr"] for r in rows})
    eqc = [m for m in ("eqclib_eq", "eqclib_co_eq", "eqclib_noneq_1e-4", "eqclib_noneq_1e-6")
           if any(r["method"] == m for r in rows)]
    lines = ["CHIMES CPU time / EQCLib CPU time (>1: EQCLib faster). "
             "'sum' = ratio of summed times, 'median' = median of per-cell ratios.", "",
             "| CHIMES | EQCLib | " + " | ".join(f"dt={dt:g} yr sum / median" for dt in years) + " |",
             "| --- | --- |" + " ---: |"*len(years)]
    for ref in ("chimes_matched", "chimes_full"):
        for m in eqc:
            cells = []
            for dt in years:
                pairs = [(t[(ref,) + k[1:]], v) for k, v in t.items()
                         if k[0] == m and k[3] == dt and (ref,) + k[1:] in t]
                a, b = map(np.array, zip(*pairs))
                cells.append(f"{a.sum()/b.sum():.1f} / {np.median(a/b):.2f}")
            lines.append(f"| {ref} | {m} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def plot_timing_accuracy(rows, out_dir):
    """Final temperature after one update relative to CHIMES (full)."""
    ref = {(r["nH"], r["T0"], r["dt_yr"]): r["T"] for r in rows if r["method"] == "chimes_full" and r["status"] == 0}
    methods = [m for m in ("chimes_matched", "eqclib_eq", "eqclib_co_eq", "eqclib_noneq_1e-4", "eqclib_noneq_1e-6")
               if any(r["method"] == m for r in rows)]
    years = sorted({r["dt_yr"] for r in rows})
    fig, axes = plt.subplots(1, len(years), figsize=(16, 4.4), sharey=True)
    for ax, dt in zip(axes, years):
        for m in methods:
            pts = [(r["T0"]*(1+0.04*(methods.index(m)-2)), math.log10(r["T"]/ref[k]))
                   for r in rows if r["method"] == m and r["dt_yr"] == dt and r["status"] == 0
                   for k in [(r["nH"], r["T0"], dt)] if k in ref and r["T"] > 0]
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, "o", ms=3, alpha=0.6, color=STYLE[m]["color"], label=STYLE[m]["label"])
        ax.axhline(0, color=BLACK, lw=0.8)
        ax.set_xscale("log")
        ax.set_title(f"dt = {dt:g} yr")
        ax.set_xlabel(r"initial $T$ [K]")
    axes[0].set_ylabel(r"$\log_{10}(T_{\rm final}/T_{\rm final}^{\rm CHIMES\,full})$")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "timing_final_temperature.png"), dpi=130)
    plt.close(fig)


def plot_cooling(data_dir, out_dir):
    rows = []
    for name in ("eqclib", "chimes_matched", "chimes_full"):
        path = os.path.join(data_dir, f"cooling_{name}.csv")
        if os.path.exists(path):
            rows += read(path)
    if not rows:
        return None
    cases = sorted({(r["T0"], r["nH"]) for r in rows}, key=lambda c: (-c[0], c[1]))
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), sharex=True)
    stats = []
    for ax, (T0, nH) in zip(axes.flat, cases):
        for m in STYLE:
            pts = sorted((r["t_yr"], r["T"]) for r in rows if r["method"] == m and r["T0"] == T0 and r["nH"] == nH
                         and r["status"] == 0 and r["t_yr"] > 0)
            if not pts:
                continue
            t, T = zip(*pts)
            ax.plot(t, T, color=STYLE[m]["color"], ls=STYLE[m]["ls"], label=STYLE[m]["label"])
            cpu = sum(r["cpu_us"] for r in rows if r["method"] == m and r["T0"] == T0 and r["nH"] == nH)
            nfail = sum(1 for r in rows if r["method"] == m and r["T0"] == T0 and r["nH"] == nH and r["status"] != 0)
            stats.append((T0, nH, m, cpu, len(pts), nfail))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(rf"$T_0 = {T0:.0e}$ K, $n_{{\rm H}} = {nH:g}$ cm$^{{-3}}$")
        ax.set_xlim(1, 1e9)
    for ax in axes[-1]:
        ax.set_xlabel("time [yr]")
    for ax in axes[:, 0]:
        ax.set_ylabel("T [K]")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Isochoric cooling from each code's own fixed-T equilibrium (181 log-spaced steps, 1 yr - 1 Gyr)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "cooling_histories.png"), dpi=130)
    plt.close(fig)
    return stats


def equilibrium_stats(eq):
    """Median |log10 ratio| and fraction within a factor 2, by temperature regime."""
    ref = eq["chimes_full"]
    regimes = (("T < 1e4 K", lambda t: t < 4.0), ("1e4 <= T < 1e5 K", lambda t: 4.0 <= t < 5.0),
               ("T >= 1e5 K", lambda t: t >= 5.0))
    lines = ["| model | quantity | regime | N | median abs(log10 ratio) | within x2 |",
             "| --- | --- | --- | ---: | ---: | ---: |"]
    for name in ("eqclib", "eqclib_co", "chimes_matched", "chimes_full_euv"):
        if name not in eq:
            continue
        for col in ("xe", "fH2", "cool", "heat"):
            for label, test in regimes:
                d = []
                for k, r in eq[name].items():
                    s = ref.get(k)
                    if not test(k[1]) or s is None or r["status"] != 0 or s["status"] != 0:
                        continue
                    if col == "fH2" and max(r[col], s[col]) < 1e-3:
                        continue
                    if r[col] > 0 and s[col] > 0:
                        d.append(abs(math.log10(r[col]/s[col])))
                if d:
                    d = np.array(d)
                    lines.append(f"| {name} | {col} | {label} | {len(d)} | {np.median(d):.2f} | "
                                 f"{np.mean(d < math.log10(2)):.0%} |")
    return "\n".join(lines)


COURANT_METHODS = ("eqclib_eq", "eqclib_noneq_1e-4", "eqclib_noneq_1e-6", "chimes_matched", "chimes_full")


def load_courant(data_dir):
    rows = []
    for name in ("eqclib", "chimes_matched", "chimes_full"):
        path = os.path.join(data_dir, f"courant_{name}.csv")
        if os.path.exists(path):
            rows += read(path)
    return {(r["method"], r["mass_msun"], r["nH"], r["T0"]): r for r in rows if r["status"] == 0}


def courant_summary(table):
    """Cost and one-step final-temperature differences on cells where every method succeeded."""
    lines = []
    regimes = (("all", 0.0, 1e99), ("T0 < 8e3 K", 0.0, 8e3), ("8e3 <= T0 < 3e4 K", 8e3, 3e4), ("T0 >= 1e5 K", 1e5, 1e99))
    pairs = (("eqclib_eq", "eqclib_noneq_1e-6", "EQCLib eq vs non-eq (tol 1e-6)"),
             ("eqclib_noneq_1e-6", "chimes_matched", "EQCLib non-eq (tol 1e-6) vs CHIMES matched"),
             ("chimes_matched", "chimes_full", "CHIMES matched vs CHIMES full"),
             ("eqclib_eq", "chimes_full", "EQCLib eq vs CHIMES full"),
             ("eqclib_noneq_1e-4", "eqclib_noneq_1e-6", "EQCLib non-eq tolerance 1e-4 vs 1e-6"))
    for mass in sorted({k[1] for k in table}):
        cells = sorted({k[2:] for k in table if k[1] == mass
                        and all((m, mass) + k[2:] in table for m in COURANT_METHODS)})
        dt = np.array([table[("chimes_full", mass) + c]["dt_yr"] for c in cells])
        lines += [f"### Particle mass {mass:g} Msun", "",
                  f"{len(cells)} cells; Courant dt median {np.median(dt):.3g} yr, range {dt.min():.3g}-{dt.max():.3g} yr, "
                  f"{np.mean(dt >= 1e6):.0%} at the 1 Myr cap.", "",
                  "| method | median us | mean us | p90 us | CPU s per 1e6 updates | CHIMES full / method (sum) |",
                  "| --- | ---: | ---: | ---: | ---: | ---: |"]
        ref = np.array([table[("chimes_full", mass) + c]["cpu_us_median"] for c in cells])
        for m in COURANT_METHODS:
            v = np.array([table[(m, mass) + c]["cpu_us_median"] for c in cells])
            lines.append(f"| {m} | {np.median(v):.1f} | {v.mean():.1f} | {np.percentile(v, 90):.1f} | "
                         f"{v.mean():.0f} | {ref.sum()/v.sum():.1f} |")
        lines += ["", "abs(log10) difference of the final temperature after one Courant step "
                  "(median dex / fraction within a factor 2):", "",
                  "| comparison | " + " | ".join(r[0] for r in regimes) + " |",
                  "| --- |" + " ---: |"*len(regimes)]
        for a, b, label in pairs:
            out = []
            for _, lo, hi in regimes:
                d = np.array([abs(math.log10(table[(a, mass) + c]["T"]/table[(b, mass) + c]["T"]))
                              for c in cells if lo <= c[1] < hi])
                out.append(f"{np.median(d):.2f} / {np.mean(d < math.log10(2)):.0%}")
            lines.append(f"| {label} | " + " | ".join(out) + " |")
        lines.append("")
    return "\n".join(lines)


def plot_courant(table, out_dir):
    masses = sorted({k[1] for k in table})
    fig, axes = plt.subplots(1, len(masses), figsize=(12, 4.6), sharey=True, squeeze=False)
    for ax, mass in zip(axes[0], masses):
        for m in COURANT_METHODS:
            pts = sorted((r["dt_yr"], r["cpu_us_median"]) for k, r in table.items() if k[0] == m and k[1] == mass)
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, "o", ms=3, alpha=0.6, color=STYLE[m]["color"], label=STYLE[m]["label"])
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Courant time step [yr]")
        ax.set_title(f"particle mass {mass:g} Msun")
    axes[0, 0].set_ylabel(r"CPU time per update [$\mu$s]")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "courant_timing.png"), dpi=130)
    plt.close(fig)


def generate_report(data_dir, out_dir):
    """Write figures and summary.md from saved CSVs, and return the summary text."""
    os.makedirs(out_dir, exist_ok=True)
    report = [
        "With the current benchmark settings, EQCLib applies an Asplund et al. (2009) CIE supplement "
        "for N, Ne, Mg, S, Ca and Fe. CHIMES matched evolves H, He, C, O and Si. CHIMES full adds "
        "the six elements at Wiersma et al. (2009) abundances.", "",
    ]
    with plt.rc_context(PLOT_STYLE):
        eq = load_equilibrium(data_dir)
        if "chimes_full" in eq:
            plot_profiles(eq, out_dir)
            plot_ratio_maps(eq, out_dir)
            plot_thermal_equilibrium(eq, out_dir)
            cpu = {n: np.median([r["cpu_us"] for r in t.values() if r["status"] == 0]) for n, t in eq.items()}
            conv = {n: f"{sum(r['status'] == 0 for r in t.values())}/{len(t)}" for n, t in eq.items()}
            report += ["## Equilibrium", "", "Converged: " + ", ".join(f"{n} {c}" for n, c in conv.items()),
                       "", "Median CPU per fixed-T equilibrium [us]: " + ", ".join(f"{n} {v:.1f}" for n, v in cpu.items()),
                       "", equilibrium_stats(eq), ""]
        rows = load_timing(data_dir)
        if rows:
            plot_timing(rows, out_dir)
            plot_timing_accuracy(rows, out_dir)
            report += ["## Timing", "", timing_summary(rows), "", speed_ratios(rows), ""]
        courant = load_courant(data_dir)
        if courant:
            plot_courant(courant, out_dir)
            report += ["## Courant-step updates", "", courant_summary(courant)]
        stats = plot_cooling(data_dir, out_dir)
        if stats:
            report += ["## Cooling histories", "", "| T0 | nH | method | total CPU ms | steps | failures |",
                       "| ---: | ---: | --- | ---: | ---: | ---: |"]
            report += [f"| {T0:.0e} | {nH:g} | {m} | {cpu/1e3:.1f} | {n} | {f} |" for T0, nH, m, cpu, n, f in stats]
    text = "\n".join(report)
    with open(os.path.join(out_dir, "summary.md"), "w") as f:
        f.write(text + "\n")
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_DIR)
    args = parser.parse_args()
    plt.switch_backend("Agg")
    print(generate_report(args.data_dir, args.out_dir))


if __name__ == "__main__":
    main()
