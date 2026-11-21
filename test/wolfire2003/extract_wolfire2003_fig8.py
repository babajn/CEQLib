"""Extract the standard-model curve from Figure 8 (arXiv:astro-ph/0207098v2).

The axes and solid polyline are vector paths on page 86 of the local PDF.
Source https://arxiv.org/pdf/astro-ph/0207098v2#page=86
Regenerating the CSV requires Poppler's pdftocairo and Python 3.
"""

import argparse
import csv
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


def figure_points(pdf_path):
    svg = subprocess.check_output(
        ["pdftocairo", "-f", "86", "-l", "86", "-svg", str(pdf_path), "-"]
    )
    paths = []
    for path in ET.fromstring(svg).iter("{http://www.w3.org/2000/svg}path"):
        if "stroke:none" in path.get("style", "") or not path.get("transform"):
            continue
        segments = [
            (kind, float(x), float(y))
            for kind, x, y in re.findall(r"([ML])\s+([-\d.]+)\s+([-\d.]+)", path.get("d", ""))
        ]
        paths.append((path, segments))

    horizontal_segments = [
        (x0, x1, y0)
        for _, segments in paths
        for (_, x0, y0), (kind, x1, y1) in zip(segments, segments[1:])
        if kind == "L" and y0 == y1 and x0 < x1
    ]
    left, right, _ = max(horizontal_segments, key=lambda edge: edge[1] - edge[0])
    bottom, top = sorted(
        {y for x0, x1, y in horizontal_segments if (x0, x1) == (left, right)}
    )
    [curve] = [
        segments
        for path, segments in paths
        if segments
        and "stroke-dasharray" not in path.get("style", "")
        and (segments[0][1], segments[-1][1]) == (left, right)
        and all(
            kind == "L" and x0 < x1
            for (_, x0, _), (kind, x1, _) in zip(segments, segments[1:])
        )
    ]

    # The horizontal axis is 0.01 to 1000 cm^-3 and the vertical axis
    # is 1 to 100000 K cm^-3. Both span five decades.
    for _, x, y in curve:
        log_density = -2 + 5 * (x - left) / (right - left)
        log_pressure = 5 * (y - bottom) / (top - bottom)
        yield 10**log_density, 10**log_pressure


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pdf",
        nargs="?",
        type=Path,
        default=here.parents[1] / "reference" / "Wolfire_2003.pdf",
        help="arXiv:astro-ph/0207098v2 PDF",
    )
    points = list(figure_points(parser.parse_args().pdf))
    with (here / "wolfire2003_fig8_standard.csv").open("w", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["nH_cm3", "P_over_k_K_cm3"])
        writer.writerows((f"{density:.9g}", f"{pressure:.9g}") for density, pressure in points)
