"""Matplotlib style shared by the EQCLib plotting scripts and notebooks."""

BLACK, GRAY = "#202020", "#555555"
BLUE, ORANGE, GREEN, PURPLE = "#0072b2", "#c44e00", "#007b5e", "#8d4b9f"
SKY, GOLD, PINK = "#56b4e9", "#b8860b", "#cc79a7"

# Colour order for comparisons; the first entry marks the reference model.
SERIES_COLORS = (BLACK, BLUE, ORANGE, GREEN, PURPLE, SKY, GOLD, PINK)

PLOT_STYLE = {
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.alpha": 0.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "font.size": 10,
}
