"""
Shared look for every write-up figure: palette, typography, axes, export.

Colors come from the validated reference palette (dataviz skill, palette.md).
Two identity channels are used, and they never share a chart:
  - variant class:  neutral blue / deleterious orange (+ WT ink, Zn control muted)
  - method:         AF3 violet / MD green
Both pairs pass the all-pairs CVD and normal-vision checks on a white card.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parents[1]

# ── Ink & chrome ──────────────────────────────────────────────────────────
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
WASH = "#f5f4f0"

# ── Identity: variant class ───────────────────────────────────────────────
NEU = "#2a78d6"
DEL = "#eb6834"
WT = INK2
ZN = MUTED
CLASS_COLOR = {"deleterious": DEL, "neutral": NEU, "WT": WT, "zn_ctrl": ZN}
CLASS_MARKER = {"deleterious": "o", "neutral": "o", "WT": "s", "zn_ctrl": "D"}
CLASS_LABEL = {"deleterious": "Deleterious (EMSA)", "neutral": "Neutral (EMSA)",
               "WT": "Wild type", "zn_ctrl": "Zn-knockout control"}

# ── Identity: method ──────────────────────────────────────────────────────
AF3 = "#6250d6"
MD = "#008300"

# ── Sequential ramp (blue, palette.md steps 100-700) ──────────────────────
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]

# ── Variants ──────────────────────────────────────────────────────────────
# Labels follow Akcesme et al. 2025 (RNA Biology 22:1). C208A is the Zn-knockout control; in the
# construct numbering used here it mutates C212 (see docs/02_MD_results.md).
DELETERIOUS = ["Y154H", "R160W", "R184C", "G204D", "C206F"]
NEUTRAL = ["S185N", "T187S", "D210E"]
VARIANT_CLASS = {"WT": "WT", "C208A": "zn_ctrl",
                 **{v: "deleterious" for v in DELETERIOUS},
                 **{v: "neutral" for v in NEUTRAL}}
DISPLAY_ORDER = ["WT", *DELETERIOUS, *NEUTRAL, "C208A"]


def apply() -> None:
    """Install the shared rcParams."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        # Arial ships a real bold face; matplotlib only sees regular Helvetica (.ttc).
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size": 9.5,
        "text.color": INK,
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK2,
        "axes.labelsize": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.titlecolor": INK,
        "axes.titlelocation": "left",
        "axes.titlepad": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 0.75,
        "grid.linestyle": "-",
        "xtick.color": AXIS,
        "ytick.color": AXIS,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.labelcolor": INK2,
        "lines.linewidth": 1.5,
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "svg.fonttype": "none",
        "svg.hashsalt": "zfp36l2",
        "figure.dpi": 110,
    })


def hgrid(ax) -> None:
    """Recessive horizontal hairlines behind the data."""
    ax.yaxis.grid(True)
    ax.set_axisbelow(True)


def vgrid(ax) -> None:
    ax.xaxis.grid(True)
    ax.set_axisbelow(True)


def top(fig, inches: float) -> float:
    """Figure-fraction y that sits `inches` below the top edge.

    Layout is measured in inches from the top so the header block looks the
    same on short and tall figures.
    """
    return 1 - inches / fig.get_figheight()


def header(fig, title: str, subtitle: str | None = None) -> None:
    """Left-aligned title + takeaway; occupies the top ~0.85 in."""
    fig.text(0.012, top(fig, 0.12), title, ha="left", va="top", fontsize=13.5,
             fontweight="bold", color=INK)
    if subtitle:
        fig.text(0.012, top(fig, 0.45), subtitle, ha="left", va="top",
                 fontsize=10, color=INK2, linespacing=1.4)


def save(fig, name: str, dpi: int = 220) -> None:
    """Write figures/<name>.svg (editable text) and figures/<name>.png."""
    svg_path = FIG_DIR / f"{name}.svg"
    png_path = FIG_DIR / f"{name}.png"
    fig.savefig(svg_path, dpi=dpi, metadata={"Date": None})  # dpi sets embedded image resolution
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    print(f"  wrote {svg_path.name} + {png_path.name}")
