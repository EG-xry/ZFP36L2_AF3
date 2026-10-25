"""
Shared look for every write-up figure: palette, typography, axes, export.

Colors come from the validated reference palette (dataviz skill, palette.md).
Two identity channels are used, and they never share a chart:
  - variant class:  neutral blue / deleterious orange (+ WT ink, Zn control muted)
  - method:         AF3 violet / MD green
Both pairs pass the all-pairs CVD and normal-vision checks on a white card.
"""
from __future__ import annotations

import re
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

FIG_DIR = Path(__file__).resolve().parents[1]

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_HREF = "{http://www.w3.org/1999/xlink}href"
for _prefix, _uri in [("", SVG_NS), ("xlink", "http://www.w3.org/1999/xlink"),
                      ("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#"),
                      ("cc", "http://creativecommons.org/ns#"),
                      ("dc", "http://purl.org/dc/elements/1.1/")]:
    ET.register_namespace(_prefix, _uri)
_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

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


FONT = "Source Sans 3"  # figures/fonts, SIL Open Font License


def apply() -> None:
    """Install the shared rcParams."""
    for ttf in sorted((FIG_DIR / "fonts").glob("*.ttf")):
        font_manager.fontManager.addfont(str(ttf))
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT, "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,  # plain hyphen for negative numbers
        "font.size": 9.5,
        "text.color": INK,
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK2,
        "axes.labelsize": 10,
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
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.labelcolor": INK2,
        "lines.linewidth": 1.5,
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT,
        "mathtext.it": f"{FONT}:italic",
        "mathtext.bf": f"{FONT}:bold",
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


def header(fig, title: str) -> None:
    """Short left-aligned title; the caption in the write-up carries the explanation."""
    fig.text(0.012, top(fig, 0.12), title, ha="left", va="top", fontsize=13,
             fontweight="bold", color=INK)


def _tag(el) -> str:
    return el.tag.rsplit("}", 1)[-1]


def _style(el) -> dict:
    pairs = (p.split(":", 1) for p in (el.get("style") or "").split(";") if ":" in p)
    return {k.strip(): v.strip() for k, v in pairs}


def _bbox(el) -> tuple[float, float, float, float] | None:
    """Bounding box of a shape or group in its own user space; None if not computable."""
    if el.get("transform") and _tag(el) != "path":
        return None
    boxes = []
    for e in el.iter():
        tag = _tag(e)
        if tag == "g":
            if e.get("transform"):
                return None
            continue
        if tag == "path":
            nums = [float(n) for n in _NUM.findall(e.get("d", ""))]
            if not nums or len(nums) % 2:
                return None
            xs, ys = nums[0::2], nums[1::2]
        elif tag == "rect":
            x, y = float(e.get("x", 0)), float(e.get("y", 0))
            xs, ys = [x, x + float(e.get("width"))], [y, y + float(e.get("height"))]
        else:
            return None
        tx = ty = 0.0
        if e.get("transform"):
            m = re.fullmatch(r"translate\(([^ ,]+)[ ,]+([^ )]+)\)", e.get("transform").strip())
            if not m or e is el:
                return None
            tx, ty = float(m[1]), float(m[2])
        st = _style(e)  # strokes reach half their width past the path (round joins and caps)
        pad = float(st.get("stroke-width", "1").rstrip("px")) / 2 if st.get("stroke", "none") != "none" else 0
        boxes.append((min(xs) + tx - pad, min(ys) + ty - pad, max(xs) + tx + pad, max(ys) + ty + pad))
    if not boxes:
        return None
    return (min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes))


def make_editable(svg_path: Path) -> None:
    """Rewrite a matplotlib SVG so each shape is its own element in a vector editor.

    Inlines every <use> clone (markers, ticks) as a real <path>, drops clipping masks
    from shapes that lie inside them anyway (images keep theirs), and unwraps wrapper
    groups that hold a single element or have no name.
    """
    tree = ET.parse(svg_path)
    root = tree.getroot()
    parents = {c: p for p in root.iter() for c in p}
    defs = {el.get("id"): el for d in root.iter(f"{{{SVG_NS}}}defs") for el in d if el.get("id")}
    for use in list(root.iter(f"{{{SVG_NS}}}use")):
        ref = defs[(use.get(XLINK_HREF) or use.get("href"))[1:]]
        new = ET.Element(ref.tag, {k: v for k, v in ref.attrib.items() if k not in ("id", "style")})
        style = {**_style(ref), **_style(use)}
        if style:
            new.set("style", "; ".join(f"{k}: {v}" for k, v in style.items()))
        x, y = float(use.get("x", 0)), float(use.get("y", 0))
        if x or y:
            new.set("transform", f"translate({x:g} {y:g})")
        if use.get("clip-path"):
            new.set("clip-path", use.get("clip-path"))
        new.tail = use.tail
        parent = parents[use]
        parent[list(parent).index(use)] = new
        parents[new] = parent

    clips = {}
    for cp in root.iter(f"{{{SVG_NS}}}clipPath"):
        r = cp.find(f"{{{SVG_NS}}}rect")
        if r is not None:
            x, y = float(r.get("x")), float(r.get("y"))
            clips[cp.get("id")] = (x, y, x + float(r.get("width")), y + float(r.get("height")))
    for el in root.iter():
        box = clips.get((el.get("clip-path") or "")[5:-1])
        bb = _bbox(el) if box else None
        tol = 1.0  # pt; an overhang this small cannot be seen
        if bb and bb[0] >= box[0] - tol and bb[1] >= box[1] - tol \
                and bb[2] <= box[2] + tol and bb[3] <= box[3] + tol:
            del el.attrib["clip-path"]
        if _tag(el) == "text" and re.fullmatch(r"rotate\(-?0 \S+ \S+\)", el.get("transform", "")):
            del el.attrib["transform"]

    used = set()
    for el in root.iter():
        for v in el.attrib.values():
            used.update(re.findall(r"url\(#([^)]+)\)", v))
        if (el.get(XLINK_HREF) or "").startswith("#"):
            used.add(el.get(XLINK_HREF)[1:])
    for d in list(root.iter(f"{{{SVG_NS}}}defs")):
        for el in list(d):
            if el.get("id") and el.get("id") not in used:
                d.remove(el)

    def unwrap(el) -> None:
        out = []
        for child in list(el):
            unwrap(child)
            plain = _tag(child) == "g" and not set(child.attrib) - {"id"}
            kids = list(child)
            if _tag(child) == "defs" and not kids:
                continue
            if plain and len(kids) == 1:
                if child.get("id") and not kids[0].get("id"):
                    kids[0].set("id", child.get("id"))
                kids[0].tail = child.tail
                out.append(kids[0])
            elif plain and not child.get("id"):
                out.extend(kids)
            else:
                out.append(child)
        for child in list(el):
            el.remove(child)
        el.extend(out)

    unwrap(root)
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)


def save(fig, name: str, dpi: int = 220) -> None:
    """Write figures/<name>.svg (every element editable) and figures/<name>.png."""
    svg_path = FIG_DIR / f"{name}.svg"
    png_path = FIG_DIR / f"{name}.png"
    fig.savefig(svg_path, dpi=dpi, metadata={"Date": None})  # dpi sets embedded image resolution
    make_editable(svg_path)
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    print(f"  wrote {svg_path.name} + {png_path.name}")
