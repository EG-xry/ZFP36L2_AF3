"""
Ray-traced view of the WT AlphaFold 3 complex that seeded the MD runs.

Run from the repo root with PyMOL's bundled Python:

    /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq figures/src/render_structure.py

Writes figures/data/structure_render.png (no labels) and
figures/data/structure_anchors.csv (pixel position of every residue that
make_figures.py labels), so labels stay editable text in the final SVG.
Pass "-- preview" after the script name for quick low-resolution renders of
candidate orientations.
"""
import sys
from pathlib import Path

import numpy as np
from pymol import cmd

ROOT = Path.cwd()
PDB = ROOT / "anvil_md" / "best_models" / "WT_best_model.pdb"
OUT_PNG = ROOT / "figures" / "data" / "structure_render.png"
OUT_ANCHORS = ROOT / "figures" / "data" / "structure_anchors.csv"
WIDTH, HEIGHT = 2400, 1500

# Palette (figures/src/style.py)
COLORS = {
    "del": "0xeb6834", "neu": "0x2a78d6", "zn_ctrl": "0x898781", "t156": "0x008300",
    "protein": "0xd8d6cf", "finger1": "0xc9d3df", "finger2": "0xdcd3c4",
    "rna": "0x3b3a37", "zn": "0xdfe2e6",
}
# ZFP36L2 numbering (construct residue + 150). D210E and C208A are drawn at the
# residues their AF3 constructs actually mutate (D219, C212).
SITES = [
    ("Y154H", 154, "del"), ("R160W", 160, "del"), ("R184C", 184, "del"),
    ("G204D", 204, "del"), ("C206F", 206, "del"),
    ("S185N", 185, "neu"), ("T187S", 187, "neu"), ("D210E", 219, "neu"),
    ("C208A", 212, "zn_ctrl"), ("T156", 156, "t156"),
]
ZN_LIGANDS = [159, 168, 174, 178, 197, 206, 212, 216]


def style() -> None:
    cmd.load(str(PDB), "cx")
    cmd.alter("cx and chain A", "resv = resv + 150")
    cmd.sort()
    cmd.remove("hydro")
    for name, hexcode in COLORS.items():
        cmd.set_color(name, [int(hexcode[i:i + 2], 16) / 255 for i in (2, 4, 6)])

    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("cartoon_loop_radius", 0.34)
    cmd.show("cartoon", "chain A")
    cmd.color("protein", "chain A")
    cmd.color("finger1", "chain A and resi 159-178")
    cmd.color("finger2", "chain A and resi 197-216")

    cmd.show("cartoon", "chain D")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.show("sticks", "chain D and not name P+OP1+OP2+O5'")
    cmd.color("rna", "chain D")
    cmd.set("stick_radius", 0.18, "chain D")

    cmd.show("spheres", "chain B+C")
    cmd.color("zn", "chain B+C")
    cmd.set("sphere_scale", 0.62, "chain B+C")
    lig = "chain A and resi " + "+".join(map(str, ZN_LIGANDS)) + " and not name N+C+O"
    cmd.show("sticks", lig)
    cmd.color("grey50", lig)
    cmd.set("stick_radius", 0.14, lig)

    for _, resi, cls in SITES:
        sel = f"chain A and resi {resi} and name CA"
        cmd.show("spheres", sel)
        cmd.set("sphere_scale", 0.95, sel)
        cmd.color(cls if cls != "t156" else "t156", sel)
    cmd.show("sticks", "chain A and resi 156 and not name N+C+O")
    cmd.color("t156", "chain A and resi 156 and not name N+C+O")
    cmd.set("stick_radius", 0.22, "chain A and resi 156")

    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_color", "0x52514e")
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.42)
    cmd.set("direct", 0.55)
    cmd.set("specular", 0.15)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)


def project(selection: str) -> tuple[float, float]:
    """Pixel (x, y) of a selection's centroid in the current orthoscopic view."""
    view = cmd.get_view()
    rot = np.array(view[0:9]).reshape(3, 3)  # rows = model axes in camera space
    cam_origin = np.array(view[9:12])
    model_origin = np.array(view[12:15])
    p = np.array(cmd.get_coords(selection)).mean(axis=0)
    cam = (p - model_origin) @ rot + cam_origin
    half_h = -cam_origin[2] * np.tan(np.radians(cmd.get_setting_float("field_of_view")) / 2)
    half_w = half_h * WIDTH / HEIGHT
    return ((cam[0] / half_w + 1) / 2 * WIDTH, (1 - cam[1] / half_h) / 2 * HEIGHT)


def main() -> None:
    style()
    cmd.orient("cx")
    preview = "preview" in sys.argv
    if preview:
        for i, (axis, angle) in enumerate([("x", 0), ("x", 90), ("x", 180), ("x", 270),
                                           ("y", 180)]):
            cmd.orient("cx")
            cmd.turn(axis, angle)
            cmd.png(str(OUT_PNG.with_name(f"preview_{i}_{axis}{angle}.png")), width=900,
                    height=560, ray=1)
        return
    cmd.turn("x", float(sys.argv[sys.argv.index("--turn-x") + 1]) if "--turn-x" in sys.argv else 0)
    cmd.zoom("cx", 3.0)
    cmd.png(str(OUT_PNG), width=WIDTH, height=HEIGHT, dpi=300, ray=1)
    rows = ["label,resi,cls,x_px,y_px"]
    for label, resi, cls in SITES:
        x, y = project(f"chain A and resi {resi} and name CA")
        rows.append(f"{label},{resi},{cls},{x:.1f},{y:.1f}")
    for label, sel in [("5'", "chain D and resi 1 and name C5'"),
                       ("3'", "chain D and resi 9 and name C3'"),
                       ("Zn1", "chain B"), ("Zn2", "chain C")]:
        x, y = project(sel)
        rows.append(f"{label},,{'rna' if 'Zn' not in label else 'zn'},{x:.1f},{y:.1f}")
    OUT_ANCHORS.write_text("\n".join(rows) + "\n")
    print(f"wrote {OUT_PNG} and {OUT_ANCHORS}")


main()
