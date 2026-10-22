"""A cohort of cut surfaces as one figure, with the OPENINGS made visible.

For the tool paper. A montage of grey aortas demonstrates nothing -- the tool's output
IS the openings, planar cuts where the closed segmentation had none, so the figure has
to show those or it is a gallery. Each cut opening is filled with a contrasting disc at
the plane recorded in that case's sidecar; the wall stays neutral.

    .venv\\Scripts\\python tools\\cohort_montage.py nubs
    .venv\\Scripts\\python tools\\cohort_montage.py arch --mark-bovine
    .venv\\Scripts\\python tools\\cohort_montage.py nubs --cols 6 --order area

Every panel uses the SAME view direction, so the eye compares anatomy rather than
camera. Panels are scaled individually to fill their cell -- an aorta cohort spans
16,000 to 87,000 mm2 and a common scale would leave the small ones unreadable -- so the
area is printed on each panel instead, and the caption must say the scaling is per
panel.
"""
import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

VIEW = np.array([0.4, -1.0, 0.15])       # nub_viewer.DEFAULT_VIEW, so figures match
UP = np.array([0.0, 0.0, 1.0])
WALL = "#d9d4cf"
OPENING = "#c2410c"                      # the cuts: the thing being demonstrated
BOVINE = "#1d4ed8"

SETS = {
    "nubs": dict(folder="output/nubs_final", suffix="_nub",
                 title="Aortas cut with the nub-in preset: every branch reduced to a "
                       "short stub, aortic ends squared"),
    "arch": dict(folder="output/arch_final", suffix="_arch",
                 title="Aortic arches isolated with the arch preset, great vessels "
                       "left as nubs"),
}
# Detected structurally by the tool -- brachiocephalic feeding three vessels -- not
# labelled by hand. tools/NUBS.md records the list and the test.
BOVINE_CASES = {"KK24_1", "KK5_1", "KR27_2", "KY160_1", "KY215_1", "KY238_1",
                "KY244_1", "KY68_1"}


def load_case(stl, sidecar):
    m = trimesh.load(str(stl), process=False)
    m.merge_vertices()
    poly = pv.PolyData(np.asarray(m.vertices),
                       np.hstack([np.full((len(m.faces), 1), 3),
                                  np.asarray(m.faces)]).ravel())
    meta = {}
    if os.path.exists(sidecar):
        meta = json.load(open(sidecar, encoding="utf-8"))
    return poly, meta


def opening_discs(meta):
    """A disc at every recorded cut plane, sized from that opening's diameter.

    Drawn from the SIDECAR rather than found in the mesh, deliberately: the sidecar is
    what the tool claims it did, so a disc that does not sit in a hole is a visible
    failure of the record rather than a hidden one.
    """
    out = []
    for o in meta.get("openings", []):
        c, n = o.get("plane_origin"), o.get("plane_normal")
        if not c or not n:
            continue
        d = o.get("diameter_mm") or 6.0
        out.append(pv.Disc(center=np.asarray(c, float), inner=0.0,
                           outer=float(d) * 0.62, normal=np.asarray(n, float),
                           c_res=48))
    return out


def cut_region_scalars(poly, meta, band_mm=4.0):
    """1.0 on a collar around each cut rim, 0.0 on the rest of the wall.

    A 3 mm branch opening on a 300 mm aorta is about four pixels in a montage panel, so
    a disc at the rim is invisible however carefully it is drawn, and enlarging it would
    misstate the size of the thing being demonstrated. A collar of the surface beside
    each cut is visible and is honestly the right size, because it is the surface the
    cut passed through.

    Measured as distance to the cut PLANE, bounded laterally by the opening's own
    radius -- not distance to the opening's centre, which cannot work at both scales at
    once: a 35 mm aortic end has its rim 17 mm from its centre, so any threshold tight
    enough to be a collar there catches nothing, and any threshold loose enough catches
    half the vessel. That was the first attempt and it turned 70% of every aorta orange.
    """
    v = np.asarray(poly.points, float)
    tag = np.zeros(len(v), dtype=float)
    for o in meta.get("openings", []):
        c, nrm = o.get("plane_origin"), o.get("plane_normal")
        if not c or not nrm:
            continue
        c = np.asarray(c, float)
        nrm = np.asarray(nrm, float)
        nn = np.linalg.norm(nrm)
        if nn < 1e-9:
            continue
        nrm = nrm / nn
        d = float(o.get("diameter_mm") or 6.0)
        rel = v - c
        along = rel @ nrm                                   # signed, along the normal
        lateral = np.linalg.norm(rel - np.outer(along, nrm), axis=1)
        tag[(np.abs(along) < band_mm) & (lateral < 0.85 * d)] = 1.0
    return tag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=sorted(SETS))
    ap.add_argument("--cols", type=int, default=6)
    ap.add_argument("--order", choices=("area", "case"), default="area")
    ap.add_argument("--mark-bovine", action="store_true")
    ap.add_argument("--cell", type=int, default=460, help="pixels per panel")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    spec = SETS[args.which]
    folder = ROOT / spec["folder"]
    stls = sorted(folder.glob(f"*{spec['suffix']}.stl"))
    if not stls:
        raise SystemExit(f"no meshes in {folder}")

    cases = []
    for stl in stls:
        cid = stl.name.replace(f"{spec['suffix']}.stl", "")
        side = str(stl).replace(".stl", "_cut_planes.json")
        poly, meta = load_case(stl, side)
        cases.append((cid, poly, meta, float(meta.get("area_mm2") or 0.0)))
    if args.order == "area":
        cases.sort(key=lambda t: t[3])

    n = len(cases)
    cols = args.cols
    rows = (n + cols - 1) // cols
    print(f"{n} cases, {rows}x{cols}")

    pl = pv.Plotter(shape=(rows, cols), off_screen=True, border=False,
                    window_size=(cols * args.cell, rows * args.cell))
    pl.set_background("white")
    for i, (cid, poly, meta, area) in enumerate(cases):
        pl.subplot(i // cols, i % cols)
        bov = args.mark_bovine and cid in BOVINE_CASES
        accent = BOVINE if bov else OPENING
        # The wall in one colour and everything the cuts produced in another, as a
        # two-value scalar rather than two actors, so the shading stays continuous
        # across the join instead of showing a seam that is not in the surface.
        poly = poly.copy()
        poly["cut"] = cut_region_scalars(poly, meta)
        pl.add_mesh(poly, scalars="cut", cmap=[WALL, accent], show_scalar_bar=False,
                    smooth_shading=True, specular=0.15, diffuse=0.9, ambient=0.22)
        for d in opening_discs(meta):
            pl.add_mesh(d, color=accent, ambient=0.45)
        n_open = len(meta.get("openings", []))
        label = f"{cid}" + ("   bovine" if bov else "")
        pl.add_text(label, position="upper_left", font_size=9,
                    color="#1d4ed8" if bov else "#222222")
        pl.add_text(f"{area:,.0f} mm²   {n_open} openings",
                    position="lower_left", font_size=6, color="#777777")
        pl.view_vector(VIEW, viewup=UP)
        pl.reset_camera()
        # No zoom past the reset. 1.25 overflowed the cell on the taller aortas and
        # ran their labels into the row below.
        pl.camera.zoom(0.92)
    for j in range(n, rows * cols):                     # blank the unused cells
        pl.subplot(j // cols, j % cols)
        pl.set_background("white")

    out = args.out or str(ROOT / "docs" / "figures" / f"montage_{args.which}.png")
    pl.screenshot(out, scale=1)
    pl.close()
    size = os.path.getsize(out) / 1024
    print(f"wrote {out}  ({size:.0f} KB)")
    if args.mark_bovine:
        found = {c[0] for c in cases} & BOVINE_CASES
        print(f"bovine marked: {len(found)} of {len(BOVINE_CASES)} -> {sorted(found)}")
        missing = BOVINE_CASES - {c[0] for c in cases}
        if missing:
            print(f"  NOT IN THIS SET: {sorted(missing)}")


if __name__ == "__main__":
    main()
