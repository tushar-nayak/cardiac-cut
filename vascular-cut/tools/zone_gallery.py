"""Render every case's cut and every zone it offers, as one strip per case.

    .venv\\Scripts\\python tools\\zone_gallery.py                       # thoracic
    .venv\\Scripts\\python tools\\zone_gallery.py Y:/aorta/abdominal/MK/JK --tag JK

One PNG per case: the whole vessel as the tool cuts it, then each zone in turn, the
zone in orange against the rest of the vessel ghosted so it is obvious WHERE on the
aorta it sits. tools/zone_gallery_page.py turns the strips into a page to scroll.

Every panel of a case shares ONE camera, so the zone is the only thing that changes
between them. Panels are NOT shared across cases -- a 250 mm aorta and a 450 mm one
would leave the small one unreadable -- so each strip carries its own scale bar.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
from region_station import region_by_station                 # noqa: E402
import branch_cut as BC                                      # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
GHOST = "#e4dfd9"
ZONE = "#c2410c"
WALL = "#d9d4cf"
CUT = "#c2410c"
VIEW = np.array([0.4, -1.0, 0.15])
UP = np.array([0.0, 0.0, 1.0])


def poly(m):
    f = np.asarray(m.faces)
    return pv.PolyData(np.asarray(m.vertices),
                       np.hstack([np.full((len(f), 1), 3), f]).ravel())


def rim_scalars(mesh, meta_openings, band_mm=3.0):
    """Mark the wall beside each recorded cut plane, so the cuts are visible."""
    v = np.asarray(mesh.vertices, float)
    tag = np.zeros(len(v))
    for o in meta_openings:
        c, n = o.get("plane_origin"), o.get("plane_normal")
        if not c or not n:
            continue
        c, n = np.asarray(c, float), np.asarray(n, float)
        nn = np.linalg.norm(n)
        if nn < 1e-9:
            continue
        n = n / nn
        d = float(o.get("diameter_mm") or 6.0)
        rel = v - c
        along = rel @ n
        lat = np.linalg.norm(rel - np.outer(along, n), axis=1)
        tag[(np.abs(along) < band_mm) & (lat < 0.85 * d)] = 1.0
    return tag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--tag", default="GT")
    ap.add_argument("--cell", type=int, default=300)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-nubs", action="store_true",
                    help="leave branches at their traced length")
    args = ap.parse_args()

    outdir = ROOT / "tools" / f"zone_gallery_{args.tag}"
    outdir.mkdir(exist_ok=True)
    cases = find_cases(args.cohort)
    if args.limit:
        cases = cases[:args.limit]
    print(f"{len(cases)} case(s) -> {outdir}")

    index = []
    t0 = time.time()
    for i, (cid, cl, mp) in enumerate(cases, 1):
        try:
            s = CutSession(str(cl), str(mp))
            if not args.no_nubs:
                # The cohort convention, and what makes a zone LOCAL: every branch
                # back to a stub 0.75 of its own radius past the aortic wall. Without
                # it a zone bounded at the left subclavian drags the whole traced
                # vessel along, which is neither local geometry nor comparable
                # between cases -- one scan traces 20 mm of subclavian and the next
                # 120 mm, and the zone's area follows the tracing rather than the
                # anatomy.
                s.nub_preset()
            st = s.recut()
            whole = s.mesh
            route = BC.aorta_route(s.tree)
            presets = s.presets()
            panels = [("as cut", whole, None)]
            for p in presets:
                m, info = region_by_station(whole, s.tree, p["s0"], p["s1"],
                                            route=route)
                panels.append((p["name"], m, info))

            n = len(panels)
            pl = pv.Plotter(shape=(1, n), off_screen=True, border=False,
                            window_size=(args.cell * n, int(args.cell * 1.28)))
            pl.set_background("white")
            # One camera for the whole strip: framed on the WHOLE vessel, so a zone
            # never appears to be a different size than the aorta it came from.
            for c, (name, m, info) in enumerate(panels):
                pl.subplot(0, c)
                if c == 0:
                    q = poly(m).copy()
                    q["cut"] = rim_scalars(m, s.sidecar().get("openings", []))
                    pl.add_mesh(q, scalars="cut", cmap=[WALL, CUT],
                                show_scalar_bar=False, smooth_shading=True,
                                specular=0.15, diffuse=0.9, ambient=0.25)
                else:
                    pl.add_mesh(poly(whole), color=GHOST, opacity=0.26,
                                smooth_shading=True, specular=0.0, ambient=0.62)
                    pl.add_mesh(poly(m), color=ZONE, smooth_shading=True,
                                specular=0.15, diffuse=0.9, ambient=0.3)
                pl.view_vector(VIEW, viewup=UP)
                pl.reset_camera()
                pl.camera.zoom(1.45)
            png = outdir / f"{cid}.png"
            # TRANSPARENT, not white. The page these strips go into renders in the
            # reader's own theme, and a white plate behind every vessel glares on a
            # dark ground and pins the design to one theme. On alpha the vessel sits
            # on whatever the page is.
            img = pl.screenshot(None, transparent_background=True, return_img=True)
            pl.close()
            # Trim the empty bands top and bottom. reset_camera frames on the whole
            # vessel and an aorta is tall and narrow, so a fixed cell leaves a third
            # of the strip empty -- which in a page of 30 strips is a lot of scrolling
            # past nothing. Columns are NOT trimmed: they carry the panel alignment.
            rows = np.nonzero((np.asarray(img)[..., 3] > 8).any(axis=1))[0]
            if len(rows):
                lo_r = max(int(rows.min()) - 6, 0)
                hi_r = min(int(rows.max()) + 7, np.asarray(img).shape[0])
                img = np.asarray(img)[lo_r:hi_r]
            import imageio.v3 as iio
            iio.imwrite(str(png), np.asarray(img))

            index.append({
                "case": cid,
                "png": png.name,
                "base_ok": bool(st["ok"]),
                "base_msg": "" if st["ok"] else str(st.get("message", ""))[:160],
                "route_mm": None if route is None else round(float(route.length), 1),
                "panels": [
                    {"name": nm,
                     "area_mm2": None if inf is None else round(inf["area_mm2"], 0),
                     "loops": None if inf is None else inf["n_loops"],
                     "faces": None if inf is None else inf["faces"],
                     "islands": None if inf is None
                     else inf["faces_dropped_as_islands"],
                     "s0": None if inf is None else inf["s0"],
                     "s1": None if inf is None else inf["s1"]}
                    for nm, _m, inf in panels],
            })
            print(f"  [{i:3d}/{len(cases)}] {cid:12s} {n - 1} zone(s)  "
                  f"{time.time() - t0:6.0f}s", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  [{i:3d}/{len(cases)}] {cid:12s} ERROR {type(e).__name__}: {e}")
            index.append({"case": cid, "error": f"{type(e).__name__}: {e}"[:200]})

    (outdir / "index.json").write_text(json.dumps(index, indent=1), encoding="utf-8")
    print(f"\nwrote {outdir}/index.json  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
