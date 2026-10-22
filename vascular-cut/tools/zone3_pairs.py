"""Where zone 3 begins: placed by the station curve, then on the subclavian's mouth.

    .venv\\Scripts\\python tools\\zone3_pairs.py
    .venv\\Scripts\\python tools\\zone3_pairs.py --only KK32_0,KK20_2

Framed on the left subclavian, whose mouth is outlined in black, because that is the
whole question: zone 3 is meant to begin immediately after that vessel, and the
complaint was that a band of zone 2 sits between them.

Top row takes the subclavian's landmark station as given. Bottom row checks it against
the vessel's own mouth and, where the landmark is more than 2 mm past it, slides zones 2,
3 and 4 back by the difference. Same camera, same colours.

Only eight of thirty move at all -- the landmark and the mouth agree to a median of
0.1 mm across the cohort. These are the eight.
"""
import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import imageio.v3 as iio

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_map import ZONE_COLOURS, UNCLAIMED, UP, hex_rgb    # noqa: E402
import zone_stub_bench as ZB                                 # noqa: E402
from zone_stub_bench import label_case, ostium_zones         # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
VIEWS = [np.array([0.4, -1.0, 0.15]), np.array([-1.0, -0.2, 0.15])]
ZOOM_OUT = 26.0           # subclavian radii of arch to keep around it


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--tag", default="GT")
    ap.add_argument("--only",
                    default="KK32_0,KK20_2,KK16_1,KK24_1,KP20_1,KR47_1,KK5_1,KP35_1")
    ap.add_argument("--clear", action="store_true",
                    help="drop any earlier pair images first")
    args = ap.parse_args()

    outdir = ROOT / "tools" / f"zone_map_{args.tag}" / "pairs"
    if args.clear and outdir.is_dir():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    want = set(args.only.split(","))

    for cid, cl, mp in find_cases(args.cohort):
        if cid not in want:
            continue
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            pl = pv.Plotter(shape=(2, 2), off_screen=True, border=False,
                            window_size=(940, 900))
            focus = dist = None
            for r, fix_on in enumerate((False, True)):
                ZB.LANDMARK_LAG_FIX = fix_on
                mesh, route, zones, lab, isb, adj, _f = label_case(s)
                reg = RS.stub_regions(mesh, route, adj=adj)
                oz = ostium_zones(s, mesh, route, reg) if reg.max() >= 0 else {}
                if reg.max() >= 0:
                    lab = RS.own_stubs(mesh, lab, reg, isb, adj=adj, region_zone=oz)
                V = np.asarray(mesh.vertices, float)
                cols = np.tile(hex_rgb(UNCLAIMED), (len(V), 1))
                for z in range(int(lab.max()) + 1):
                    cols[lab == z] = hex_rgb(ZONE_COLOURS[z % len(ZONE_COLOURS)])
                lsa = [j for j, z in oz.items() if z == 2]
                e = np.asarray(mesh.edges_unique, int)
                if lsa:
                    sel = np.isin(reg, lsa)
                    cut = sel[e[:, 0]] != sel[e[:, 1]]
                    rim = np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))
                    cols[rim] = np.array([15, 15, 15], float)
                    if focus is None:
                        focus = V[sel].mean(axis=0)
                        dist = ZOOM_OUT * float(
                            np.linalg.norm(V[sel] - focus, axis=1).max())
                if focus is None:
                    focus, dist = V.mean(axis=0), 220.0
                F = np.asarray(mesh.faces)
                poly = pv.PolyData(V, np.hstack([np.full((len(F), 1), 3), F]).ravel())
                poly.point_data["rgb"] = cols.astype(np.uint8)
                for c, view in enumerate(VIEWS):
                    pl.subplot(r, c)
                    pl.add_mesh(poly.copy(), scalars="rgb", rgb=True,
                                smooth_shading=True, specular=0.12, diffuse=0.92,
                                ambient=0.3)
                    v = view / np.linalg.norm(view)
                    pl.camera.position = focus + v * dist
                    pl.camera.focal_point = focus
                    pl.camera.up = UP
            ZB.LANDMARK_LAG_FIX = True
            img = np.asarray(pl.screenshot(None, return_img=True))
            pl.close()
            iio.imwrite(str(outdir / f"{cid}.png"), img)
            print(f"  {cid}: wrote {outdir / (cid + '.png')}", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid}: ERROR {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()
