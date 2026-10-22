"""Close-ups of the arch branches, labelled before and after stub ownership.

    .venv\\Scripts\\python tools\\zone_stub_pairs.py --only DJ18_2,KK16_1,KP35_1,KR52_1

The whole-aorta maps cannot show this defect: a streak down one flank of a 6 mm vessel
is a few pixels at that zoom. These pair the same view with the stub rule off and on,
framed on the branches, so the change is the only thing in the picture.

Top row is before, bottom row after, and the camera is identical between them.
"""
import argparse
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
from zone_map import ZONE_COLOURS, UNCLAIMED, VIEWS, UP, hex_rgb   # noqa: E402
from zone_stub_bench import label_case, own_by_stub, ostium_zones  # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
ZOOM_OUT = 3.2            # how much aorta to keep around the branch cluster


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--tag", default="GT")
    ap.add_argument("--only", default="DJ18_2,KK16_1,KP35_1,KR52_1")
    args = ap.parse_args()

    outdir = ROOT / "tools" / f"zone_map_{args.tag}" / "pairs"
    outdir.mkdir(parents=True, exist_ok=True)
    want = set(args.only.split(","))
    for cid, cl, mp in find_cases(args.cohort):
        if cid not in want:
            continue
        s = CutSession(str(cl), str(mp))
        s.nub_preset()
        s.recut()
        mesh, route, _zones, lab, isb, adj, _flat = label_case(s)
        reg = RS.stub_regions(mesh, route, adj=adj)
        if reg.max() < 0:
            print(f"  {cid}: no stub region")
            continue
        new = own_by_stub(mesh, lab, reg, isb, adj,
                          region_zone=ostium_zones(s, mesh, route, reg))
        V = np.asarray(mesh.vertices, float)
        focus = V[reg >= 0].mean(axis=0)
        dist = ZOOM_OUT * float(np.linalg.norm(V[reg >= 0] - focus, axis=1).max())
        F = np.asarray(mesh.faces)
        poly = pv.PolyData(V, np.hstack([np.full((len(F), 1), 3), F]).ravel())
        pl = pv.Plotter(shape=(2, 2), off_screen=True, border=False,
                        window_size=(920, 860))
        for r, L in enumerate((lab, new)):
            cols = np.tile(hex_rgb(UNCLAIMED), (len(V), 1))
            for zi in range(int(L.max()) + 1):
                cols[L == zi] = hex_rgb(ZONE_COLOURS[zi % len(ZONE_COLOURS)])
            for c, view in enumerate(VIEWS):
                pl.subplot(r, c)
                q = poly.copy()
                q.point_data["rgb"] = cols.astype(np.uint8)
                pl.add_mesh(q, scalars="rgb", rgb=True, smooth_shading=True,
                            specular=0.12, diffuse=0.92, ambient=0.3)
                v = view / np.linalg.norm(view)
                pl.camera.position = focus + v * dist
                pl.camera.focal_point = focus
                pl.camera.up = UP
        img = np.asarray(pl.screenshot(None, return_img=True))
        pl.close()
        iio.imwrite(str(outdir / f"{cid}.png"), img)
        print(f"  {cid}: {int((new != lab).sum())} label(s) repainted "
              f"-> {outdir / (cid + '.png')}", flush=True)


if __name__ == "__main__":
    main()
