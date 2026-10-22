"""Your marks, against each candidate rule's boundary. The score, made visible.

    .venv\\Scripts\\python tools\\zone_score_shots.py
    .venv\\Scripts\\python tools\\zone_score_shots.py --rules current,nudge3

tools/zone_score.py reports one number per rule: mean millimetres of wall between the
marks and the boundary. That number decides which rule ships, so it has to be lookable
at -- a rule can post a better mean while putting the line somewhere obviously wrong,
and only the picture shows that.

Red spheres are the marks. The white ring is that rule's zone 2 / zone 3 boundary. The
black outline is the mouth of the vessel zone 2 owns. One row per rule, same camera.
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

import json                                                  # noqa: E402
from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
from zone_map import ZONE_COLOURS, UNCLAIMED, UP, hex_rgb    # noqa: E402
import zone_stub_bench as ZB                                 # noqa: E402
from zone_stub_bench import label_case, ostium_zones         # noqa: E402
from zone_score import RULES, boundary_distance              # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
MARKS = ROOT / "tools" / "zone_marks.json"
VIEWS = [np.array([0.4, -1.0, 0.15]), np.array([-1.0, -0.2, 0.15])]
ZOOM = 7.0                # marks-spread radii of arch to keep around them


def label_for(s, rule):
    """Labels under one rule, exactly as zone_score scores them."""
    nudge = float(RULES[rule].get("_nudge_mm", 0.0))
    for k, v in RULES[rule].items():
        if not k.startswith("_"):
            setattr(ZB, k, v)
    mesh, route, zones, lab, isb, adj, _f = label_case(s)
    reg = RS.stub_regions(mesh, route, adj=adj)
    oz = ostium_zones(s, mesh, route, reg) if reg.max() >= 0 else {}
    if reg.max() >= 0:
        lab = RS.own_stubs(mesh, lab, reg, isb, adj=adj, region_zone=oz)
    if nudge > 0:
        frozen = reg >= 0
        zs = sorted(set(lab[lab >= 0].tolist()))
        if 2 in zs and 3 in zs:
            lab, _ = RS.nudge_boundary(mesh, lab, 2, 3, nudge, adj, frozen)
            if 4 in zs:
                lab, _ = RS.nudge_boundary(mesh, lab, 3, 4, nudge, adj, frozen)
    return mesh, lab, reg, oz, adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--rules", default="current,nudge3")
    ap.add_argument("--tag", default="GT")
    args = ap.parse_args()
    if not MARKS.exists():
        raise SystemExit("no marks yet -- run tools/zone_mark.py first")
    marks = json.loads(MARKS.read_text())
    rules = [r for r in args.rules.split(",") if r in RULES]
    outdir = ROOT / "tools" / f"zone_marks_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    for cid, cl, mp in find_cases(args.cohort):
        if cid not in marks:
            continue
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            pts = np.array([m["point"] for m in marks[cid]], float)
            pl = pv.Plotter(shape=(len(rules), 2), off_screen=True, border=False,
                            window_size=(940, 430 * len(rules)))
            focus = pts.mean(axis=0)
            dist = None
            scores = []
            for r, rule in enumerate(rules):
                mesh, lab, reg, oz, adj = label_for(s, rule)
                V = np.asarray(mesh.vertices, float)
                e = np.asarray(mesh.edges_unique, int)
                cols = np.tile(hex_rgb(UNCLAIMED), (len(V), 1))
                for z in range(int(lab.max()) + 1):
                    cols[lab == z] = hex_rgb(ZONE_COLOURS[z % len(ZONE_COLOURS)])
                on = (((lab[e[:, 0]] == 2) & (lab[e[:, 1]] == 3))
                      | ((lab[e[:, 0]] == 3) & (lab[e[:, 1]] == 2)))
                if on.any():
                    cols[np.unique(e[on])] = np.array([255, 255, 255], float)
                lsa = [j for j, z in oz.items() if z == 2]
                if lsa:
                    sel = np.isin(reg, lsa)
                    cut = sel[e[:, 0]] != sel[e[:, 1]]
                    cols[np.unique(np.concatenate([e[cut, 0], e[cut, 1]]))] = np.zeros(3)
                dd = boundary_distance(mesh, adj, reg, lab, 2, marks[cid])
                scores.append(np.mean(dd) if dd else float("nan"))
                if dist is None:
                    # Framed on the MARKS, not on the vessel: the whole aorta is 200 mm
                    # long and the question lives in the 30 mm around the subclavian.
                    spread = float(np.linalg.norm(pts - focus, axis=1).max())
                    dist = max(ZOOM * max(spread, 6.0), 90.0)
                F = np.asarray(mesh.faces)
                poly = pv.PolyData(V, np.hstack([np.full((len(F), 1), 3), F]).ravel())
                poly.point_data["rgb"] = cols.astype(np.uint8)
                for c, view in enumerate(VIEWS):
                    pl.subplot(r, c)
                    pl.add_mesh(poly.copy(), scalars="rgb", rgb=True,
                                smooth_shading=True, specular=0.1, diffuse=0.92,
                                ambient=0.35)
                    pl.add_points(pts, color="#e11d48", point_size=11,
                                  render_points_as_spheres=True)
                    v = view / np.linalg.norm(view)
                    pl.camera.position = focus + v * dist
                    pl.camera.focal_point = focus
                    pl.camera.up = UP
            img = np.asarray(pl.screenshot(None, return_img=True))
            pl.close()
            iio.imwrite(str(outdir / f"{cid}.png"), img)
            print(f"  {cid:10s} " + "  ".join(
                f"{r} {v:.1f} mm" for r, v in zip(rules, scores))
                + f"  -> {outdir / (cid + '.png')}", flush=True)
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}", flush=True)
    for k, v in RULES["current"].items():
        setattr(ZB, k, v)


if __name__ == "__main__":
    main()
