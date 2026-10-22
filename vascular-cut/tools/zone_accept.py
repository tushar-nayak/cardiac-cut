"""Acceptance test for zone geometry, calibrated on cases a human approved.

    .venv\\Scripts\\python tools\\zone_accept.py
    .venv\\Scripts\\python tools\\zone_accept.py --ref KK15_1,KR6_1,KR47_1

WHY THIS EXISTS. Every property stated as a requirement -- boundaries planar, zones
tiling once, one piece each, branches whole, zone 3 owning nothing -- has a check, and
all of them pass on versions that were then rejected on sight. The thing being judged
was PROPORTION: whether a zone is the right size for the stretch of aorta it covers.
Nothing measured that, so there was no way to tell an improvement from a regression
without asking, and no way to work unattended.

The reference cases supply it. On a bend a zone covers more outer wall than inner, so
area does not scale with span -- but the RATIO of a zone's area share to its span share
is scale-free, near 1 on a straight tube, and above 1 where the arch is tight. The
approved cases fix the band that ratio may occupy. Anything outside it is a zone the
wrong size for its span, which is what "zone 2 is infringing on zone 3" looks like as a
number.

Zone 4 is excluded from the shares: it runs the whole descending aorta and would swamp
everything else.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
import branch_cut as BC                                      # noqa: E402
from anat_bench import find_cases                            # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
TILING = ("Ishimaru zone 0", "Ishimaru zone 1", "Ishimaru zone 2",
          "Ishimaru zone 3", "Ishimaru zone 4")
SHARE_ZONES = 4                       # zones 0..3; zone 4 is the whole descending


def label_case(s):
    """Exactly what tools/zone_map.py produces, so the test judges what ships."""
    mesh = s.mesh
    route = BC.aorta_route(s.tree)
    L = float(route.length)
    station, _own, onr = RS.station_field(mesh, s.tree, route)
    is_branch = ~onr
    u = RS.harmonic_coordinate_capped(mesh, route, is_branch)
    if u is None:
        u = RS.harmonic_coordinate(mesh, route)
    ss, uu = RS.calibrate_u_to_station(mesh, route, u)
    zones = sorted((p for p in s.presets() if p["name"] in TILING),
                   key=lambda p: p["s0"])
    lab = -np.ones(len(u), int)
    for zi, p in enumerate(zones):
        lo = -1.0 if zi == 0 else float(np.interp(p["s0"], ss, uu))
        hi = (2.0 if (p["s1"] is None or zi == len(zones) - 1)
              else float(np.interp(p["s1"], ss, uu)))
        m = (~is_branch) & (u >= lo) & (u < hi)
        inb = RS.branch_inside(np.where(is_branch, station, -1e9),
                               p["s0"], p["s1"], L)
        lab[m | (is_branch & inb)] = zi
    lab, _ = RS.enforce_connected_labels(mesh, lab)
    adj = RS.surface_graph(mesh)
    for zi in range(len(zones) - 1):
        lab, _ = RS.flatten_boundary(mesh, lab, zi, zi + 1, adj)
    return mesh, route, zones, lab, is_branch


def shape_ratios(mesh, zones, lab, is_branch, L):
    """Area share / span share, per zone, over zones 0..3.

    Wall only: a branch's area belongs to the vessel, not to the stretch of aorta the
    zone covers, and including it would make a zone look large because a big vessel
    happens to leave there.
    """
    F = np.asarray(mesh.faces, int)
    wall_f = ~is_branch[F].any(axis=1)
    fa = np.asarray(mesh.area_faces, float)
    n = min(SHARE_ZONES, len(zones))
    areas, spans = [], []
    for zi in range(n):
        p = zones[zi]
        sel = wall_f & (lab[F] == zi).all(axis=1)
        areas.append(float(fa[sel].sum()))
        hi = L if p["s1"] is None else float(p["s1"])
        spans.append(max(hi - float(p["s0"]), 1e-6))
    areas, spans = np.array(areas), np.array(spans)
    if areas.sum() <= 0:
        return None
    return (areas / areas.sum()) / (spans / spans.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--ref", default="KK15_1,KR6_1,KR47_1")
    ap.add_argument("--tol", type=float, default=1.5,
                    help="how many reference SDs a zone may sit from the reference "
                         "mean before it is flagged")
    args = ap.parse_args()
    refs = set(args.ref.split(","))

    rows = {}
    for cid, cl, mp in find_cases(args.cohort):
        try:
            s = CutSession(str(cl), str(mp))
            s.nub_preset()
            s.recut()
            mesh, route, zones, lab, isb = label_case(s)
            r = shape_ratios(mesh, zones, lab, isb, float(route.length))
            if r is None:
                continue
            rows[cid] = r
            print(f"  {cid:10s} " + "  ".join(f"z{i}:{v:5.2f}"
                                              for i, v in enumerate(r)))
        except Exception as e:                               # noqa: BLE001
            print(f"  {cid:10s} ERROR {type(e).__name__}: {e}")

    have = [c for c in refs if c in rows]
    if not have:
        raise SystemExit("none of the reference cases produced a result")
    n = min(len(rows[c]) for c in rows)
    R = np.array([rows[c][:n] for c in have])
    mu, sd = R.mean(axis=0), np.maximum(R.std(axis=0), 0.05)
    print(f"\n  reference ({', '.join(sorted(have))}): "
          + "  ".join(f"z{i} {mu[i]:.2f}+-{sd[i]:.2f}" for i in range(n)))

    print(f"\n  cases outside {args.tol} reference SD on any zone")
    bad = []
    for cid, r in sorted(rows.items()):
        z = np.abs(r[:n] - mu) / sd
        if (z > args.tol).any():
            worst = int(np.argmax(z))
            bad.append((float(z.max()), cid, worst, r[worst]))
    for zmax, cid, wi, val in sorted(bad, reverse=True):
        print(f"    {cid:10s} worst z{wi}: ratio {val:5.2f} "
              f"({zmax:4.1f} SD from reference)")
    print(f"\n  {len(rows) - len(bad)} of {len(rows)} within tolerance")
    out = ROOT / "tools" / "zone_accept.json"
    out.write_text(json.dumps(
        {"reference": sorted(have), "mean": mu.tolist(), "sd": sd.tolist(),
         "ratios": {k: v.tolist() for k, v in rows.items()}}, indent=1),
        encoding="utf-8")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
