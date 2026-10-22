"""Read a zone_bench run and say whether the zones are SANE, not just whether they cut.

    .venv\\Scripts\\python tools\\zone_bench_report.py tools/zone_bench_GT.jsonl

"It cut" is a low bar for a level-set boundary: it cuts by construction. So the checks
here are the ones that could still be wrong -- a zone that came out in pieces, one that
kept wall from outside its own span, one whose surface is unusable for measuring, and
one whose length bears no relation to the span asked for.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def pct(n, d):
    return f"{100.0 * n / d:.0f}%" if d else "--"


def main():
    path = Path(sys.argv[1] if len(sys.argv) > 1
                else ROOT / "tools" / "zone_bench_GT.jsonl")
    recs = [json.loads(x) for x in path.open(encoding="utf-8")]
    print(f"{path.name}: {len(recs)} case(s)")

    errs = [r for r in recs if "error" in r]
    if errs:
        print(f"\n{len(errs)} case(s) failed outright:")
        for r in errs:
            print(f"  {r['case']:12s} {r['error']}")
    recs = [r for r in recs if "error" not in r]

    base_bad = [r for r in recs if not r.get("base_ok")]
    if base_bad:
        print(f"\n{len(base_bad)} case(s) whose BRANCH cut already fails "
              f"(zone results there are on a mesh the tool would not hand over):")
        for r in base_bad:
            print(f"  {r['case']:12s} {r['base_msg']}")

    zones = [(r["case"], z) for r in recs for z in r.get("zones", [])]
    print(f"\n{len(zones)} zone preset(s) over {len(recs)} case(s)")

    # ---- availability -------------------------------------------------------
    by_name = defaultdict(list)
    for cid, z in zones:
        by_name[z["name"]].append((cid, z))
    print(f"\n{'zone':26s} {'offered':>8s} {'station':>9s} {'planar':>9s}")
    for name in sorted(by_name):
        zz = [z for _, z in by_name[name]]
        s = sum(1 for z in zz if z.get("station_ok"))
        p = sum(1 for z in zz if z.get("planar_ok"))
        has_p = sum(1 for z in zz if "planar_ok" in z)
        print(f"  {name:24s} {len(zz):8d} {s:5d} {pct(s, len(zz)):>4s} "
              f"{p:5d} {pct(p, has_p):>4s}")

    # ---- what the planar path refused --------------------------------------
    ref = [(c, z) for c, z in zones if "planar_ok" in z and not z["planar_ok"]]
    if ref:
        print(f"\n{len(ref)} planar refusal(s), all of which the station cut made:")
        for c, z in ref:
            print(f"  {c:12s} {z['name']:24s} {z.get('planar_msg', '')[:70]}")

    # ---- is the station zone SANE ------------------------------------------
    print("\nstation-cut sanity")
    bad = defaultdict(list)
    lens = []
    for c, z in zones:
        if not z.get("station_ok"):
            bad["cut raised"].append((c, z["name"], z.get("station_err", "")))
            continue
        st = z["station"]
        if not st["winding_ok"]:
            bad["winding inconsistent"].append((c, z["name"], ""))
        if st["degenerate"]:
            bad["zero-area faces"].append((c, z["name"], str(st["degenerate"])))
        if st["n_loops"] < 1:
            bad["no boundary loop"].append((c, z["name"], ""))
        if st["faces"] < 50:
            bad["almost nothing kept"].append((c, z["name"],
                                               f"{st['faces']} faces"))
        if st["faces_dropped_as_islands"] > 0.02 * max(st["faces"], 1):
            bad["came out in pieces"].append(
                (c, z["name"], f"{st['faces_dropped_as_islands']} face(s) discarded"))
        for k, lab in (("overrun_prox", "proximal"), ("overrun_dist", "distal")):
            v = st.get(k)
            if v is not None and v > 1.0:
                bad[f"kept wall past its {lab} bound"].append(
                    (c, z["name"], f"{v:.1f} mm"))
        if z["s1"] is not None and st["wall_min"] is not None:
            asked = z["s1"] - z["s0"]
            got = st["wall_max"] - st["wall_min"]
            lens.append((c, z["name"], asked, got))
    if not bad:
        print("  no defects found")
    for k in sorted(bad):
        v = bad[k]
        print(f"  {k}: {len(v)}")
        for c, n, extra in v[:8]:
            print(f"      {c:12s} {n:24s} {extra}")
        if len(v) > 8:
            print(f"      ... and {len(v) - 8} more")

    # ---- does the zone span what it says it spans ---------------------------
    if lens:
        r = np.array([g / a for _, _, a, g in lens if a > 1e-6])
        print(f"\nzone length, wall span / span asked for  (n={len(r)})")
        print(f"  median {np.median(r):.3f}   "
              f"5-95% {np.percentile(r, 5):.3f}-{np.percentile(r, 95):.3f}")
        worst = sorted(lens, key=lambda t: abs(t[3] / max(t[2], 1e-6) - 1.0))[-6:]
        for c, n, a, g in reversed(worst):
            print(f"    {c:12s} {n:24s} asked {a:6.1f} mm, wall spans {g:6.1f} mm")

    # ---- what the planar path costs -----------------------------------------
    dr = [(c, z["name"], z.get("planar_drift_prox"), z.get("planar_drift_dist"))
          for c, z in zones if "planar_drift_prox" in z]
    flat = [abs(x) for _, _, a, b in dr for x in (a, b) if x is not None]
    if flat:
        f = np.array(flat)
        print(f"\nplanar boundary drift from the station anatomy asks for "
              f"(n={len(f)} boundaries)")
        print(f"  median {np.median(f):.1f} mm   90th {np.percentile(f, 90):.1f} mm   "
              f"max {f.max():.1f} mm")
        print(f"  over 5 mm: {int((f > 5).sum())} ({pct(int((f > 5).sum()), len(f))})"
              f"   over 15 mm: {int((f > 15).sum())} "
              f"({pct(int((f > 15).sum()), len(f))})")
        worst = sorted(((max(abs(a or 0), abs(b or 0)), c, n)
                        for c, n, a, b in dr), reverse=True)[:8]
        for v, c, n in worst:
            print(f"    {c:12s} {n:24s} {v:6.1f} mm")


if __name__ == "__main__":
    main()
