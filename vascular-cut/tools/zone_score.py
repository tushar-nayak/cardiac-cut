"""Score a boundary rule against where a person put the line.

    .venv\\Scripts\\python tools\\zone_score.py
    .venv\\Scripts\\python tools\\zone_score.py --rules current,landmark_lag,rim

Reads the marks made in tools/zone_mark.py and, for each candidate rule, measures how
far that rule's boundary lands from them -- across the wall, in millimetres, which is
the distance a person is looking at.

THIS IS THE ONLY SCORE IN THIS DIRECTORY THAT IS NOT CHECKED AGAINST ITSELF. Every other
number here was derived from the same geometry the rule under test was derived from,
and three of them in a row endorsed a change that was rejected on sight. A mark is
independent by construction: it is where somebody said the line goes.

A rule is better only if it moves the boundary TOWARD the marks. Boundary flatness,
zone islands and the tiling checks still have to hold -- they are necessary, and what
they were never able to say is whether the zone is in the right place.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from vascular_api import CutSession                          # noqa: E402
import region_station as RS                                  # noqa: E402
from anat_bench import find_cases                            # noqa: E402
import zone_stub_bench as ZB                                 # noqa: E402
from zone_stub_bench import final_labels                     # noqa: E402

THORACIC = r"Z:/Martino/branches_CAMP/GT"
MARKS = ROOT / "tools" / "zone_marks.json"

# Each rule is a set of flags on zone_stub_bench, so a rule here IS the shipped code
# path with its switches thrown, not a reimplementation of it.
RULES = {
    "current":      dict(CROSS_SECTION_PLACEMENT=False, RIM_PLACEMENT=False,
                         LANDMARK_LAG_FIX=False),
    "cross_section": dict(CROSS_SECTION_PLACEMENT=True, RIM_PLACEMENT=False,
                          LANDMARK_LAG_FIX=False),
    "rim":          dict(CROSS_SECTION_PLACEMENT=False, RIM_PLACEMENT=True,
                         LANDMARK_LAG_FIX=False),
    "landmark_lag": dict(CROSS_SECTION_PLACEMENT=False, RIM_PLACEMENT=False,
                         LANDMARK_LAG_FIX=True),
}
# A nudge is not a flag, it is a distance, so these are generated: the 2/3 boundary is
# moved back by that many millimetres of wall and 3/4 with it, so zone 3 keeps its
# width. Marks say the offset is case-dependent, so a constant is not expected to win
# outright -- what is worth knowing is whether ANY constant beats leaving it alone.
for _mm in (2, 3, 5, 8, 12):
    RULES[f"nudge{_mm}"] = dict(CROSS_SECTION_PLACEMENT=False, RIM_PLACEMENT=False,
                                LANDMARK_LAG_FIX=False, _nudge_mm=float(_mm),
                                _plane=False)
# THESE TWO KNOBS ARE NOW INERT. The shipped anchor is RS.mouth_clearing_plane, which
# takes no percentile and no offset, so every "_plane" rule below scores the SAME
# labelling whatever pct and offset say. The sweep is kept only because the numbers it
# produced are quoted in tools/ZONES.md; do not read a fresh run of it as a comparison
# between anchors. A lip plus a constant was measured against eight cases judged by eye
# and failed in both directions at once -- see ZONES.md, "Neither lip".
#
# THE PLANE RULE, swept over its two knobs rather than guessed. Which rim vertex counts
# as the distal lip (pct) and how far past it the cut sits (offset) are the only choices
# it makes, and the marks can settle both.
for _pct in (70, 80, 90, 95):
    for _off in (0.0, 2.0, 4.0):
        RULES[f"plane{_pct}+{_off:g}"] = dict(
            CROSS_SECTION_PLACEMENT=False, RIM_PLACEMENT=False,
            LANDMARK_LAG_FIX=False, _plane=True, _pct=float(_pct),
            _offset=float(_off))
for _k in ("current", "cross_section", "rim", "landmark_lag"):
    RULES[_k]["_plane"] = False


def boundary_distance(mesh, adj, reg, lab, b, marks):
    """Millimetres of wall from each mark to the b/(b+1) boundary."""
    e = np.asarray(mesh.edges_unique, int)
    on = (((lab[e[:, 0]] == b) & (lab[e[:, 1]] == b + 1))
          | ((lab[e[:, 0]] == b + 1) & (lab[e[:, 1]] == b)))
    ring = np.unique(e[on])
    ring = ring[reg[ring] < 0]
    if len(ring) < 6:
        return None
    keep = np.flatnonzero(reg < 0)
    pos = -np.ones(len(reg), int)
    pos[keep] = np.arange(len(keep))
    seed = pos[ring]
    seed = seed[seed >= 0]
    if not len(seed):
        return None
    d = dijkstra(adj[keep][:, keep], directed=False, indices=seed, min_only=True)
    out = []
    for m in marks:
        v = pos[int(m["vertex"])]
        if v >= 0 and np.isfinite(d[v]):
            out.append(float(d[v]))
    return out or None


def shape_guards(mesh, lab, reg, adj):
    """(worst boundary flatness in mm, zone islands) -- what a rule must not break.

    Distance to the marks says whether the boundary is in the right PLACE. These say
    whether it is still a boundary: a rule that drags it toward the marks by shredding
    it has not improved anything, and the marks alone cannot see that.
    """
    from zone_islands import islands
    V = np.asarray(mesh.vertices, float)
    e = np.asarray(mesh.edges_unique, int)
    worst = 0.0
    zs = sorted(set(lab[lab >= 0].tolist()))
    for a, b in zip(zs, zs[1:]):
        on = (((lab[e[:, 0]] == a) & (lab[e[:, 1]] == b))
              | ((lab[e[:, 0]] == b) & (lab[e[:, 1]] == a)))
        ring = np.unique(e[on])
        ring = ring[reg[ring] < 0]
        if len(ring) < 8:
            continue
        P = V[ring]
        c = P.mean(axis=0)
        _u, _s, vt = np.linalg.svd(P - c, full_matrices=False)
        worst = max(worst, float(np.abs((P - c) @ vt[2]).max()))
    return worst, len(islands(mesh, lab))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohort", nargs="?", default=THORACIC)
    ap.add_argument("--rules", default="current")
    args = ap.parse_args()
    if not MARKS.exists():
        raise SystemExit(f"no marks yet -- run tools/zone_mark.py first ({MARKS})")
    marks = json.loads(MARKS.read_text())
    if not marks:
        raise SystemExit("the marks file is empty")
    rules = [r for r in args.rules.split(",") if r in RULES]
    if not rules:
        raise SystemExit(f"unknown rule; pick from {', '.join(RULES)}")

    cases = [c for c in find_cases(args.cohort) if c[0] in marks]
    print(f"{len(cases)} marked case(s), {sum(len(v) for v in marks.values())} mark(s)")
    per = {r: [] for r in rules}
    shape = {r: [] for r in rules}
    for cid, cl, mp in cases:
        s = CutSession(str(cl), str(mp))
        s.nub_preset()
        s.recut()
        line = f"  {cid:10s}"
        for r in rules:
            nudge = float(RULES[r].get("_nudge_mm", 0.0))
            for k, v in RULES[r].items():
                if not k.startswith("_"):
                    setattr(RS, k, v)   # the flags live in region_station now
            RS.PLANE_BOUNDARY = bool(RULES[r].get("_plane", False))
            RS.MOUTH_PCT = float(RULES[r].get("_pct", 95.0))
            RS.MOUTH_OFFSET_MM = float(RULES[r].get("_offset", 2.0))
            try:
                mesh, route, zones, lab, isb, adj, _f, reg, oz = final_labels(s)
                if nudge > 0:
                    frozen = reg >= 0
                    zs = sorted(set(lab[lab >= 0].tolist()))
                    if 2 in zs and 3 in zs:
                        lab, _ = RS.nudge_boundary(mesh, lab, 2, 3, nudge, adj, frozen)
                        if 4 in zs:
                            lab, _ = RS.nudge_boundary(mesh, lab, 3, 4, nudge, adj,
                                                       frozen)
                shape[r].append(shape_guards(mesh, lab, reg, adj))
                by_b = {}
                for m in marks[cid]:
                    by_b.setdefault(int(str(m["boundary"]).split("/")[0]),
                                    []).append(m)
                got = []
                for b, ms in by_b.items():
                    dd = boundary_distance(mesh, adj, reg, lab, b, ms)
                    if dd:
                        got += dd
                if got:
                    per[r] += got
                    line += f"   {r} {np.mean(got):5.1f} mm"
                else:
                    line += f"   {r}   n/a"
            except Exception as ex:                          # noqa: BLE001
                line += f"   {r} ERR {type(ex).__name__}"
        print(line, flush=True)
    for k, v in RULES["current"].items():
        setattr(RS, k, v)

    print()
    base = None
    for r in rules:
        d = np.array(per[r]) if per[r] else np.array([np.nan])
        # TRIMMED, because a stray click should not decide anything. KK32_0 carries one
        # mark 44 mm from every other, sitting on a different zone entirely.
        trim = (np.array(sorted(per[r])[:max(1, int(0.9 * len(per[r])))])
                if per[r] else d)
        fl = [x[0] for x in shape[r]] or [np.nan]
        isl = sum(x[1] for x in shape[r])
        line = (f"  {r:14s} marks: mean {np.nanmean(d):5.1f} mm  "
                f"trimmed {np.nanmean(trim):5.1f}  median {np.nanmedian(d):5.1f}  "
                f"| flatness worst {np.nanmax(fl):5.2f} mm  islands {isl}")
        if base is None:
            base = np.nanmean(trim)
        elif np.nanmean(trim) < base - 0.05:
            line += "   BETTER"
        print(line)
    print("\n  lower is better, and only this decides. Flatness, islands and the "
          "tiling checks still have to hold.")


if __name__ == "__main__":
    main()
