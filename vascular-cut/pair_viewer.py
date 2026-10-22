"""The same patient branched and unbranched, side by side, one case at a time.

    .venv\\Scripts\\python pair_viewer.py

    right / space / n     next case          g   makes sense
    left / p              previous case      b   does not make sense
    f   next unjudged                        u   undecided (clears a verdict)
    r   back to the starting orientation     s   save a screenshot
    q   quit

The two views are LINKED: rotating either rotates both. The meshes come from the same
segmentation in the same coordinates, so a feature stays in the same place on screen in
both panels and the eye can compare them directly rather than by memory.

WHY THIS EXISTS

Paired over the same 30 patients, delta K comes out

    branched = 0.90 x unbranched + 0.00384        R2 0.828

-- slope near one, and an intercept twice the size of a normal aorta's entire delta K.
So the three branch nubs add a near-constant lump that has almost nothing to do with
the aorta's shape (the offset correlates with the unbranched value at rho -0.23), and
it compresses the normal-to-diseased separation from 3.6x to 1.7x.

That is a claim about geometry, and the numbers cannot show whether it is the nubs
themselves, the rims where they were cut, or something about how they were trimmed.
Looking at the pairs can. Each header carries both values and the offset, so a case
whose offset is unusual can be found and inspected rather than assumed about.

Verdicts go to pair_review.json beside the branched meshes.
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import nub_viewer as NV                                            # noqa: E402

BRANCHED = "#c2410c"
UNBRANCHED = "#1d4ed8"
MESH_RGB = (0.86, 0.86, 0.88)
SEEDS = list(range(42, 50))


def read_pairs(branched_glob, unbranched_glob):
    """Both runs as {case: (mean, se)}, keyed the same way."""
    def load(pat):
        by = {}
        for f in glob.glob(str(ROOT / pat)):
            for line in open(f, encoding="utf-8"):
                if not line.strip():
                    continue
                r = json.loads(line)
                if "error" in r:
                    continue
                by.setdefault(r["case"], {})[r["seed"]] = (r["dk"], r["area"])
        out = {}
        for c, v in by.items():
            if len(v) < 2:
                continue
            x = np.array([v[s][0] for s in sorted(v)])
            out[c] = (float(x.mean()), float(x.std(ddof=1) / np.sqrt(len(x))),
                      float(list(v.values())[0][1]))
        return out
    return load(branched_glob), load(unbranched_glob)


class PairViewer:
    def __init__(self, cases, review, folder):
        self.cases = cases
        self.review = review
        self.folder = folder
        self.i = 0
        self.p = pv.Plotter(shape=(1, 2), window_size=(1800, 950), border=False)
        self.p.set_background("white")
        for k, d in (("Right", +1), ("space", +1), ("n", +1), ("Left", -1), ("p", -1)):
            self.p.add_key_event(k, (lambda dd: (lambda: self.step(dd)))(d))
        self.p.add_key_event("g", lambda: self.verdict("good"))
        self.p.add_key_event("b", lambda: self.verdict("bad"))
        self.p.add_key_event("u", lambda: self.verdict(None))
        self.p.add_key_event("f", self.next_unjudged)
        self.p.add_key_event("r", self.home)
        self.p.add_key_event("s", self.shot)

    def key(self):
        return self.cases[self.i][0]

    def step(self, d):
        self.i = (self.i + d) % len(self.cases)
        self.draw(home=True)

    def verdict(self, v):
        self.review.set(self.key(), v)
        self.step(+1) if v is not None else self.draw()

    def next_unjudged(self):
        for d in range(1, len(self.cases) + 1):
            j = (self.i + d) % len(self.cases)
            if self.review.get(self.cases[j][0]) is None:
                self.i = j
                self.draw(home=True)
                return
        print("  every case has a verdict")

    def home(self):
        self.draw(home=True)

    def shot(self):
        out = self.folder / f"{self.key()}_pair.png"
        self.p.screenshot(str(out))
        print(f"  saved {out}")

    def draw(self, home=False):
        cid, bstl, ustl, bd, ud = self.cases[self.i]
        # ONCE, before the loop. Plotter.clear_actors() clears EVERY renderer, not
        # the active subplot, so calling it inside the loop drew the left panel and
        # then wiped it while setting up the right one.
        self.p.clear_actors()
        for col, (stl, colour, label, stat) in enumerate(
                ((bstl, BRANCHED, "BRANCHED (nubs)", bd),
                 (ustl, UNBRANCHED, "UNBRANCHED", ud))):
            self.p.subplot(0, col)
            m = trimesh.load(str(stl), process=False)
            m.merge_vertices()
            poly = pv.PolyData(np.asarray(m.vertices),
                               np.hstack([np.full((len(m.faces), 1), 3),
                                          np.asarray(m.faces)]).ravel())
            self.p.add_mesh(poly, color=MESH_RGB, smooth_shading=True, specular=0.2)
            txt = f"{label}\n"
            if stat:
                txt += f"dK {stat[0]:.5f} +/- {stat[1]:.5f}    area {stat[2]:,.0f} mm2"
            self.p.add_text(txt, position="upper_left", font_size=11, color=colour)
        self.p.subplot(0, 0)
        good, bad = self.review.counts()
        v = self.review.get(cid)
        mark = {"good": "MAKES SENSE", "bad": "DOES NOT", None: "-- not judged --"}[v]
        off = (bd[0] - ud[0]) if (bd and ud) else float("nan")
        ratio = (bd[0] / ud[0]) if (bd and ud and ud[0]) else float("nan")
        self.p.add_text(
            f"[{self.i + 1}/{len(self.cases)}]  {cid}   "
            f"{'normal' if cid.startswith(('KY', 'SA')) else 'diseased'}\n"
            f"{mark}     judged: {good} sense, {bad} not\n"
            f"branch offset {off:+.5f}   ({ratio:.1f}x)",
            position="lower_left", font_size=11, color="black")
        self.p.subplot(0, 1)
        self.p.add_text("g makes sense   b does not   u clear   f next unjudged\n"
                        "<- ->  move    r reset view    s screenshot    q quit",
                        position="lower_left", font_size=9, color="#555555")
        self.p.link_views()
        if home:
            for col in (0, 1):
                self.p.subplot(0, col)
                self.p.view_vector(NV.DEFAULT_VIEW, viewup=NV.DEFAULT_UP)
                self.p.reset_camera()
        self.p.render()

    def run(self):
        self.draw(home=True)
        self.p.show()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--branched", default="output/nubs_final")
    ap.add_argument("--unbranched", default="output/unbranched_src")
    ap.add_argument("--order", choices=("area", "offset", "case"), default="area",
                    help="area (default), offset (biggest branch contribution first), case")
    args = ap.parse_args()

    bdir, udir = ROOT / args.branched, ROOT / args.unbranched
    bstats, ustats = read_pairs("tools/cohort_seeds_[0-9]*.jsonl",
                                "tools/cohort_seeds_unbranched_*.jsonl")
    cases = []
    for stl in sorted(bdir.glob("*_nub.stl")):
        cid = stl.name.replace("_nub.stl", "")
        u = udir / f"{cid}_camp_input.stl"
        if not u.exists():
            continue
        cases.append((cid, stl, u, bstats.get(cid), ustats.get(cid)))
    if not cases:
        raise SystemExit(f"no pairs found between {bdir} and {udir}")
    if args.order == "area":
        cases.sort(key=lambda t: t[3][2] if t[3] else 0)
    elif args.order == "offset":
        cases.sort(key=lambda t: -((t[3][0] - t[4][0]) if (t[3] and t[4]) else 0))

    missing = [c[0] for c in cases if not (c[3] and c[4])]
    print(f"{len(cases)} pairs, ordered by {args.order}")
    if missing:
        print(f"  no delta K for: {', '.join(missing)} (mesh still shown)")
    review = NV.Review(bdir / "pair_review.json")
    PairViewer(cases, review, bdir).run()


if __name__ == "__main__":
    raise SystemExit(main())
