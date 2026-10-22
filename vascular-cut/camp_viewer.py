"""Scroll through the cohort in area order and see each geometry beside its delta K.

    .venv\\Scripts\\python camp_viewer.py --patches      # plasma integrated-K surface
    .venv\\Scripts\\python camp_viewer.py                # plain mesh with rims
    .venv\\Scripts\\python camp_viewer.py --patches --compare tools/camp_nub_constpatch.jsonl

Same window and the same keys as nub_viewer, which this reuses outright, so the mesh,
the rims and their labels are drawn by the code that already draws them:

    right / space / n     next case          g   makes sense
    left / p              previous case      b   does not make sense
    f   next unjudged                        u   undecided (clears a verdict)
    r   back to the starting orientation     s   save a screenshot
    q   quit

--patches colours the surface by INTEGRATED K, the way campsrc's own figure does it:
one value per face, taken from the patch it belongs to, on Plasma, with the patch
seams drawn. It reads a cache built by tools/camp_patch_cache.py, because CAMP takes
most of a minute per case and that is unusable behind an arrow key. Two more keys:

    c   colour limits: this case 2-98 pct (default), this case full range,
        or one scale shared across the whole cohort -- all centred on zero
    e   patch seams on / off (off by default; ~17,000 edges per aorta)

Every scale is symmetric about zero, so the midpoint colour always means K = 0 and
the same colour means the same thing on every case.

Per-case limits show each aorta's own structure; the cohort scale is the only way to
compare two of them, since a case whose K never leaves +/-0.05 and one that reaches
-0.86 look identical when each is stretched to fill its own scale. The default is the
2nd-98th percentile of the case's own patches, because on the full range a handful of
outliers own the scale and everything else saturates into one colour.

WHAT THIS IS FOR

delta K rises with surface area across this cohort, and the question is whether that
is anatomy or arithmetic. A plot cannot answer it -- a point at 3.9 looks the same
whether the aorta is genuinely contorted or the number is an artefact. So the cases are
ordered by ASCENDING AREA and shown one at a time with the number attached: walking the
order is walking up the trend, and the geometry either gets visibly more contorted as
the number climbs or it does not.

Normalised the way the cohort plots are normalised -- to the mean of the KY and SA
group -- so 1.00 is a typical normal and the figure on screen is the figure on the
plot. The bar is that value drawn to scale, for reading the run at a glance without
comparing decimals.

--compare puts a second run on the same geometry, which is how a change to CAMP gets
judged: same aorta, same view, both numbers, and the rank each run gives it.

Verdicts go to camp_review.json in the geometry folder, written on every keystroke.
It is kept separate from review.json, which holds verdicts on the CUTS -- a geometry
can be cut correctly and still carry a delta K that makes no sense.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import nub_viewer as NV                                            # noqa: E402

NORMAL_PREFIXES = ("KY", "SA")
BAR_WIDTH = 34             # characters at the widest normalised value on screen
PATCH_DIR = "output/camp_intk"
SEAM = "#20202a"
CLIM_MODES = ("full", "robust", "cohort")
CLIM_LABEL = {"full": "this case, full range",
              "robust": "this case, 2-98 pct",
              "cohort": "COHORT scale"}


def _symmetric(lo, hi):
    """A range centred on zero, so the midpoint colour always means K = 0.

    Integrated K is signed -- negative at a saddle, positive at a cap -- and its
    range is lopsided, typically reaching -2.5 while barely passing +0.5. Fitted to
    that raw range the colour bar puts zero at four fifths of its length, so the
    colour standing for flat has to be looked up on every case and the same colour
    means different things on two of them. Widening to +/-max(|lo|,|hi|) costs some
    contrast on the short side and buys a scale that reads the same everywhere.
    """
    m = max(abs(float(lo)), abs(float(hi)))
    return (-m, m) if m > 0 else (-1.0, 1.0)


def read_run(path, trim):
    """One CAMP run as {case: {...}}, normalised to the KY+SA mean.

    Ranks are computed here rather than at draw time so a case carries its own place
    in the run, which is the thing --compare is actually comparing: two runs can agree
    on every value to two decimals and still order the cohort differently.
    """
    recs = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        if "error" in r:
            continue
        dk = r["dk"].get(str(trim))
        if dk is None:
            raise SystemExit(f"{path} has no trim {trim}; it has "
                             f"{', '.join(sorted(r['dk']))}")
        recs[r["case"]] = dict(r, dk_at=float(dk))
    if not recs:
        raise SystemExit(f"no usable records in {path}")
    norm = [v["dk_at"] for k, v in recs.items() if k.startswith(NORMAL_PREFIXES)]
    base = float(np.mean(norm)) if norm else float(np.mean([v["dk_at"] for v in recs.values()]))
    order = sorted(recs, key=lambda c: recs[c]["dk_at"])
    for c, v in recs.items():
        v["norm"] = v["dk_at"] / base if base else float("nan")
        v["rank"] = order.index(c) + 1
    return recs


class CampViewer(NV.Viewer):
    def __init__(self, cases, review, folder, runs, labels, trim, patch_dir=None,
                 off_screen=False):
        super().__init__(cases, review, folder, off_screen=off_screen)
        self.runs = runs
        self.labels = labels
        self.trim = trim
        self.span = max((v["norm"] for r in runs for v in r.values()), default=1.0)
        self.patch_dir = patch_dir
        # Off by default. At ~850 patches an aorta carries ~17,000 seam edges, which
        # at this scale reads as a wireframe laid over the surface and hides the very
        # colours it is there to delimit. The colouring alone shows the patches.
        self.seams = False
        self.clim_mode = 1          # robust by default; see CLIM_MODES
        self._cache = {}
        self._cohort_clim = None
        if patch_dir is not None and not off_screen:
            self.p.add_key_event("c", self.toggle_clim)
            self.p.add_key_event("e", self.toggle_seams)

    # -- integrated-K colouring ---------------------------------------------
    def toggle_clim(self):
        self.clim_mode = (self.clim_mode + 1) % len(CLIM_MODES)
        self.draw(reset_camera=False)

    def clim(self, d):
        """Colour limits for the case in hand, under the current mode.

        FULL is what campsrc does and is honest but usually unreadable: a handful of
        patches reach K = -0.9 while the bulk sits near zero, so the bulk saturates
        into one colour and the surface looks uniform when it is not. ROBUST clips to
        the 2nd-98th percentile of PATCH values, which puts the contrast where the
        patches actually are and pushes the outliers to the ends of the bar -- they
        stay visible as extremes, they just stop owning the scale.
        """
        K = d["patch_K"]
        mode = CLIM_MODES[self.clim_mode]
        if mode == "cohort":
            lo, hi = self.cohort_clim()
        elif mode == "robust":
            lo, hi = (float(x) for x in np.percentile(K, [2, 98]))
            if hi <= lo:
                lo, hi = float(K.min()), float(K.max())
        else:
            lo, hi = float(K.min()), float(K.max())
        return _symmetric(lo, hi)

    def toggle_seams(self):
        self.seams = not self.seams
        self.draw(reset_camera=False)

    def patches(self, cid):
        """The cached partition for one case, kept in memory once loaded.

        Scrolling back and forth over thirty cases is the normal way this gets used,
        so re-reading and re-decompressing an npz on every keystroke is worth avoiding
        even though each one is only a couple of megabytes.
        """
        if cid not in self._cache:
            f = Path(self.patch_dir) / f"{cid}.npz"
            if not f.exists():
                self._cache[cid] = None
            else:
                with np.load(f) as z:
                    self._cache[cid] = {k: z[k] for k in z.files}
        return self._cache[cid]

    def cohort_clim(self):
        """One colour range covering every case, so two of them can be compared.

        Taken over PATCH values rather than face values: they are the same numbers,
        but there are ~850 of them per case instead of 100,000, and every case has to
        be opened to find the range.
        """
        if self._cohort_clim is None:
            lo, hi = np.inf, -np.inf
            for _, rec in self.cases:
                d = self.patches(rec["case"])
                if d is None:
                    continue
                lo = min(lo, float(d["patch_K"].min()))
                hi = max(hi, float(d["patch_K"].max()))
            self._cohort_clim = (lo, hi) if lo < hi else (-1.0, 1.0)
        return self._cohort_clim

    def draw(self, reset_camera=True, home=False):
        """The plain mesh, or the integrated-K surface when a cache is loaded.

        Falls back to nub_viewer's own draw for any case with no cached partition, so
        a half-built cache still scrolls rather than crashing on the gap.
        """
        if self.patch_dir is None:
            return super().draw(reset_camera=reset_camera, home=home)
        _, rec = self.cases[self.i]
        d = self.patches(rec["case"])
        if d is None:
            return super().draw(reset_camera=reset_camera, home=home)

        # Same camera contract as nub_viewer: keep the direction, re-frame the
        # distance, because each case sits at its own scanner coordinates.
        look, up = (NV.DEFAULT_VIEW, NV.DEFAULT_UP) if home else self.look_direction()
        cam = None if reset_camera else self.p.camera_position
        self.p.clear_actors()

        verts = d["verts"].astype(float)
        faces = d["faces"]
        poly = pv.PolyData(verts, np.hstack([np.full((len(faces), 1), 3),
                                             faces]).ravel())
        poly.cell_data["integrated K"] = d["face_K"]
        clim = self.clim(d)
        # Flat shading, not smooth: the value is per PATCH, and interpolating it
        # across a patch boundary draws a gradient where the partition has a step.
        self.p.add_mesh(poly, scalars="integrated K", cmap="plasma", clim=clim,
                        smooth_shading=False, show_scalar_bar=True,
                        scalar_bar_args=dict(title="integrated K", color="black",
                                             n_labels=5, vertical=True,
                                             position_x=0.89, position_y=0.28,
                                             height=0.46, width=0.04))
        if self.seams and len(d["seam_edges"]):
            e = d["seam_edges"]
            lines = pv.PolyData(verts)
            lines.lines = np.hstack([np.full((len(e), 1), 2), e]).ravel()
            self.p.add_mesh(lines, color=SEAM, line_width=1, opacity=0.25)

        self.p.add_text(self._header(rec), position="upper_left", font_size=11,
                        color="black")
        self.p.add_text(self._footer(rec), position="lower_left", font_size=9,
                        color="#555555")
        if cam is not None:
            self.p.camera_position = cam
        else:
            self.p.view_vector(look, viewup=up)
            self.p.reset_camera()
        self.p.render()

    def _bar(self, x):
        n = int(round(BAR_WIDTH * x / self.span)) if self.span > 0 else 0
        return "#" * max(n, 1 if x > 0 else 0)

    def _header(self, rec):
        cid = rec["case"]
        good, bad = self.review.counts()
        v = self.review.get(cid)
        mark = {"good": "MAKES SENSE", "bad": "DOES NOT",
                None: "-- not judged --"}[v]
        primary = self.runs[0].get(cid)
        area = primary["area"] if primary else float("nan")
        group = "normal" if cid.startswith(NORMAL_PREFIXES) else "diseased"
        lines = [f"[{self.i + 1}/{len(self.cases)}]  {cid}   {group}   "
                 f"{area:,.0f} mm2   (area rank {self.i + 1})",
                 f"{mark}          judged: {good} sense, {bad} not"]
        for run, label in zip(self.runs, self.labels):
            r = run.get(cid)
            if r is None:
                lines.append(f"{label:>14s}  -- not in this run --")
                continue
            lines.append(f"{label:>14s}  dK {r['norm']:5.2f} x normal   "
                         f"(dK rank {r['rank']}/{len(run)}, "
                         f"{r['n_patches']} patches, abar {r['area']/r['n_patches']:.0f} mm2)")
            lines.append(f"{'':>14s}  {self._bar(r['norm'])}")
        return "\n".join(lines)

    def _footer(self, rec):
        lines = [f"trim {self.trim:g}   normalised to the KY+SA mean   "
                 f"ordered by ascending area"]
        d = self.patches(rec["case"]) if self.patch_dir else None
        if d is not None:
            lo, hi = self.clim(d)
            top = int(np.argmax(d["patch_share"]))
            lines.append(
                f"integrated K  [{lo:+.3f}, {hi:+.3f}]  "
                f"{CLIM_LABEL[CLIM_MODES[self.clim_mode]]}   |   "
                f"biggest single patch carries {100 * d['patch_share'][top]:.0f}% of "
                f"this case's dK  (K {d['patch_K'][top]:+.3f}, "
                f"{d['patch_area'][top]:.0f} mm2)")
        for n in rec.get("notes", []):
            lines.append(n)
        keys = ("g makes sense   b does not   u clear   f next unjudged   "
                "<- ->  move   r reset view   s screenshot   q quit")
        lines.append(keys + ("   |   c scale   e seams" if self.patch_dir else ""))
        return "\n".join(lines[-5:])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="output/nubs_final")
    ap.add_argument("--results", default="tools/camp_nub_bothfix_100k.jsonl")
    ap.add_argument("--compare", default=None,
                    help="a second run to show on the same geometry")
    ap.add_argument("--trim", type=float, default=0.0)
    ap.add_argument("--patches", action="store_true",
                    help="colour by integrated K (needs tools/camp_patch_cache.py)")
    ap.add_argument("--patch-dir", default=PATCH_DIR)
    args = ap.parse_args()

    folder = Path(args.dir)
    if not folder.exists():
        raise SystemExit(f"no such folder: {folder}")
    runs = [read_run(args.results, args.trim)]
    labels = [Path(args.results).stem.replace("camp_nub_", "")]
    if args.compare:
        runs.append(read_run(args.compare, args.trim))
        labels.append(Path(args.compare).stem.replace("camp_nub_", ""))

    cases = NV.load_cases(folder)
    # Ascending area, and only what the primary run actually measured -- a geometry
    # with no number would be a blank stop on a walk up the trend.
    keyed = [(stl, rec) for stl, rec in cases if rec["case"] in runs[0]]
    missing = sorted({rec["case"] for _, rec in cases} - set(runs[0]))
    keyed.sort(key=lambda t: runs[0][t[1]["case"]]["area"])
    if not keyed:
        raise SystemExit(f"no geometry in {folder} matches {args.results}")

    print(f"{len(keyed)} geometries, area "
          f"{runs[0][keyed[0][1]['case']]['area']:,.0f} to "
          f"{runs[0][keyed[-1][1]['case']]['area']:,.0f} mm2, trim {args.trim:g}")
    for label, run in zip(labels, runs):
        vals = [v["norm"] for c, v in run.items() if c.startswith(NORMAL_PREFIXES)]
        if vals:
            print(f"  {label:>14s}  normals {min(vals):.2f}-{max(vals):.2f}  "
                  f"CV {np.std(vals, ddof=1) / np.mean(vals):.2f}")
    if missing:
        print(f"  no number for: {', '.join(missing)}")

    pdir = None
    if args.patches:
        pdir = Path(args.patch_dir)
        have = {f.stem for f in pdir.glob("*.npz")} if pdir.exists() else set()
        want = {rec["case"] for _, rec in keyed}
        if not have:
            raise SystemExit(f"no cached partitions in {pdir}; run\n"
                             f"    .venv\\Scripts\\python tools\\camp_patch_cache.py")
        if want - have:
            print(f"  no cached partition for: {', '.join(sorted(want - have))} "
                  f"(those cases fall back to the plain mesh)")

    review = NV.Review(folder / "camp_review.json")
    CampViewer(keyed, review, folder, runs, labels, args.trim, pdir).run()


if __name__ == "__main__":
    raise SystemExit(main())
