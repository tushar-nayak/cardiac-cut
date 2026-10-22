"""delta K in the web app: the job, its result, and both ways it can fail.

The computation runs in a process of its own, so almost everything that can go wrong
here is plumbing rather than mathematics: a job that never starts, a failure that
arrives as a blank panel, arrays that do not line up with the faces they colour. None
of that is visible from the browser, and a 30-second run is not something you debug by
clicking. Drive it headlessly.
"""
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "webapp"))

from shape import ShapeJob, ShapeError                     # noqa: E402
import shape_metrics as SM                                  # noqa: E402

SURFACE = ROOT / "output" / "KY160_1_open.stl"
fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


def wait(job, limit=240):
    """Poll the way the app does, and report the phases it passed through."""
    seen, t0 = [], time.monotonic()
    while True:
        st = job.poll()
        if st["state"] == "running":
            if st["message"] not in seen:
                seen.append(st["message"])
            if time.monotonic() - t0 > limit:
                job.cancel()
                return {"state": "failed", "message": "timed out"}, seen
            time.sleep(0.25)
            continue
        return st, seen


print("=" * 78)
print("a surface that already exists on this machine")
print("=" * 78)
if not SURFACE.exists():
    raise SystemExit(f"missing {SURFACE}; run the cut tool once to produce it")

job = ShapeJob()
job.adopt_path(SURFACE)
info = job.source_info()
print(f"  {info}")
check(info["euler"] == -3, "topology read from the file, not from the cut that made it")
check(info["openings"] == 5, "five openings")
check(not info["closed"], "not reported as closed")

try:
    job.adopt_path(ROOT / "no_such_file.stl")
    check(False, "a missing file is refused")
except ShapeError as e:
    check("does not exist" in str(e), f"a missing file is refused: {e}")

print()
print("=" * 78)
print("the run")
print("=" * 78)
t0 = time.monotonic()
job.start(50_000, 0.01)
st, phases = wait(job)
print(f"  {time.monotonic() - t0:.1f}s, state={st['state']}")
for p in phases:
    print(f"    phase: {p}")
check(st["state"] == "done", f"the run finished: {st['message']}")
# What the poller happened to CATCH is a race: CAMP is silent for most of a run and
# the phases move faster than a 0.25 s poll on a warm machine. Asserting it saw two
# made this suite fail about half the time. Check the phases it did see are real
# ones, then test the phase logic itself directly, where there is no timing at all.
KNOWN = {"Remeshing to a uniform triangulation",
         "Partitioning: splitting patches until every Monge fit converges",
         "Trimming and integrating curvature"}
check(bool(phases) and set(phases) <= KNOWN,
      f"every phase reported is a real one: {phases}")


class _FakeLog:
    """A ShapeJob whose log says whatever we want, so phase() can be checked
    against each marker state instead of against the clock."""

    def __init__(self, text):
        self.text = text

    def _log_text(self):
        return self.text


PHASE_CASES = [
    ("", "Remeshing to a uniform triangulation"),
    ("[input] loaded",
     "Partitioning: splitting patches until every Monge fit converges"),
    ("[input] loaded\n[done] ok", "Trimming and integrating curvature"),
]
for text, want in PHASE_CASES:
    got = ShapeJob.phase(_FakeLog(text))
    check(got == want, f"log {text.splitlines()[-1] if text else '(empty)'!r} "
                       f"-> {got!r}")

s = job.summary
print(f"  delta K = {s['k_fluct']:.6g}   A = {s['A_mm2']:,.0f} mm2   "
      f"patches = {s['n_patches']}   faces = {s['faces']:,}")
# The value moved on 2026-08-19 when two partition defects were fixed in
# camp/campsrc.py, and it moved because the defects were real, not because the
# measurement drifted. Attribution on this exact surface at these exact settings:
#
#     both fixes off (the old golden)  0.01598563576
#     connectivity split only          0.00954629042
#     fan gate only                    0.01135047087
#     both on (this golden)            0.00937026880
#
# Setting fan_thresh=None and split_disconnected=False in shape_metrics.PARAMS
# reproduces the old number exactly, which is the check that this is a stale
# golden rather than a regression. Re-derive it that way before changing it again.
check(abs(s["k_fluct"] - 0.00937026880) < 1e-9,
      "delta K matches the value shape_metrics produces directly")
check(s["euler_number"] == -3, "remeshing did not change the topology")
check(s["campsrc"] == "vendored (camp/)", "the vendored CAMP answered")

print()
print("=" * 78)
print("the arrays that colour the patches")
print("=" * 78)
verts, faces = job.polydata_arrays()
n_f = len(faces)
check(len(job.arrays["face_labels"]) == n_f, "one label per face")
check(int(faces.max()) < len(verts), "every face indexes a vertex that exists")

# The colour bar's range is the only thing with choices in it now: one metric,
# three ways to scale it. Robust must be tighter than Full or it is not doing
# anything, and Manual must be obeyed exactly or the bar lies about itself.
full = job.face_values(scale="full")
rob = job.face_values(scale="robust")
tight = job.face_values(scale="robust", clip_pct=15.0)
man = job.face_values(scale="manual", lo=-0.25, hi=0.25)
for tag, fv in (("full", full), ("robust 5%", rob), ("robust 15%", tight),
                ("manual", man)):
    lo, hi = fv["clim"]
    sat = float((fv["values"] < lo).mean() + (fv["values"] > hi).mean())
    print(f"  {tag:<11} clim=({lo:10.4g}, {hi:9.4g})  saturating {100*sat:5.1f}%  "
          f"clipped={fv['clipped']}")
    check(len(fv["values"]) == n_f and np.isfinite(fv["values"]).all(),
          f"{tag}: one finite value per face")
    check(hi > lo, f"{tag}: the bar spans something")

check(full["range"] == (float(job.arrays["integrated_K"].min()),
                        float(job.arrays["integrated_K"].max())),
      "the full range is the actual data range")
check(not full["clipped"], "full is reported as unclipped")
check(rob["clipped"] and rob["clim"][0] > full["clim"][0]
      and rob["clim"][1] < full["clim"][1],
      "robust is strictly inside full, and says it is clipped")
check(tight["clim"][1] - tight["clim"][0] < rob["clim"][1] - rob["clim"][0],
      "a bigger clip gives a tighter bar")
check(man["clim"] == (-0.25, 0.25), "manual bounds are used exactly")
check(job.face_values(scale="manual", lo=1.0, hi=1.0)["clim"] != (1.0, 1.0),
      "an inverted or empty manual range falls back rather than drawing nothing")
check(fv["title"] == "Integrated K" and fv["cmap"] == "plasma",
      "the colouring is integrated K on CAMP's plasma map")

mk = job.mean_K()
w = job.arrays["patch_area"][job.arrays["keep"]]
q = job.arrays["integrated_K"][job.arrays["keep"]]
by_hand = float(np.sum(q * w / w.sum()))
check(abs(mk - by_hand) < 1e-12, f"mean K is area weighted over kept patches ({mk:.5g})")
# delta K is the variance about exactly that mean, which is the whole reason the
# deviation view is worth drawing.
var = float(np.sum((w / w.sum()) * (q - mk) ** 2))
check(abs(var - s["k_fluct"]) < 1e-12, "delta K is the variance about that mean")

keep_f = job.kept_face_mask()
check(len(keep_f) == n_f and keep_f.sum() < n_f,
      f"the trim excludes some faces ({n_f - int(keep_f.sum()):,} of {n_f:,})")

edges = job.patch_edges()
check(edges.ndim == 2 and edges.shape[1] == 2 and len(edges) > 0,
      f"patch outlines built ({len(edges):,} segments)")
check(int(edges.max()) < len(verts), "every outline segment indexes a real vertex")
check(job.patch_edges() is edges, "outlines are built once, not per frame")

print()
print("=" * 78)
print("what the panel says about comparability")
print("=" * 78)
notes = job.notes()
for level, text in notes:
    print(f"  [{level}] {text}")
check(any("vertex floor" in t for _, t in notes),
      "88.9% at the vertex floor is flagged as outside the reference band")
check(any("50,000 faces" in t for _, t in notes),
      "running below the reference's 100,000 faces is flagged")
# The partition note comes first because it is the only one that changes how every
# other number on the panel should be read: a delta K measures the partition as well
# as the anatomy, and a difference smaller than the band is not a difference.
#
# It replaced a "Repeatability" note on 2026-08-20. That note quoted the REPEATABILITY
# table, which was measured before the connectivity split and the fan gate and so
# mixed partition noise with defects since fixed -- it warned of a 2.6x band at trim 0
# where the pipeline now carries about 16% for one partition and 6% for eight.
check(notes and ("ONE partition" in notes[0][1]
                 or "independent partitions" in notes[0][1]),
      "the partition-uncertainty note leads, before any comparability warning")
check(any("fixed mesh density" in t for _, t in notes),
      "the panel says delta K is only comparable at a fixed face count")
(med, worst), at = job.repeatability()
check(at == 0.01 and med < 1.5,
      f"the band is quoted at the trim actually used ({at:g}: {med}x, {worst}x)")
# KY160_1 carries 8 folded edges inside 5 mm. That is a segmentation defect, and it
# has to be said before someone reads a number off the surface, not after.
check(any("folded edge" in t for _, t in notes),
      "folds in the source surface are reported")
check(not any("fold" in t and "correct" in t.lower() and "NOT" not in t
              for _, t in notes),
      "and are never presented as a correction applied to the value")
check(any("both tails" in t.lower() for _, t in notes),
      "the trim is stated as taking both tails, not one")
print(f"  headline: {job.headline()}")

print()
print("=" * 78)
print("export")
print("=" * 78)
out = Path(tempfile.mkdtemp(prefix="vcshape_test_"))
table, side = job.export(out, "check")
head = table.read_text(encoding="utf-8").splitlines()
print(f"  {table.name}: {len(head) - 1} rows")
print(f"  {head[0]}")
check(len(head) - 1 == s["n_patches"], "one CSV row per patch")
check("integrated_K" in head[0] and "kept_after_trim" in head[0],
      "the table carries the per-patch integrated K and the trim flag")
check("target_faces" in side.read_text(encoding="utf-8"),
      "the sidecar records the conditions the number is only valid under")
shutil.rmtree(out, ignore_errors=True)

print()
print("=" * 78)
print("a failure reaches the panel instead of hanging it")
print("=" * 78)
bad = ShapeJob()
bad.adopt_path(SURFACE)
bad.source_euler = 99            # the worker asserts the topology it was handed
bad.start(25_000, 0.01)
st, _ = wait(bad, limit=180)
print(f"  state={st['state']}")
print(f"  {st['message']}")
check(st["state"] == "failed", "a worker that refuses the job reports failed")
check("euler" in st["message"], "the reason is in the message, not only in the log")
check(not bad.done, "no stale result is left behind for the view to draw")

print()
print("=" * 78)
print("straight from a cut in the same session")
print("=" * 78)
from vascular_api import CutSession                         # noqa: E402

cut = CutSession(str(ROOT / "KY160_1_CentrelinePoints.csv"),
                 str(ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"))
cut.recut()
hand = ShapeJob()
hand.adopt_mesh(cut.mesh, "KY160_1_open", expect_euler=int(cut.mesh.euler_number))
print(f"  {hand.source_info()}")
check(hand.source_kind == "cut", "the handoff is marked as coming from the cut")
check(hand.source_euler == int(cut.mesh.euler_number),
      "the exported STL welds back to the topology the cut reported")

try:
    ShapeJob().adopt_mesh(cut.mesh, "x", expect_euler=0)
    check(False, "a topology change during export is caught")
except ShapeError as e:
    check("euler" in str(e), f"a topology change during export is caught: {e}")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
