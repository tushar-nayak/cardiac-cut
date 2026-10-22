"""The nub-in cut, made by the batch and then restored into the GUI, must agree.

The workbench opens each case in the desktop tool showing what the preset produced,
so that a person starts from the batch cut rather than from a default offset on every
opening. If the restore is even slightly off, the geometry someone signs off in the
window is not the geometry the batch scored, and nothing downstream would say so.

Self-contained: cuts the repository's own KY160 rather than a cohort case, so it runs
anywhere the other suites do.
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets                                         # noqa: E402
import branch_cut as BC                                            # noqa: E402
import nub_cut as NC                                               # noqa: E402
import nub_session                                                 # noqa: E402
from branch_gui import CutWindow                                   # noqa: E402
from branch_tree import load_tree                                  # noqa: E402

CSV = ROOT / "KY160_1_CentrelinePoints.csv"
MAT = ROOT / "Geometries & Meshes/mesh/KY160_1.mat"

app = QtWidgets.QApplication(sys.argv)
fails = []

print("=" * 78)
print("BATCH: cut the nubs headlessly")
print("=" * 78)
mesh = BC.load_mat_mesh(str(MAT))
tree = load_tree(str(CSV))
BC.classify_terminals(tree, mesh, verbose=False)
_tn, vnames, _ = BC.name_anatomy(tree, mesh)
for seg in tree.segments.values():
    if getattr(seg, "_radius_mm", None) is None:
        seg._radius_mm = BC.segment_radius(mesh, seg)

out, plan, notes, info = NC.cut_nubs(mesh, tree, k=NC.NUB_K, anat=vnames)
cuts = []
for t in info["terminals"]:
    if t.label not in info["offsets"]:
        continue
    cuts.append({"label": t.label, "segment": t.segment.name, "name": t.name,
                 "offset_mm": float(info["offsets"][t.label]), "mode": t.mode,
                 "interior": bool(getattr(t, "is_interior", False)),
                 "trunk": bool(getattr(t, "is_trunk", False)),
                 "subtree": sorted(getattr(t, "subtree_segments", None) or [])})
print(f"  {len(out.faces)} faces, {len(BC.boundary_loops(out))} openings, "
      f"euler {int(out.euler_number)}")
for c in cuts:
    print(f"    {c['name'][:34]:36s} {c['offset_mm']:7.2f} mm  mode={c['mode']:6s}"
          f"{'  TRIM' if c['interior'] else ''}")

side = Path(tempfile.mkdtemp()) / "ky160_nub.json"
side.write_text(json.dumps({"case": "local/KY160_1", "k": NC.NUB_K, "cuts": cuts}))

print()
print("=" * 78)
print("GUI: open the case and restore that cut")
print("=" * 78)
w = CutWindow(str(CSV), str(MAT))
w.show()
for _ in range(10):
    app.processEvents()
print(f"  before restore: {len(w.rows)} rows, {len(w.trims)} trims, "
      f"{len(BC.boundary_loops(w.cut))} openings")

for note in nub_session.apply_from_file(w, side):
    print(f"    {note}")
for _ in range(10):
    app.processEvents()

print(f"  after restore:  {len(w.rows)} rows, {len(w.trims)} trims, "
      f"{len(BC.boundary_loops(w.cut))} openings")
print(f"  status: {w.status.text()[:70]}")

# -- the checks -------------------------------------------------------------
want_trims = sorted(c["segment"] for c in cuts if c["interior"])
got_trims = sorted(tr["cand"]["segment"] for tr in w.trims)
print(f"\n  trims wanted {want_trims}  got {got_trims}")
if want_trims != got_trims:
    fails.append(("trims", f"wanted {want_trims}, got {got_trims}"))

by_seg = {}
for c in cuts:
    by_seg.setdefault((c["segment"], c["interior"]), c)
print("  offsets:")
for r in w.rows:
    t = r.terminal
    c = by_seg.get((t.segment.name, bool(getattr(t, "is_interior", False))))
    if c is None:
        fails.append(("unmatched-row", t.name))
        print(f"    {t.name[:34]:36s} NOT IN THE BATCH CUT")
        continue
    d = abs(float(r.value()) - c["offset_mm"])
    print(f"    {t.name[:34]:36s} gui {r.value():7.2f}  batch {c['offset_mm']:7.2f}"
          f"  diff {d:.3f}  mode {t.mode}")
    if d > 0.05:
        fails.append(("offset", f"{t.name}: gui {r.value():.2f} vs {c['offset_mm']:.2f}"))
    if t.mode != c["mode"]:
        fails.append(("mode", f"{t.name}: gui {t.mode} vs {c['mode']}"))

print("\n  geometry:")
print(f"    batch {len(out.faces)} faces / {len(BC.boundary_loops(out))} openings"
      f"   gui {len(w.cut.faces)} faces / {len(BC.boundary_loops(w.cut))} openings")
if len(w.cut.faces) != len(out.faces):
    fails.append(("faces", f"batch {len(out.faces)}, gui {len(w.cut.faces)}"))
if len(BC.boundary_loops(w.cut)) != len(BC.boundary_loops(out)):
    fails.append(("openings", "count differs"))
if w.problem is not None:
    fails.append(("gui-invalid", str(w.problem)[:70]))

# Same faces, not merely the same number of them.
a = np.sort(np.asarray(out.triangles_center).round(4), axis=0)
b = np.sort(np.asarray(w.cut.triangles_center).round(4), axis=0)
same = a.shape == b.shape and bool(np.allclose(a, b, atol=1e-3))
print(f"    identical face set: {same}")
if not same:
    fails.append(("geometry", "the restored cut is not the batch cut"))

print()
print("=" * 78)
if fails:
    print(f"{len(fails)} FAILURE(S)")
    for k, v in fails:
        print(f"  {k}: {v}")
else:
    print("ALL PASS")
print("=" * 78)
