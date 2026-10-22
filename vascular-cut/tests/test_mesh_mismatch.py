"""A mesh that does not match the centerline must be obvious BEFORE the anatomy is.

preflight already refuses to cut on this and explains it well. But it runs at cut
time, and the family call, every vessel name and every region preset are computed
from the same radii long before that -- so the wrong anatomy is on screen for as long
as the user has been looking at it.

Measured on the thoracic cohort: pairing a branched centerline with the trunk-only
mesh sitting beside it put 4 to 5 of 6 terminals OUTSIDE the surface, one at -68 mm.
outflow_split then finds no dividing junction, the family quietly falls back to the
2nd/1st radius rule, and KK40_3, KR47_1 and KR52_1 were all called abdominal. All
three read correctly with the matching mesh.

Reproduced here on the repository's own KY160, by pairing its branched centerline
with a mesh of a different extent -- no cohort access needed.

    python tests/test_mesh_mismatch.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import trimesh                                                      # noqa: E402
from branch_tree import load_tree                                   # noqa: E402
from branch_cut import (load_mat_mesh, classify_terminals,          # noqa: E402
                        preflight, name_anatomy,
                        MIN_CENTERLINE_COVERAGE)

CSV = ROOT / "KY160_1_CentrelinePoints.csv"
MAT = ROOT / "Geometries & Meshes/mesh/KY160_1.mat"

fails = []


def check(cond, what):
    print(f"  {'ok  ' if cond else 'FAIL'}  {what}")
    if not cond:
        fails.append(what)


print("=" * 78)
print("the matching pair: nothing is flagged")
print("=" * 78)
tree = load_tree(str(CSV))
mesh = load_mat_mesh(str(MAT))
classify_terminals(tree, mesh, verbose=False)
outside = list(getattr(tree, "radius_outside", None) or [])
n_trunk = sum(1 for t in tree.terminals if getattr(t, "is_trunk", False))
print(f"  {len(tree.terminals)} terminals, {n_trunk} trunk end(s), "
      f"{len(outside)} outside the mesh")
check(not outside, "no centerline runs outside the mesh it belongs to")
check(getattr(tree, "family_trusted", None) is True,
      "and the family call is marked trustworthy")
check(n_trunk == 2, "KY160 is read as thoracic")

print()
print("=" * 78)
print("a mesh of the wrong extent: flagged, and the call marked untrustworthy")
print("=" * 78)
# Keep only the faces in the lower half of z. That is what a trunk-only export looks
# like to this code: a real surface, in the right frame and units, that simply does
# not contain the branches the centerline describes.
half = mesh.copy()
zmid = float(np.median(half.vertices[:, 2]))
keep = half.triangles_center[:, 2] < zmid
half.update_faces(keep)
half.remove_unreferenced_vertices()
print(f"  cropped mesh: {len(half.faces):,} of {len(mesh.faces):,} faces")

t2 = load_tree(str(CSV))
classify_terminals(t2, half, verbose=False)
out2 = list(getattr(t2, "radius_outside", None) or [])
rad = {t.label: float(getattr(t.segment, "_radius_mm", None) or 0.0)
       for t in t2.terminals}
print(f"  radii: {[f'{v:.1f}' for v in rad.values()]}")
print(f"  coverage: {100 * t2.coverage:.1f}%   outside: {out2}")

# Coverage is the check that carries here. Cropping leaves an OPEN surface, and on an
# open surface segment_radius falls back to unsigned distance and returns cheerfully
# positive nonsense -- 79 mm on a 22 mm aorta, flagged by nothing. That is the real
# second failure mode, not an artefact of the test.
check(t2.coverage < MIN_CENTERLINE_COVERAGE,
      f"the centerline is measured as outside this mesh ({100 * t2.coverage:.0f}%)")
check(t2.mesh_matches is False, "so the pair is not accepted as the same case")
check(getattr(t2, "family_trusted", True) is False,
      "and the family call is marked NOT trustworthy")
check(max(rad.values()) > 3 * max(v for v in
      {t.label: float(getattr(t.segment, "_radius_mm", None) or 0.0)
       for t in tree.terminals}.values()),
      "the unsigned fallback really does return nonsense here, positive and huge")

print()
print("=" * 78)
print("naming still runs, but the caller has been told not to believe it")
print("=" * 78)
names, _v, _n = name_anatomy(t2, half)
print(f"  {len(names)} opening(s) named on a mesh that does not match")
check(True, "name_anatomy does not crash on a mismatched pair")

print()
print("=" * 78)
print("the desktop caption leads with the mismatch")
print("=" * 78)
from qtpy import QtWidgets                                          # noqa: E402
app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
from branch_gui import CutWindow                                    # noqa: E402

w = CutWindow(None, None)
w.tree, w.base = t2, half
w.anat_names, w.anat_vessels, w.anat_notes = names, {}, []
w._update_anat_label()
text = w.anat_lbl.text()
print(f"  {text[:110]}")
check("do not look like the same export" in text,
      "the caption says the two are not the same export")
check("read from the anatomy" in text,
      "and that nothing shown was read from the anatomy")

w.tree, w.base = tree, mesh
w.anat_names, w.anat_vessels, w.anat_notes = name_anatomy(tree, mesh)
w._update_anat_label()
good = w.anat_lbl.text()
print(f"  matching pair: {good[:90]}")
check("do not look like the same export" not in good,
      "and a matching pair gets the ordinary caption back")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
