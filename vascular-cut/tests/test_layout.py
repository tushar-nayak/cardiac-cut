"""The new panel: does it fit a small screen, and is everything reachable?

argv[1] = jk5|ky160
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from qtpy import QtWidgets, QtCore
from branch_gui import CutWindow, Section
from branch_cut import region_presets

MODELS = {"jk5": ("JK5_1_CentrelinePoints_branches.csv", "JK5_1_M10.mat"),
          "ky160": ("KY160_1_CentrelinePoints.csv",
                    "Geometries & Meshes/mesh/KY160_1.mat")}
which = sys.argv[1] if len(sys.argv) > 1 else "jk5"
csv, mat = MODELS[which]

app = QtWidgets.QApplication(sys.argv[:1])
# a fresh settings scope, so a previous run's fold state cannot skew the numbers
QtCore.QSettings("branches", "cuttool_test_layout").clear()
fails = []
w = CutWindow(str(ROOT / csv), str(ROOT / mat))
w.show()
for _ in range(15):
    app.processEvents()

secs = [("Input", w.sec_input), ("Centerline", w.sec_cl),
        ("View & measure", w.sec_view), ("Output", w.sec_out)]
print("=" * 78)
print(f"{which.upper()} ,  new collapsible panel")
print("=" * 78)
print("\n  sections after load (setup should be folded):")
for name, s in secs:
    print(f"    {name:16s} open={s.header.isChecked()!s:5s}  header={s.header.text()!r}")
folded = [n for n, s in secs if n != "View & measure" and s.header.isChecked()]
if folded:
    fails.append(("not-folded", folded))
# a folded section must still say what it is set to
for name, s in secs:
    if not s.header.isChecked() and s.header.text().strip() == name:
        fails.append(("no-summary", name))

print(f"\n  openings box height: {w.openings_box.height()} px "
      f"(was 54 in the old layout)")
if w.openings_box.height() < 200:
    fails.append(("openings-cramped", w.openings_box.height()))

# how short can the window go now?
w.resize(1200, 200)
for _ in range(10):
    app.processEvents()
floor_h = w.height()
print(f"  window minimum height: {floor_h} px  (was 851, 941 with a region)")
for target, label in ((728, "1366x768 laptop"), (672, "1080p at 150% scaling")):
    ok = floor_h <= target
    print(f"    fits {label:24s} ({target} px): {'yes' if ok else 'NO'}")
    if not ok:
        fails.append((f"too-tall-for-{label}", floor_h))

# with the window that short, is the footer still on screen and is the rest scrollable?
for label, wid in (("status", w.status), ("Export", w.export_btn), ("Reset", w.reset_btn)):
    y = wid.mapTo(w, wid.rect().topLeft()).y()
    on = 0 <= y and y + wid.height() <= w.height() and wid.isVisible()
    print(f"    footer {label:8s} y={y}..{y + wid.height()} on-screen={on}")
    if not on:
        fails.append(("footer-off-screen", label))
sb = w.panel_scroll.verticalScrollBar()
print(f"    panel scrolls instead of clipping: range 0..{sb.maximum()} "
      f"({'scrollable' if sb.maximum() > 0 else 'fits without scrolling'})")

# region on: the banner must be reachable at this tiny size
presets = region_presets(w.tree, w.base, w.anat_names, w.anat_vessels)
if presets:
    w._set_region((presets[0]["s0"], presets[0]["s1"]), presets[0]["name"])
    for _ in range(10):
        app.processEvents()
    print(f"\n  region on at {w.width()}x{w.height()}: "
          f"window min now {w.minimumSizeHint().height()} px")
    banner = w.openings_lay.itemAt(0).widget()
    btn = banner.findChild(QtWidgets.QPushButton) if banner else None
    print(f"    banner button present: {btn.text() if btn else None}")
    if btn is None:
        fails.append("no-banner-small")
    else:
        # reachable = on screen, or scrollable into view
        w.panel_scroll.ensureWidgetVisible(btn)
        for _ in range(6):
            app.processEvents()
        y = btn.mapTo(w, btn.rect().topLeft()).y()
        on = 0 <= y and y + btn.height() <= w.height()
        print(f"    after scrolling to it: y={y}..{y + btn.height()}, on-screen={on}")
        if not on:
            fails.append("banner-unreachable-small")
    w._set_region(None)

# the splitter must actually move
w.resize(1400, 900)
for _ in range(8):
    app.processEvents()
before = w.split.sizes()[0]
w.split.setSizes([620, 780])
for _ in range(8):
    app.processEvents()
after = w.split.sizes()[0]
print(f"\n  splitter: panel {before} -> {after} px "
      f"({'resizable' if after != before else 'STUCK'})")
if after == before:
    fails.append("splitter-stuck")

# reopening a folded section must restore its controls
w.sec_input.set_open(True)
for _ in range(6):
    app.processEvents()
print(f"  reopened Input: body visible={w.sec_input.body.isVisible()}, "
      f"Load button visible={w.load_btn.isVisible()}")
if not w.load_btn.isVisible():
    fails.append("reopen-broken")

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
sys.stdout.flush()

# Report the verdict to the SHELL, not only to a reader. Without this the script
# exits 0 while printing FAILURES, so every automated runner reads it as green -- which
# is how a stale colour-label assertion sat broken without anyone noticing.
raise SystemExit(1 if fails else 0)
