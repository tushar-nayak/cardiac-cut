"""When a visitor leaves, the hosted app must hold nothing of theirs.

The web app in public mode admits one visitor at a time (webapp/visit_gate.py) and
wipes everything the moment they leave, so the next person cannot see their case,
their submission id or their files. This drives that wipe on the REAL app, with no
browser: import cut_app in public mode, load a real case, dirty every place a visitor
can leave a trace, end the visit, and look at each of those places again.

The browser-level proof -- a second browser refused, the next visitor's page clean,
a load that finishes after its visitor left -- is tests/web/check_privacy.py. This is
the fast one, and it looks where a browser cannot: module globals, the scene graph,
the camera, the render buffer.

WHAT IT PINS
  1. A fresh public process has no upload token, and its state is the blank page.
  2. Admitting a visitor mints a token; ending the visit withdraws it, and the old
     one is refused -- including the trap where "no token" equalled "no header".
  3. After the visit ends: every state key equals a fresh process's, and no trace of
     the visitor (their case name, typed export name, submission id) is anywhere in
     state; the session, pending uploads, download links, zone map, delta K job and
     framing memory are empty; the scene holds no surface, pyvista holds no dataset,
     the camera is back where it started, and the rendered frame is blank.
  4. Work that finishes after its visitor left is dropped: a session still building,
     and a zone map still being computed.
  5. The next visitor gets a new token.

    .venv\\Scripts\\python tests\\test_visit_wipe.py
"""
import asyncio
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT))

# Public mode needs a store. A new folder each run, under the gitignored tests/_out;
# this run stores one small submission there (KY160, ~2 MB).
STORE_DIR = ROOT / "tests" / "_out" / f"visit_wipe_{time.strftime('%Y%m%d_%H%M%S')}"
assert not STORE_DIR.exists(), STORE_DIR
os.environ["VC_PUBLIC"] = "1"
os.environ["VC_STORE_DIR"] = str(STORE_DIR)

import cut_app as APP                                              # noqa: E402

state = APP.state
CSV = ROOT / "KY160_1_CentrelinePoints.csv"
MAT = ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"
fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def state_diff():
    """Keys whose value is not a fresh process's."""
    return [k for k, v in APP.VISIT_DEFAULTS.items() if state[k] != v]


def ink():
    """Fraction of the rendered frame that is not background, outside the corner the
    L/P/S axes live in (bottom-left). 0 for an empty scene; several % with a model."""
    SCENE = APP.SCENE
    SCENE.render()
    img = np.asarray(SCENE.plotter.screenshot(return_img=True))[..., :3].astype(int)
    h, w = img.shape[:2]
    bg = img[2, w - 3]                                   # top-right: never drawn on
    lit = np.abs(img - bg).max(axis=2) > 24
    lit[int(0.70 * h):, :int(0.30 * w)] = False          # the axes triad
    return float(lit.mean())


print("=" * 78)
print("1. a fresh public process")
print("=" * 78)
check("the gate is on in public mode", APP.GATE is not None)
check("no upload token before anyone is admitted",
      APP.HTTP_TOKEN is None and state.vc_token == "")
check("state is the blank page", not state_diff(), f"{state_diff()[:6]}")

print("\n" + "=" * 78)
print("2. admit: a token for this visit")
print("=" * 78)
APP._admit_visit()
tok = APP.HTTP_TOKEN
check("admitting mints a token and hands it to the page",
      bool(tok) and state.vc_token == tok)
check("it is accepted", APP._token_ok(tok))
check("a wrong one is not", not APP._token_ok("x" * len(tok)))
check("no header is not", not APP._token_ok(None) and not APP._token_ok(""))

print("\n" + "=" * 78)
print("3. a visitor leaves traces everywhere, then leaves")
print("=" * 78)
APP._open_session(str(CSV), str(MAT))
s = APP.SESSION["s"]
check("control: the case loaded", s is not None and state.loaded and state.openings,
      state.status_text)
SECRET = "wfSECRETname"
SUB = "20260925T000000Z-" + "ab" * 16
state.export_name = SECRET
state.submission_id = SUB
APP.SESSION["sub"] = SUB
state.marks = [{"label": "Mark 1", "mm": 12.3, "vessel": "Left subclavian"}]
state.marks_note = "12.3 mm from the Left subclavian node"
state.up_mesh_note = "KY160_1.mat  (1.7 MB, 37,786 faces)"
state.status_text = "Loading KY160_1.mat"
checked = APP.STORE.check_upload("surface", "KY160_1.mat", MAT.read_bytes())
APP.PENDING["surface"] = checked
APP.HTTP_DL["tok"] = (f"{SECRET}.zip", b"zip bytes", "application/zip")
APP.ZONES.update(rgb=np.zeros((10, 3)), legend=[{"zone": 0}], for_cut=(1, 2, True, None))
APP._FRAMED["src"] = f"/tmp/{SECRET}.stl"
# What a hostile page could have written, had it got past the gate's write filter:
# trame's own page template -- which the NEXT visitor's page is built from -- and a
# key no fresh page has. The wipe must undo both on its own.
PLANT = "wfPLANTED"
state["trame__template_main"] = f"<div>{PLANT}</div>"
state["wfPLANTkey"] = PLANT
APP.JOB.adopt_mesh(s.mesh, SECRET, expect_euler=int(s.mesh.euler_number))
APP.SCENE.plotter.camera.azimuth = 70
before = ink()
check("control: the loaded case is drawn", before > 0.02, f"{100 * before:.1f}% lit")
check("control: the camera moved", APP.SCENE.plotter.camera_position != APP.INITIAL_CAMERA)
check("control: the delta K job holds their surface",
      APP.JOB.source is not None and APP.JOB.source_mesh is not None)

APP._end_visit()

diff = state_diff()
check("every state key is a fresh process's", not diff, f"{diff[:8]}")
# Every key trame holds, pushed or pending. NOT state.to_dict(), which returns only
# what has been pushed to a page -- nothing, in a process that is not serving -- and
# made this check pass vacuously on its first version.
now = APP._state_now()
check("control: the whole state is being looked at",
      len(now) >= len(APP.VISIT_DEFAULTS) and "trame__template_main" in now, f"{len(now)}")
blob = json.dumps(now, default=str)
leaks = [n for n in (SECRET, SUB, "KY160", "subclavian", "37,786", PLANT) if n in blob]
check("nothing of theirs anywhere in state", not leaks, f"{leaks}")
full = APP._SNAPSHOT["full"]
check("control: the snapshot is the whole first state, trame's keys included",
      full is not None and len(full) >= len(APP.VISIT_DEFAULTS)
      and "trame__template_main" in full, f"{len(full or {})}")
bad = [k for k in full if state[k] != full[k]]
check("every key -- trame's own included -- is back to what the first visitor got",
      not bad, f"{bad[:6]}")
check("  ...the page template among them",
      state["trame__template_main"] == full["trame__template_main"])
check("a key no fresh page has is blanked", state["wfPLANTkey"] is None)
check("no token, and theirs is refused",
      APP.HTTP_TOKEN is None and state.vc_token == "" and not APP._token_ok(tok))
check("no session, no submission", APP.SESSION == {"s": None, "sub": None})
check("no pending uploads", APP.PENDING == {"surface": None, "centerline": None})
check("no download links", not APP.HTTP_DL)
check("no zone map",
      APP.ZONES["rgb"] is None and not APP.ZONES["legend"] and APP.ZONES["for_cut"] is None)
check("the delta K job holds nothing of theirs",
      APP.JOB.source is None and APP.JOB.source_mesh is None and APP.JOB.source_name == ""
      and APP.JOB.summary is None and APP.JOB.proc is None)
check("framing memory cleared", APP._FRAMED["src"] is None)
check("pyvista holds no dataset",
      APP.SCENE.plotter.mesh is None and APP.SCENE.plotter.mapper is None)
check("the camera is back where it started",
      APP.SCENE.plotter.camera_position == APP.INITIAL_CAMERA)
after = ink()
check("the rendered frame is blank", after < 0.0005, f"{100 * after:.3f}% lit")

print("\n" + "=" * 78)
print("4. work that finishes after its visitor left is dropped")
print("=" * 78)
APP._admit_visit()
surf = APP.STORE.check_upload("surface", "KY160_1.mat", MAT.read_bytes())
cl = APP.STORE.check_upload("centerline", CSV.name, CSV.read_bytes())


async def leave_during(coro, delay):
    """Run coro, and 'leave' -- change the visit -- while it is inside an await."""
    task = asyncio.ensure_future(coro)
    await asyncio.sleep(delay)
    check("  (control: it was still running when they left)", not task.done())
    APP.GATE.gen += 1                    # what the gate does when the visitor leaves
    await task

APP.LOAD["seq"] += 1
asyncio.run(leave_during(APP._load_uploads(APP._gen(), APP.LOAD["seq"], surf, cl), 0.4))
check("a session built after its visitor left is never installed",
      APP.SESSION["s"] is None and not state.loaded and not state.openings)

APP._admit_visit()


async def replaced_during(delay):
    """A second pair completes while the first is still building (same visitor)."""
    APP.LOAD["seq"] += 1
    task = asyncio.ensure_future(
        APP._load_uploads(APP._gen(), APP.LOAD["seq"], surf, cl))
    await asyncio.sleep(delay)
    check("  (control: the first load was still building when it was replaced)",
          not task.done())
    APP.LOAD["seq"] += 1                 # what a newer complete pair does
    await task

asyncio.run(replaced_during(0.4))
check("a load replaced by a newer one is never installed, and no submission is "
      "paired with it", APP.SESSION["s"] is None and APP.SESSION["sub"] is None)
APP.LOAD["seq"] += 1
asyncio.run(APP._load_uploads(APP._gen(), APP.LOAD["seq"], surf, cl))
check("the newest load installs, with its own submission beside it",
      APP.SESSION["s"] is not None and APP.SESSION["sub"]
      and APP.SESSION["sub"] == state.submission_id, state.status_text)
APP._end_visit()

APP._admit_visit()
APP._open_session(str(CSV), str(MAT))
check("control: a case to find zones on", APP.SESSION["s"] is not None
      and APP.SESSION["s"].ok, APP.SESSION["s"] and APP.SESSION["s"].problem)
state.zones_busy = True
asyncio.run(leave_during(APP._zones_job(APP._gen(), APP.SESSION["s"]), 0.05))
check("a zone map computed after its visitor left is never installed",
      APP.ZONES["rgb"] is None and not state.zone_legend)
APP._end_visit()

print("\n" + "=" * 78)
print("5. the next visitor")
print("=" * 78)
APP._admit_visit()
check("gets a new token", APP.HTTP_TOKEN and APP.HTTP_TOKEN != tok
      and state.vc_token == APP.HTTP_TOKEN)
check("and a blank page", state_diff() == ["vc_token"], f"{state_diff()[:6]}")
APP._end_visit()

print()
if fails:
    print(f"{len(fails)} FAILURE(S):")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
