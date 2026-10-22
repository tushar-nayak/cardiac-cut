"""Hosted: "Apply to this case" rebuilds an upload with a corrected anatomy call.

Every upload is read with the anatomy detected (_load_uploads forces auto), so on the
hosted site a thoracic aorta misread as abdominal could never reach nub-in or the
Ishimaru zones: changing the Anatomy box did nothing until the next load, and there
was no next load. reload_case rebuilds the case from the STORED upload with the boxes
as they are now. It must keep the submission (the ID the visitor noted is still the
one their files are under), never read a path from state, refuse when there is
nothing to rebuild, and drop a rebuild whose visitor has left.

This drives the real app in public mode with no browser, as test_visit_wipe does.

    .venv\\Scripts\\python tests\\test_reload_case.py
"""
import asyncio
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT))

# Public mode needs a store. A new folder each run, under the gitignored tests/_out;
# this run stores one small submission there (JK5, ~2.3 MB).
STORE_DIR = ROOT / "tests" / "_out" / f"reload_case_{time.strftime('%Y%m%d_%H%M%S')}"
assert not STORE_DIR.exists(), STORE_DIR
os.environ["VC_PUBLIC"] = "1"
os.environ["VC_STORE_DIR"] = str(STORE_DIR)

import cut_app as APP                                              # noqa: E402

state = APP.state
CSV = ROOT / "JK5_1_CentrelinePoints_branches.csv"
MAT = ROOT / "JK5_1_M10.mat"
fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        fails.append(name)


async def until_idle(limit=300.0):
    t0 = time.time()
    while state.busy and time.time() - t0 < limit:
        await asyncio.sleep(0.1)


print("=" * 78)
print("1. an abdominal upload, then Apply with Anatomy = Thoracic")
print("=" * 78)
APP._admit_visit()
surf = APP.STORE.check_upload("surface", MAT.name, MAT.read_bytes())
cl = APP.STORE.check_upload("centerline", CSV.name, CSV.read_bytes())
APP.LOAD["seq"] += 1
asyncio.run(APP._load_uploads(APP._gen(), APP.LOAD["seq"], surf, cl))
sub = APP.SESSION["sub"]
check("control: the upload loaded, with a submission", APP.SESSION["s"] is not None
      and state.loaded and sub and state.submission_id == sub, state.status_text)
check("control: read as abdominal", state.anat_kind == "abdominal"
      and not state.is_thoracic and state.loaded_family == "auto")
check("the case is named after the upload", state.case_name == "JK5_1_M10",
      state.case_name)
seq0, s0 = state.load_seq, APP.SESSION["s"]

state.family = "thoracic"
state.mat_path = r"C:\somewhere\else.mat"      # a path in state is never read


async def apply():
    APP.reload_case()
    check("  (control: the rebuild runs off the event loop)", state.busy)
    await until_idle()

asyncio.run(apply())
check("the case was rebuilt", APP.SESSION["s"] is not s0 and state.loaded
      and not state.status_bad, state.status_text)
check("now thoracic, set by hand", state.is_thoracic is True
      and state.anat_kind == "thoracic" and state.anat_how == "set by hand",
      f"{state.anat_kind} / {state.anat_how}")
check("nub-in is offered", state.can_nub is True)
check("the Apply offer is satisfied", state.loaded_family == "thoracic")
check("the submission is kept", APP.SESSION["sub"] == sub
      and state.submission_id == sub)
check("so is the case name", state.case_name == "JK5_1_M10", state.case_name)
check("the load counter moved (the page reopens Cuts)", state.load_seq == seq0 + 1,
      f"{seq0} -> {state.load_seq}")

print("\n" + "=" * 78)
print("2. a rebuild whose visitor left is dropped")
print("=" * 78)
state.family = "abdominal"
s1 = APP.SESSION["s"]


async def leave_during(delay):
    APP.reload_case()
    await asyncio.sleep(delay)
    check("  (control: still building when they left)", state.busy)
    APP.GATE.gen += 1
    t0 = time.time()
    while APP.SESSION["s"] is s1 and time.time() - t0 < 30:
        await asyncio.sleep(0.2)
        if all(t.done() for t in asyncio.all_tasks() if t is not asyncio.current_task()):
            break

asyncio.run(leave_during(0.3))
check("the old visitor's rebuild is never installed", APP.SESSION["s"] is s1
      and state.is_thoracic is True)
APP._end_visit()

print("\n" + "=" * 78)
print("3. nothing to rebuild")
print("=" * 78)
APP._admit_visit()
APP.reload_case()
check("refused, saying so", state.status_bad
      and "nothing to read again" in state.status_text.lower(), state.status_text)
check("  ...and nothing started", not state.busy)
APP.SESSION["sub"] = "20260101T000000Z-" + "0" * 32      # a well-formed id, not stored
APP.SESSION["s"] = object()
APP.reload_case()
check("an upload that cannot be read again is refused",
      state.status_bad and "could not be read again" in state.status_text,
      state.status_text)
APP._end_visit()

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
