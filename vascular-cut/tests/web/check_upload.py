"""Upload a surface and a centerline, download the cut, and look at what was kept.

The acceptance test for tools/HOSTING.md job 1 (uploads instead of paths) and job 7
(keep what is uploaded, under generated ids). It was written BEFORE either was built,
and run against the app as it stood, so that it is known to FAIL when the feature is
missing. A check that passes against code that does not exist measures nothing.

Unlike the checks in tests/web/, this one starts the app itself, because what it
tests only exists in one configuration: hosted mode (VC_PUBLIC=1), on loopback, with
its storage and its temp files pointed at folders this script can inspect afterwards.
Against a server somebody else started, the storage check would be looking in the
wrong place and would say nothing.

What it drives, which is therefore the contract the implementation has to meet:

  VC_STORE_DIR       the storage root, read from the environment. Every upload lands
                     under <root>/<generated id>/..., and the original filename is
                     recorded as metadata (a JSON file under that id), never as a path.
  file inputs        real <input type=file> elements in the cut section. Either one
                     labelled Mesh / Surface and one labelled Centerline (the visible
                     label, aria-label or accept list is enough to tell them apart),
                     or a single input with `multiple` that takes both.
  loading            either automatic once both files are in, or a button whose name
                     starts with Load or Upload. Whichever, the openings appear as
                     .op-card rows and the footer status line changes.
  download           a button or link whose name contains "Download". It produces a
                     browser download: a .zip holding the cut mesh and its sidecar, or
                     the mesh and the sidecar as two downloads.

The checks, each with its own verdict:

  2  SECURITY   no enabled text box that takes a filesystem path, in either section
  2b DISCLOSURE no server filesystem path reaching the browser at all: in an input's
                value (a DISABLED box still shows it), in the text, or on the websocket
  2c REFUSAL    acting as a hostile client -- server paths written into state from
                page JS, the Load / Load surface triggers fired directly -- nothing
                loads. The widget being disabled is not the control; this is
  3  UPLOAD     KY160_1.mat + KY160_1_CentrelinePoints.csv through the real file
                inputs; the openings and the status line are read off the page and
                compared with the same pair cut headlessly
  3L LARGE      the JK5 geometry as an ASCII STL over WSLINK_MAX_MSG_SIZE, under a
                filename carrying a fake MRN. Proves the websocket carries it, the
                bytes arrive intact, and the STL branch of the loader cuts it
  4  DOWNLOAD   captured by Playwright's download event, then checked without the
                app: trimesh loads it, its open boundary loops match the openings,
                and the sidecar parses and names every cut plane
  5  STORAGE    stored copies are byte-identical, sit under a generated id, the
                original filename is kept as metadata and appears in NO path written
                anywhere this script can see; the cut and sidecar are kept with them
  E  PAGE       no page errors, and no server-side traceback in the server's log

Run:
    $env:PYTHONIOENCODING="utf-8"
    .venv\\Scripts\\python tests\\web\\check_upload.py
    .venv\\Scripts\\python tests\\web\\check_upload.py --skip-large     # no JK5 files
    .venv\\Scripts\\python tests\\web\\check_upload.py --app <other app.py>

Every run writes into a new timestamped folder under <out root>/runs/: the server log,
the store the app was given, the downloads, screenshots. The out root is
$PROJECTS_SHARE/branches/_cache/check_upload unless --out-root is given, and the script
refuses to run without one; the app's TEMP is a new local folder under tests/_out/.
Nothing is overwritten and nothing is deleted. Exits 1 on any failure. Windows only
for now: it starts the app with .venv\\Scripts\\python.exe and stops it with taskkill.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path

import numpy as np
import trimesh
from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                                   # the repo
PY = ROOT / ".venv" / "Scripts" / "python.exe"
sys.path.insert(0, str(ROOT))

SMALL_MESH = ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"
SMALL_CL = ROOT / "KY160_1_CentrelinePoints.csv"
# The large case is JK5, not KY160 again: it has 11 openings to KY160's 5, so the page
# visibly changes when it loads. Re-uploading the same geometry would leave the first
# upload's rows on screen and a silent failure would read as a pass.
LARGE_SRC_MESH = ROOT / "JK5_1_M10.mat"
LARGE_SRC_CL = ROOT / "JK5_1_CentrelinePoints_branches.csv"
# Rebound in main() to <out root>/fixtures -- the large fixture is ~10 MB, and the
# project rules send regenerable heavy intermediates to the share, not C:.
FIXTURES = HERE / "fixtures"   # placeholder; main() rebinds to the share
SENTINEL = "MRN0012345"                  # a fake identifier, in the filename only
LARGE_MESH = FIXTURES / f"wf{SENTINEL}_JK5_1.stl"
LARGE_CL = FIXTURES / f"wf{SENTINEL}_JK5_1_centerline.csv"

STORE_ENV = "VC_STORE_DIR"
# Anything that could be a generated id: uuid4 hex, a canonical uuid, or a
# secrets.token_urlsafe(16) or longer. Checked separately for containing a name.
ID_RE = re.compile(r"^(?:[0-9a-f]{32}|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}"
                   r"-[0-9a-f]{12}|[A-Za-z0-9_-]{22,})$")
# Fragments of the uploaded filenames. Case IDs, not just full names: a stored
# "KY160_open.stl" leaks the identifier as surely as "KY160_1.mat" does.
NEEDLES = ["KY160", "CentrelinePoints", SENTINEL, "JK5_1"]
MESH_EXT = {".stl", ".ply", ".obj", ".off", ".vtp", ".vtk", ".glb", ".gltf"}

BUILD = re.search(r'^BUILD = "([^"]+)"', (ROOT / "webapp" / "cut_app.py")
                  .read_text(encoding="utf-8"), re.M).group(1)

# ── verdicts ─────────────────────────────────────────────────────────────────
CHECKS = {}          # id -> {"title", "notes": [(ok, what)], "skipped": reason}
ORDER = []
_cur = {"id": None}


def check(cid, title):
    CHECKS[cid] = {"title": title, "notes": []}
    ORDER.append(cid)
    _cur["id"] = cid
    print()
    print("=" * 74)
    print(f"[{cid}] {title}")
    print("=" * 74)


def note(ok, what):
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}")
    CHECKS[_cur["id"]]["notes"].append((bool(ok), what))
    return ok


def info(what):
    print(f"        {what}")


def verdict(cid):
    n = CHECKS[cid]["notes"]
    if not n:
        return "FAIL"            # a check that asserted nothing has not passed
    return "PASS" if all(ok for ok, _ in n) else "FAIL"


# ── helpers ──────────────────────────────────────────────────────────────────
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.5)
        return s.connect_ex(("127.0.0.1", port)) == 0


def listing(folder, recursive=False):
    if not folder.exists():
        return set()
    it = folder.rglob("*") if recursive else folder.iterdir()
    return {str(p.relative_to(folder)) for p in it}


def open_loops(mesh):
    """Open boundary loops, counted from the edges and nothing else."""
    e = np.sort(mesh.edges, axis=1)
    uniq, counts = np.unique(e, axis=0, return_counts=True)
    b = uniq[counts == 1]
    if len(b) == 0:
        return 0
    return len(trimesh.graph.connected_components(b, min_len=1))


def reference(csv, mesh):
    """The same pair cut headlessly, with the settings the page starts with."""
    from vascular_api import CutSession
    s = CutSession(str(csv), str(mesh), family="auto", smooth_mm=None,
                   clean_rims=False)
    st = s.status()
    return {"ok": st["ok"], "n_openings": st.get("n_openings"),
            "euler": st.get("euler_number"), "n_faces": st.get("n_faces"),
            "names": sorted(o["name"] for o in s.openings()),
            "message": st["message"]}


def make_large_fixtures():
    """JK5 as ASCII STL (about 9.6 MB) and its centerline, under a leaky name."""
    FIXTURES.mkdir(parents=True, exist_ok=True)
    if not LARGE_MESH.exists():
        from branch_cut import load_mat_mesh
        load_mat_mesh(LARGE_SRC_MESH).export(str(LARGE_MESH), file_type="stl_ascii")
    if not LARGE_CL.exists():
        shutil.copyfile(LARGE_SRC_CL, LARGE_CL)


def start_app(app, port, run, store, tmp):
    env = dict(os.environ)
    env.update({
        "VC_PUBLIC": "1", "VC_CUT_PORT": str(port), "VC_CUT_HOST": "127.0.0.1",
        STORE_ENV: str(store), "PYTHONIOENCODING": "utf-8",
        # Temp files land where this script can look at them afterwards. An upload
        # handler that spools to tempfile under the original name is a leak too.
        "TMP": str(tmp), "TEMP": str(tmp), "TMPDIR": str(tmp),
    })
    log = open(run / "server.log", "w", encoding="utf-8")
    proc = subprocess.Popen([str(PY), str(app)], cwd=str(ROOT), env=env,
                            stdout=log, stderr=subprocess.STDOUT)
    return proc, log


def stop_app(proc, port):
    """Kill the whole tree. The venv python.exe is a launcher whose CHILD is the real
    interpreter holding the port, so terminating the launcher alone leaves the server
    listening -- the 'two ids for one process' trap in tests/web/README.md."""
    if proc.poll() is None:
        subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                       capture_output=True)
    try:
        proc.wait(timeout=20)
    except subprocess.TimeoutExpired:
        proc.kill()
    for _ in range(20):
        if not port_in_use(port):
            return True
        time.sleep(0.5)
    return False


def wait_listening(proc, port, run, timeout=240):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if proc.poll() is not None:
            return False
        if port_in_use(port):
            return True
        time.sleep(1.0)
    return False


def open_panel(page, title):
    """Open a section by its title, checking its state rather than toggling blind."""
    panels = page.locator(".v-expansion-panel",
                          has=page.locator(".v-expansion-panel-title",
                                           has_text=re.compile(rf"^\s*{re.escape(title)}\s*$", re.I)))
    if panels.count() == 0:
        return False
    p = panels.first
    if "v-expansion-panel--active" not in (p.get_attribute("class") or ""):
        p.locator(".v-expansion-panel-title").first.click()
        page.wait_for_timeout(700)
    return True


def open_all_panels(page):
    titles = page.evaluate("""() => [...document.querySelectorAll(
        '.v-navigation-drawer .v-expansion-panel-title')].map(t => t.textContent.trim())""")
    for t in titles:
        open_panel(page, t)
    return titles


def set_section(page, which):
    tog = page.locator(".top-tog button", has_text=re.compile(which, re.I))
    if tog.count():
        tog.first.click()
        page.wait_for_timeout(1500)


INPUTS_JS = """() => [...document.querySelectorAll('input, textarea, [contenteditable]')]
  .map(el => {
    const fld = el.closest('.fld'), vin = el.closest('.v-input');
    const lab = [];
    if (el.labels) for (const l of el.labels) lab.push(l.textContent);
    if (fld) { const l = fld.querySelector('label'); if (l) lab.push(l.textContent); }
    if (vin) { const l = vin.querySelector('.v-label'); if (l) lab.push(l.textContent); }
    const panel = el.closest('.v-expansion-panel');
    const pt = panel ? panel.querySelector('.v-expansion-panel-title') : null;
    return {tag: el.tagName.toLowerCase(), type: (el.type || '').toLowerCase(),
            disabled: !!el.disabled || !!el.closest('fieldset[disabled]')
                      || !!(vin && vin.classList.contains('v-input--disabled')),
            readOnly: !!el.readOnly, placeholder: el.placeholder || '',
            value: (el.value || '').slice(0, 400),
            aria: el.getAttribute('aria-label') || '',
            label: lab.join(' | ').trim(), accept: el.accept || '',
            multiple: !!el.multiple,
            panel: pt ? pt.textContent.trim() : '',
            visible: !!(el.offsetWidth || el.offsetHeight || el.getClientRects().length)};
  })"""

PATHY_LABEL = re.compile(r"\b(mesh|cent(?:er|re)[ -]?line|surface|path|file|folder|director(y|ies))\b",
                         re.I)
PATHY_PLACEHOLDER = re.compile(r"\.(mat|stl|ply|obj|csv|vtp|vtk|vtu|json)\b", re.I)
PATHY_VALUE = re.compile(r"^[A-Za-z]:[\\/]|^/|^\\\\|[\\/]\.\.[\\/]|^\.\.[\\/]")


def path_inputs(rows):
    """Enabled, editable text boxes that look like they take a path."""
    out = []
    for r in rows:
        if r["tag"] == "input" and r["type"] not in {"text", "search", "url", ""}:
            continue
        if r["disabled"] or r["readOnly"]:
            continue
        blob = " ".join([r["label"], r["aria"], r["placeholder"]])
        if (PATHY_LABEL.search(blob) or PATHY_PLACEHOLDER.search(r["placeholder"])
                or PATHY_VALUE.search(r["value"])):
            out.append(r)
    return out


def page_state(page):
    names = page.evaluate("""() => [...document.querySelectorAll('.op-card .op-name input')]
        .map(e => e.value.trim())""")
    status = ""
    sb = page.locator(".statusbar")
    if sb.count():
        status = sb.first.inner_text().strip()
    cls = sb.first.get_attribute("class") if sb.count() else ""
    chips = page.evaluate("""() => [...document.querySelectorAll('.v-toolbar .v-chip')]
        .map(c => c.textContent.trim())""")
    return {"cards": page.locator(".op-card").count(), "names": names,
            "status": status, "status_bad": bool(re.search(r"\b(text|bg)-error\b", cls or "")),
            "chips": chips}


def chip_int(chips, pattern):
    for c in chips:
        m = re.search(pattern, c, re.I)
        if m:
            return int(m.group(1))
    return None


def classify_file_inputs(rows):
    """(mesh index, centerline index, multiple index) among the file inputs."""
    files = [r for r in rows if r["tag"] == "input" and r["type"] == "file"]
    mesh = cl = multi = None
    for i, r in enumerate(files):
        blob = " ".join([r["label"], r["aria"]])
        acc = r["accept"].lower()
        is_cl = bool(re.search(r"cent(?:er|re)[ -]?line", blob, re.I))
        is_mesh = bool(re.search(r"\b(mesh|surface|geometry)\b", blob, re.I))
        if not (is_cl or is_mesh):          # fall back on what it accepts
            is_mesh = any(x in acc for x in (".mat", ".stl", ".ply", ".obj"))
            is_cl = (not is_mesh) and any(x in acc for x in (".csv", ".vtp", ".json"))
        if is_mesh and not is_cl and mesh is None:
            mesh = i
        elif is_cl and not is_mesh and cl is None:
            cl = i
        if r["multiple"] and multi is None:
            multi = i
    return files, mesh, cl, multi


def do_upload(page, mesh_file, cl_file, expect, label, timeout_s=180):
    """Set the real file inputs and wait for the page to show the model. Returns the
    page state, or None if there was no way to upload at all."""
    rows = page.evaluate(INPUTS_JS)
    files, mi, ci, multi = classify_file_inputs(rows)
    note(len(files) > 0, f"{label}: an <input type=file> is on the page "
                         f"({len(files)} found)")
    if not files:
        return None
    for i, r in enumerate(files):
        info(f"file input {i}: label={r['label']!r} accept={r['accept']!r} "
             f"multiple={r['multiple']} disabled={r['disabled']}")
    loc = page.locator("input[type=file]")
    before = page_state(page)
    # The rows and the status line are text; the 3-D view is a picture pushed by the
    # server. The first positive run passed every text assertion with the PREVIOUS
    # case still in the viewport, so hash the viewport's pixels (not the canvas via
    # toDataURL, which reports "unchanged" for the remote view -- tests/web/README).
    view_before = vp_hash(page)
    if mi is not None and ci is not None:
        loc.nth(mi).set_input_files(str(mesh_file))
        page.wait_for_timeout(500)
        loc.nth(ci).set_input_files(str(cl_file))
    elif multi is not None:
        loc.nth(multi).set_input_files([str(mesh_file), str(cl_file)])
    else:
        note(False, f"{label}: cannot tell the mesh input from the centerline input")
        return None

    # The right rows are not enough if they were already on screen: the same case
    # loaded some other way (the 2c probe, on an app that fails it) would satisfy
    # this before the upload was even read. Require the status to have moved.
    seen = {"change": before["cards"] == 0}

    def arrived(st):
        if st["status"] != before["status"]:
            seen["change"] = True
        return (seen["change"] and st["cards"] == expect["n_openings"]
                and sorted(st["names"]) == expect["names"]
                and re.match(r"ready", st["status"], re.I) is not None)

    # An error status is not the end on its own: an app that loads as each file lands
    # can briefly pair the new mesh with the old centerline. Stop early only once an
    # error has sat unchanged for 20 s.
    t0, clicked = time.time(), False
    bad_since, bad_text = None, None
    st = page_state(page)
    while time.time() - t0 < timeout_s:
        page.wait_for_timeout(1000)
        st = page_state(page)
        if arrived(st):
            break
        if st["status_bad"]:
            if st["status"] != bad_text:
                bad_since, bad_text = time.time(), st["status"]
            elif time.time() - bad_since > 20:
                break
        else:
            bad_since = bad_text = None
        if not clicked and time.time() - t0 > 10:
            # Nothing happened on its own, so look for the button that does it.
            btn = page.get_by_role("button", name=re.compile(r"^\s*(load|upload)\b", re.I))
            for k in range(btn.count()):
                b = btn.nth(k)
                if b.is_visible() and b.is_enabled():
                    info(f"clicking {b.inner_text().strip()!r}")
                    b.click()
                    break
            clicked = True
    info(f"waited {time.time() - t0:.0f}s")
    info(f"status before : {before['status'][:110]!r}")
    info(f"status after  : {st['status'][:110]!r}")
    info(f"toolbar chips : {st['chips']}")
    info(f"opening rows  : {st['cards']}  {st['names']}")
    note(st["cards"] > 0, f"{label}: the openings table has rows ({st['cards']})")
    note(seen["change"] or st["status"] != before["status"],
         f"{label}: the status line moved during the upload")
    note(re.match(r"ready, \d+ openings", st["status"], re.I) is not None
         and not st["status_bad"],
         f"{label}: the status line says the cut is ready ({st['status'][:60]!r})")
    note(st["cards"] == expect["n_openings"],
         f"{label}: {expect['n_openings']} opening rows, as the same pair cut headlessly "
         f"({st['cards']})")
    note(sorted(st["names"]) == expect["names"],
         f"{label}: the rows are the uploaded case's openings, by name")
    note(chip_int(st["chips"], r"(-?\d+)\s+openings") == expect["n_openings"]
         and chip_int(st["chips"], r"euler\s+(-?\d+)") == expect["euler"],
         f"{label}: toolbar reads {expect['n_openings']} openings, euler {expect['euler']}")
    t1, view_after = time.time(), vp_hash(page)
    while view_after == view_before and time.time() - t1 < 20:
        page.wait_for_timeout(1000)
        view_after = vp_hash(page)
    info(f"viewport {view_before} -> {view_after} "
         f"({time.time() - t1:.0f}s after the rows arrived)")
    note(view_after != view_before,
         f"{label}: the 3-D view repainted with the uploaded model")
    return st


def vp_hash(page):
    # The canvas only: the empty-view card disappearing on load would otherwise make
    # "the 3-D view repainted" pass without the model ever being drawn.
    hide = page.add_style_tag(content='[data-test="view-overlay"]{visibility:hidden!important}')
    try:
        loc = page.locator(".v-main")
        png = loc.first.screenshot() if loc.count() else page.screenshot()
    finally:
        hide.evaluate("e => e.remove()")
    return hashlib.sha1(png).hexdigest()[:12]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--app", default=str(ROOT / "webapp" / "cut_app.py"))
    ap.add_argument("--port", type=int, default=8092)
    ap.add_argument("--skip-large", action="store_true")
    # Where run folders and fixtures go. Each run holds the store, screenshots and the
    # download -- about 10 MB -- and they accumulate. Projects/CLAUDE.md rule 4: a heavy
    # run root comes from PROJECTS_SHARE, and the script refuses without it unless
    # given an explicit root.
    ap.add_argument("--out-root", default=None,
                    help="where runs/ and fixtures/ go (default: the share cache)")
    args = ap.parse_args()
    port = args.port
    if port == 8080:
        sys.exit("Not on 8080: the user may have the app open there.")

    share = os.environ.get("PROJECTS_SHARE")
    if args.out_root:
        base = Path(args.out_root)
    elif share:
        base = Path(share) / ROOT.name / "_cache" / "check_upload"
    else:
        sys.exit("PROJECTS_SHARE is unset, so there is nowhere agreed to put the run "
                 "folders. Set it, or pass --out-root.")
    print(f"[outputs] {base}", flush=True)
    global FIXTURES, LARGE_MESH, LARGE_CL
    FIXTURES = base / "fixtures"
    LARGE_MESH = FIXTURES / LARGE_MESH.name
    LARGE_CL = FIXTURES / LARGE_CL.name
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run = base / "runs" / stamp
    store, dl = run / "store", run / "download"
    # The app's TEMP stays on this disk. It holds a few KB of job files, and reading a
    # folder the just-stopped server had open back off the SMB share failed with
    # WinError 59 on two runs out of two -- the check then died before its verdict.
    # tests/_out/ is local and gitignored, and each run gets its own folder.
    tmp = ROOT / "tests" / "_out" / "check_upload_tmp" / stamp
    for d in (store, tmp, dl):
        d.mkdir(parents=True, exist_ok=False)
    url = f"http://127.0.0.1:{port}/"
    export_name = f"wf_upload_{stamp}"
    print(f"app      : {args.app}")
    print(f"url      : {url}   (VC_PUBLIC=1)")
    print(f"run dir  : {run}")
    print(f"build    : {BUILD} expected on screen")

    # ── expectations, computed without the app ───────────────────────────────
    print("\ncutting the samples headlessly for the expected result...")
    exp_small = reference(SMALL_CL, SMALL_MESH)
    print(f"  KY160_1 : {exp_small['message']}  ({exp_small['n_faces']} faces)")
    exp_large = None
    if not args.skip_large:
        make_large_fixtures()
        # Expected from the .mat: the fixture is the same geometry written as STL.
        exp_large = reference(LARGE_SRC_CL, LARGE_SRC_MESH)
        from wslink.protocol import MAX_MSG_SIZE
        print(f"  JK5_1   : {exp_large['message']}  ({exp_large['n_faces']} faces)")
        print(f"  large fixture {LARGE_MESH.name}: {LARGE_MESH.stat().st_size:,} bytes, "
              f"WSLINK_MAX_MSG_SIZE {MAX_MSG_SIZE:,}")
        assert LARGE_MESH.stat().st_size > MAX_MSG_SIZE, "fixture is not over the cap"

    watched = {"output": ROOT / "output", "repo root": ROOT, "webapp": ROOT / "webapp"}
    before = {k: listing(v) for k, v in watched.items()}
    if (ROOT / "output" / f"{export_name}.stl").exists():
        sys.exit(f"output/{export_name}.stl already exists; refusing to run")

    if port_in_use(port):
        sys.exit(f"Something is already listening on {port}. Not touching it.")

    proc, log = start_app(Path(args.app), port, run, store, tmp)
    page_errs, frame_hits = [], []
    shot_n = [0]
    try:
        if not wait_listening(proc, port, run):
            log.flush()
            print((run / "server.log").read_text(encoding="utf-8", errors="replace")[-3000:])
            sys.exit("The app never started listening.")
        print(f"server pid {proc.pid} listening on {port}")

        needles_wire = [str(ROOT), str(ROOT).replace("\\", "/"),
                        str(ROOT).replace("\\", "\\\\"), str(store)]

        def on_ws(ws):
            def scan(payload):
                data = payload if isinstance(payload, bytes) else payload.encode("utf-8", "replace")
                for n in needles_wire:
                    if n.encode("utf-8") in data:
                        i = data.find(n.encode("utf-8"))
                        frame_hits.append(data[max(0, i - 60): i + len(n) + 60]
                                          .decode("utf-8", "replace"))
            ws.on("framereceived", scan)

        with sync_playwright() as p:
            b = p.chromium.launch()
            page = b.new_page(viewport={"width": 1500, "height": 1000},
                              accept_downloads=True)
            page.on("pageerror", lambda e: page_errs.append(str(e)))
            page.on("websocket", on_ws)
            downloads = []
            page.on("download", lambda d: downloads.append(d))

            def shot(tag):
                shot_n[0] += 1
                path = run / f"{shot_n[0]:02d}_{tag}.png"
                page.screenshot(path=str(path))
                return path

            page.goto(url, wait_until="domcontentloaded", timeout=60000)
            page.wait_for_selector(".v-application", timeout=60000)
            page.wait_for_timeout(2500)
            body = page.inner_text("body")
            print(f"\n  build stamp {BUILD!r} on screen: {BUILD in body}")
            set_section(page, "cut")
            shot("start")

            # ── 2 security ───────────────────────────────────────────────────
            check("2", "SECURITY: no text box takes a filesystem path in public mode")
            rows_all = []
            for section in ("cut", "shape"):
                set_section(page, section)
                titles = open_all_panels(page)
                rows = page.evaluate(INPUTS_JS)
                for r in rows:
                    r["section"] = section
                rows_all += rows
                info(f"{section} section: {len(rows)} inputs over panels {titles}")
                for r in rows:
                    if r["tag"] == "input" and r["type"] in {"text", "search", "url", ""}:
                        info(f"  text  {r['panel'][:18]:<18} label={r['label'][:20]!r:<22} "
                             f"disabled={r['disabled']!s:<5} value={r['value'][:48]!r}")
            set_section(page, "cut")
            bad = path_inputs(rows_all)
            for r in bad:
                info(f"PATH BOX: {r['section']}/{r['panel']} label={r['label']!r} "
                     f"placeholder={r['placeholder']!r}")
            note(not bad, f"no enabled text input takes a path ({len(bad)} found)")
            note(len(rows_all) > 5, f"and the scan saw the panels' inputs at all "
                                    f"({len(rows_all)}), so the zero is not an empty page")

            # 2b: disabling a path box does not stop its VALUE reaching the browser.
            check("2b", "DISCLOSURE: no server filesystem path reaches the browser")
            root_s, root_f = str(ROOT), str(ROOT).replace("\\", "/")
            leaked = [r for r in rows_all
                      if root_s.lower() in r["value"].lower()
                      or root_f.lower() in r["value"].lower()]
            for r in leaked:
                info(f"server path in {r['section']}/{r['panel']} {r['label']!r}: "
                     f"{r['value']!r}")
            note(not leaked, f"no input's value is a server filesystem path "
                             f"({len(leaked)} found)")
            vis = page.inner_text("body")
            note(root_s not in vis and root_f not in vis,
                 "no server filesystem path in the visible text")
            for h in dict.fromkeys(frame_hits):
                info(f"on the wire: ...{h!r}...")
            note(not frame_hits, f"no server filesystem path sent over the websocket "
                                 f"({len(frame_hits)} frame hit(s))")

            # ── 2c the refusal behind the widget ─────────────────────────────
            # A disabled box is a hint to a person. What stops a hostile client is the
            # server refusing a path it was handed, so be that client: write paths into
            # state from the page and fire the Load triggers directly. Run before any
            # upload, so a model on screen can only have come from the path.
            check("2c", "SECURITY: a client that bypasses the page cannot load by path")
            probe = page.evaluate("""() => {
              const t = window.trame;
              if (!t || !t.state) return {err: 'window.trame is not there'};
              const all = t.state.get();
              const keys = Object.keys(all).filter(k => k.startsWith('trame__template'));
              const tpl = keys.map(k => String(all[k] || '')).join('\\n');
              const found = {};
              let parsed = 0;
              // trame writes <VBtn>, not <v-btn>. Accept both.
              const re = /<(?:v-btn|VBtn)\\b([^>]*)>\\s*([^<]*?)\\s*<\\/(?:v-btn|VBtn)>/g;
              let m;
              while ((m = re.exec(tpl))) {
                const trg = /trigger\\('([^']+)'/.exec(m[1]);
                if (trg) parsed += 1;
                // Any button whose label starts with Load, so a rename ("Load case")
                // does not quietly take it out of the probe.
                if (trg && /^load\\b/i.test(m[2].trim())) found[m[2].trim()] = trg[1];
              }
              return {keys, found, parsed};
            }""")
            info(f"templates {probe.get('keys')}, {probe.get('parsed')} triggered buttons "
                 f"parsed, load triggers {probe.get('found')}")
            note("err" not in probe, f"the page's trame client is reachable "
                                     f"({probe.get('err', 'yes')})")
            # Guard against this probe going blind: it once matched <v-btn> while the
            # template says <VBtn>, found nothing, fired nothing and passed.
            note((probe.get("parsed") or 0) >= 5,
                 f"the template parser sees the page's buttons ({probe.get('parsed')}), so "
                 f"a missing Load trigger means it is really not on the page")
            if not probe.get("found"):
                info("no Load trigger in the template: unreachable from a client")
            fire = """async ([name, kv]) => {
              await window.trame.state.update(kv);
              try { return String(await window.trame.trigger(name, [], {})); }
              catch (e) { return 'threw: ' + (e && (e.message || e)); }
            }"""
            kv = {"mat_path": str(SMALL_MESH), "csv_path": str(SMALL_CL),
                  "sm_path": str(SMALL_MESH)}
            for label, trig in (probe.get("found") or {}).items():
                r = page.evaluate(fire, [trig, kv])
                info(f"fired {trig} ({label}) with server paths in state -> {r[:60]}")
            page.wait_for_timeout(6000)
            st0 = page_state(page)
            sm_info = page.evaluate("() => window.trame && window.trame.state.get('sm_info')")
            info(f"after: {st0['cards']} rows, status {st0['status'][:90]!r}, "
                 f"sm_info {str(sm_info)[:60]}")
            note(st0["cards"] == 0 and not re.match(r"ready", st0["status"], re.I),
                 "the cut section did not load the case named by path")
            note(not sm_info, "the shape section did not load the surface named by path")
            # Put the boxes back, so nothing later reads the probe's paths.
            page.evaluate("() => window.trame.state.update("
                          "{mat_path: '', csv_path: '', sm_path: ''})")
            page.wait_for_timeout(500)

            # ── 3 upload ─────────────────────────────────────────────────────
            check("3", "UPLOAD: KY160_1.mat + KY160_1_CentrelinePoints.csv via file inputs")
            open_panel(page, "Input")
            st_small = do_upload(page, SMALL_MESH, SMALL_CL, exp_small, "KY160")
            shot("after_upload")

            # ── 4 download ───────────────────────────────────────────────────
            check("4", "DOWNLOAD: the cut arrives in the browser and stands up on its own")
            open_panel(page, "Output")
            name_box = page.locator(".v-expansion-panel", has=page.locator(
                ".v-expansion-panel-title", has_text=re.compile(r"^\s*Output\s*$", re.I))
            ).locator("input[type=text]")
            if name_box.count() and name_box.first.is_enabled():
                name_box.first.fill(export_name)
                name_box.first.press("Tab")
                page.wait_for_timeout(400)
                info(f"export name set to {export_name}")
            ctrl = None
            for role in ("button", "link"):
                c = page.get_by_role(role, name=re.compile(r"download", re.I))
                for k in range(c.count()):
                    if c.nth(k).is_visible():
                        ctrl = c.nth(k)
                        break
                if ctrl is not None:
                    break
            exp_btn = page.get_by_role("button", name=re.compile(r"export open mesh", re.I))
            if exp_btn.count():
                info(f"'Export open mesh' present, enabled={exp_btn.first.is_enabled()} "
                     f"(writes to the server's output/, not a download)")
            note(ctrl is not None, "a Download control is on the page")
            got = []
            if ctrl is not None:
                enabled = ctrl.is_enabled()
                note(enabled, f"and it is enabled ({ctrl.inner_text().strip()!r})")
                if enabled:
                    n0 = len(downloads)
                    ctrl.click()
                    t0 = time.time()
                    while time.time() - t0 < 60 and len(downloads) == n0:
                        page.wait_for_timeout(500)
                    page.wait_for_timeout(3000)          # a second file, if two
                    got = downloads[n0:]
            note(len(got) > 0, f"the browser received a download ({len(got)} file(s))")
            saved = []
            for d in got:
                fail = d.failure()
                sname = re.sub(r"[^A-Za-z0-9._-]", "_", Path(d.suggested_filename).name) or "dl"
                dest = dl / sname
                k = 1
                while dest.exists():
                    dest = dl / f"{k}_{sname}"
                    k += 1
                d.save_as(str(dest))
                info(f"download {d.suggested_filename!r} -> {dest.name} "
                     f"({dest.stat().st_size:,} bytes){' FAILED ' + fail if fail else ''}")
                saved.append(dest)
            files = []
            for f in saved:
                if zipfile.is_zipfile(f):
                    out = dl / (f.stem + "_unzipped")
                    out.mkdir()
                    with zipfile.ZipFile(f) as z:
                        for m in z.infolist():
                            if m.is_dir():
                                continue
                            tgt = (out / m.filename).resolve()
                            if not str(tgt).startswith(str(out.resolve())):
                                note(False, f"zip member escapes the folder: {m.filename!r}")
                                continue
                            tgt.parent.mkdir(parents=True, exist_ok=True)
                            tgt.write_bytes(z.read(m))
                            files.append(tgt)
                            info(f"  unzipped {m.filename} ({m.file_size:,} bytes)")
                else:
                    files.append(f)
            meshes = [f for f in files if f.suffix.lower() in MESH_EXT]
            sides = [f for f in files if f.suffix.lower() == ".json"]
            dl_mesh = None
            if got:
                note(len(meshes) == 1, f"exactly one mesh in the download ({[m.name for m in meshes]})")
                note(len(sides) >= 1, f"a sidecar JSON in the download ({[s.name for s in sides]})")
            loops = None
            if meshes:
                try:
                    dl_mesh = trimesh.load(str(meshes[0]), force="mesh")
                    loops = open_loops(dl_mesh)
                    info(f"mesh: {len(dl_mesh.faces):,} faces, {len(dl_mesh.vertices):,} "
                         f"vertices, euler {dl_mesh.euler_number}, {loops} open loops")
                    ok_geom = (len(dl_mesh.faces) > 1000
                               and np.isfinite(dl_mesh.vertices).all())
                    note(ok_geom, "trimesh loads it as a real triangle mesh")
                    note(loops == exp_small["n_openings"],
                         f"{exp_small['n_openings']} open boundary loops, one per opening "
                         f"({loops})")
                    note(st_small is not None and loops == st_small["cards"],
                         "the same number as opening rows on the page")
                    note(dl_mesh.euler_number == exp_small["euler"],
                         f"euler {exp_small['euler']} as reported ({dl_mesh.euler_number})")
                    note(abs(len(dl_mesh.faces) - exp_small["n_faces"])
                         <= 0.005 * exp_small["n_faces"],
                         f"face count matches the headless cut ({len(dl_mesh.faces):,} "
                         f"vs {exp_small['n_faces']:,})")
                except Exception as e:                           # noqa: BLE001
                    note(False, f"trimesh could not load it: {type(e).__name__}: {e}")
            if sides:
                try:
                    sc = json.loads(sides[0].read_text(encoding="utf-8"))
                    ops = sc.get("openings") or []
                    names = [str(o.get("name") or "") for o in ops]
                    info(f"sidecar {sides[0].name}: {len(ops)} openings: {names}")
                    note(isinstance(ops, list) and len(ops) == exp_small["n_openings"],
                         f"the sidecar lists {exp_small['n_openings']} cut planes ({len(ops)})")
                    note(all(names), "every plane is named")
                    note(sorted(names) == exp_small["names"],
                         "with the names shown on the page")

                    def vec3(v):
                        return (isinstance(v, list) and len(v) == 3
                                and all(isinstance(x, (int, float)) and np.isfinite(x)
                                        for x in v))
                    note(all(vec3(o.get("plane_origin")) and vec3(o.get("plane_normal"))
                             for o in ops), "every plane has an origin and a normal")
                    txt = sides[0].read_text(encoding="utf-8")
                    note(str(ROOT) not in txt and str(ROOT).replace("\\", "\\\\") not in txt
                         and str(store).replace("\\", "\\\\") not in txt,
                         "the sidecar names no server path")
                except Exception as e:                           # noqa: BLE001
                    note(False, f"sidecar did not parse: {type(e).__name__}: {e}")
            note(not (ROOT / "output" / f"{export_name}.stl").exists(),
                 "nothing was written to the server's output/ under the export name")
            shot("after_download")

            # ── 3L large upload ──────────────────────────────────────────────
            if exp_large is not None:
                check("3L", f"LARGE UPLOAD: {LARGE_MESH.stat().st_size / 1e6:.1f} MB ASCII "
                            f"STL, over WSLINK_MAX_MSG_SIZE, under a leaky filename")
                open_panel(page, "Input")
                do_upload(page, LARGE_MESH, LARGE_CL, exp_large, "JK5 large")
                shot("after_large")

            b.close()

        # ── 5 storage (after the uploads and the download) ───────────────────
        check("5", "STORAGE: stored under a generated id, the filename only as metadata")
        stored = sorted(p for p in store.rglob("*") if p.is_file())
        for p in stored[:40]:
            info(f"stored: {p.relative_to(store)}  ({p.stat().st_size:,} bytes)")
        note(len(stored) > 0, f"something was stored under {STORE_ENV} ({len(stored)} files)")
        by_sha = {}
        for p in stored:
            by_sha.setdefault(sha(p), []).append(p)
        uploads = [("KY160_1.mat", SMALL_MESH), ("KY160_1_CentrelinePoints.csv", SMALL_CL)]
        if exp_large is not None:
            uploads += [(LARGE_MESH.name, LARGE_MESH), (LARGE_CL.name, LARGE_CL)]
        sub_dirs = set()
        for orig, src in uploads:
            hits = by_sha.get(sha(src), [])
            note(bool(hits), f"{orig}: a byte-identical copy is stored")
            for h in hits:
                parts = h.relative_to(store).parts
                ids = [x for x in parts[:-1] if ID_RE.match(x)
                       and not any(n.lower() in x.lower() for n in NEEDLES)]
                note(bool(ids), f"{orig}: stored under a generated id "
                                f"({'/'.join(parts)})")
                if ids:
                    sub_dirs.add(store.joinpath(*parts[:parts.index(ids[0]) + 1]))
        # The original filename, kept -- but as data, not as a path.
        for orig, _src in uploads:
            meta_hit = []
            for d in sub_dirs:
                for j in d.rglob("*.json"):
                    try:
                        if orig in j.read_text(encoding="utf-8", errors="replace"):
                            meta_hit.append(j)
                    except OSError:
                        pass
            note(bool(meta_hit), f"{orig}: original filename recorded as metadata "
                                 f"({[str(m.relative_to(store)) for m in meta_hit][:2]})")
        # Job 7: the cut and its sidecar are kept with the uploads.
        kept_mesh = kept_side = False
        for d in sub_dirs:
            for f in d.rglob("*"):
                if f.suffix.lower() in MESH_EXT and dl_mesh is not None:
                    try:
                        m = trimesh.load(str(f), force="mesh")
                        if len(m.faces) == len(dl_mesh.faces):
                            kept_mesh = True
                    except Exception:                           # noqa: BLE001
                        pass
                if f.suffix.lower() == ".json":
                    try:
                        if "openings" in json.loads(f.read_text(encoding="utf-8")):
                            kept_side = True
                    except Exception:                           # noqa: BLE001
                        pass
        note(kept_mesh, "the cut mesh that was downloaded is kept with the uploads")
        note(kept_side, "and so is its sidecar")

        # Every path this script can see being written, searched for the filenames.
        written = [f"store/{x}" for x in listing(store, recursive=True)]
        written += [f"tmp/{x}" for x in listing(tmp, recursive=True)]
        for k, v in watched.items():
            new = sorted(listing(v) - before[k])
            written += [f"{k}/{x}" for x in new]
            for x in new:
                info(f"new in {k}/: {x}  (possibly another process)")
        leaks = [w for w in written if any(n.lower() in w.lower() for n in NEEDLES)]
        for w in leaks:
            info(f"LEAK: {w}")
        note(not leaks, f"no path written contains an uploaded filename or case id "
                        f"({len(written)} paths searched for {NEEDLES})")
        left = listing(tmp, recursive=True)
        info(f"temp files left behind: {len(left)}")

    finally:
        stopped = stop_app(proc, port)
        log.close()
        print(f"\nserver stopped, port {port} free: {stopped}")

    # ── E page and server errors ─────────────────────────────────────────────
    check("E", "PAGE: no page errors, no server-side traceback")
    for e in page_errs[:6]:
        info(f"page error: {e[:140]}")
    note(not page_errs, f"no page errors ({len(page_errs)})")
    text = (run / "server.log").read_text(encoding="utf-8", errors="replace")
    # One block per traceback: from "Traceback" to its unindented exception line.
    # Windows' proactor loop logs a ConnectionResetError from _call_connection_lost
    # when the browser goes away; that is the socket closing, not the app failing,
    # and it is reported rather than hidden.
    blocks, cur = [], None
    for ln in text.splitlines():
        if ln.startswith("Traceback"):
            cur = [ln]
        elif cur is not None:
            cur.append(ln)
            if ln and not ln[0].isspace():
                blocks.append(cur)
                cur = None
    benign = [b for b in blocks if "_call_connection_lost" in "\n".join(b)
              and b[-1].startswith("ConnectionResetError")]
    real = [b for b in blocks if b not in benign]
    for b in benign:
        info(f"server log (socket closed, not counted): {b[-1][:110]}")
    for b in real:
        info(f"server log: {b[-1][:140]}")
        for ln in b[-5:-1]:
            info(f"    {ln.strip()[:130]}")
    note(not real, f"no traceback in the server log ({len(real)})")
    note(stopped, f"the server was stopped and port {port} released")

    print()
    print("=" * 74)
    print("verdicts")
    print("=" * 74)
    fails = []
    for cid in ORDER:
        v = verdict(cid)
        bad = [w for ok, w in CHECKS[cid]["notes"] if not ok]
        print(f"  {cid:<3} {v}  {CHECKS[cid]['title']}")
        if bad:
            print(f"        first failure: {bad[0]}")
            fails.append(cid)
    print(f"\n  run folder -> {run}")
    print()
    print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): {', '.join(fails)}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
