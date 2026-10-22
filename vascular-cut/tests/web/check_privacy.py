"""Two visitors, one process: nothing of one may reach the other.

The acceptance test for visitor privacy in hosted mode (VC_PUBLIC=1). The fix it
holds in place -- webapp/visit_gate.py and the wipe in webapp/cut_app.py -- admits one
browser per process, rotates the upload token per admission, wipes every trace of a
visitor the moment their socket closes, and uses a generation counter so work still
running for a visitor who has gone cannot install its result. Written before the fix
and run against the app as it stood, where A, A2, A3, C, D, E0-E3 and F failed; a
check that passes against code that does not exist measures nothing.

Like check_upload.py it starts the app itself: hosted mode (VC_PUBLIC=1), loopback,
store and temp folders pointed where this script can look.

What it drives, which is therefore the contract:

  admission  ONE browser at a time. While a visitor is admitted, any other connection
             -- a real page, or a raw websocket that never loads the page -- receives
             no state, no published update and no picture. Refusing the wslink hello
             gives exactly that: an unauthenticated client is sent nothing
             (wslink/protocol.py:289-297, 377-395, 472-481). The hello secret is
             trame's public default, so it is no barrier on its own.
  told       the refused page says why, in words REFUSED_RE matches.      (A2, policy)
  token      vc_token in state: new for each admitted visitor, unchanged by refused
             attempts, refused (not 200) once its visitor has gone.
  wipe       when the admitted socket closes, however it closes (tab closed without
             beforeunload, reload, navigation), the trame state goes back key for key
             to what a fresh process sends (vc_token excepted), the 3-D view holds
             nothing but the axes triad in either section, and nothing pending --
             half a pair, an upload still streaming, a download link, a load, a zone
             map, a delta K job -- reaches the next visitor.
  reload     a reload is a new visitor: admitted, not refused as "someone else", and
             blank.

Checks, each with its own verdict:

  0   SETUP      a fresh visitor is admitted; the probes are shown to see (needles
                 absent from a fresh page, fresh view empty, loaded view not)
  A   REFUSED    browser B while A is admitted: no state, publish or picture on B's
                 socket; nothing of A's in B's DOM, client state or HTTP responses
  A2  TOLD       B's page says the tool is in use                         (policy)
  A3  HOSTILE    raw websocket: hello, trame.state.get, then listen while A works
  B   UNHARMED   A's session, token, rows, view and Download survive B and A3
  C   CLEAN      A closes its tab; the next visitor starts from nothing
  D   STALE      A's token is refused and cannot complete C's half pair; A's unfollowed
                 download link is gone
  E0  HALF-PAIR  a visitor leaves with one file pending; it does not wait for the next
  E1  IN-FLIGHT  an upload still streaming when its visitor leaves reaches nobody
  E2  MID-LOAD   visitor leaves during a ~19 s build; the next visitor is already on the
                 page when it finishes, and never sees it
  E3  LATE-LOAD  a ~1 s build finishes with nobody connected; the next visitor is clean
  E4  JOB        (--with-shape) a delta K job running at disconnect never reports
  F   RELOAD     a loaded visitor reloads: admitted, blank, new token
  ERR PAGE       no page errors on admitted pages, no traceback, server still up

Run:
    $env:PYTHONIOENCODING="utf-8"
    .venv\\Scripts\\python tests\\web\\check_privacy.py
    .venv\\Scripts\\python tests\\web\\check_privacy.py --with-shape

Every run writes a new timestamped folder under <out root>/runs/. The out root is
$PROJECTS_SHARE/branches/_cache/check_privacy unless --out-root is given, and the script
refuses to run without one. The app's TEMP is a new local folder under tests/_out/
(gitignored; reading it back off SMB failed in check_upload). Nothing is overwritten
and nothing is deleted. Exits 1 on any failure.
"""
import argparse
import asyncio
import hashlib
import http.client
import io
import json
import os
import re
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
# The verdict book and the page helpers are check_upload's, so both checks read alike
# and a fix to a helper reaches both.
from check_upload import (                                        # noqa: E402
    check, note, info, verdict, CHECKS, ORDER, ROOT, BUILD, SMALL_MESH, SMALL_CL,
    port_in_use, listing, start_app, stop_app, wait_listening, open_panel,
    open_all_panels, set_section, INPUTS_JS, page_state, classify_file_inputs,
    do_upload, reference,
)

DEFAULT_PORT = 8093     # 8080 the user's; 8087 README; 8092 check_upload; 8099 HOSTING
SENT_A = "wfPRIVA"      # visitor A's export name carries it, and nothing else does
SENT_E = "wfPRIVE"      # in the race fixture's filename
RACE_NAME = f"{SENT_E}_KY160_sub2.stl"
RACE_FACES = 636_832    # KY160_1.mat subdivided twice; built headlessly in 19.2 s
REFUSED_RE = re.compile(r"\b(in use|busy|someone else|another "
                        r"(visitor|person|browser|session))\b", re.I)

# The 3-D view. Calibrated on the app's own Scene, rendered off-screen at 1070x900 and
# JPEG'd: fresh 0 px of ink outside the triad corner; KY160 loaded 71,446 px (7.9%);
# a lone leftover centerline 5,352 px, lone boundary loops 4,454 px. The triad fills
# 0.8-1.0% of its corner in every frame, which is how an image is known to have come.
AXES = 0.25             # pyvista add_axes viewport is (0, 0, 0.2, 0.2); plus margin
INK = 24                # this far from white, in any channel, is ink
EMPTY_MAX = 0.0005      # <= 0.05% ink outside the corner is empty (~480 px)
LOADED_MIN = 0.02       # the loaded-view calibration must clear this
RENDERED_MIN = 0.003    # the triad's share of its corner, proving a frame arrived
ADMIT_S = 15            # the slot must be free this soon after a visitor leaves

VIEW_BOX_JS = """() => {
  const main = document.querySelector('.v-main');
  if (!main) return null;
  let best = null;
  for (const c of main.querySelectorAll('canvas, img')) {
    const r = c.getBoundingClientRect();
    if (r.width > 50 && r.height > 50 && (!best || r.width * r.height > best.width * best.height))
      best = r;
  }
  if (!best) {   // v-main is PADDED under the drawer, toolbar and footer: take its content box
    const r = main.getBoundingClientRect(), s = getComputedStyle(main);
    const l = parseFloat(s.paddingLeft), t = parseFloat(s.paddingTop);
    best = {x: r.x + l, y: r.y + t, width: r.width - l - parseFloat(s.paddingRight),
            height: r.height - t - parseFloat(s.paddingBottom)};
  }
  return {x: Math.round(best.x), y: Math.round(best.y),
          width: Math.round(best.width), height: Math.round(best.height)};
}"""

STATE_JS = """() => { try {
  const s = window.trame.state.get();
  return JSON.parse(JSON.stringify(Object.fromEntries(
    Object.entries(s).filter(([k]) => !k.startsWith('trame__')))));
} catch (e) { return null; } }"""

ADMITTED_JS = """() => { try {
  return !!(window.trame && window.trame.state && window.trame.state.get('vc_token')
            && document.querySelector('.v-navigation-drawer'));
} catch (e) { return false; } }"""


class Stop(Exception):
    """A step the rest depends on did not happen; go straight to the verdicts."""


# ── one person ───────────────────────────────────────────────────────────────
class Visitor:
    """A browser context of its own, its page, and everything its socket received."""

    def __init__(self, browser, label, url):
        self.label, self.url = label, url
        self.ctx = browser.new_context(viewport={"width": 1500, "height": 1000},
                                       accept_downloads=True)
        self.page = self.ctx.new_page()
        self.frames, self.errors, self.responses, self.downloads = [], [], [], []
        self.t_admitted = self.t_closed = None
        p = self.page
        p.on("pageerror", lambda e: self.errors.append(str(e)))
        p.on("download", lambda d: self.downloads.append(d))
        p.on("response", lambda r: self.responses.append(r))

        def on_ws(ws):
            ws.on("framereceived", lambda b: self.frames.append(
                b if isinstance(b, bytes) else b.encode("utf-8", "replace")))
        p.on("websocket", on_ws)

    def go(self, wait_s=30):
        self.page.goto(self.url, wait_until="domcontentloaded", timeout=60000)
        return self.wait_admitted(wait_s)

    def wait_admitted(self, wait_s):
        t0 = time.time()
        while time.time() - t0 < wait_s:
            if self.page.evaluate(ADMITTED_JS):
                self.t_admitted = time.time()
                self.page.wait_for_timeout(2000)            # let the first picture land
                return True
            self.page.wait_for_timeout(500)
        return False

    def state(self):
        return self.page.evaluate(STATE_JS)

    def token(self):
        return (self.state() or {}).get("vc_token") or ""

    def http_bodies(self):
        out = []
        for r in self.responses:
            try:
                out.append(r.body())
            except Exception:                                   # noqa: BLE001
                pass
        return out

    def close(self, before_unload=False):
        # run_before_unload=False is Playwright's default and the hard case: no
        # beforeunload, so trame's client never sends its client_exited lifecycle
        # message. Only the socket closes -- as for a crashed tab or a dropped network.
        self.page.close(run_before_unload=before_unload)
        self.t_closed = time.time()
        self.ctx.close()


def admit(browser, label, url, within_s=ADMIT_S):
    """The next person arriving, retrying as a person would. (visitor, attempts, s)."""
    t0, n = time.time(), 0
    while time.time() - t0 < within_s:
        n += 1
        v = Visitor(browser, label if n == 1 else f"{label}#{n}", url)
        if v.go(wait_s=max(3.0, within_s - (time.time() - t0))):
            return v, n, time.time() - t0
        v.close()
    return None, n, time.time() - t0


def must_admit(browser, label, url, after=None):
    v, n, s = admit(browser, label, url)
    gap = f", {v.t_admitted - after.t_closed:.1f}s after {after.label} left" if (
        v and after and after.t_closed) else ""
    note(v is not None, f"{label} is admitted within {ADMIT_S}s ({n} attempt(s){gap})")
    if v is None:
        raise Stop(f"{label} never got in")
    return v


# ── the websocket ────────────────────────────────────────────────────────────
def decode(frames):
    """Reassemble wslink's chunked msgpack frames into messages."""
    from wslink.chunking import UnChunker
    u = UnChunker()
    u.set_max_message_size(1 << 31)
    out = []
    for f in frames:
        try:
            m = u.process_chunk(f)
        except Exception as e:                                  # noqa: BLE001
            out.append({"id": "undecodable", "error": str(e)})
            continue
        if m is not None:
            out.append(m)
    return out


def kind(m):
    i = str(m.get("id", ""))
    if i.startswith("publish:"):
        return "publish:" + i.split(":")[1]
    if "error" in m:
        return "error"
    r = m.get("result")
    if isinstance(r, dict) and "state" in r:
        return "STATE"                  # the reply to trame.state.get
    if isinstance(r, dict) and "clientID" in r:
        return "hello-ok"
    return "result"


def needle_hits(blobs, needles):
    return sorted({n for b in blobs for n in needles if n and n.encode() in b})


def hostile(port, listen_s=12):
    """A client that never loads the page: hello with trame's public default secret,
    ask for the whole state, then listen. Runs in its own thread and event loop."""
    import aiohttp
    import msgpack
    from wslink.chunking import generate_chunks, UnChunker
    out = {"raw": [], "msgs": [], "state_get": "no answer", "upgrade": None}

    def chunks(m):
        return list(generate_chunks(msgpack.packb(m), 4 * 2 ** 20))

    async def go():
        async with aiohttp.ClientSession() as s:
            try:
                ws = await s.ws_connect(f"http://127.0.0.1:{port}/ws",
                                        timeout=aiohttp.ClientWSTimeout(ws_close=10))
            except aiohttp.WSServerHandshakeError as e:
                out["upgrade"] = f"refused {e.status}"
                return
            out["upgrade"] = "accepted"
            for c in chunks({"wslink": "1.0", "id": "system:c0:0",
                             "method": "wslink.hello",
                             "args": [{"secret": "wslink-secret"}], "kwargs": {}}):
                await ws.send_bytes(c)
            await asyncio.sleep(0.5)
            for c in chunks({"wslink": "1.0", "id": "rpc:cX:1",
                             "method": "trame.state.get", "args": [], "kwargs": {}}):
                await ws.send_bytes(c)
            u = UnChunker()
            u.set_max_message_size(1 << 31)
            t_end = time.time() + listen_s
            while time.time() < t_end:
                try:
                    msg = await ws.receive(timeout=max(0.1, t_end - time.time()))
                except asyncio.TimeoutError:
                    break
                if msg.type == aiohttp.WSMsgType.BINARY:
                    out["raw"].append(msg.data)
                    m = u.process_chunk(msg.data)
                    if m is not None:
                        out["msgs"].append(m)
                        if m.get("id") == "rpc:cX:1":
                            out["state_get"] = "error" if "error" in m else "ANSWERED"
                elif msg.type in (aiohttp.WSMsgType.CLOSE, aiohttp.WSMsgType.CLOSING,
                                  aiohttp.WSMsgType.CLOSED):
                    out["closed"] = True
                    break
            await ws.close()

    asyncio.run(go())
    return out


# ── HTTP, as a client that is not the page ───────────────────────────────────
def post(req, kind_, token, path, name):
    return req.post(f"vc/upload/{kind_}", data=Path(path).read_bytes(), headers={
        "X-VC-Token": token, "X-VC-Name": quote(name),
        "Content-Type": "application/octet-stream"})


def slow_post(port, kind_, token, path, name, rate_bps, box):
    """One upload, streamed slowly, so it is still arriving when its visitor leaves.
    Real, not contrived: a reverse proxy that buffers request bodies (nginx does by
    default) hands the app a complete upload after the uploader's socket is gone."""
    data = Path(path).read_bytes()
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=180)
    c.putrequest("POST", f"/vc/upload/{kind_}")
    for k, v in (("X-VC-Token", token), ("X-VC-Name", quote(name)),
                 ("Content-Type", "application/octet-stream"),
                 ("Content-Length", str(len(data)))):
        c.putheader(k, v)
    c.endheaders()
    box["total"], t0, step = len(data), time.time(), 16384
    # The server may answer, and drop the connection, before the body is all sent:
    # that is the refusal this checks for, not a failure of the check. Either way the
    # outcome is recorded, never raised, so the thread always reports.
    try:
        for i in range(0, len(data), step):
            c.send(data[i:i + step])
            box["sent"] = min(len(data), i + step)
            time.sleep(max(0.0, box["sent"] / rate_bps - (time.time() - t0)))
    except OSError as e:
        box["send_error"] = type(e).__name__
    try:
        r = c.getresponse()
        box["status"], box["body"] = r.status, r.read()[:300]
    except (OSError, http.client.HTTPException) as e:
        box["status"], box["body"] = f"dropped ({type(e).__name__})", b""


# ── the page ─────────────────────────────────────────────────────────────────
def upload_one(page, kind_, path):
    files, mi, ci, _multi = classify_file_inputs(page.evaluate(INPUTS_JS))
    idx = mi if kind_ == "surface" else ci
    if idx is None:
        return note(False, f"cannot find the {kind_} input")
    page.locator("input[type=file]").nth(idx).set_input_files(str(path))
    return True


def wait_status(page, pattern, timeout_s):
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        st = page_state(page)
        if re.search(pattern, st["status"], re.I):
            return st
        page.wait_for_timeout(150)
    return None


def dom_blob(page):
    open_all_panels(page)          # a collapsed panel's content is not in the DOM
    vals = [r["value"] for r in page.evaluate(INPUTS_JS)]
    return page.inner_text("body") + "\n" + "\n".join(vals)


def view(page):
    box = page.evaluate(VIEW_BOX_JS)
    if not box or box["width"] < 50 or box["height"] < 50:
        return None
    # The canvas only: the page's own cards and buttons over the view are not what
    # a visitor could have left behind, and they differ between states by design.
    hide = page.add_style_tag(content='[data-test="view-overlay"]{visibility:hidden!important}')
    try:
        a = np.asarray(Image.open(io.BytesIO(page.screenshot(clip=box)))
                       .convert("RGB")).astype(np.int16)
    finally:
        hide.evaluate("e => e.remove()")
    h, w, _ = a.shape
    corner = np.zeros((h, w), bool)
    corner[int(h * (1 - AXES)):, :int(w * AXES)] = True
    ink = (255 - a).max(axis=2) > INK
    return {"out": float(ink[~corner].mean()), "px": int(ink[~corner].sum()),
            "n": int((~corner).sum()), "corner": float(ink[corner].mean()),
            "sha": hashlib.sha1(a.tobytes()).hexdigest()[:12], "box": box}


# What "empty" looks like in THIS browser: the fresh process's own view, measured in
# check 0. The off-screen calibration read 0 px for an empty scene, but the page's view
# carries a few hundred px of static chrome (the drawer's edge) that is there in every
# frame, fresh or not -- 804 px on the first run. So empty is judged against the fresh
# view, not against zero.
BASELINE = {"px": None}


def settled_view(page, timeout_s=25):
    """The server's picture, once one has arrived (the triad is drawn in every frame)
    and it has stopped changing. A blank canvas is 'not rendered', never 'empty'."""
    prev, t0 = None, time.time()
    while time.time() - t0 < timeout_s:
        v = view(page)
        if (v and prev and v["corner"] >= RENDERED_MIN and v["sha"] == prev["sha"]):
            v["settled"] = True
            return v
        prev = v
        page.wait_for_timeout(1000)
    if prev:
        prev["settled"] = False
    return prev


def note_view_empty(page, label):
    v = settled_view(page)
    if v is None:
        return note(False, f"{label}: there is no 3-D view to look at")
    info(f"{label}: {v['px']} px of ink outside the triad ({100 * v['out']:.3f}%), "
         f"triad {100 * v['corner']:.2f}% of its corner")
    note(v["settled"], f"{label}: a picture arrived from the server, so 'empty' means empty")
    if BASELINE["px"] is None:
        # The fresh process defines empty. Sanity-bound it, so a baseline that already
        # held a model could not make everything after it look clean.
        BASELINE["px"] = v["px"]
        return note(v["settled"] and v["out"] <= 0.005,
                    f"{label}: nearly empty ({v['px']} px), and taken as the baseline")
    extra = v["px"] - BASELINE["px"]
    return note(v["settled"] and extra <= EMPTY_MAX * v["n"],
                f"{label}: nothing in the view a fresh process does not show "
                f"({extra:+d} px against it)")


def state_diff(fresh, now):
    keys = sorted((set(fresh) | set(now)) - {"vc_token"})

    def js(x):
        return json.dumps(x, sort_keys=True, default=str)
    out = []
    for k in keys:
        if k not in now:
            out.append((k, "(fresh)", "(missing)"))
        elif k not in fresh:
            out.append((k, "(absent when fresh)", now[k]))
        elif js(fresh[k]) != js(now[k]):
            out.append((k, fresh[k], now[k]))
    return out


def assert_clean(v, fresh, needles, old_tokens, frames_from=0):
    """Everything C must show; reused for every visitor who follows another."""
    lab = v.label
    st = v.state() or {}
    diff = state_diff(fresh["state"], st)
    for k, a, b in diff[:15]:
        info(f"{lab}: {k}: fresh {str(a)[:40]!r} -> now {str(b)[:70]!r}")
    note(not diff, f"{lab}: state equals a fresh process's, key for key "
                   f"({len(diff)} of {len(fresh['state'])} differ)")
    tok = st.get("vc_token") or ""
    note(bool(tok) and tok not in old_tokens,
         f"{lab}: an upload token of its own, not any earlier visitor's")
    ps = page_state(v.page)
    note(ps["cards"] == 0, f"{lab}: no opening rows ({ps['cards']})")
    note(ps["status"] == fresh["status"],
         f"{lab}: status line as a fresh page has it ({ps['status'][:60]!r})")
    blob = dom_blob(v.page)
    hit = [n for n in needles if n and n.lower() in blob.lower()]
    note(not hit, f"{lab}: nothing of an earlier visitor in the page text or any input {hit}")
    fhit = needle_hits(v.frames[frames_from:], needles)
    note(not fhit, f"{lab}: nothing of an earlier visitor on its websocket {fhit}")
    set_section(v.page, "cut")
    note_view_empty(v.page, f"{lab} cut view")
    set_section(v.page, "shape")
    note_view_empty(v.page, f"{lab} shape view")
    s2 = v.state() or {}
    note(not s2.get("sm_info") and not s2.get("sm_result") and not s2.get("sm_busy"),
         f"{lab}: the shape section holds no surface, result or running job")
    set_section(v.page, "cut")
    return tok


def log_tail(run, since_bytes):
    t = (run / "server.log").read_text(encoding="utf-8", errors="replace")
    return t[since_bytes:], len(t)


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--app", default=str(ROOT / "webapp" / "cut_app.py"))
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--out-root", default=None,
                    help="where runs/ and fixtures/ go (default: the share cache)")
    ap.add_argument("--with-shape", action="store_true",
                    help="E4: a delta K job running at disconnect (about 2 min more)")
    args = ap.parse_args()
    port = args.port
    if port == 8080:
        sys.exit("Not on 8080: the user may have the app open there.")

    # Projects/CLAUDE.md rule 4: a heavy run root comes from PROJECTS_SHARE, and the
    # script refuses without it unless given one. A run stores ~40 MB of uploads.
    share = os.environ.get("PROJECTS_SHARE")
    if args.out_root:
        base = Path(args.out_root)
    elif share:
        if not Path(share).is_dir():
            sys.exit(f"PROJECTS_SHARE ({share}) is not reachable. Not writing elsewhere.")
        base = Path(share) / ROOT.name / "_cache" / "check_privacy"
    else:
        sys.exit("PROJECTS_SHARE is unset, so there is nowhere agreed to put the run "
                 "folders. Set it, or pass --out-root.")
    print(f"[outputs] {base}", flush=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run = base / "runs" / stamp
    store = run / "store"
    fixtures = base / "fixtures"
    tmp = ROOT / "tests" / "_out" / "check_privacy_tmp" / stamp   # local, gitignored
    for d in (store, tmp):
        d.mkdir(parents=True, exist_ok=False)
    fixtures.mkdir(parents=True, exist_ok=True)
    url = f"http://127.0.0.1:{port}/"

    # The race fixture: KY160 twice subdivided, so the same 5 openings and a build
    # slow enough (19 s) for the next visitor to be on the page when it finishes.
    race = fixtures / RACE_NAME
    if not race.exists():
        from branch_cut import load_mat_mesh
        load_mat_mesh(SMALL_MESH).subdivide().subdivide().export(str(race))
    if race.stat().st_size != 84 + 50 * RACE_FACES:
        sys.exit(f"{race} is not the {RACE_FACES:,}-face binary STL this check expects "
                 f"(a partial write?). Not overwriting it; move it aside by hand.")

    print("\ncutting the samples headlessly for the expected results and timings...")
    t = time.time()
    exp_small = reference(SMALL_CL, SMALL_MESH)
    t_small = time.time() - t
    t = time.time()
    exp_race = reference(SMALL_CL, race)
    t_race = time.time() - t
    print(f"  KY160_1 : {exp_small['message']}  built in {t_small:.1f}s")
    print(f"  race    : {exp_race['message']}  built in {t_race:.1f}s")

    if port_in_use(port):
        sys.exit(f"Something is already listening on {port}. Not touching it.")
    # The E3 gap has nobody connected on purpose; a slow run must not meet the idle exit.
    os.environ["VC_IDLE_TIMEOUT"] = "900"
    proc, log = start_app(Path(args.app), port, run, store, tmp)
    admitted_pages, alive, stopped = [], None, False
    try:
        if not wait_listening(proc, port, run):
            log.flush()
            print((run / "server.log").read_text(encoding="utf-8", errors="replace")[-3000:])
            sys.exit("The app never started listening.")
        with sync_playwright() as p:
            b = p.chromium.launch()
            req = p.request.new_context(base_url=url)
            try:
                flow(b, req, url, port, run, store, exp_small, t_small, t_race, race,
                     admitted_pages, args.with_shape)
            except Stop as e:
                print(f"\n  STOPPED: {e}. The checks after this point did not run.")
            alive = proc.poll() is None
            b.close()
    finally:
        stopped = stop_app(proc, port)
        log.close()
        print(f"\nserver stopped, port {port} free: {stopped}")

    check("ERR", "PAGE: no page errors, no traceback, the server outlived every visitor")
    errs = [f"{v.label}: {e}" for v in admitted_pages for e in v.errors]
    for e in errs[:6]:
        info(f"page error: {e[:140]}")
    note(not errs, f"no page errors on any admitted page ({len(errs)})")
    text = (run / "server.log").read_text(encoding="utf-8", errors="replace")
    blocks, cur = [], None                     # as check_upload.py:907-925
    for ln in text.splitlines():
        if ln.startswith("Traceback"):
            cur = [ln]
        elif cur is not None:
            cur.append(ln)
            if ln and not ln[0].isspace():
                blocks.append(cur)
                cur = None
    real = [b_ for b_ in blocks if not ("_call_connection_lost" in "\n".join(b_)
                                        and b_[-1].startswith("ConnectionResetError"))]
    for b_ in real:
        info(f"server log: {b_[-1][:140]}")
    note(not real, f"no traceback in the server log ({len(real)})")
    note(alive is True, "the server was still running after the last visitor left")
    note(stopped, f"the server was stopped and port {port} released")

    print("\n" + "=" * 74 + "\nverdicts\n" + "=" * 74)
    fails = []
    for cid in ORDER:
        bad = [w for ok, w in CHECKS[cid]["notes"] if not ok]
        print(f"  {cid:<3} {verdict(cid)}  {CHECKS[cid]['title']}")
        if bad:
            print(f"        first failure: {bad[0]}")
            fails.append(cid)
    print(f"\n  run folder -> {run}\n")
    print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): {', '.join(fails)}")
    sys.exit(1 if fails else 0)


def flow(b, req, url, port, run, store, exp_small, t_small, t_race, race,
         admitted_pages, with_shape):
    def shot(v, tag):
        v.page.screenshot(path=str(run / f"{v.label}_{tag}.png"))

    # ── 0 the first visitor to a fresh process defines "fresh" ──────────────
    check("0", "SETUP: a fresh visitor, and probes shown to see")
    A = Visitor(b, "A", url)
    note(A.go(), "A is admitted to a fresh process")
    if A.t_admitted is None:
        raise Stop("the first visitor was not admitted")
    admitted_pages.append(A)
    info(f"build stamp {BUILD!r} on screen: {BUILD in A.page.inner_text('body')}")
    set_section(A.page, "cut")
    fresh = {"state": A.state() or {}, "status": page_state(A.page)["status"]}
    fresh_blob = dom_blob(A.page)
    tokA = fresh["state"].get("vc_token") or ""
    note(bool(tokA) and len(fresh["state"]) > 50,
         f"A holds an upload token, and {len(fresh['state'])} state keys were read")
    note_view_empty(A.page, "a fresh process's view")
    shot(A, "fresh")
    open_panel(A.page, "Input")
    do_upload(A.page, SMALL_MESH, SMALL_CL, exp_small, "A: KY160")
    sidA = (A.state() or {}).get("submission_id") or ""
    note(bool(sidA), f"A has a submission id ({sidA[:26]})")
    # A hostile page. trame accepts any key from an admitted page, and the next page
    # is BUILT from trame's own keys -- trame__template_main, trame__scripts,
    # trame__vue_use run in the next visitor's browser. Written here: trame__favicon
    # (the same kind of key, harmless to A's own page) and a key A invents. The server
    # must refuse both; anything of either on a later visitor's socket is a leak.
    log_at, err_at = log_tail(run, 0)[1], len(A.errors)
    A.page.evaluate("() => { window.trame.state.set('trame__favicon', 'wfPLANTfav');"
                    " window.trame.state.set('wfPLANTkey', 'wfPLANTval'); }")
    A.page.wait_for_timeout(1500)
    note("refused a page write" in log_tail(run, log_at)[0],
         "the server refused A's writes to a trame key and an invented one")
    # trame's client acts on its own trame__ keys locally, and setting one by hand
    # can throw in A's OWN page ("n[i[e]] is not a function" on the first run). That
    # is the hostile page breaking itself, not the app, so it is reported here rather
    # than counted by ERR.
    for e in A.errors[err_at:]:
        info(f"A's page after its own hostile write: {e[:100]}")
    del A.errors[err_at:]
    vA = settled_view(A.page)
    note(vA is not None and vA["out"] >= LOADED_MIN,
         f"the view probe sees A's model ({100 * (vA or {}).get('out', 0):.2f}% ink; "
         f"empty is <= {100 * EMPTY_MAX:.2f}%)")
    shot(A, "loaded")
    # Things of A's that must not outlive A: a typed name, zone colours, a delta K
    # source, and a minted download link that is never followed.
    open_panel(A.page, "Output")
    box = A.page.locator(".v-expansion-panel", has=A.page.locator(
        ".v-expansion-panel-title", has_text=re.compile(r"^\s*Output\s*$", re.I))
    ).locator("input[type=text]")
    if box.count():
        box.first.fill(f"{SENT_A}_export")
        box.first.press("Tab")
        A.page.wait_for_timeout(500)
    note((A.state() or {}).get("export_name") == f"{SENT_A}_export",
         "A's export name carries A's sentinel")
    mint = "async () => String(await window.trame.trigger('download_cut'))"
    href1, href2 = A.page.evaluate(mint), A.page.evaluate(mint)
    r = req.get(href1)
    note(href1.startswith("vc/download/") and r.status == 200,
         f"a minted link works while A is here ({r.status}); a second is kept unfollowed")
    zb = A.page.get_by_role("button", name=re.compile(r"^\s*Find zones", re.I))
    if zb.count():
        zb.first.click()
        for _ in range(40):
            if (A.state() or {}).get("zone_legend"):
                break
            A.page.wait_for_timeout(500)
    info(f"A's zone legend: {len((A.state() or {}).get('zone_legend') or [])} zones")
    needles = [sidA, tokA, "KY160", "CentrelinePoints", SENT_A, "wfPLANT"]
    blind = [n for n in needles if n.lower() in fresh_blob.lower()]
    note(not blind, f"none of A's needles is on a fresh page, so finding one is a leak {blind}")

    # ── A: a second browser while A is here ─────────────────────────────────
    check("A", "REFUSED: a second browser gets nothing while A is admitted")
    Bs = []
    for k in range(3):                      # a second person, reloading impatiently
        B = Visitor(b, f"B{k}", url)
        adm = B.go(wait_s=12)
        Bs.append(B)
        msgs = decode(B.frames)
        kinds = [kind(m) for m in msgs]
        info(f"B{k}: {len(B.frames)} frames, {sum(map(len, B.frames)):,} bytes: {kinds[:6]}")
        (run / f"frames_B{k}.json").write_text(json.dumps(kinds), encoding="utf-8")
        note(not adm, f"B{k}: not admitted")
        note("STATE" not in kinds, f"B{k}: no state on its socket")
        note(not any(x.startswith("publish:") for x in kinds),
             f"B{k}: no published update and no picture on its socket")
        note(not needle_hits(B.frames, needles), f"B{k}: none of A's needles on its socket")
        blob = B.page.inner_text("body") + json.dumps(B.state() or {}, default=str)
        note(not [n for n in needles if n.lower() in blob.lower()],
             f"B{k}: nothing of A's in its page or its client state")
        note(not needle_hits(B.http_bodies(), needles),
             f"B{k}: nothing of A's in any HTTP response it got")
        vb = view(B.page)
        note(vb is None or vb["out"] <= EMPTY_MAX, f"B{k}: no picture of A's model")
        shot(B, "refused")
    check("A2", "TOLD: the refused browser is told the tool is in use (policy)")
    for B in Bs:
        body = B.page.inner_text("body")
        note(REFUSED_RE.search(body), f"{B.label}: says so ({body.strip()[:70]!r})")
        B.close()

    check("A3", "HOSTILE: a raw websocket that skips the page is refused too")
    hbox = {}
    th = threading.Thread(target=lambda: hbox.update(hostile(port)), daemon=True)
    th.start()
    A.page.wait_for_timeout(3000)
    # Make A's session publish while the raw socket listens: a state push and a render.
    A.page.evaluate("() => window.trame.state.set('show_planes', false)")
    A.page.wait_for_timeout(1500)
    A.page.evaluate("() => window.trame.state.set('show_planes', true)")
    th.join(40)
    hk = [kind(m) for m in hbox.get("msgs", [])]
    info(f"upgrade {hbox.get('upgrade')}, {len(hbox.get('raw', []))} frames: {hk[:6]}")
    note(hbox.get("upgrade") is not None, "the probe reached the server")
    note(hbox.get("state_get") != "ANSWERED", f"trame.state.get was not answered "
                                              f"({hbox.get('state_get')})")
    note("STATE" not in hk and not any(x.startswith("publish:") for x in hk),
         "no state and no publish reached it while A was working")
    note(not needle_hits(hbox.get("raw", []), needles), "none of A's needles reached it")

    check("B", "UNHARMED: A's page survives the refused browsers and the raw socket")
    A.page.wait_for_timeout(5000)     # a wipe keyed to ANY socket closing has run by now
    s = A.state() or {}
    note(s.get("submission_id") == sidA, "A's submission id is unchanged")
    note(s.get("vc_token") == tokA, "A's token is unchanged: a refusal must not rotate it")
    ps = page_state(A.page)
    note(ps["cards"] == exp_small["n_openings"] and sorted(ps["names"]) == exp_small["names"],
         f"A's {exp_small['n_openings']} opening rows are all there")
    note(re.match(r"ready", ps["status"], re.I), f"A's status still reads ready")
    note(not re.search(r"connection (closed|error)", A.page.inner_text("body"), re.I),
         "A's page was not replaced by trame's connection-error template")
    v = settled_view(A.page)
    note(v is not None and v["out"] >= LOADED_MIN, "A's view still shows A's model")
    n0 = len(A.downloads)
    open_panel(A.page, "Output")
    A.page.get_by_role("button", name=re.compile(r"download", re.I)).first.click()
    for _ in range(60):
        if len(A.downloads) > n0:
            break
        A.page.wait_for_timeout(500)
    note(len(A.downloads) > n0, "A can still download its cut (session, submission, "
                                "routes all intact)")
    shot(A, "after_B")

    # ── C: A closes the tab, someone else arrives ───────────────────────────
    check("C", "CLEAN: after A closes its tab, the next visitor starts from nothing")
    A.close()
    C = must_admit(b, "C", url, after=A)
    admitted_pages.append(C)
    tokens = {tokA}
    tokC = assert_clean(C, fresh, needles, tokens)
    tokens.add(tokC)
    shot(C, "clean")

    # ── D: A's leftovers are dead ────────────────────────────────────────────
    check("D", "STALE: A's token and A's link are dead once A has gone")
    before = listing(store, recursive=True)
    r1 = post(req, "surface", tokA, SMALL_MESH, f"{SENT_A}_stale.mat")
    r2 = post(req, "centerline", tokA, SMALL_CL, f"{SENT_A}_stale.csv")
    note(r1.status != 200 and r2.status != 200,
         f"A's token is refused for both kinds ({r1.status}, {r2.status}; 403 today "
         f"means a stale token)")
    note(post(req, "surface", "", SMALL_MESH, "x.mat").status != 200, "no token: refused")
    C.page.wait_for_timeout(1500)
    s = C.state() or {}
    note(not s.get("up_mesh_note") and not s.get("up_cl_note"), "C's upload notes untouched")
    note(listing(store, recursive=True) == before, "nothing was stored")
    # The pairing attack: C sends half a pair with ITS token; A's token sends the rest.
    upload_one(C.page, "centerline", SMALL_CL)
    got = wait_status(C.page, r"now add the surface", 30)
    note(got is not None, "C's own token works: its centerline arrived, surface asked for")
    r3 = post(req, "surface", tokA, SMALL_MESH, f"{SENT_A}_pair.mat")
    note(r3.status != 200, f"A's token cannot send C the other half ({r3.status})")
    C.page.wait_for_timeout(8000)
    s, ps = C.state() or {}, page_state(C.page)
    note(not s.get("loaded") and not s.get("busy") and not s.get("submission_id")
         and not re.match(r"loading", ps["status"], re.I),
         f"C's half pair was not completed: nothing loads ({ps['status'][:50]!r})")
    note(len(listing(store) - {x.split(os.sep)[0] for x in before}) == 0,
         "no submission was created")
    r4 = req.get(href2)
    note(r4.status == 404, f"A's unfollowed download link is gone ({r4.status})")

    # ── E0: C leaves with a centerline pending ───────────────────────────────
    check("E0", "HALF-PAIR: a file left pending does not wait for the next visitor")
    C.close()
    E = must_admit(b, "E", url, after=C)
    admitted_pages.append(E)
    tokE = assert_clean(E, fresh, needles + [tokC], tokens)
    tokens.add(tokE)
    upload_one(E.page, "surface", race)       # 30 MB; E2 needs it loading anyway
    got = wait_status(E.page, r"now add the centerline|loading", 90)
    note(got is not None and re.search(r"now add the centerline", got["status"], re.I),
         f"E's surface arrived and E is asked for a centerline: C's did not carry over "
         f"({(got or {}).get('status', '')[:60]!r})")

    # ── E2: E leaves while its 19 s build runs; W is there when it finishes ─
    check("E2", "MID-LOAD: the next visitor is present when a departed visitor's "
                "build finishes, and never sees it")
    upload_one(E.page, "centerline", SMALL_CL)
    t0, sidE = time.time(), ""
    while time.time() - t0 < 120 and not sidE:  # stored, so the build itself is running
        sidE = (E.state() or {}).get("submission_id") or ""
        E.page.wait_for_timeout(200)
    loaded_at_close = bool((E.state() or {}).get("loaded"))
    note(bool(sidE) and not loaded_at_close,
         "E left while its build was running (submission stored, case not yet shown)")
    log_mark = log_tail(run, 0)[1]
    E.close()
    W = must_admit(b, "W", url, after=E)
    admitted_pages.append(W)
    lead = W.t_admitted - E.t_closed
    note(lead < 0.6 * t_race, f"W was on the page {lead:.1f}s after E left, well inside "
                              f"E's {t_race:.0f}s build, so the race was exercised")
    needles += [SENT_E, sidE, tokE]
    dirty, n = [], 0
    while time.time() < E.t_closed + 1.5 * t_race + 20:
        s = W.state() or {}
        if (s.get("loaded") or s.get("n_openings") or s.get("submission_id")
                or s.get("openings") or re.search(r"loading|ready", str(s.get("status_text")), re.I)):
            dirty.append((round(time.time() - E.t_closed, 1), s.get("status_text")))
        n += 1
        W.page.wait_for_timeout(1000)
    for d in dirty[:3]:
        info(f"W at +{d[0]}s: {str(d[1])[:80]!r}")
    note(not dirty, f"W's state stayed clean through the build's finish ({n} samples)")
    tail, _ = log_tail(run, log_mark)
    for ln in [x for x in tail.splitlines() if x.startswith(("[visit]", "[uploads]"))][:6]:
        info(f"server: {ln[:120]}")
    tokW = assert_clean(W, fresh, needles, tokens)
    tokens.add(tokW)
    info(f"E's submission kept in the store: {any(sidE in x for x in listing(store))}")

    # ── E1: an upload still streaming when its visitor leaves ───────────────
    check("E1", "IN-FLIGHT: an upload arriving after its visitor left reaches nobody")
    pbox = {}
    th = threading.Thread(target=slow_post, daemon=True, args=(
        port, "surface", tokW, SMALL_MESH, f"{SENT_E}_inflight.mat", 150_000, pbox))
    th.start()
    while pbox.get("sent", 0) < 0.25 * SMALL_MESH.stat().st_size:
        time.sleep(0.1)
    W.close()                                  # the body is still arriving
    L = must_admit(b, "L", url, after=W)
    admitted_pages.append(L)
    note(pbox.get("sent", 0) < pbox.get("total", 1),
         "L was admitted while W's upload was still arriving")
    th.join(60)
    note(not th.is_alive() and pbox.get("status") not in (None, 200),
         f"the upload that finished after W left was refused ({pbox.get('status')}, "
         f"{pbox.get('body', b'')[:60]!r})")
    tokL = assert_clean(L, fresh, needles, tokens)
    tokens.add(tokL)
    upload_one(L.page, "centerline", SMALL_CL)
    got = wait_status(L.page, r"now add the surface|loading", 30)
    note(got is not None and re.search(r"now add the surface", got["status"], re.I),
         "L's centerline alone does not start a load: W's surface is not pending")

    # ── E3: a short build finishes while nobody is connected ─────────────────
    check("E3", "LATE-LOAD: a build that finishes with nobody connected is not kept "
                "for the next visitor")
    upload_one(L.page, "surface", SMALL_MESH)
    got = wait_status(L.page, r"^loading", 30)
    loaded_at_close = bool((L.state() or {}).get("loaded"))
    L.close()
    note(got is not None and not loaded_at_close,
         "L left while its case was loading")
    time.sleep(3 * t_small + 8)                # the build finishes to an empty room
    M = must_admit(b, "M", url, after=L)
    admitted_pages.append(M)
    tokM = assert_clean(M, fresh, needles, tokens)
    tokens.add(tokM)

    # ── F: reload ────────────────────────────────────────────────────────────
    check("F", "RELOAD: a visitor who reloads is admitted and starts blank (documented)")
    open_panel(M.page, "Input")
    do_upload(M.page, SMALL_MESH, SMALL_CL, exp_small, "M: KY160")
    sidM = (M.state() or {}).get("submission_id") or ""
    mark = len(M.frames)
    M.page.reload(wait_until="domcontentloaded", timeout=60000)   # beforeunload runs
    ok = M.wait_admitted(ADMIT_S + 10)
    body = M.page.inner_text("body")
    note(ok and not REFUSED_RE.search(body),
         "the reloaded page is admitted, not refused as someone else")
    if not ok:
        raise Stop("the reloading visitor was shut out")
    assert_clean(M, fresh, needles + [sidM, tokM], tokens, frames_from=mark)
    shot(M, "after_reload")

    if with_shape:
        check("E4", "JOB: a delta K job running at disconnect never reports")
        do_upload(M.page, SMALL_MESH, SMALL_CL, exp_small, "M: KY160 again")
        M.page.evaluate("() => window.trame.state.set('sm_faces', '25000')")
        open_panel(M.page, "Output")
        # Not /measure/ alone: the cut section's "Measure" button (cut_app.py:2195)
        # would match first. The handoff (2306-2309) switches to the shape section.
        M.page.get_by_role("button", name=re.compile(r"measure\s+δk\s+on\s+this\s+cut",
                                                      re.I)).first.click()
        M.page.wait_for_timeout(1500)
        M.page.get_by_role("button", name=re.compile(r"compute\s+δk", re.I)).first.click()
        for _ in range(60):
            if (M.state() or {}).get("sm_busy"):
                break
            M.page.wait_for_timeout(250)
        note(bool((M.state() or {}).get("sm_busy")), "M's delta K job is running")
        M.close()
        N = must_admit(b, "N", url, after=M)
        admitted_pages.append(N)
        late = []
        while time.time() < M.t_closed + 120:
            s = N.state() or {}
            if s.get("sm_result") or s.get("sm_busy") or s.get("sm_info"):
                late.append(round(time.time() - M.t_closed, 1))
            N.page.wait_for_timeout(2000)
        note(not late, f"nothing of M's job reached N in two minutes {late[:3]}")
        assert_clean(N, fresh, needles, tokens)


if __name__ == "__main__":
    main()
