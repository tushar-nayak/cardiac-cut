# Getting the web app ready to host publicly

Instructions for work on `webapp/` and `vascular_api.py`. The cutting geometry
(`branch_cut.py`) is being worked on elsewhere -- **do not change it here**, or the two
efforts will collide. If something in this task seems to need a change to the cutter,
write down what and why and leave it.

## What this is

A tool that turns a closed vascular surface mesh into an OPEN one with clean planar
inlets and outlets, for geometric measurement. Two front ends over the same geometry:

| | |
|---|---|
| `branch_gui.py` | the desktop tool (Qt + PyVista). The reference implementation |
| `webapp/cut_app.py` | the web tool (trame). Same features, drives the same cutter |
| `vascular_api.py` | the middle layer with no Qt in it, which the web tool uses |
| `branch_cut.py` | the geometry. **Off limits for this task** |

The desktop tool is what the lab uses. The web tool is for a paper: readers upload a
surface and a centerline, get a cut mesh back. It is at feature parity with the desktop
as of `8f9f354` -- cutting, trimming a branch group, isolating a region, the measure
plane, branch colouring, click-to-pin marks, the pull-back button.

## Run it

```powershell
.venv\Scripts\python branch_gui.py KY160_1_CentrelinePoints.csv "Geometries & Meshes/mesh/KY160_1.mat"
$env:VC_CUT_PORT=8099; .venv\Scripts\python webapp\cut_app.py --server   # then :8099
```

`VC_CUT_PORT` sets the port -- there is no `--port` flag, and 8080 is often already
taken. trame holds ONE SHARED STATE PER PROCESS -- every tab on it sees the same case --
and in public mode the app admits one visitor per process (job 5), so a second copy
needs a second port.

## The goal

Public hosting on a small CPU-only container. No GPU: the scene is ~37k faces and a few
discs, which software rendering handles.

> **δK is NOT dropped, and job 2 below was never built.** Corrected 28 Aug 2026, after
> a parity audit found this document contradicting the code. `VC_PUBLIC` exists
> (`webapp/cut_app.py`) but gates only path loading and the path widgets; the shape
> panel, `sm_use_cut` and `sm_run` are ungated, and CAMP is imported by every session.
> What was done instead: the partition count now defaults to **1** for every session
> of `cut_app`, with the resulting **±16%** printed beside the value everywhere it
> appears. That makes δK affordable rather than absent. Job 2 stands as an option, not
> as a description of the build -- decide which before hosting, but do not read the
> rest of this file as if δK were already gone.

## The jobs, in the order they matter

### 1. Uploads instead of file paths -- DONE 24 Sep 2026, with job 7

**Status: built and verified in a real browser.** `tests/web/check_upload.py` starts the
app in public mode and drives it: 62 checks, ALL PASS. Two file pickers (surface,
centerline) replace the path boxes in public mode, the case loads once both arrive, and
a Download button hands back one zip of the STL and its `_cut_planes.json`. The typed
path stays on loopback, where it is how a researcher points at a lab share; public mode
does not render it, and `load_model` still refuses it server-side. Storage is
`webapp/submissions.py` (`tests/test_submissions.py`, 127 checks). Things the work found
that were not on this list, all fixed:

- **One 13-byte message could make the server reserve 4 GB.** After the handshake
  wslink raises every client's reassembly limit to 4 GB and allocates the claimed size
  from the first chunk's header. cut_app now clamps that limit to the largest accepted
  upload; the page refuses an oversize file before sending, so only a client that
  ignored the page is cut off. The clamp patches a wslink internal: **pin wslink**.
- **trimesh opened files named inside an uploaded mesh.** Given a path it resolves an
  OBJ's `mtllib` and a PLY's `TextureFile` against the mesh's folder with no traversal
  check. `vascular_api._load` now reads non-.mat meshes from bytes, which opens nothing;
  public mode still refuses .obj and .ply.
- **STL could not be cut at all.** Read unprocessed, every triangle's corners are
  separate vertices -- KY160 came out as 39,802 "openings". Merged now, on the web path.
  **`branch_gui.py` has the same bug** and is frozen, so the desktop still has it.
- **Public mode sent this machine's absolute paths** (account name included) to every
  visitor in the first state message. Blanked when public.
- Re-picking the same file was silently ignored, and `accept=` on VFileInput is dropped
  without a word unless passed through `__properties` -- both handled.

**The freeze is fixed (24 Sep 2026).** Uploads used to travel the websocket, which
carries everything in order, so a large file froze the whole page for as long as it took
to send; and the case was then built on the event loop, freezing it again. Now each file
goes over HTTP beside the socket (`webapp/static/vc_upload.js`, and the routes cut_app
adds through trame's aiohttp-only `on_server_bind` hook), with progress shown, and the
session is built in a worker thread while only the drawing stays on the loop. Downloads
are a one-time link rather than bytes over the socket. Measured by toggling a view
setting and timing the 3-D view's response while an 8 MB, 159k-face mesh uploaded at
1 MiB/s and loaded: **19.4 s** worst before, **0.36 s** after. With uploads off the
socket, its reassembly clamp is now 8 MB rather than the upload size.

Two things this does not change. A `.vtp/.vtk/.vtu` CENTERLINE is read through VTK, and
VTK work off the render thread has crashed this process before, so that case is still
built on the loop. And every recut after the first -- each slider release -- still runs
on the loop, so a large mesh is slow to work with once loaded: about 6-11 s per 100k
faces. The delta K panel still takes its surface only from the cut section when public.

The original brief, kept for reference:

The app currently takes a PATH typed into a text field and the server reads it
(`state.csv_path`, `state.mat_path`, `load_model`). Public, that is arbitrary file read
on the host. It must become a file upload before the app is reachable from anywhere.

Uploads are kept -- see job 7 -- so build this and the storage together rather than
bolting retention on afterwards: the upload handler is where a submission gets its id,
its stored copy and its recorded original filename.

Accepted centerline formats are wider than they look -- `load_spline_csv` reads
delimited tables with coordinates under any usual header spelling (or none, taking the
first three numeric columns), `.vtp/.vtk/.vtu` polylines as VMTK writes them, 3D Slicer
markup JSON, and plain JSON point lists. Say so on the page: a reader needs a surface
AND a centerline, and the tool cannot derive one from the other.

### 2. A flag that turns δK off — NOT BUILT, and possibly not wanted

**Status as of 28 Aug 2026: not implemented.** `VC_PUBLIC` exists but does not do this.
The cost was addressed a different way -- one partition by default, with the ±16% band
printed beside the value -- so a hosted δK is roughly a tenth of what the figures below
assume for an 8-partition run. Keep this job if the import cost alone is unacceptable;
drop it if a cheap, honestly-labelled δK is worth having on the page. Do not assume it
is done.

Were it built: something like `VC_PUBLIC=1` that hides the shape section, refuses its
triggers, and **never imports CAMP**, so the import cost is not paid by every session.
Not deleting it -- the desktop keeps it. Measured on KY160:

| | time | peak memory |
|---|---|---|
| cutting (load + cut, 37,786 faces) | 4.8 s | 287 MB |
| δK (remesh to 100k + CAMP partition) | 15.4 s | **989 MB** |

δK is also 15 seconds of pegged CPU per press, which is the obvious way to knock a
public box over.

### 3. Bound the re-cut cache

`CutSession.recut` keeps the mesh after each cut so a slider move resumes rather than
re-cutting everything (`_cut_steps`). It is worth ~5x on an 11-opening model and it
costs a lot of memory, because it holds live `Trimesh` objects whose cached edge and
adjacency tables balloon 0.7 MB of geometry into far more:

    python baseline                 25 MB
    + vascular_api imports         149 MB
    + one model loaded and cut     290 MB
    + re-cut cache warm            645 MB
    the running web app            854 MB   (with rendering, one session)

Cache the raw arrays and rebuild on demand, or bound it to the last few steps. Measure
before and after; keep as much of the speedup as you can. This decides whether the
instance is 2 GB or 4 GB.

### 4. Limits

An upload size cap and a face-count ceiling. Cutting a two-million-face mesh is still a
way to hurt the box even with δK gone.

**Done (24 Sep 2026).** Two limits, and they protect different things.

The **face-count limit** is 1,000,000 (`VC_MAX_FACES`), checked in
`SubmissionStore.check_upload` BEFORE the mesh is built, by `webapp/meshcount.py`, which
reads the count from the file's own structure for every accepted format. It was set from
evidence rather than guessed. This lab's cohort runs 14k-74k faces; the Vascular Model
Repository's 131 aortic surfaces run median 91k, 95th percentile 267k, 99th 529k, max
666k. A limit near this lab's data would have turned away 43% of those; a million turns
away none, and 2 of all 254 of its models. What it costs, measured on 25 Sep 2026 on
the development workstation after that day's speed fix (commit 0276d58): 147k-159k
faces load in 4-6 s and take 8-12 s for the nub-in cut; 589k-637k faces, 18-30 s and
31-48 s. Nothing has been timed at the 1M limit, and an earlier straight-line guess of
"35 s to load, 2 min to cut" predates the fix. Memory is unsettled: 0.6 GB per 100k
faces was estimated here, 0.93 GB was measured at 637k on the ui-redesign branch.
**Measure on the real host, and size it for the limit.**

`scipy.io.whosmat` turned out NOT to be safe for this: it decompresses each variable to
report its shape, 122 MB for a 0.18 MB file. meshcount decompresses in a stream, keeps
only each variable's header, and stops at 1 GB -- 1.9 MB of memory for the same file.
Headers that lie are checked against the file: a binary STL must be exactly 84 + 50n
bytes, and a PLY's declared faces must fit in its bytes. `tests/test_meshcount.py`
confirms the count equals the real loader on every format and all 32 .mat files.

The **byte cap** is still 50 MB for a surface and 10 MB for a centerline
(`VC_MAX_MESH_MB`, `VC_MAX_CENTERLINE_MB`). It bounds the transfer, not the load. Note it
binds before the face limit for ASCII STL, which is about five times larger than binary:
50 MB of ASCII STL is roughly 190k faces. The refusal says to save as binary STL, which
reaches the face limit at 50 MB. Raise the byte cap if ASCII uploads matter; HTTP upload
makes a larger one practical.

### 5. Session isolation

**Privacy between visitors: built 25 Sep 2026. Serving several at once: still open.**

trame holds one state per process and pushes every change to every browser on it, and
a browser that connects is handed the whole current state first. So two visitors on
one process shared a case, a cut, a camera and a submission id; one arriving within
the idle shutdown after another had left was handed the last person's case; and two
uploads at once could be stored as one person's surface with another's centerline.

In public mode the app now fixes that itself, whatever the hosting does
(`webapp/visit_gate.py`, the wipe in `webapp/cut_app.py`):

- **One visitor per process.** Everyone else gets a static "This tool is in use"
  page and no websocket -- refused before they receive a byte of state. A raw client
  that skips the page and presents trame's default secret is refused at the hello.
- **Wiped the moment they leave** -- tab closed, reload, crash, or a dead network once
  the heartbeat gives up (~45 s at trame's default 30 s heartbeat): every state key back
  to a fresh process's, the session, pending uploads, download links, zone map and
  delta K job emptied, the scene cleared and re-rendered, the camera reset.
- **Work still running for a visitor who has gone is dropped**, not installed for the
  next: an upload still streaming, a session still building in a thread, a zone map, a
  delta K run. A submission stored after its visitor left is logged by id only.
- **The upload token is per visit**, withdrawn when they leave.
- **An idle tab is closed** after `VC_VISIT_IDLE_MIN` minutes with no message from it
  (default 30; 0 turns it off), since an open tab never lets the idle shutdown run.
- **wslink's `application.exit` RPCs are removed** in public mode: any admitted
  client could otherwise stop the process.
- **A page may write only the app's own state keys.** trame accepts ANY key from an
  admitted page, including its own `trame__template_main`, `trame__scripts` and
  `trame__vue_use` -- and the next visitor's page is built from those and runs them.
  Found by an adversarial review before this shipped: one visitor could have planted
  script in the page the next one uploads through. Writes to any other key are
  dropped and logged, and the wipe restores EVERY key, trame's included, from the
  state as it stood before the first visitor, so a write that got past the filter
  still would not outlive its visit. Lifecycle calls are limited to the three trame's
  page sends, and an RPC that raises no longer sends the page the server's traceback.
- **Two uploads in one visit cannot cross.** Only the newest pair's load can reach
  the screen, and the submission id a download is filed under changes only together
  with the case it belongs to.
- **delta K run files are named per run**, not per setting, so a crashed run can
  never pick up an earlier visitor's result or error text; and a worker that will not
  stop makes the wipe fail closed.

Consequences to know: a **reload starts over** (the page says to note the submission
id first), and a person who opens a second tab is refused like anyone else. Local
runs are unchanged -- no gate, tabs share a state, as before.

**One visitor per process can be held on purpose.** A script that claims (POST
`paraview/`) again each time its 10 s ticket lapses, or that is admitted and keeps
sending messages, keeps everyone else out. That is inherent to one-visitor-per-process
with no login. Put a **rate limit on POST `paraview/` per client address at the
proxy**, and serve more than one process.

Tested by `tests/test_visit_gate.py` (the gate's rules), `tests/test_visit_wipe.py`
(the wipe, on the real app, no browser) and `tests/web/check_privacy.py` (browsers:
refused, told, a raw socket, the next visitor clean, stale token and link dead, an
upload and a load finishing after their visitor left, reload).

**Still open: capacity.** One visitor per process means concurrent visitors = number
of processes. Serving more than one at a time needs something that starts a process
per visitor and routes them to it (trame's launcher, or a container per session;
Kitware publish deployment guidance and images for this). The gate makes a reused
process safe, so a small pool of long-lived processes also works. The claim URL the
page opens is `ws://USE_HOST/ws`; under a path prefix set `VC_SESSION_URL` (lower case,
e.g. `ws://USE_HOST/vc/ws`) -- untested behind a proxy.

### 6. A Dockerfile

Headless OpenGL (OSMesa or EGL, or xvfb). A plain `pip install` in a bare container
renders nothing. Preload a sample geometry so the tool demonstrates itself to someone
who arrives without files.

### 7. Keep what is uploaded, and say so before it is

**Built 24 Sep 2026; the wording is not decided.** Uploads are stored under
`VC_STORE_DIR` as `<id>/surface.<ext>`, `<id>/centerline.<ext>` and a `submission.json`,
with each downloaded cut kept beside them. The id is a UTC time plus 128 random bits;
the original filename is metadata only and never a path; nothing is overwritten,
moved or deleted, and a collision fails rather than clobbers. Public mode refuses to
start without `VC_STORE_DIR`. The page shows the submission id after upload, because
it is the only thing linking a person to their files for a removal request.

The notice above the pickers is written, in `webapp/submissions.py`. Of its six
policy blanks, four were **decided on 25 Sep 2026**: the files are used to run and
improve the tool **and for research** on aortic geometry; they are kept **until the end
of the project, and no longer than five years**; **the holding group and its
collaborators** can see them; and uploading is **a condition of use**, not an opt-in.
**Two are still placeholders:** which institution or group holds the files, and a
monitored address for removal requests. Until they are filled, **the app refuses to
bind any non-loopback address** -- placeholder text on a public page would read as
policy. Two consequences of the decisions: research use makes institutional sign-off
(below) a real prerequisite, and nothing deletes on a schedule yet, so the five-year cap
is kept by hand until a clean-up job is approved and built. Each submission records a fingerprint of the notice it was
shown and whether that was still a draft.

Two things the notice now also tells people, found while building it: identifiers sit
INSIDE centerline files as well as in their names -- a Scan-IP export's first line can
record the project's full folder path and the export date, and its Spline name column
names each curve -- and the store keeps the bytes exactly as sent, so neither is
scrubbed. Whether to scrub on arrival is open; doing it means the stored file no longer
matches what was uploaded.

**Decided: uploads are retained.** Both files -- the surface and the centerline -- and
keep the cut mesh and its sidecar JSON with them, since the sidecar records the offsets
and modes that produced the cut and without it the pair is much less useful.

Container filesystems are ephemeral, so this needs a persistent volume or an
object store; writing into the image is not storage. Budget for it growing without
bound and decide a retention period rather than discovering one.

Store under generated IDs, with the original filename kept as metadata rather than as
the path. **Filenames are the likeliest place identifiers leak** -- people name exports
after the patient or the MRN far more often than they put anything in an STL -- and the
`Spline name` column of a Scan-IP centerline export is a second such place. Neither is
something the tool can reliably scrub, which is why the notice has to do some work.

**The page must tell users before they upload, not after**:

- that both files are kept, along with the cut produced from them
- what they will be used for
- how long they are held, and who can access them
- how to ask for removal, and a real address to ask at
- an instruction not to put names, MRNs or dates in filenames

Whether retention is a condition of use or an opt-in tick box is the PI's call. Opt-in
gets better data and is easier to get approved; a condition of use is simpler to build.
Either way the wording must be plain, above the upload control, and not buried in a
link nobody opens.

**Still a prerequisite, not an implementation detail:** collecting patient-derived data
from users at other institutions needs institutional sign-off, and that is the PI's to
obtain. This file records what to build, not that it may be switched on.

## Traps this codebase has already sprung

Every one of these was silent -- no console error, nothing in the log.

- **`classes=`, never `class_=`.** trame maps `classes` to the HTML class attribute and
  ignores `class_`. Whole stylesheets have quietly not applied.
- **A bare `>` inside a `v-if` closes the tag.** `v_if="openings.length > 1"` made the
  entire element vanish from the DOM. Compute the condition in Python and pass a flag.
- **Nested quotes in a Vue expression, or an apostrophe in a neighbouring `title=`,
  can take down the WHOLE panel body** -- every control in it, including ones that had
  worked for months. Compose text in Python; keep logic out of the template.
- **State read back inside the same handler is not reliably what you just wrote.**
  Pass the value to the helper instead of re-reading `state.x`.
- **Playwright's `inner_text()` returns RENDERED text**, and Vuetify uppercases button
  labels in CSS. Match case-insensitively or your check will report a missing button
  that is on screen.
- **Panels default closed.** A check that does not open the panel will report every
  control in it as missing. Match the title case-insensitively too ("View & measure" in
  the DOM, "VIEW & MEASURE" on screen).

## How to verify

**The 23 desktop suites must pass**, one per process -- two `QtInteractor` windows in
one interpreter segfault VTK:

```powershell
foreach ($t in Get-ChildItem tests\test_*.py) { .venv\Scripts\python $t.FullName }
```

They print `ALL PASS` or a list of failures, except `test_empty_region`, which prints
`NO BUG`. Some take `jk5` or `ky160`; see `tests/README.md`. These cover the desktop,
but they exercise `vascular_api`'s geometry too, so run them after any change here.

**Drive the browser, do not infer.** `tests/web/` has the pattern -- Playwright,
`VC_URL` overridable. Assert what a person would see, and take a screenshot and LOOK at
it; several bugs above passed every assertion.

**If you ever do touch `branch_cut.py`** (you should not), the cohort benchmark is
mandatory: `python tools/cut_bench.py run` over 66 real aortas, and `tools/CUTTING.md`
has the accept rules. The last full run at `86b0e7c` was byte-identical across 3,384
cuts, which is the standard for a change that is meant to be performance-only.

## Rules

- **Never commit to `main` directly.** Branch, then merge. (This was written after
  breaking it twice, which is how it got written down.)
- **Do not read anything outside this project folder** except the two cohort roots in
  `tools/cohorts.local.json`, and only meshes and centerlines within them. Imaging
  directories are off limits. If a case looks incomplete, say what is missing and stop.
- Never force-push, never `reset --hard` on `main`, never move the
  `cutter-shipped-2026-08-04` or `cutter-merged-2026-08-05` tags.

## Where performance stands

Desktop against web, same model, 5 openings, measured end to end:

| | desktop | web |
|---|---|---|
| commit a slider | 0.62 s | 0.63 s |
| Reset | 0.61 s | 0.71 s |
| live plane preview while dragging | 36 ms | 80 ms |

The 80 ms is the server round trip and is the floor for a server-rendered view. Moving
rendering into the browser (`VtkLocalView`, or `VtkRemoteLocalView`) would make camera
work instant and would cut hosting cost, since the server would stop rendering
entirely. It would NOT fix the 80 ms preview, which is dominated by computing the plane
pose rather than drawing it. Worth doing for the bill, not for the latency, and only
after the jobs above.
