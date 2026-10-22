"""Storage for uploads and the cuts made from them (submissions.py).

A plain script like the rest of tests/: prints each check, then ALL PASS or FAILURES,
and exits nonzero on failure.

What it proves, by running rather than by reading:
  * a hostile original filename cannot influence where anything is written -- every
    filesystem write the process makes is captured with an audit hook, not inferred
    from the directory afterwards;
  * nothing is ever overwritten, moved or deleted -- a forced ID collision fails and
    leaves the first submission byte-identical;
  * the extension allow-lists are exactly what the loaders dispatch on;
  * size limits refuse rather than truncate, and a refused upload writes nothing;
  * a real case (KY160) stored this way loads, cuts, and keeps its cut and sidecar
    beside it, with no original filename reaching the sidecar.

Writes ONLY under <this folder>/_out/submissions_<stamp>/, a new folder every run.
Deletes nothing.

    python test_submissions.py            everything, including the KY160 cut (~10 s)
    python test_submissions.py --no-cut   skip the cut (no geometry imports)
"""
import hashlib
import io
import json
import os
import re
import secrets
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in [HERE, *HERE.parents] if (p / "vascular_api.py").is_file())
for p in (HERE, ROOT, ROOT / "webapp"):             # the module ships in webapp/
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import submissions as S                                    # noqa: E402

DO_CUT = "--no-cut" not in sys.argv[1:]
RUN = HERE / "_out" / f"submissions_{time.strftime('%Y%m%dT%H%M%S')}_{secrets.token_hex(3)}"
RUN.mkdir(parents=True)                  # a NEW folder: fails rather than reuse one
print(f"writing under {RUN}")

fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}"
          + (f"\n         {detail}" if detail and not ok else ""))
    if not ok:
        fails.append((name, detail))


# ── every filesystem write, as the interpreter performs it ───────────────────
WRITE_FLAGS = (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND)
MUTATING = {"os.rename", "os.remove", "os.rmdir", "shutil.rmtree", "shutil.move",
            "os.truncate", "os.symlink", "os.link", "os.chmod", "os.utime",
            "shutil.copyfile", "shutil.copytree"}
WATCH = {"on": False, "events": []}


def _hook(event, args):
    if not WATCH["on"]:
        return
    if event == "open":
        path, mode, flags = args[0], args[1], args[2]
        if path is None or isinstance(path, int):
            return
        writes = (isinstance(mode, str) and any(c in mode for c in "wax+")) \
            or (isinstance(flags, int) and flags & WRITE_FLAGS)
        if writes:
            WATCH["events"].append(("write", os.path.abspath(os.fsdecode(path))))
    elif event == "os.mkdir":
        WATCH["events"].append(("mkdir", os.path.abspath(os.fsdecode(args[0]))))
    elif event in MUTATING:
        WATCH["events"].append((event, repr(args)[:200]))


sys.addaudithook(_hook)


def watched(fn, *a, **k):
    """(result, exception, events) of one call, with every write it made."""
    WATCH["events"] = []
    WATCH["on"] = True
    try:
        return fn(*a, **k), None, list(WATCH["events"])
    except Exception as e:                        # noqa: BLE001 -- recorded
        return None, e, list(WATCH["events"])
    finally:
        WATCH["on"] = False


def under(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def tree(root):
    root = Path(root)
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())


MESH = b"solid t\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\n" \
       b"vertex 0 1 0\nendloop\nendfacet\nendsolid t\n"
CL = b"x,y,z\n0,0,0\n0,0,1\n0,0,2\n"


def _tiny_mat():
    """A real one-triangle Simpleware-style .mat, built in memory."""
    import numpy as np
    import scipy.io as sio
    buf = io.BytesIO()
    sio.savemat(buf, {"X_surface_vertices": np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                                                     float),
                      "X_surface_triangles": np.array([[1, 2, 3]], dtype=np.uint64)})
    return buf.getvalue()


# A surface is now face-counted before it is accepted, so a filename check has to send
# a real mesh of the type its name claims. The type comes from the STORE's own reading
# of the name, so 'x<NUL>.mat' and '.mat' are paired with a .mat like the store expects.
TINY = {".mat": _tiny_mat(), ".stl": MESH,
        ".ply": (b"ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\n"
                 b"property float y\nproperty float z\nelement face 1\n"
                 b"property list uchar int vertex_indices\nend_header\n"
                 b"0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n"),
        ".obj": b"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n"}


def mesh_for(name):
    return TINY.get(S._clean_original_name(name)[2], MESH)

# ── 1. the allow-lists are what the loaders actually dispatch on ─────────────
print("\n1. extension allow-lists")
import branch_tree                                          # noqa: E402
import trimesh                                              # noqa: E402

reader = (set(branch_tree.TABLE_SUFFIXES) - {""}) | set(branch_tree.POLYDATA_SUFFIXES) \
    | {".json"}
check("centerline list == the suffixes load_spline_csv names (less '')",
      set(S.CENTERLINE_EXTENSIONS) == reader,
      f"list={sorted(S.CENTERLINE_EXTENSIONS)} reader={sorted(reader)}")
check("load_spline_csv still dispatches .json itself (source check)",
      'suf == ".json"' in (ROOT / "branch_tree.py").read_text(encoding="utf-8"))
gui = (ROOT / "branch_gui.py").read_text(encoding="utf-8")
m = re.search(r"Meshes \(([^)]*)\)", gui)
desk = {s.lstrip("*") for s in m.group(1).split()} if m else set()
check("mesh list == the desktop Open dialog's Meshes filter",
      set(S.MESH_EXTENSIONS) == desk, f"list={S.MESH_EXTENSIONS} desktop={desk}")
fmts = set(trimesh.available_formats())
check("every non-.mat mesh extension is one trimesh.load reads",
      all(e[1:] in fmts for e in S.MESH_EXTENSIONS if e != ".mat"))
va = (ROOT / "vascular_api.py").read_text(encoding="utf-8")
# The loader dispatches on suffix, which is why the store keeps it. It also now reads
# non-.mat meshes FROM BYTES rather than from the path: handed a path, trimesh opens
# files named inside an OBJ or PLY with no traversal check. So this checks both --
# that the dispatch is intact, and that the path form has not crept back in.
check("vascular_api._load still sends .mat to load_mat_mesh, the rest to trimesh",
      'self.base = load_mat_mesh(self.mesh_path)' in va
      and 'if suffix == ".mat":' in va
      and "trimesh.load(io.BytesIO(self.mesh_path.read_bytes())" in va)
check("vascular_api._load never hands trimesh a path (it follows names inside files)",
      "trimesh.load(str(self.mesh_path)" not in va)
never = [".zip", ".gz", ".bz2", ".tar.gz", ".gltf", ".glb", ".svg", ".dxf", ".xml",
         ".html", ".htm", ".exe", ".py", ".sip", ".dcm", ".mhd", ".raw", ".nii", ""]
bad = [e for e in never if e in S.MESH_EXTENSIONS + S.CENTERLINE_EXTENSIONS]
check("no archive, scene, script, imaging or empty extension on either list", not bad,
      f"found {bad}")
check("public lists drop the formats that name other files (.obj, .ply)",
      not set(S.PUBLIC_MESH_EXTENSIONS + S.PUBLIC_CENTERLINE_EXTENSIONS)
      & {".obj", ".ply"}
      and set(S.PUBLIC_MESH_EXTENSIONS) == {".mat", ".stl"})

# ── 2. hostile filenames cannot choose where anything is written ─────────────
print("\n2. hostile original filenames")
HOSTILE_ROOT = RUN / "store_hostile"
st = S.SubmissionStore(HOSTILE_ROOT, public=False)
LONG = "A" * 296 + ".mat"                                  # 300 characters
FAKE_ID = "20260101T000000Z-" + "0" * 32
surface_cases = [
    # (sent name, accepted?, recorded name if accepted)
    ("../../x.mat", True, "x.mat"),
    ("..\\..\\x.mat", True, "x.mat"),
    ("C:\\Windows\\x.mat", True, "x.mat"),
    ("/etc/x.mat", True, "x.mat"),
    ("a/b/c.stl", True, "c.stl"),
    ("x\x00.mat", True, "x?.mat"),
    ("evil.mat\x00.exe", False, None),
    (LONG, True, None),                                     # checked separately
    ("John_Smith_MRN123456.mat", True, "John_Smith_MRN123456.mat"),
    ("\\\\server\\share\\x.stl", True, "x.stl"),
    ("file:///etc/x.mat", True, "x.mat"),
    ("~/x.mat", True, "x.mat"),
    ("CON.stl", True, "CON.stl"),
    ("NUL.mat", True, "NUL.mat"),
    ("a\uff0fb.stl", True, "a\uff0fb.stl"),                 # fullwidth solidus
    ("\u202eltsx.mat", True, "\u202eltsx.mat"),             # right-to-left override
    (f"{FAKE_ID}.mat", True, f"{FAKE_ID}.mat"),             # a name posing as an ID
    ("X.STL", True, "X.STL"),
    (".mat", True, ".mat"),
    ("x.stl:stream", False, None),                          # NTFS alternate stream
    ("x.stl.", False, None),
    ("x.stl ", False, None),
    ("x.tar.gz", False, None),
    ("x.gltf", False, None),
    ("submission.json", False, None),
    ("", False, None),
    (None, False, None),
]
centerline_cases = [
    ("../../x.csv", True, "x.csv"),
    ("C:\\Windows\\x.csv", True, "x.csv"),
    ("a/b/c.csv", True, "c.csv"),
    ("x\x00.csv", True, "x?.csv"),
    ("C" * 296 + ".csv", True, None),
    ("John_Smith_MRN123456.csv", True, "John_Smith_MRN123456.csv"),
    ("submission.json", True, "submission.json"),           # must not clobber it
    ("surface.mat", False, None),
    ("case.mrk.json", True, "case.mrk.json"),
    ("x", False, None),
]

ids_made = []
for kind, cases in (("surface", surface_cases), ("centerline", centerline_cases)):
    for sent, accept, recorded in cases:
        surf = (sent, mesh_for(sent)) if kind == "surface" else ("s.stl", MESH)
        cl = (sent, CL) if kind == "centerline" else ("c.csv", CL)
        res, exc, ev = watched(st.create_submission, surf, cl)
        tag = f"{kind} {sent!r:.60}"
        if not accept:
            check(f"{tag}: refused, and nothing written",
                  isinstance(exc, S.UploadRejected) and not ev,
                  f"exc={exc!r} events={ev}")
            continue
        if exc is not None:
            check(f"{tag}: accepted", False, repr(exc))
            continue
        ids_made.append(res.id)
        folder = HOSTILE_ROOT / res.id
        ext_s = res.surface.extension
        ext_c = res.centerline.extension
        expect = {f"{res.id}/surface{ext_s}", f"{res.id}/centerline{ext_c}",
                  f"{res.id}/submission.json"}
        writes = {Path(p).resolve().relative_to(HOSTILE_ROOT.resolve()).as_posix()
                  for kind_, p in ev if kind_ == "write" and under(p, HOSTILE_ROOT)}
        outside = [e for e in ev if e[0] in ("write", "mkdir")
                   and not under(e[1], HOSTILE_ROOT)]
        mutating = [e for e in ev if e[0] not in ("write", "mkdir")]
        mkdirs = {Path(p).resolve() for k_, p in ev if k_ == "mkdir"}
        on_disk = {f"{res.id}/{p}" for p in tree(folder)}
        got_rec = (res.record["files"][kind]["original_name"])
        ok = (not outside and not mutating and writes == expect and on_disk == expect
              and mkdirs <= {HOSTILE_ROOT.resolve(), folder.resolve()}
              and (folder / f"surface{ext_s}").read_bytes() == surf[1]   # what was sent
              and (folder / f"centerline{ext_c}").read_bytes() == CL
              and (recorded is None or got_rec == recorded))
        check(f"{tag}: written only as {res.id[:8]}../{{surface,centerline,"
              f"submission.json}}, recorded as metadata", ok,
              f"outside={outside} mutating={mutating} writes={sorted(writes)} "
              f"disk={sorted(on_disk)} recorded={got_rec!r}")

# The 300-character names specifically.
recs = [json.loads((HOSTILE_ROOT / i / "submission.json").read_text("ascii"))
        for i in ids_made]
longs = [r["files"]["surface"]["original_name"] for r in recs
         if r["files"]["surface"]["original_name"].startswith("AAAA")]
check("300-char name: accepted, recorded at <= 255 chars with its .mat kept",
      len(longs) == 1 and len(longs[0]) <= 255 and longs[0].endswith("~.mat"),
      f"{[len(x) for x in longs]} {[x[-8:] for x in longs]}")
# A name posing as an ID never becomes the folder.
check("a name shaped like an ID does not become the folder",
      not (HOSTILE_ROOT / FAKE_ID).exists())
# Identifying fragments appear in no path anywhere under the run.
frag = ["John", "Smith", "MRN", "123456", "Windows", "etc", "server", "evil", "AAAA",
        "CCCC", "CON", "NUL", "x.mat", "x.csv"]
leaks = [p for p in (x.relative_to(RUN).as_posix() for x in RUN.rglob("*"))
         if any(f in p for f in frag)]
check("no identifying or hostile fragment in any path under the run", not leaks,
      f"{leaks[:5]}")
check("every submission.json is plain ASCII JSON",
      all((HOSTILE_ROOT / i / "submission.json").read_bytes().isascii()
          for i in ids_made))
check("sent names with a directory part are recorded as the bare name, flagged",
      all(r["files"]["surface"]["original_name_altered"] for r in recs
          if r["files"]["surface"]["original_name"] in ("x.mat", "c.stl")))

# ── 3. never overwrite: a forced ID collision fails and changes nothing ──────
print("\n3. exclusive create")
COLL = RUN / "store_collision"
st2 = S.SubmissionStore(COLL, public=False)
fixed = "20260924T000000Z-" + "ab" * 16
real_new_id = S._new_id
S._new_id = lambda ts: fixed
try:
    first, e1, _ = watched(st2.create_submission, ("a.stl", MESH), ("a.csv", CL))
    before = {p: (COLL / p).read_bytes() for p in tree(COLL)}
    mt = {p: os.stat(COLL / p).st_mtime_ns for p in before}
    second, e2, ev2 = watched(st2.create_submission,
                              ("b.stl", MESH + b"different"), ("b.csv", CL + b"9,9,9\n"))
finally:
    S._new_id = real_new_id
after = {p: (COLL / p).read_bytes() for p in tree(COLL)}
check("first submission stored", e1 is None and first.id == fixed, repr(e1))
check("second with the same ID raises StorageError",
      isinstance(e2, S.StorageError), repr(e2))
check("... and opened no file for writing",
      not [e for e in ev2 if e[0] == "write"], f"{ev2}")
check("... and the first is byte-identical, same mtimes",
      before == after and all(os.stat(COLL / p).st_mtime_ns == mt[p] for p in mt))
probe = COLL / fixed / "submission.json"
_, e3, _ = watched(S._write_new, probe, b"clobber")
check("_write_new onto an existing file raises FileExistsError, content kept",
      isinstance(e3, FileExistsError) and probe.read_bytes() == before[f"{fixed}/submission.json"],
      repr(e3))

# ── 4. cuts are kept with their submission, numbered, never replaced ─────────
print("\n4. cuts")
CUTS = RUN / "store_cuts"
st3 = S.SubmissionStore(CUTS, public=False, max_cuts=5)
sub = st3.create_submission(("k.stl", MESH), ("k.csv", CL))
side = {"mesh": "surface.stl", "openings": []}
c1, e, ev = watched(st3.add_cut, sub.id, b"STL-ONE", side)
c2, e_, ev_ = watched(st3.add_cut, sub.id, b"STL-TWO", json.dumps(side))
check("two cuts -> 001 and 002",
      e is None and e_ is None and (c1.cut_no, c2.cut_no) == (1, 2), f"{e!r} {e_!r}")
check("cut writes stay inside the submission folder, nothing mutating",
      all(under(p, CUTS / sub.id) for k_, p in ev + ev_ if k_ in ("write", "mkdir"))
      and not [x for x in ev + ev_ if x[0] not in ("write", "mkdir")])
# A folder already sitting at 003 (say, a crashed earlier write) is skipped, not reused.
squat = CUTS / sub.id / "cuts" / "003"
os.mkdir(squat)
with open(squat / "sentinel.txt", "xb") as f:
    f.write(b"leave me")
c4, e4, _ = watched(st3.add_cut, sub.id, b"STL-FOUR", side)
check("a taken number is skipped (003 -> 004), its contents untouched",
      e4 is None and c4.cut_no == 4 and (squat / "sentinel.txt").read_bytes() == b"leave me"
      and tree(squat) == ["sentinel.txt"], repr(e4))
rec = st3.load(sub.id)
check("load() lists complete cuts only, with sidecar paths relative to the folder",
      [c["cut_no"] for c in rec["cuts"]] == [1, 2, 4]
      and rec["cuts"][0]["sidecar"]["stored_as"] == "cuts/001/open_cut_planes.json"
      and rec["cuts"][0]["mesh"]["stored_as"] == "cuts/001/open.stl")
stl_p, side_p = st3.cut_paths(sub.id, 2)
check("cut_paths returns the stored bytes, sha256 matching the record",
      stl_p.read_bytes() == b"STL-TWO"
      and hashlib.sha256(side_p.read_bytes()).hexdigest()
      == rec["cuts"][1]["sidecar"]["sha256"])
_, e5, _ = watched(st3.add_cut, sub.id, b"STL", side)       # 005, the last allowed
_, e6, ev6 = watched(st3.add_cut, sub.id, b"STL", side)     # over max_cuts=5
check("past max_cuts: UploadRejected, no file written",
      e5 is None and isinstance(e6, S.UploadRejected)
      and not [x for x in ev6 if x[0] == "write"], f"{e5!r} {e6!r}")
_, e7, ev7 = watched(st3.add_cut, sub.id, b"STL", "{not json")
check("a sidecar that is not JSON is refused before anything is written",
      isinstance(e7, S.StorageError) and not ev7, f"{e7!r} {ev7}")
_, e8, _ = watched(st3.add_cut, "20260924T000000Z-" + "cd" * 16, b"STL", side)
check("a cut for an unknown submission: SubmissionNotFound",
      isinstance(e8, S.SubmissionNotFound), repr(e8))

# ── 4b. the face-count limit, enforced before anything is built ──────────────
print("\n4b. face-count limit")
FACES = RUN / "store_faces"
two = (b"solid t\n" + b"facet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\n"
       b"vertex 0 1 0\nendloop\nendfacet\n" * 2 + b"endsolid t\n")
over = S.SubmissionStore(FACES, public=False, max_faces=1)
_, e, ev = watched(over.create_submission, ("m.stl", two), ("c.csv", CL))
check("a 2-face mesh over a 1-face limit: refused, says both counts, nothing written",
      isinstance(e, S.UploadRejected) and "2 faces" in str(e) and "up to 1" in str(e)
      and "Nothing was stored" in str(e) and not ev, f"{e!r}")
print(f"         message: {e}")
at = S.SubmissionStore(FACES, public=False, max_faces=2)
r, e, _ = watched(at.create_submission, ("m.stl", two), ("c.csv", CL))
check("the same mesh at a 2-face limit: accepted", e is None, repr(e))
check("the default limit is a million, and VC_MAX_FACES overrides it",
      S.SubmissionStore(FACES, public=False).max_faces == 1_000_000
      and S.SubmissionStore.from_env(ROOT, public=False,
                                     environ={"VC_MAX_FACES": "250,000"}).max_faces
      == 250_000)
try:
    S.SubmissionStore.from_env(ROOT, public=False, environ={"VC_MAX_FACES": "lots"})
    check("a non-numeric VC_MAX_FACES stops the app at startup", False)
except S.StorageConfigError:
    check("a non-numeric VC_MAX_FACES stops the app at startup", True)

# ── 5. size limits refuse, never truncate ─────────────────────────────────────
print("\n5. size limits")
SIZE = RUN / "store_size"
st4 = S.SubmissionStore(SIZE, public=False, max_mesh_bytes=1000,
                        max_centerline_bytes=100)
_, e, ev = watched(st4.create_submission, ("m.stl", b"x" * 1001), ("c.csv", CL))
check("1001-byte mesh over a 1000-byte limit: refused, clear message, nothing written",
      isinstance(e, S.UploadRejected) and "limit" in str(e)
      and "Nothing was stored" in str(e) and not ev and not SIZE.exists(),
      f"{e!r} {ev}")
print(f"         message: {e}")
AT_LIMIT = MESH + b"\n" * (1000 - len(MESH))              # a valid STL, 1000 bytes
r, e, _ = watched(st4.create_submission, ("m.stl", AT_LIMIT), ("c.csv", CL))
check("exactly at the limit: accepted, all 1000 bytes kept",
      e is None and r.surface.path.read_bytes() == AT_LIMIT, repr(e))
stream = io.BytesIO(b"y" * 5000)
_, e, ev = watched(st4.create_submission, ("m.stl", stream), ("c.csv", CL))
check("a stream over the limit: refused, read stopped at limit+1, nothing written",
      isinstance(e, S.UploadRejected) and stream.tell() <= 1001 and not ev,
      f"{e!r} tell={stream.tell()}")
_, e, ev = watched(st4.create_submission, ("m.stl", b""), ("c.csv", CL))
check("an empty file: refused, nothing written",
      isinstance(e, S.UploadRejected) and "empty" in str(e) and not ev, repr(e))
_, e, ev = watched(st4.create_submission, ("m.stl", MESH), ("c.csv", b"1" * 101))
check("the centerline has its own limit", isinstance(e, S.UploadRejected)
      and "centerline" in str(e) and not ev, repr(e))

# ── 6. configuration from the environment ────────────────────────────────────
print("\n6. environment")
fake_repo = RUN / "fake_repo"
s_loc = S.SubmissionStore.from_env(fake_repo, public=False, environ={})
check("unset, local: root is <repo>/output/submissions, and nothing is created",
      s_loc.root == fake_repo / "output" / "submissions" and not fake_repo.exists())
check("defaults: 50 MB mesh, 10 MB centerline",
      s_loc.limits == {"surface": 50 * 2**20, "centerline": 10 * 2**20})
try:
    S.SubmissionStore.from_env(fake_repo, public=True, environ={})
    check("public with no VC_STORE_DIR is refused", False)
except S.StorageConfigError as ex:
    check("public with no VC_STORE_DIR is refused", "VC_STORE_DIR" in str(ex))
vol = RUN / "vol"
s_pub = S.SubmissionStore.from_env(fake_repo, public=True,
                                   environ={"VC_STORE_DIR": str(vol),
                                            "VC_MAX_MESH_MB": "0.5"})
check("VC_STORE_DIR and VC_MAX_MESH_MB are honoured",
      s_pub.root == vol.resolve() and s_pub.limits["surface"] == 2**19
      and s_pub.public)
for badv in ("abc", "0", "-1", "nan", "inf", "5000", "1e9"):
    try:
        S.SubmissionStore.from_env(fake_repo, public=False,
                                   environ={"VC_MAX_MESH_MB": badv})
        check(f"VC_MAX_MESH_MB={badv!r} refused at startup", False)
    except S.StorageConfigError:
        check(f"VC_MAX_MESH_MB={badv!r} refused at startup", True)
_, e, ev = watched(s_pub.create_submission, ("m.obj", MESH), ("c.csv", CL))
check("public store refuses an .obj surface, with the reason, writing nothing",
      isinstance(e, S.UploadRejected) and "refer to other files" in str(e) and not ev,
      repr(e))
_, e, ev = watched(s_pub.create_submission, ("m.stl", MESH), ("c.ply", CL))
check("public store refuses a .ply centerline", isinstance(e, S.UploadRejected)
      and not ev, repr(e))
check("accept attribute for the file input",
      s_pub.accept_attribute("surface") == ".mat,.stl")

# ── 7. the record ────────────────────────────────────────────────────────────
print("\n7. submission.json")
REC = RUN / "store_record"
st5 = S.SubmissionStore(REC, public=True)
t0 = datetime.now(timezone.utc)
sub5 = st5.create_submission(("case.stl", MESH), ("case.csv", CL), app_version="ui-17")
disk = json.loads((REC / sub5.id / "submission.json").read_text("ascii"))
check("keys are exactly the documented ones",
      set(disk) == {"schema", "id", "created_utc", "mode", "app_version", "notice",
                    "files"}, sorted(disk))
check("id matches the folder and the format", disk["id"] == sub5.id
      and S.ID_RE.fullmatch(sub5.id) and (REC / sub5.id).is_dir())
ts = datetime.strptime(disk["created_utc"], "%Y-%m-%dT%H:%M:%SZ").replace(
    tzinfo=timezone.utc)
check("created_utc is UTC and now", abs((ts - t0).total_seconds()) < 60,
      disk["created_utc"])
check("mode is public", disk["mode"] == "public")
f_ = disk["files"]
check("sizes and SHA-256 match the bytes",
      f_["surface"]["sha256"] == hashlib.sha256(MESH).hexdigest()
      and f_["centerline"]["sha256"] == hashlib.sha256(CL).hexdigest()
      and f_["surface"]["size_bytes"] == len(MESH)
      and f_["surface"]["stored_as"] == "surface.stl"
      and f_["surface"]["original_name"] == "case.stl")
check("notice fingerprint recorded, marked draft while placeholders remain",
      disk["notice"] == {"version": S.notice_version(), "draft": True})


def all_keys(o):
    if isinstance(o, dict):
        for k, v in o.items():
            yield k
            yield from all_keys(v)
    elif isinstance(o, list):
        for v in o:
            yield from all_keys(v)


who = {"ip", "ip_address", "remote_addr", "addr", "host", "hostname", "user",
       "user_id", "username", "user_agent", "email", "client", "session", "cookie"}
check("no IP, user, host or session field anywhere",
      not (set(k.lower() for k in all_keys(disk)) & who))
check("upload_paths gives what to hand CutSession",
      st5.upload_paths(sub5.id) == {"surface": REC / sub5.id / "surface.stl",
                                    "centerline": REC / sub5.id / "centerline.csv"})

# ── 8. IDs coming back from a client are validated before any path is built ──
print("\n8. IDs from outside")
good = sub5.id
arabic = good.translate(str.maketrans("0123456789", "\u0660\u0661\u0662\u0663\u0664"
                                                    "\u0665\u0666\u0667\u0668\u0669"))
bad_ids = ["../../x", "..\\x", "C:\\x", good + "/../x", good + "\\..", good.upper(),
           arabic, good + "\x00", " " + good, good[:-1], "", None, 5, [good]]
for b in bad_ids:
    for fn in (st5.load, st5.upload_paths, lambda i: st5.cut_paths(i, 1),
               lambda i: st5.add_cut(i, b"x", {})):
        _, e, ev = watched(fn, b)
        if not (isinstance(e, S.UploadRejected) and not ev):
            check(f"bad id {b!r:.40} refused by {getattr(fn, '__name__', 'fn')}", False,
                  repr(e))
            break
    else:
        check(f"bad id {b!r:.40} refused by load/upload_paths/cut_paths/add_cut", True)
for n in ("1/../..", True, 0, -1, 1000, 1.0):
    _, e, _ = watched(st5.cut_paths, good, n)
    check(f"cut number {n!r} refused", isinstance(e, S.UploadRejected), repr(e))

# ── 9. the retention notice ──────────────────────────────────────────────────
print("\n9. retention notice")
T = S.NOTICE_TEXT
low = T.lower()
check("says both files and the cut are kept",
      "both files" in low and "surface mesh and the centerline" in low
      and "every cut mesh you export" in low)
check("says what they are used for: the tool, and research (decided 25 Sep 2026)",
      "run this tool on them" in low and "research on aortic geometry" in low)
check("says how long: the project, capped at five years (decided 25 Sep 2026)",
      "until the end of this project" in low and "no longer than five years" in low)
check("says who can see them: the group and its collaborators (decided 25 Sep 2026)",
      "members of that group" in low and "collaborators" in low)
check("says who holds them (placeholder)",
      "[PI TO SUPPLY: institution or group responsible for the stored files]" in T)
check("says uploading is agreeing: a condition of use, not opt-in (decided 25 Sep 2026)",
      "by uploading you agree" in low and "do not upload" in low
      and "checkbox" not in low)
check("says how to ask for removal, at a real address (placeholder), quoting the ID",
      "removed" in low and "[PI TO SUPPLY: a monitored contact address]" in T
      and "submission id" in low)
check("tells people to keep names, MRNs and dates out of filenames",
      all(w in low for w in ("names", "medical record numbers", "dates of birth",
                             "scan dates", "file names")))
check("names the two places inside a Scan-IP centerline where the same turn up",
      "spline name column" in low and "first line" in low)
ph = S.notice_placeholders()
check("two placeholders left, the holder and the contact address",
      len(ph) == 2 and any("institution" in p for p in ph)
      and any("contact address" in p for p in ph), f"{ph}")
check("no apostrophe, quote, backtick or angle bracket (HOSTING.md template trap)",
      not any(c in T for c in "'\"`<>"))
check("fingerprint is stable and changes with the wording",
      S.notice_version() == S.notice_version(T) and S.notice_version(T + ".") != S.notice_version())

# ── 10. a real case, end to end ──────────────────────────────────────────────
if DO_CUT:
    print("\n10. KY160 stored, loaded from the store, cut, and the cut kept")
    from vascular_api import CutSession                      # noqa: E402
    mat_src = ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat"
    csv_src = ROOT / "KY160_1_CentrelinePoints.csv"
    REAL = RUN / "store_real"
    st6 = S.SubmissionStore(REAL, public=True)
    sub6 = st6.create_submission((mat_src.name, mat_src.read_bytes()),
                                 (csv_src.name, csv_src.read_bytes()))
    p = st6.upload_paths(sub6.id)
    t = time.perf_counter()
    s = CutSession(str(p["centerline"]), str(p["surface"]))
    dt = time.perf_counter() - t
    check("CutSession loads and cuts from the stored copies", bool(s.ok),
          str(s.problem))
    print(f"         {len(s.openings())} openings, {len(s.mesh.faces)} faces, "
          f"load+cut {dt:.1f} s")
    stl = s.mesh.export(file_type="stl")
    sc = s.sidecar()
    cut, e, ev = watched(st6.add_cut, sub6.id, stl, sc)
    check("the cut and its sidecar are kept with the upload",
          e is None and cut.stl_path.parent.parent.parent == REAL / sub6.id
          and cut.stl_path.read_bytes() == stl, repr(e))
    kept = json.loads(cut.sidecar_path.read_text("utf-8"))
    check("sidecar names the stored files, not the uploader's",
          kept["mesh"] == "surface.mat" and kept["centerline"] == "centerline.csv",
          f"{kept['mesh']} {kept['centerline']}")
    check("the original stem appears nowhere in the stored sidecar",
          b"KY160" not in cut.sidecar_path.read_bytes())
    # Against what the app's own export writes today, into a NEW folder.
    ref_dir = RUN / "export_reference"
    assert not ref_dir.exists()
    ref_stl, ref_side = s.export(ref_dir, "wf_ref")
    check("stored STL is byte-identical to CutSession.export's",
          ref_stl.read_bytes() == stl)
    check("stored sidecar parses identical to CutSession.export's",
          json.loads(ref_side.read_text()) == kept)
    rb, kb = ref_side.read_bytes(), cut.sidecar_path.read_bytes()
    crlf_r, crlf_k = rb.count(b"\r\n"), kb.count(b"\r\n")
    print(f"         sidecar bytes: export {len(rb)} ({crlf_r} CRLF), "
          f"stored {len(kb)} ({crlf_k} CRLF); equal once CRLF->LF: "
          f"{rb.replace(bytes([13, 10]), bytes([10])) == kb}")
    rec6 = st6.load(sub6.id)
    check("record: original names as metadata, one cut, sidecar path",
          rec6["files"]["surface"]["original_name"] == "KY160_1.mat"
          and rec6["files"]["centerline"]["original_name"] == "KY160_1_CentrelinePoints.csv"
          and [c["sidecar"]["stored_as"] for c in rec6["cuts"]]
          == ["cuts/001/open_cut_planes.json"])
else:
    print("\n10. skipped (--no-cut)")

# ── the whole run: only the expected shapes of path, and nothing deleted ─────
print("\nwhole run")
shape = re.compile(
    r"store_[a-z]+/[0-9]{8}T[0-9]{6}Z-[0-9a-f]{32}/("
    r"(surface|centerline)\.[a-z0-9]+|submission\.json|"
    r"cuts/[0-9]{3}/(open\.stl|open_cut_planes\.json|cut\.json|sentinel\.txt))"
    r"|export_reference/wf_ref(\.stl|_cut_planes\.json)")
odd = [p for p in tree(RUN) if not shape.fullmatch(p)]
check("every file under the run has a generated-ID path of a known shape", not odd,
      f"{odd[:6]}")

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): {[f[0] for f in fails]}")
raise SystemExit(1 if fails else 0)
