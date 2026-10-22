# Where uploads, and the cuts made from them, are kept.
#
# tools/HOSTING.md job 7: uploads are RETAINED -- the surface and the centerline -- and
# every cut exported from them is kept beside them with its sidecar JSON. This module
# is only the storage. It knows nothing about trame, VTK or the geometry, so it can be
# tested on its own and so a change to the app cannot quietly change what is kept.
#
# The layout, one folder per submission, nothing in it ever rewritten:
#
#   <root>/<id>/surface.<ext>              the mesh exactly as uploaded
#   <root>/<id>/centerline.<ext>           the centerline exactly as uploaded
#   <root>/<id>/submission.json            written LAST: its presence means complete
#   <root>/<id>/cuts/001/open.stl          a cut, as the app would have exported it
#   <root>/<id>/cuts/001/open_cut_planes.json   its sidecar
#   <root>/<id>/cuts/001/cut.json          written LAST for that cut
#
# Four rules, and why:
#
# * The ID is generated here and is the only thing from which a path is built. The
#   uploader's filename is recorded in submission.json and NEVER becomes part of a
#   path: filenames are where patient names and MRNs turn up, and a name is also the
#   classic way to walk out of a directory. The stored names are fixed ("surface",
#   "centerline"), which also keeps an identifying stem out of the sidecar -- the
#   sidecar records mesh_path.name, and branch_tree names an unlabelled curve after
#   the file's stem.
# * The extension is kept, lower-cased, from an allow-list, because the loaders
#   dispatch on it (vascular_api._load: .mat -> load_mat_mesh, anything else ->
#   trimesh.load; branch_tree.load_spline_csv: polyline suffixes -> pyvista, .json ->
#   JSON, anything else -> a table). A file whose extension is not on the list is
#   refused, not guessed at.
# * Nothing is ever overwritten, moved or deleted. Folders are made with os.mkdir and
#   files opened O_CREAT|O_EXCL, so a collision FAILS instead of clobbering. There is
#   no delete function: removal on request is an operator's job, by hand, against the
#   ID the uploader was shown.
# * Everything is checked in memory before the first byte touches the disk, so a
#   refused upload leaves nothing behind -- this module cannot tidy up after itself,
#   by design.

import hashlib
import json
import os
import re
import secrets
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from meshcount import MeshCountError, count_faces

SCHEMA_SUBMISSION = "vc-submission/1"
SCHEMA_CUT = "vc-cut/1"

# ── configuration, in the VC_* style of VC_PUBLIC / VC_CUT_PORT ──────────────
ENV_STORE_DIR = "VC_STORE_DIR"
ENV_MAX_MESH_MB = "VC_MAX_MESH_MB"
ENV_MAX_CENTERLINE_MB = "VC_MAX_CENTERLINE_MB"
ENV_MAX_FACES = "VC_MAX_FACES"

# Real meshes run 10-20 MB (the samples are 1.4-1.7 MB); centerlines under 1 MB.
# MB here means 2**20 bytes. A byte cap is NOT a memory cap: a compressed .mat
# expands roughly a thousandfold on load, so job 4's face-count ceiling is still
# needed on top of this.
DEFAULT_MAX_MESH_MB = 50.0
DEFAULT_MAX_CENTERLINE_MB = 10.0

# THE FACE-COUNT LIMIT, and why it is a million. Set from two measurements, not a guess.
# What people have: this lab's thoracic cohort runs 14k-74k faces, but the Vascular Model
# Repository -- the main public collection -- has aortic surfaces with median 91k, 95th
# percentile 267k, 99th 529k and maximum 666k faces (131 aortic models, headers read
# 24 Sep 2026). A limit near this lab's own data would have turned away 43% of those; a
# million turns away none of them, and 2 of all 254 of the repository's models.
# What it costs, measured 25 Sep 2026 on the development workstation: 589k-637k faces
# load in 18-30 s and take 31-48 s for the first full cut; nothing has been timed at
# the limit itself (tools/HOSTING.md job 4). Generous by design; the host has to be
# sized for it.
DEFAULT_MAX_FACES = 1_000_000
MAX_CONFIGURABLE_FACES = 20_000_000
MAX_CONFIGURABLE_MB = 2048.0

# Each press of Export on the same upload is a new cut folder. Bounded so one visitor
# cannot fill the volume from one submission; a fresh upload starts a fresh count.
DEFAULT_MAX_CUTS = 50

# ── what may be stored ───────────────────────────────────────────────────────
# Surface: .mat goes to load_mat_mesh; the rest are what the desktop's Open dialog
# offers (branch_gui.py "Meshes (*.mat *.stl *.ply *.obj)"). trimesh.load itself
# takes ~50 formats, among them .zip/.tar.gz archives, .gltf (which references
# other files) and .svg/.dxf (which are not surfaces): none of those belong here.
MESH_EXTENSIONS = (".mat", ".stl", ".ply", ".obj")

# Centerline: exactly the suffixes branch_tree.load_spline_csv names --
# TABLE_SUFFIXES (less the empty one), POLYDATA_SUFFIXES, and .json (which covers
# Slicer's .mrk.json, whose last suffix is .json). A file with NO extension is read as
# a table by the loader, but refused here: that is guessing.
CENTERLINE_EXTENSIONS = (".csv", ".tsv", ".txt", ".dat", ".asc", ".pts", ".xyz",
                         ".json", ".vtp", ".vtk", ".vtu", ".obj", ".ply")

# Formats that can NAME ANOTHER FILE, which the loader then opens: an OBJ's `mtllib`
# line and a PLY's `comment TextureFile` line. Measured with an audit hook
# (scratch_wf/storage-retention/probe_external_refs.py): trimesh.load(path) on an
# uploaded OBJ opened ../../outside/not_yours.mtl and an absolute path, and a PLY
# opened ../../outside/not_yours.png. On a public host that is a file-read primitive,
# and `mtllib /dev/zero` would be a memory exhaustion. Loading from a stream
# (trimesh.load(io.BytesIO(data), file_type=...)) opened nothing outside, so these can
# come back to the public lists once vascular_api._load does that. The centerline
# route (pyvista -> vtkOBJReader/vtkPLYReader) opens files in C++, which the probe
# cannot see, so it is excluded on the same grounds until someone checks it.
NAMES_OTHER_FILES = (".obj", ".ply")
PUBLIC_MESH_EXTENSIONS = tuple(e for e in MESH_EXTENSIONS
                               if e not in NAMES_OTHER_FILES)
PUBLIC_CENTERLINE_EXTENSIONS = tuple(e for e in CENTERLINE_EXTENSIONS
                                     if e not in NAMES_OTHER_FILES)

KINDS = {
    # kind: (stored stem, label for messages)
    "surface": ("surface", "surface mesh"),
    "centerline": ("centerline", "centerline"),
}

CUT_STEM = "open"                       # open.stl + open_cut_planes.json, the same
SIDECAR_SUFFIX = "_cut_planes.json"     # <name>.stl / <name>_cut_planes.json pairing
                                        # CutSession.export and the viewers use

MAX_ORIGINAL_NAME = 255

# The ID: UTC time for a sortable listing (the upload time, not anything about the
# patient), then 128 random bits so it cannot be guessed -- it is what a download
# link or a removal request will quote. [0-9], not \d: \d matches any Unicode digit.
ID_RE = re.compile(r"[0-9]{8}T[0-9]{6}Z-[0-9a-f]{32}")


# ── errors ───────────────────────────────────────────────────────────────────
class UploadRejected(ValueError):
    """The upload or request is not acceptable. The message is safe to show the
    person who sent it: it names no server path."""


class SubmissionNotFound(LookupError):
    """A well-formed ID with no complete submission behind it."""


class StorageError(RuntimeError):
    """A server-side failure. Log the message; show the user something generic."""


class StorageConfigError(StorageError):
    """The environment asks for something this store will not do."""


# ── the retention notice (HOSTING.md job 7) ──────────────────────────────────
# Goes ABOVE the upload control, in full, not behind a link. Every [PI TO SUPPLY: ...]
# is a decision that is not engineering's to make. Nothing here may be filled in with
# a guess: an invented period or address would be read as policy.
#
# No apostrophes, quotes or angle brackets anywhere in it, on purpose: HOSTING.md
# records an apostrophe near a Vue template taking down a whole panel. Render each
# paragraph as plain text in its own element, composed in Python, not as HTML.
NOTICE_HEADING = "Your files are kept. Please read this before uploading."

NOTICE_PARAGRAPHS = (
    # 1. that both files are kept, with the cut produced from them
    "What is kept: both files you upload, the surface mesh and the centerline, "
    "exactly as you send them, and every cut mesh you export from them, together "
    "with the record of the cut settings that produced it. The name of each file is "
    "recorded as you send it.",

    # 2. what they will be used for -- decided 25 Sep 2026: the tool, and research.
    # Research use of uploads is what makes institutional sign-off a prerequisite
    # (HOSTING.md), not a formality.
    "What they are used for: to run this tool on them for you, to improve the tool, "
    "and in research on aortic geometry by the group that holds them.",

    # 3. how long they are held, and who can access them -- period and access
    # decided 25 Sep 2026. Nothing deletes on a schedule yet: the five-year cap is
    # kept by hand until a clean-up job is approved and built.
    "How long they are kept: until the end of this project, and no longer than five "
    "years from upload. "
    "They are held by [PI TO SUPPLY: institution or group responsible for the "
    "stored files]. "
    "Who can see them: members of that group, and collaborators working with them "
    "on this tool.",

    # 4. how to ask for removal, and a real address
    "To have your files removed, write to [PI TO SUPPLY: a monitored contact "
    "address] and quote the submission ID this page shows you after uploading. "
    "Keep that ID: it is how your files are found.",

    # 5. no names, MRNs or dates in filenames -- and the two places inside a
    # centerline export where the same things turn up
    "Do not put names, medical record numbers, dates of birth, scan dates or any "
    "other identifying detail in your file names. Rename the files before you "
    "upload them. Check inside the centerline as well: a Scan-IP export names each "
    "curve in its Spline name column, and its first line can record the folder the "
    "project was saved in and the date it was exported. Upload only a surface and a "
    "centerline, never scans, DICOM files or reports.",

    # the basis for keeping them -- decided 25 Sep 2026: a condition of use, so
    # there is no opt-in tick box and no path that uploads without keeping
    "By uploading you agree to your files being kept and used as described above. "
    "If you do not agree, do not upload.",
)

NOTICE_TEXT = NOTICE_HEADING + "\n\n" + "\n\n".join(NOTICE_PARAGRAPHS)
PLACEHOLDER_RE = re.compile(r"\[PI TO SUPPLY: [^\]]+\]")


def notice_version(text=NOTICE_TEXT):
    """A short fingerprint of the notice wording, recorded with each submission so
    it is known afterwards which wording the uploader was shown."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def notice_placeholders(text=NOTICE_TEXT):
    """Every [PI TO SUPPLY: ...] still in the text. Non-empty means draft."""
    return PLACEHOLDER_RE.findall(text)


# ── small helpers ────────────────────────────────────────────────────────────
def _utc_now():
    return datetime.now(timezone.utc)


def _iso(ts):
    return ts.isoformat(timespec="seconds").replace("+00:00", "Z")


def _faces_from_env(env):
    raw = (env.get(ENV_MAX_FACES) or "").strip().replace("_", "").replace(",", "")
    if not raw:
        return DEFAULT_MAX_FACES
    try:
        n = int(raw)
    except ValueError:
        raise StorageConfigError(f"{ENV_MAX_FACES}={raw!r} is not a whole number")
    if not 1 <= n <= MAX_CONFIGURABLE_FACES:
        raise StorageConfigError(
            f"{ENV_MAX_FACES} must be between 1 and {MAX_CONFIGURABLE_FACES:,}")
    return n


def _new_id(ts):
    return f"{ts:%Y%m%dT%H%M%SZ}-{secrets.token_hex(16)}"


def _mb(n_bytes):
    return n_bytes / (1024 * 1024)


def _size(n_bytes):
    """A size a person can read: MB (2**20), KB or bytes."""
    if n_bytes >= 1024 * 1024:
        return f"{_mb(n_bytes):.1f} MB"
    if n_bytes >= 1024:
        return f"{n_bytes / 1024:.1f} KB"
    return f"{n_bytes} bytes"


def _clean_original_name(name):
    """(basename as recorded, changed?, extension) -- METADATA ONLY, never a path.

    Only the last component is kept: a browser sends a bare name, and anything more
    (C:\\Users\\<someone>\\Desktop\\...) says who the uploader is, which is not wanted.
    Control characters, the NUL byte among them, become '?' so the JSON reads cleanly
    in any tool, and the result is capped at MAX_ORIGINAL_NAME characters.
    """
    raw = "" if name is None else str(name)
    base = re.split(r"[\\/]", raw)[-1]
    clean = "".join("?" if (ord(c) < 32 or ord(c) == 127) else c for c in base)
    ext = _extension(clean)                 # from the WHOLE name, before any cut
    if len(clean) > MAX_ORIGINAL_NAME:
        # Shorten before the extension, which decides the format and should stay
        # visible; '~' marks the cut.
        keep = ext if len(ext) <= 16 else ""
        clean = clean[:MAX_ORIGINAL_NAME - 1 - len(keep)] + "~" + keep
    return clean, clean != raw, ext


def _extension(clean_name):
    """The last '.suffix', lower-cased, or '' if there is none. ASCII only."""
    i = clean_name.rfind(".")
    if i < 0 or i == len(clean_name) - 1:
        return ""
    ext = clean_name[i:].lower()
    return ext if ext.isascii() else ""


def _as_bytes(data, limit, label):
    """The upload as bytes, refusing anything over `limit` without reading past it."""
    if isinstance(data, (bytes, bytearray, memoryview)):
        n = memoryview(data).nbytes
        if n > limit:
            raise UploadRejected(
                f"The {label} file is {_size(n)}, over the {_size(limit)} limit for "
                f"a {label} on this site. Nothing was stored.")
        return bytes(data)
    if hasattr(data, "read"):
        chunks, total = [], 0
        while True:
            chunk = data.read(min(1 << 20, limit + 1 - total))
            if not chunk:
                break
            chunks.append(bytes(chunk))
            total += len(chunk)
            if total > limit:
                raise UploadRejected(
                    f"The {label} file is over the {_size(limit)} limit for a "
                    f"{label} on this site. Nothing was stored.")
        return b"".join(chunks)
    raise TypeError(f"expected bytes or a binary file object, got {type(data).__name__}")


def _write_new(path, data):
    """Create `path` and write `data`; fail if anything is already there."""
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    fd = os.open(path, flags, 0o640)
    try:
        view = memoryview(data)
        while view:
            n = os.write(fd, view)
            view = view[n:]
        os.fsync(fd)
    finally:
        os.close(fd)


def _json_bytes(obj):
    return (json.dumps(obj, indent=2, ensure_ascii=True) + "\n").encode("ascii")


@dataclass
class CheckedUpload:
    name: str               # the recorded original name -- metadata, never a path
    name_altered: bool      # directory part dropped, control characters or length
    extension: str
    data: bytes = field(repr=False)
    faces: int | None = None    # surfaces only, counted before loading


@dataclass
class StoredFile:
    kind: str
    path: Path
    original_name: str
    extension: str
    size_bytes: int
    sha256: str


@dataclass
class Submission:
    id: str
    folder: Path
    surface: StoredFile
    centerline: StoredFile
    record: dict = field(repr=False)


@dataclass
class CutRecord:
    submission_id: str
    cut_no: int
    stl_path: Path
    sidecar_path: Path
    record: dict = field(repr=False)


# ── the store ────────────────────────────────────────────────────────────────
class SubmissionStore:
    """Uploads and cuts under generated IDs. Adds only; never rewrites or removes."""

    def __init__(self, root, *, public,
                 max_mesh_bytes=int(DEFAULT_MAX_MESH_MB * 1024 * 1024),
                 max_centerline_bytes=int(DEFAULT_MAX_CENTERLINE_MB * 1024 * 1024),
                 max_cuts=DEFAULT_MAX_CUTS, max_faces=DEFAULT_MAX_FACES):
        self.root = Path(root)
        self.public = bool(public)
        self.limits = {"surface": int(max_mesh_bytes),
                       "centerline": int(max_centerline_bytes)}
        if min(self.limits.values()) <= 0:
            raise StorageConfigError("upload size limits must be positive")
        if not 1 <= int(max_cuts) <= 999:
            raise StorageConfigError("max_cuts must be between 1 and 999")
        self.max_cuts = int(max_cuts)
        if not 1 <= int(max_faces) <= MAX_CONFIGURABLE_FACES:
            raise StorageConfigError(
                f"max_faces must be between 1 and {MAX_CONFIGURABLE_FACES:,}")
        self.max_faces = int(max_faces)
        self._allowed = {
            "surface": PUBLIC_MESH_EXTENSIONS if self.public else MESH_EXTENSIONS,
            "centerline": (PUBLIC_CENTERLINE_EXTENSIONS if self.public
                           else CENTERLINE_EXTENSIONS),
        }

    @classmethod
    def from_env(cls, repo_root, *, public, environ=None, **kw):
        """The store the app should use. Nothing is created until the first upload.

        VC_STORE_DIR   where submissions go. Unset: <repo>/output/submissions, which
                       is gitignored -- but only when NOT public. A container's own
                       filesystem is thrown away with the container (HOSTING.md job
                       7), so the hosted configuration must name a persistent volume
                       rather than fall back to one silently.
        VC_MAX_MESH_MB, VC_MAX_CENTERLINE_MB   per-file limits, in MB (2**20 bytes).
        """
        env = os.environ if environ is None else environ
        raw = (env.get(ENV_STORE_DIR) or "").strip()
        if raw:
            root = Path(raw).expanduser().resolve()
        elif public:
            raise StorageConfigError(
                f"VC_PUBLIC is set but {ENV_STORE_DIR} is not. A hosted instance must "
                f"keep uploads on persistent storage: set {ENV_STORE_DIR} to a mounted "
                f"volume. The default, output/submissions inside the app folder, is "
                f"lost when the container is replaced.")
        else:
            root = Path(repo_root) / "output" / "submissions"
        return cls(root, public=public,
                   max_mesh_bytes=_limit_from_env(env, ENV_MAX_MESH_MB,
                                                  DEFAULT_MAX_MESH_MB),
                   max_centerline_bytes=_limit_from_env(env, ENV_MAX_CENTERLINE_MB,
                                                        DEFAULT_MAX_CENTERLINE_MB),
                   max_faces=_faces_from_env(env),
                   **kw)

    # ── what the page should say ─────────────────────────────────────────────
    def accepted_extensions(self, kind):
        return tuple(self._allowed[kind])

    def accept_attribute(self, kind):
        """For the file input's `accept=`. A hint to the browser, not a control."""
        return ",".join(self._allowed[kind])

    def limit_mb(self, kind):
        return _mb(self.limits[kind])

    # ── checking an upload, in memory ────────────────────────────────────────
    def check_upload(self, kind, original_name, data):
        """CheckedUpload(name, name_altered, extension, data), or UploadRejected.

        Writes nothing. `name` is the recorded original name (metadata only).
        """
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {sorted(KINDS)}")
        _stem, label = KINDS[kind]
        clean, altered, ext = _clean_original_name(original_name)
        allowed = self._allowed[kind]
        if not ext:
            raise UploadRejected(
                f"The {label} file has no extension, so its format cannot be told. "
                f"Rename it to end in one of: {', '.join(allowed)}.")
        if ext not in allowed:
            extra = ""
            if self.public and ext in NAMES_OTHER_FILES:
                extra = (f" {ext.upper()[1:]} files can refer to other files on the "
                         f"computer that reads them, so this site does not take them.")
            shown = (f"ending in {ext}" if ext[1:].isalnum() and len(ext) <= 12
                     else "with that ending")
            raise UploadRejected(
                f"A {label} file {shown} is not accepted here.{extra} "
                f"Accepted: {', '.join(allowed)}.")
        blob = _as_bytes(data, self.limits[kind], label)
        if not blob:
            raise UploadRejected(f"The {label} file is empty. Nothing was stored.")
        faces = None
        if kind == "surface":
            # Counted from the file's own structure BEFORE anything builds the mesh. A
            # byte cap does not bound the load: a compressed .mat expands about
            # 1,600-fold, and a file can simply lie about its size in its header.
            try:
                faces = count_faces(blob, ext).faces
            except MeshCountError as e:
                raise UploadRejected(str(e)) from None
            if faces > self.max_faces:
                raise UploadRejected(
                    f"This mesh has {faces:,} faces, and this site takes up to "
                    f"{self.max_faces:,}. Reduce it first, for example with Decimate in "
                    f"3D Slicer or MeshLab. A cut does not need more, and the curvature "
                    f"measurement resamples to 100,000 faces anyway. Nothing was stored.")
        return CheckedUpload(clean, altered, ext, blob, faces)

    # ── storing ──────────────────────────────────────────────────────────────
    def create_submission(self, surface, centerline, *, notice_text=NOTICE_TEXT,
                          app_version=None):
        """Store one surface + centerline pair as a new submission.

        surface, centerline: (original filename, bytes or binary file object).
        notice_text: the notice as the page displayed it; its fingerprint, and
            whether it still held placeholders, go in the record.

        Both files are checked before anything is written. Returns a Submission
        whose .surface.path and .centerline.path are what to hand to CutSession.
        """
        checked = {}
        for kind, item in (("surface", surface), ("centerline", centerline)):
            if not (isinstance(item, (tuple, list)) and len(item) == 2):
                raise TypeError(f"{kind} must be (original_name, data)")
            checked[kind] = self.check_upload(kind, item[0], item[1])

        ts = _utc_now()
        sub_id = _new_id(ts)
        if not ID_RE.fullmatch(sub_id):
            raise StorageError("generated ID failed its own format check")
        try:
            self.root.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise StorageError(f"cannot create the store root: {e}") from e
        folder = self._folder(sub_id)
        try:
            os.mkdir(folder, 0o750)
        except FileExistsError:
            raise StorageError(f"submission {sub_id} already exists; nothing was "
                               f"written") from None

        files, stored = {}, {}
        try:
            for kind, c in checked.items():
                stem, _label = KINDS[kind]
                target = self._inside(folder / f"{stem}{c.extension}")
                _write_new(target, c.data)
                digest = hashlib.sha256(c.data).hexdigest()
                files[kind] = {
                    "stored_as": target.name,
                    "original_name": c.name,
                    "original_name_altered": c.name_altered,
                    "extension": c.extension,
                    "size_bytes": len(c.data),
                    "sha256": digest,
                }
                stored[kind] = StoredFile(kind, target, c.name, c.extension,
                                          len(c.data), digest)
        except OSError as e:
            raise StorageError(
                f"writing submission {sub_id} failed ({e}); it has no "
                f"submission.json, so it reads as incomplete, and it was left as "
                f"it is") from e

        record = {
            "schema": SCHEMA_SUBMISSION,
            "id": sub_id,
            "created_utc": _iso(ts),
            "mode": "public" if self.public else "local",
            "app_version": app_version,
            "notice": {"version": notice_version(notice_text),
                       "draft": bool(notice_placeholders(notice_text))},
            "files": files,
        }
        # Last, so a submission.json means every file above is whole.
        _write_new(self._inside(folder / "submission.json"), _json_bytes(record))
        return Submission(sub_id, folder, stored["surface"], stored["centerline"],
                          record)

    def add_cut(self, submission_id, stl_bytes, sidecar, *, app_version=None):
        """Keep one exported cut with its submission, as the next numbered cut.

        stl_bytes: the cut mesh as STL, e.g. session.mesh.export(file_type="stl").
        sidecar:   session.sidecar() (a dict), or its JSON already serialised.
        """
        folder = self._complete_folder(submission_id)
        if not isinstance(stl_bytes, (bytes, bytearray, memoryview)) \
                or not memoryview(stl_bytes).nbytes:
            raise StorageError("add_cut needs the cut mesh as non-empty bytes")
        if isinstance(sidecar, dict):
            # The same serialisation CutSession.export uses (json.dumps, indent=2,
            # ASCII). It writes through write_text, so on Windows its file has CRLF
            # line ends; this is LF. Parsed, the two are identical.
            side_blob = json.dumps(sidecar, indent=2).encode("utf-8")
        elif isinstance(sidecar, (bytes, bytearray)):
            side_blob = bytes(sidecar)
        elif isinstance(sidecar, str):
            side_blob = sidecar.encode("utf-8")
        else:
            raise TypeError("sidecar must be a dict, bytes or str")
        try:
            json.loads(side_blob)
        except ValueError as e:
            raise StorageError(f"the sidecar is not valid JSON: {e}") from e

        cuts = folder / "cuts"
        try:
            os.mkdir(cuts, 0o750)
        except FileExistsError:
            pass                                    # the folder of earlier cuts
        cut_no, cut_dir = None, None
        for n in range(1, self.max_cuts + 1):       # claim the next free number;
            d = cuts / f"{n:03d}"                   # mkdir is the lock, so two
            try:                                    # writers cannot share one
                os.mkdir(d, 0o750)
            except FileExistsError:
                continue
            cut_no, cut_dir = n, d
            break
        if cut_no is None:
            raise UploadRejected(
                f"This upload already has {self.max_cuts} exported cuts, the most one "
                f"upload may keep. Upload the files again to start a new one.")

        stl_path = self._inside(cut_dir / f"{CUT_STEM}.stl")
        side_path = self._inside(cut_dir / f"{CUT_STEM}{SIDECAR_SUFFIX}")
        stl_blob = bytes(stl_bytes)
        try:
            _write_new(stl_path, stl_blob)
            _write_new(side_path, side_blob)
        except OSError as e:
            raise StorageError(
                f"writing cut {cut_no:03d} of {submission_id} failed ({e}); it has "
                f"no cut.json, so it reads as incomplete, and it was left as it is"
            ) from e
        rel = lambda p: p.relative_to(folder).as_posix()       # noqa: E731
        record = {
            "schema": SCHEMA_CUT,
            "submission_id": submission_id,
            "cut_no": cut_no,
            "created_utc": _iso(_utc_now()),
            "app_version": app_version,
            "mesh": {"stored_as": rel(stl_path), "size_bytes": len(stl_blob),
                     "sha256": hashlib.sha256(stl_blob).hexdigest()},
            "sidecar": {"stored_as": rel(side_path), "size_bytes": len(side_blob),
                        "sha256": hashlib.sha256(side_blob).hexdigest()},
        }
        _write_new(self._inside(cut_dir / "cut.json"), _json_bytes(record))
        return CutRecord(submission_id, cut_no, stl_path, side_path, record)

    # ── reading back (downloads, operators) ──────────────────────────────────
    def load(self, submission_id):
        """submission.json plus a "cuts" list of every complete cut, in order."""
        folder = self._complete_folder(submission_id)
        rec = json.loads((folder / "submission.json").read_text(encoding="ascii"))
        cuts = []
        cdir = folder / "cuts"
        if cdir.is_dir():
            for d in sorted(cdir.iterdir()):
                f = d / "cut.json"
                if re.fullmatch(r"[0-9]{3}", d.name) and f.is_file():
                    cuts.append(json.loads(f.read_text(encoding="ascii")))
        rec["cuts"] = cuts
        return rec

    def upload_paths(self, submission_id):
        """{"surface": Path, "centerline": Path} for a complete submission."""
        folder = self._complete_folder(submission_id)
        rec = json.loads((folder / "submission.json").read_text(encoding="ascii"))
        return {k: self._inside(folder / v["stored_as"])
                for k, v in rec["files"].items()}

    def cut_paths(self, submission_id, cut_no):
        """(stl, sidecar) of one complete cut -- what a download serves."""
        if isinstance(cut_no, bool) or not isinstance(cut_no, int) \
                or not 1 <= cut_no <= 999:
            raise UploadRejected("That is not a cut number.")
        folder = self._complete_folder(submission_id)
        d = folder / "cuts" / f"{cut_no:03d}"
        if not (d / "cut.json").is_file():
            raise SubmissionNotFound("No such cut for this submission.")
        return (self._inside(d / f"{CUT_STEM}.stl"),
                self._inside(d / f"{CUT_STEM}{SIDECAR_SUFFIX}"))

    # ── paths: only ever from a validated ID and fixed names ─────────────────
    def _folder(self, submission_id):
        if not isinstance(submission_id, str) or not ID_RE.fullmatch(submission_id):
            raise UploadRejected("That is not a submission ID.")
        return self._inside(self.root / submission_id)

    def _complete_folder(self, submission_id):
        folder = self._folder(submission_id)
        if not (folder / "submission.json").is_file():
            raise SubmissionNotFound("No submission with that ID.")
        return folder

    def _inside(self, path):
        """Belt and braces: every path built here must resolve inside the root."""
        root = self.root.resolve()
        p = Path(path)
        try:
            p.resolve().relative_to(root)
        except ValueError:
            raise StorageError(f"refusing a path outside the store: {p}") from None
        return p


def _limit_from_env(env, name, default_mb):
    raw = (env.get(name) or "").strip()
    if not raw:
        mb = default_mb
    else:
        try:
            mb = float(raw)
        except ValueError:
            raise StorageConfigError(f"{name}={raw!r} is not a number of MB") from None
        if not (0 < mb <= MAX_CONFIGURABLE_MB):
            raise StorageConfigError(
                f"{name}={raw!r} must be above 0 and at most {MAX_CONFIGURABLE_MB:g}")
    return int(mb * 1024 * 1024)
