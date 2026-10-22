"""A mesh's face count, read without loading it (webapp/meshcount.py).

A plain script like the rest of tests/: prints each check, then ALL PASS or FAILURES,
and exits nonzero on failure.

What it proves, by running rather than by reading:
  * the count equals what the real loader builds, for every accepted format -- binary
    and ASCII STL, ASCII and binary PLY, OBJ, and MATLAB .mat from the repo and the
    thoracic cohort when it is mounted;
  * a file that lies about itself is refused: an STL header claiming four billion
    triangles, a PLY claiming a trillion faces, a truncated file, a MATLAB v7.3 file;
  * a compressed .mat that expands enormously is counted in FLAT memory -- the thing
    scipy.io.whosmat does not do -- and one past the decompression ceiling is refused;
  * counting a large file is cheap.

Writes nothing to disk: every test mesh is built in memory.
"""
import io
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import trimesh

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for p in (ROOT, ROOT / "webapp", ROOT / "tools"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import meshcount as MC                                      # noqa: E402
from branch_cut import load_mat_mesh                         # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def refused(data, ext):
    try:
        MC.count_faces(data, ext)
    except MC.MeshCountError as e:
        return str(e)
    return None


# ── every format against the real loader ─────────────────────────────────────
print("\nthe count equals what the loader builds")
base = load_mat_mesh(ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat")
quads = trimesh.Trimesh(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                        faces=[[0, 1, 2], [0, 2, 3]])
built = {
    "binary STL": (base.export(file_type="stl"), ".stl"),
    "ASCII STL": (base.export(file_type="stl_ascii"), ".stl"),
    "ASCII PLY": (base.export(file_type="ply", encoding="ascii"), ".ply"),
    "binary PLY": (base.export(file_type="ply", encoding="binary"), ".ply"),
    "OBJ": (base.export(file_type="obj"), ".obj"),
}
for label, (blob, ext) in built.items():
    blob = blob.encode() if isinstance(blob, str) else blob
    got = MC.count_faces(blob, ext).faces
    loaded = trimesh.load(io.BytesIO(blob), file_type=ext[1:], process=False)
    check(f"{label}: {got:,} == loader's {len(loaded.faces):,}", got == len(loaded.faces))

# A binary STL whose 80-byte header happens to start with "solid" -- the case that
# fools the obvious ASCII test.
sneaky = bytearray(built["binary STL"][0])
sneaky[:5] = b"solid"
check("binary STL whose header begins 'solid' is still read as binary",
      MC.count_faces(bytes(sneaky), ".stl").kind == "binary STL")

# An OBJ of quads: each becomes two triangles when read.
quad_obj = b"v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\n"
got = MC.count_faces(quad_obj, ".obj").faces
check(f"OBJ quad counts as {got} triangles, as the loader builds",
      got == len(trimesh.load(io.BytesIO(quad_obj), file_type="obj", process=False).faces))

mats = [ROOT / "Geometries & Meshes" / "mesh" / "KY160_1.mat", ROOT / "JK5_1_M10.mat"]
try:
    from anat_bench import find_cases
    mats += [Path(mp) for _cid, _cl, mp in find_cases(r"Z:/Martino/branches_CAMP/GT")]
except Exception:                                            # noqa: BLE001
    print("  (thoracic cohort not mounted; .mat checked on the repo samples only)")
bad = []
for m in mats:
    if not m.exists():
        continue
    got = MC.count_faces(m.read_bytes(), ".mat")
    real = len(load_mat_mesh(m).faces)
    if got.faces != real:
        bad.append(f"{m.name}: counted {got.faces:,}, loaded {real:,}")
check(f".mat: count == loader on all {len([m for m in mats if m.exists()])} files", not bad,
      "; ".join(bad[:3]))

# ── files that lie about themselves ──────────────────────────────────────────
print("\nfiles that lie about themselves are refused")
liar = bytearray(84 + 50 * 10)
liar[80:84] = (4_000_000_000).to_bytes(4, "little")
msg = refused(bytes(liar), ".stl")
check("binary STL claiming 4,000,000,000 triangles", msg is not None, msg or "")
ply = (b"ply\nformat binary_little_endian 1.0\nelement vertex 3\nproperty float x\n"
       b"property float y\nproperty float z\nelement face 1000000000000\n"
       b"property list uchar int vertex_indices\nend_header\n" + b"\0" * 200)
msg = refused(ply, ".ply")
check("PLY claiming a trillion faces", msg is not None, msg or "")
msg = refused(built["binary STL"][0][:5000], ".stl")
check("truncated binary STL", msg is not None, msg or "")
msg = refused(b"MATLAB 7.3 MAT-file, Platform: x" + b"\0" * 600, ".mat")
check("MATLAB v7.3 (HDF5) refused with a way forward", msg is not None and "-v7" in msg,
      msg or "")
msg = refused(b"this is not a mesh at all\n" * 40, ".stl")
check("text that is not an STL", msg is not None, msg or "")
msg = refused(b"garbage" * 30, ".mat")
check("garbage .mat", msg is not None, msg or "")

# ── the compressed .mat, in flat memory ──────────────────────────────────────
print("\na compressed .mat is counted in flat memory")
CHILD = r'''
import io, sys, threading, time, psutil
import numpy as np, scipy.io as sio
sys.path.insert(0, sys.argv[1])
import meshcount as MC
buf = io.BytesIO()
n = 4_000_000
sio.savemat(buf, {"X_surface_vertices": np.zeros((n, 3)), "X_surface_triangles":
                  np.zeros((n, 3), dtype=np.uint64)}, do_compression=True)
data = buf.getvalue(); del buf
p = psutil.Process(); peak = [p.memory_info().rss]; stop = [False]
def w():
    while not stop[0]:
        peak[0] = max(peak[0], p.memory_info().rss); time.sleep(0.003)
threading.Thread(target=w, daemon=True).start()
base = p.memory_info().rss; t0 = time.time()
c = MC.count_faces(data, ".mat")
dt = time.time() - t0; stop[0] = True
MC.MAT_DECOMPRESS_CAP = 64 << 20
try:
    MC.count_faces(data, ".mat"); capped = "not refused"
except MC.MeshCountError as e:
    capped = "refused"
print(len(data), c.faces, round((peak[0] - base) / 2**20, 1), round(dt, 3), capped)
'''
r = subprocess.run([sys.executable, "-c", CHILD, str(ROOT / "webapp")],
                   capture_output=True, text=True, timeout=600)
try:
    size, faces, mb, secs, capped = r.stdout.split()
    check(f"a {int(size) / 2**20:.2f} MB file declaring {int(faces):,} faces counted "
          f"with {mb} MB extra memory", float(mb) < 32, f"{secs} s")
    check("past the decompression ceiling it is refused, not loaded", capped == "refused")
except ValueError:
    check("pathological .mat subprocess ran", False, (r.stderr or r.stdout)[-400:])

# ── cheap on a large file ────────────────────────────────────────────────────
print("\ncounting is cheap")
big = trimesh.remesh.subdivide(*trimesh.remesh.subdivide(base.vertices, base.faces))
bigmesh = trimesh.Trimesh(*big, process=False)
blob = bigmesh.export(file_type="stl")
t0 = time.time()
got = MC.count_faces(blob, ".stl").faces
dt = time.time() - t0
check(f"binary STL of {got:,} faces ({len(blob) / 2**20:.0f} MB) counted in {dt * 1000:.1f} ms",
      got == len(bigmesh.faces) and dt < 0.5)
ab = bigmesh.export(file_type="stl_ascii")
ab = ab.encode() if isinstance(ab, str) else ab
t0 = time.time()
got = MC.count_faces(ab, ".stl").faces
dt = time.time() - t0
check(f"ASCII STL of {got:,} faces ({len(ab) / 2**20:.0f} MB) counted in {dt:.2f} s",
      got == len(bigmesh.faces) and dt < 5)

print()
if fails:
    print(f"FAILURES ({len(fails)}): {fails}")
    sys.exit(1)
print("ALL PASS")
