"""How many faces a mesh has, read WITHOUT loading it.

A face-count limit is only a protection if it is checked before the mesh is built, and
a byte cap is no substitute for one. A compressed MATLAB file of 0.14 MB was measured
taking 225 MB of memory to load, about 1,600 times its size, because zlib squeezes a
large array of repeated values to almost nothing. So every upload is counted here first,
from the file's own structure, in memory that does not grow with what the file claims.

scipy.io.whosmat looks like the obvious tool for .mat and is not one: it decompresses
each variable in full to report its shape. On a 0.18 MB file holding two 4-million-row
arrays it used 122 MB -- and within a 50 MB upload cap, a crafted file could ask for tens
of gigabytes that way. The reader below decompresses in a stream, keeps only the first
few hundred bytes of each variable (its header: name, class, dimensions), discards the
rest as it goes, and stops at a hard ceiling.

Every count is checked against the real loader in tests/test_meshcount.py.

    count_faces(data, ".stl")  ->  MeshCount(faces=..., vertices=..., kind=...)
    raises MeshCountError, whose message is safe to show to the person who uploaded.
"""
import re
import struct
import zlib
from dataclasses import dataclass

# How much a single .mat file may decompress to, in total, before it is refused. A mesh
# at any sensible face limit needs a small fraction of this: a million-face surface is
# about 36 MB of arrays. It exists to bound the TIME and memory a lying file can cost.
MAT_DECOMPRESS_CAP = 1 << 30            # 1 GiB
MAT_HEADER_PEEK = 512                   # bytes of each variable kept for its header
_CHUNK = 1 << 20


class MeshCountError(ValueError):
    """Said to the person, verbatim. Never carries a server path."""


@dataclass
class MeshCount:
    faces: int              # triangles the loader will build (an upper bound for FE files)
    vertices: int | None
    kind: str               # e.g. "binary STL", "MATLAB surface mesh"
    alloc_bytes: int | None = None   # arrays the loader will allocate (.mat only)


# ── STL ──────────────────────────────────────────────────────────────────────
def _stl(data):
    n = len(data)
    if n >= 84:
        claimed = struct.unpack_from("<I", data, 80)[0]
        # A binary STL is EXACTLY 84 + 50 bytes per triangle, so the header can be
        # checked against the file rather than believed. This test is also how binary
        # is told from ASCII: a binary STL may begin with the word "solid" in its
        # 80-byte header, which is what fools the naive check.
        if n == 84 + 50 * claimed:
            if claimed == 0:
                raise MeshCountError("This STL file holds no triangles.")
            return MeshCount(claimed, None, "binary STL")
    head = data[:512].lstrip().lower()
    if head.startswith(b"solid"):
        # "endfacet" rather than "facet": the one word closes every triangle and cannot
        # be confused with anything else in the file.
        count = sum(1 for _ in re.finditer(rb"(?i)\bendfacet\b", data))
        if count:
            return MeshCount(count, None, "ASCII STL")
        raise MeshCountError("This STL file holds no triangles.")
    # Only quote the header's claim if the file could plausibly BE a binary STL. Text
    # read as one yields a nonsense count, and "its header claims 1,763,734,377
    # triangles" about a text file is true and useless.
    printable = sum(32 <= b < 127 or b in (9, 10, 13) for b in data[:512])
    if n >= 84 and printable < 0.9 * min(n, 512):
        claimed = struct.unpack_from("<I", data, 80)[0]
        raise MeshCountError(
            f"This does not read as an STL file. Its header claims {claimed:,} "
            f"triangles, which would need {84 + 50 * claimed:,} bytes, and the file "
            f"is {n:,}. It may be truncated or not an STL at all.")
    if n >= 84:
        raise MeshCountError("This does not read as an STL file, binary or ASCII.")
    raise MeshCountError("This file is too short to be an STL.")


# ── PLY ──────────────────────────────────────────────────────────────────────
_PLY_SIZES = {"char": 1, "uchar": 1, "int8": 1, "uint8": 1, "short": 2, "ushort": 2,
              "int16": 2, "uint16": 2, "int": 4, "uint": 4, "int32": 4, "uint32": 4,
              "float": 4, "float32": 4, "double": 8, "float64": 8}


def _ply(data):
    end = data.find(b"end_header")
    if not data[:3].lower() == b"ply" or end < 0 or end > 1 << 16:
        raise MeshCountError("This is not a PLY file: it has no PLY header.")
    header = data[:end].decode("ascii", errors="replace").splitlines()
    fmt, faces, verts = None, None, None
    for line in header:
        tok = line.split()
        if len(tok) >= 2 and tok[0] == "format":
            fmt = tok[1]
        elif len(tok) >= 3 and tok[0] == "element" and tok[1] == "face":
            faces = int(tok[2])
        elif len(tok) >= 3 and tok[0] == "element" and tok[1] == "vertex":
            verts = int(tok[2])
    if faces is None:
        raise MeshCountError("This PLY file declares no faces.")
    if faces <= 0:
        raise MeshCountError("This PLY file holds no faces.")
    body = len(data) - end
    # Check the claim against the bytes that are actually there. The smallest possible
    # face is a one-byte vertex count and three one-byte indices in binary, and at
    # least "3 0 0 0" plus a line break in ASCII -- so a header claiming more faces than
    # that could fit is a lie, and is refused before anything is allocated for it.
    smallest = 4 if (fmt or "").startswith("binary") else 8
    if faces * smallest > body:
        raise MeshCountError(
            f"This PLY file claims {faces:,} faces, more than its {len(data):,} bytes "
            f"could hold. It may be truncated or damaged.")
    return MeshCount(faces, verts, f"PLY ({fmt or 'unknown format'})")


# ── OBJ ──────────────────────────────────────────────────────────────────────
def _obj(data):
    # A polygon of k vertices becomes k - 2 triangles when it is read, so count those
    # rather than face lines: a quad mesh would otherwise be counted at half its cost.
    tri, verts = 0, 0
    for line in data.splitlines():
        if line[:2] == b"f " or line[:2] == b"f\t":
            k = len(line.split()) - 1
            if k >= 3:
                tri += k - 2
        elif line[:2] == b"v " or line[:2] == b"v\t":
            verts += 1
    if not tri:
        raise MeshCountError("This OBJ file holds no faces.")
    return MeshCount(tri, verts, "OBJ")


# ── MATLAB v5 / v7 ───────────────────────────────────────────────────────────
# The format, as MathWorks documents it: a 128-byte header, then a run of data
# elements, each an 8-byte tag (type, byte count) and its data. A variable is a
# miMATRIX element, or a miCOMPRESSED element whose zlib stream holds one. Inside a
# miMATRIX the first three sub-elements are the array flags (class in the low byte),
# the dimensions, and the name -- all in the first few hundred bytes.
_MI_INT32, _MI_UINT32, _MI_MATRIX, _MI_COMPRESSED = 5, 6, 14, 15
_CLASS_ITEMSIZE = {6: 8, 7: 4, 8: 1, 9: 1, 10: 2, 11: 2, 12: 4, 13: 4, 14: 8, 15: 8,
                   4: 2}           # double, single, int8..uint64, char
_CELL, _STRUCT = 1, 2


def _tag(buf, pos, endian):
    """(type, nbytes, header_len) for the data element tag at pos."""
    first = struct.unpack_from(endian + "I", buf, pos)[0]
    if first >> 16:                      # the small data element format
        return first & 0xFFFF, first >> 16, 4
    return first, struct.unpack_from(endian + "I", buf, pos + 4)[0], 8


def _matrix_header(buf, endian):
    """(name, class, dims) from the start of a miMATRIX element's data."""
    pos = 0
    t, n, h = _tag(buf, pos, endian)            # array flags
    if t != _MI_UINT32:
        raise MeshCountError("This MATLAB file is not laid out as the format requires.")
    cls = struct.unpack_from(endian + "I", buf, pos + h)[0] & 0xFF
    pos += h + ((n + 7) // 8) * 8 if h == 8 else 8
    t, n, h = _tag(buf, pos, endian)            # dimensions
    if t != _MI_INT32 or n % 4 or n > 64:
        raise MeshCountError("This MATLAB file is not laid out as the format requires.")
    dims = list(struct.unpack_from(endian + f"{n // 4}i", buf, pos + h))
    pos += h + ((n + 7) // 8) * 8 if h == 8 else 8
    t, n, h = _tag(buf, pos, endian)            # name
    name = bytes(buf[pos + h: pos + h + n]).decode("ascii", errors="replace")
    return name, cls, dims


def _mat(data):
    if len(data) < 128:
        raise MeshCountError("This file is too short to be a MATLAB file.")
    if data[:19].lower().startswith(b"matlab 7.3") or data[512:516] == b"\x89HDF":
        raise MeshCountError("MATLAB v7.3 files are not supported. Save the mesh with "
                             "save(..., '-v7') and upload that.")
    ind = bytes(data[126:128])
    endian = "<" if ind == b"IM" else ">" if ind == b"MI" else None
    if endian is None:
        raise MeshCountError("This does not read as a MATLAB .mat file.")
    variables, total_out, pos = [], 0, 128
    while pos + 8 <= len(data):
        t, n, h = _tag(data, pos, endian)
        body = memoryview(data)[pos + h: pos + h + n]
        pos += h + (n if t == _MI_COMPRESSED else ((n + 7) // 8) * 8)
        if t == _MI_COMPRESSED:
            d = zlib.decompressobj()
            peek, out = b"", 0
            for i in range(0, len(body), _CHUNK):
                chunk = d.decompress(body[i: i + _CHUNK], _CHUNK)
                while True:
                    if len(peek) < MAT_HEADER_PEEK:
                        peek += chunk[: MAT_HEADER_PEEK - len(peek)]
                    out += len(chunk)
                    if total_out + out > MAT_DECOMPRESS_CAP:
                        raise MeshCountError(
                            "This MATLAB file expands to more than "
                            f"{MAT_DECOMPRESS_CAP / 2 ** 30:.0f} GB when decompressed, "
                            "far beyond any vascular surface. It was not loaded.")
                    if not d.unconsumed_tail:
                        break
                    chunk = d.decompress(d.unconsumed_tail, _CHUNK)
            total_out += out
            if len(peek) < 8:
                continue
            t2, n2, h2 = _tag(peek, 0, endian)
            if t2 != _MI_MATRIX:
                continue
            variables.append(_matrix_header(peek[h2:], endian))
        elif t == _MI_MATRIX:
            variables.append(_matrix_header(bytes(body[:MAT_HEADER_PEEK]), endian))
        if len(variables) > 256:
            raise MeshCountError("This MATLAB file holds too many variables to be a mesh.")
    return _mat_mesh(variables)


def _prod(dims):
    p = 1
    for d in dims:
        p *= max(int(d), 0)
    return p


def _mat_mesh(variables):
    """Which variables are the mesh, and what loading it will cost.

    Mirrors branch_cut.load_mat_mesh: a Simpleware SURFACE export holds
    <name>_surface_vertices and <name>_surface_triangles; a SOLID finite-element export
    holds <name>_FE_vertices and tetrahedra or hexahedra, from which the loader
    extracts the wall. Its face count is bounded by every element's faces, since the
    loader builds all of them before keeping the boundary.
    """
    by_name = {name: (cls, dims) for name, cls, dims in variables}
    alloc = sum(_prod(d) * _CLASS_ITEMSIZE.get(c, 8) for c, d in by_name.values()
                if c not in (_CELL, _STRUCT))
    for name, (cls, dims) in by_name.items():
        if name.endswith("_surface_triangles") and len(dims) == 2:
            v = by_name.get(name[: -len("_surface_triangles")] + "_surface_vertices")
            return MeshCount(dims[0], v[1][0] if v else None,
                             "MATLAB surface mesh", alloc)
    for suffix, per in (("_FE_tetrahedra", 4), ("_FE_hexahedra", 12)):
        for name, (cls, dims) in by_name.items():
            if name.endswith(suffix) and len(dims) == 2:
                v = by_name.get(name[: -len(suffix)] + "_FE_vertices")
                return MeshCount(dims[0] * per, v[1][0] if v else None,
                                 "MATLAB solid mesh (upper bound)", alloc)
    raise MeshCountError(
        "This MATLAB file holds no mesh this tool can read: it expects the "
        "_surface_vertices and _surface_triangles arrays that Simpleware exports.")


_READERS = {".stl": _stl, ".ply": _ply, ".obj": _obj, ".mat": _mat}


def count_faces(data, ext):
    """MeshCount for a mesh file's bytes, or MeshCountError. Loads nothing."""
    reader = _READERS.get(str(ext).lower())
    if reader is None:
        raise MeshCountError(f"No face count is available for {ext} files.")
    try:
        return reader(bytes(data) if not isinstance(data, (bytes, bytearray)) else data)
    except MeshCountError:
        raise
    except (struct.error, zlib.error, ValueError, IndexError, UnicodeError):
        raise MeshCountError(f"This {ext[1:].upper()} file is damaged or truncated, "
                             "so it could not be read.")
