"""Freeze a folder of cut geometries, so a later run cannot quietly replace them.

    python tools/arch_lock.py lock                     # freeze output/arch_final
    python tools/arch_lock.py verify                   # has anything moved?
    python tools/arch_lock.py unlock                   # deliberate, and says so

WHY A LOCK AND NOT A BACKUP

A backup answers "what did it used to be". It does not answer the question that
actually bites, which is "is the file I am measuring the one I reviewed". These
meshes are reviewed by eye, one at a time, and then fed to a curvature measure whose
numbers are meaningless across cut extents -- so a single silently re-cut case does
not announce itself as an error, it announces itself as an outlier, months later, in
a plot.

So the lock records a SHA-256 per file and the provenance of the set: the commit that
produced it, the settings it was cut at, and what each case's two region boundaries
came out at. verify re-hashes and reports anything added, missing or altered.

The files are also made read-only. That stops the accident -- a batch re-run pointed
at the wrong --export -- rather than a determined overwrite, which is what `unlock` is
for. output/ is gitignored, so the MANIFEST is the part that travels with the repo and
is what a reader checks a folder against.
"""
import argparse
import hashlib
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIR = ROOT / "output" / "arch_final"
MANIFEST = "MANIFEST.sha256"
PROVENANCE = "PROVENANCE.json"


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def members(folder):
    """Every delivered file, in a stable order. The manifest is not a member."""
    return sorted(p for p in folder.iterdir()
                  if p.is_file() and p.name not in (MANIFEST, PROVENANCE)
                  and p.name != "review.json")


def read_manifest(folder):
    f = folder / MANIFEST
    if not f.exists():
        return None
    out = {}
    for line in f.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        h, _, name = line.partition("  ")
        out[name.strip()] = h.strip()
    return out


def git(*args):
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:                                               # noqa: BLE001
        return ""


def bounds_of(folder):
    """Each case's two region boundaries, as the export recorded them.

    Carried in the provenance because it is the cheapest way for a reader to tell
    WHICH cut a file is without loading the mesh, and because these are the numbers
    the set was accepted on.
    """
    out = {}
    for js in sorted(folder.glob("*_cut_planes.json")):
        try:
            reg = (json.loads(js.read_text(encoding="utf-8")).get("region") or {})
        except Exception:                                           # noqa: BLE001
            continue
        if reg.get("start_mm_along_aorta") is None:
            continue
        out[js.name.replace("_cut_planes.json", "")] = {
            "start_mm_along_aorta": reg["start_mm_along_aorta"],
            "end_mm_along_aorta": reg.get("end_mm_along_aorta"),
            "length_mm": reg.get("length_mm"),
        }
    return out


def set_readonly(paths, on=True):
    n = 0
    for p in paths:
        try:
            mode = os.stat(p).st_mode
            os.chmod(p, (mode & ~stat.S_IWRITE) if on
                     else (mode | stat.S_IWRITE))
            n += 1
        except OSError:
            pass
    return n


def cmd_lock(args):
    folder = Path(args.dir)
    files = members(folder)
    if not files:
        print(f"nothing to lock in {folder}")
        return 1
    lines = [f"{digest(p)}  {p.name}" for p in files]
    (folder / MANIFEST).write_text("\n".join(lines) + "\n", encoding="utf-8")
    prov = {
        "folder": folder.name,
        "n_files": len(files),
        "git_commit": git("rev-parse", "HEAD"),
        "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain")),
        "note": args.note or "",
        "bounds": bounds_of(folder),
    }
    (folder / PROVENANCE).write_text(json.dumps(prov, indent=2), encoding="utf-8")
    n = set_readonly(files, on=True)
    print(f"  locked {len(files)} file(s) in {folder}")
    print(f"  {MANIFEST} + {PROVENANCE} written, {n} made read-only")
    if prov["git_dirty"]:
        print("  NOTE: the working tree is dirty, so git_commit alone does not "
              "reproduce this set")
    return 0


def cmd_verify(args):
    folder = Path(args.dir)
    want = read_manifest(folder)
    if want is None:
        print(f"no {MANIFEST} in {folder} -- nothing is locked here")
        return 1
    have = {p.name: p for p in members(folder)}
    changed, missing, added = [], [], []
    for name, h in want.items():
        p = have.get(name)
        if p is None:
            missing.append(name)
        elif digest(p) != h:
            changed.append(name)
    for name in have:
        if name not in want:
            added.append(name)
    print(f"  {len(want)} file(s) in the manifest, {len(have)} on disk")
    for tag, names in (("CHANGED", changed), ("MISSING", missing),
                       ("ADDED (not in the manifest)", added)):
        for n in sorted(names):
            print(f"    {tag}: {n}")
    if not (changed or missing or added):
        print("  intact -- every file matches the hash it was locked at")
        return 0
    return 1


def cmd_unlock(args):
    folder = Path(args.dir)
    n = set_readonly(members(folder), on=False)
    print(f"  {n} file(s) writable again in {folder}")
    print(f"  {MANIFEST} left in place: verify still reports what changed")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("lock", "verify", "unlock"):
        p = sub.add_parser(nm)
        p.add_argument("--dir", default=str(DEFAULT_DIR))
        if nm == "lock":
            p.add_argument("--note", default=None,
                           help="what this set is, in one line")
    args = ap.parse_args()
    return {"lock": cmd_lock, "verify": cmd_verify, "unlock": cmd_unlock}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
