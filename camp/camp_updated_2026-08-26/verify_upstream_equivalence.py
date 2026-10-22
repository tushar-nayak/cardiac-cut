"""Prove that the additions change nothing when they are switched off.

Every claim in CHANGES.md rests on this: with the new options disabled, campsrc.py here
must return exactly what the July 2026 upstream copy returns. If that ever stops being
true, the additions are not additions and the comparison to prior work is broken.

    python verify_upstream_equivalence.py <surface.mat>

Runs the vendored upstream copy in reference/ and this one, on the same surface, with
    split_disconnected=False, fan_thresh=None, orphan_adjacency=False, seed=42
and compares delta K exactly.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OFF = dict(split_disconnected=False, fan_thresh=None, seed=42)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dk(mod, mat, **over):
    kw = dict(trim_frac=0.0, make_figure=False)
    kw.update(over)
    # Upstream has no orphan_adjacency; pass it only where it exists.
    import inspect
    if "orphan_adjacency" in inspect.signature(mod.camp).parameters:
        kw.setdefault("orphan_adjacency", False)
    if "seed" not in inspect.signature(mod.camp).parameters:
        kw.pop("seed", None)
    if "split_disconnected" not in inspect.signature(mod.camp).parameters:
        kw.pop("split_disconnected", None)
    if "fan_thresh" not in inspect.signature(mod.camp).parameters:
        kw.pop("fan_thresh", None)
    out = mod.camp(str(mat), **kw)
    K = np.asarray(out["K_integrated_patch"], float)
    A = np.asarray(out["patch_area_kept"] if "patch_area_kept" in out
                   else out["patch_area"], float)
    w = A / A.sum()
    return float(np.sum(w * (K - np.sum(w * K)) ** 2))


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    mat = Path(sys.argv[1])
    up = load(HERE / "reference" / "upstream_campsrc.py", "upstream_campsrc")
    new = load(HERE / "campsrc.py", "campsrc_updated")

    a = dk(up, mat)
    b = dk(new, mat, **OFF)
    print(f"  upstream                     {a:.12f}")
    print(f"  this copy, additions OFF     {b:.12f}")
    same = abs(a - b) < 1e-12
    print(f"  {'IDENTICAL' if same else 'DIFFERENT -- investigate before using'}")

    c = dk(new, mat, split_disconnected=True, fan_thresh=90.0, seed=42)
    print(f"\n  this copy, additions ON      {c:.12f}"
          f"   ({100 * (c - a) / a:+.1f}% vs upstream)")
    return 0 if same else 1


if __name__ == "__main__":
    raise SystemExit(main())
