"""Build the collaborator distribution: the desktop tool, and nothing else.

This is what goes on the share drive. It is NOT a snapshot of the repository. A
collaborator who wants to cut a geometry should not have to look at the web front
end, the vendored CAMP source, the browser checks or the integration notes, and
should not be asked to install trame, Streamlit, plotly or scikit-learn to open a
window.

    python make_dist.py            -> dist/vascular-cut-<date>.zip

Use `git archive` instead if you want the whole repository.
"""
import re
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "dist"
PREFIX = "vascular-cut"

# The program, what it imports, how to launch it, and something to open.
FILES = [
    "branch_gui.py",            # entry point
    "branch_cut.py",            # cutting, measurement, validation, mesh I/O
    "branch_tree.py",           # centerline parsing and cut-plane geometry
    "cutclean.py",              # remeshing, imported lazily by branch_cut
    "cl.py",                    # centerline generation from a volume
    "run.bat",
    "run.sh",
    "environment.yml",
    "pyproject.toml",
    "README.md",
    "docs/User-Guide.pdf",
    "docs/Whats-New.pdf",
    "docs/screenshot.png",
    "docs/annotated.png",
    "Geometries & Meshes/mesh/KY160_1.mat",
    "KY160_1_CentrelinePoints.csv",
    "KK2_5.mat",
    "KK2_5.csv",
]

# Everything after this marker in requirements.txt is for the web front end and the
# curvature metrics. Shipping it makes a desktop install several hundred megabytes
# larger and gives a collaborator four more ways for pip to fail.
WEB_MARKER = "# Web front end."


def desktop_requirements():
    text = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    core = text.split(WEB_MARKER)[0].rstrip()
    return (core + "\n\n"
            "# The web front end and the curvature metrics need more than this.\n"
            "# They are not part of this distribution; see the project repository.\n")


def build():
    missing = [f for f in FILES if not (ROOT / f).exists()]
    if missing:
        raise SystemExit("missing from the checkout:\n  " + "\n  ".join(missing))

    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / f"{PREFIX}-{date.today():%Y-%m-%d}.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for f in FILES:
            z.write(ROOT / f, f"{PREFIX}/{f}")
        z.writestr(f"{PREFIX}/requirements.txt", desktop_requirements())

    with zipfile.ZipFile(out) as z:
        names = z.namelist()
        total = sum(i.file_size for i in z.infolist())
    print(f"wrote {out}")
    print(f"  {len(names)} files, {out.stat().st_size / 1e6:.2f} MB "
          f"({total / 1e6:.2f} MB uncompressed)")
    for n in sorted(names):
        print(f"    {n}")
    return out


if __name__ == "__main__":
    build()
