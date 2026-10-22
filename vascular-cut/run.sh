#!/usr/bin/env bash
# ===================================================================
#  Vascular cut tool - one-click launcher (macOS / Linux)
#
#    ./run.sh              # first run sets up the env, then launches
#
#  Prefers conda if present, otherwise builds a local venv.
# ===================================================================
set -e
cd "$(dirname "$0")"

# ---- Route 1: conda (preferred on managed machines) ---------------
if command -v conda >/dev/null 2>&1; then
    echo "[setup] conda found."
    if ! conda env list | grep -q '/vascular-cut$\|^vascular-cut '; then
        echo "[setup] Creating conda environment \"vascular-cut\" (first run only)..."
        conda env create -f environment.yml
    fi
    echo "[run] Launching the GUI..."
    exec conda run --no-capture-output -n vascular-cut python branch_gui.py "$@"
fi

# ---- Route 2: local venv + pip -----------------------------------
echo "[setup] conda not found - using a local Python virtual environment."
VENV=".venv"

if [ -x "$VENV/bin/python" ] && [ -f "$VENV/.installed" ]; then
    echo "[run] Launching the GUI..."
    exec "$VENV/bin/python" branch_gui.py "$@"
fi

# Find a compatible interpreter: 3.12, then 3.11, then 3.10.
PYEXE=""
for cand in python3.12 python3.11 python3.10; do
    if command -v "$cand" >/dev/null 2>&1; then PYEXE="$cand"; break; fi
done
# Fall back to python3 only if it is in the supported range.
if [ -z "$PYEXE" ] && command -v python3 >/dev/null 2>&1; then
    if python3 -c 'import sys; raise SystemExit(0 if (3,10) <= sys.version_info < (3,13) else 1)'; then
        PYEXE="python3"
    fi
fi
if [ -z "$PYEXE" ]; then
    echo "[error] No compatible Python found (need 3.10, 3.11, or 3.12)."
    echo "        Install one from https://www.python.org/downloads/  -- or"
    echo "        install Miniforge/Anaconda and re-run this launcher, which"
    echo "        will then use conda automatically."
    exit 1
fi

echo "[setup] Using $PYEXE to create the environment..."
"$PYEXE" -m venv "$VENV"
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -r requirements.txt
touch "$VENV/.installed"
echo "[run] Launching the GUI..."
exec "$VENV/bin/python" branch_gui.py "$@"
