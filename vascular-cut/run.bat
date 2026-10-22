@echo off
REM ===================================================================
REM  Vascular cut tool - one-click launcher (Windows)
REM
REM  Double-click this file, or run it from a terminal.
REM  First run sets up the environment (a few minutes); later runs are
REM  instant. Prefers conda if present, otherwise builds a local venv.
REM ===================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

REM ---- Find conda even when it is not on PATH -----------------------
REM
REM  The Miniforge and Miniconda installers deliberately leave conda OFF the
REM  system PATH -- they expect you to launch their own prompt shortcut -- so
REM  "where conda" finds nothing on a machine that plainly has conda, and this
REM  launcher fell through to building a venv, or failed. Two people hit it; one
REM  fixed it by hard-coding his own install path into this file, which then
REM  worked only on his account. So look for it in the places it is installed.
set "CONDA_ACT="
where conda >nul 2>nul
if %errorlevel% neq 0 (
    for %%R in (
        "%CONDA_ROOT%"
        "%USERPROFILE%\miniforge3"        "%LOCALAPPDATA%\miniforge3"
        "%USERPROFILE%\mambaforge"        "%LOCALAPPDATA%\mambaforge"
        "%USERPROFILE%\miniconda3"        "%LOCALAPPDATA%\miniconda3"
        "%USERPROFILE%\anaconda3"         "%LOCALAPPDATA%\anaconda3"
        "%USERPROFILE%\AppData\Local\Continuum\anaconda3"
        "%ProgramData%\miniforge3"        "%ProgramData%\Miniconda3"
        "%ProgramData%\Anaconda3"         "%ProgramFiles%\Miniforge3"
        "%ProgramFiles%\Miniconda3"       "%ProgramFiles%\Anaconda3"
    ) do (
        if not defined CONDA_ACT if exist "%%~R\condabin\activate.bat" (
            set "CONDA_ACT=%%~R\condabin\activate.bat"
        )
    )
    REM CONDA_EXE is set on any machine where conda init has ever run, and points
    REM into the install, so it finds one in a place this list has never heard of.
    if not defined CONDA_ACT if defined CONDA_EXE (
        for %%E in ("%CONDA_EXE%") do (
            if exist "%%~dpE..\condabin\activate.bat" (
                set "CONDA_ACT=%%~dpE..\condabin\activate.bat"
            )
        )
    )
    if defined CONDA_ACT (
        echo [setup] conda is installed but not on PATH.
        echo         Initialising from "!CONDA_ACT!"
        call "!CONDA_ACT!" >nul 2>nul
    )
)

REM ---- Route 1: conda (preferred on managed machines) ---------------
where conda >nul 2>nul
if %errorlevel%==0 (
    echo [setup] conda found.
    call conda env list | findstr /C:"vascular-cut" >nul 2>nul
    if !errorlevel! neq 0 (
        echo [setup] Creating conda environment "vascular-cut"...
        echo         ^(first run only - this downloads Python + packages^)
        call conda env create -f environment.yml
        if !errorlevel! neq 0 (
            echo [error] conda env create failed. See messages above.
            pause
            exit /b 1
        )
    )
    echo [run] Launching the GUI...
    call conda run --no-capture-output -n vascular-cut python branch_gui.py %*
    goto :eof
)

REM ---- Route 2: local venv + pip -----------------------------------
echo [setup] conda not found on PATH, and not in any of the usual install
echo         locations - using a local Python virtual environment instead.
echo         ^(If you DO have conda, open its own prompt and run this from
echo         there, or set CONDA_ROOT to the install folder.^)

set "VENV=.venv"
if exist "%VENV%\Scripts\python.exe" if exist "%VENV%\.installed" (
    echo [run] Launching the GUI...
    "%VENV%\Scripts\python.exe" branch_gui.py %*
    goto :eof
)

REM Find a compatible interpreter: 3.12, then 3.11, then 3.10.
set "PYEXE="
for %%V in (3.12 3.11 3.10) do (
    if not defined PYEXE (
        py -%%V -c "import sys" >nul 2>nul && set "PYEXE=py -%%V"
    )
)
if not defined PYEXE (
    echo [error] No compatible Python found ^(need 3.10, 3.11, or 3.12^).
    echo         Install one from https://www.python.org/downloads/  -- or
    echo         install Miniforge/Anaconda and re-run this launcher, which
    echo         will then use conda automatically.
    pause
    exit /b 1
)
echo [setup] Using %PYEXE% to create the environment...
%PYEXE% -m venv "%VENV%" || (echo [error] venv creation failed & pause & exit /b 1)
"%VENV%\Scripts\python.exe" -m pip install --upgrade pip
"%VENV%\Scripts\python.exe" -m pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo [error] Dependency install failed. See messages above.
    pause
    exit /b 1
)
echo installed> "%VENV%\.installed"
echo [run] Launching the GUI...
"%VENV%\Scripts\python.exe" branch_gui.py %*
