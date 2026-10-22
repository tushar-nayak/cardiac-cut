@echo off
REM Start the site. This is the one to run.
REM
REM Both sections live here: cut a geometry, and measure delta K over the CAMP
REM partition. Switch between them on the toolbar. run.bat starts a different,
REM plainer app that predates this one.
REM
REM Binds to 127.0.0.1 only. Do not change that without putting authentication in
REM front of it -- there is none, and this machine has a routable address.
REM
REM Set VC_CUT_PORT to run a second copy beside this one; trame serves one session
REM per process, so that is the only way to have two.

cd /d "%~dp0.."

if not exist ".venv\Scripts\python.exe" (
    echo Could not find .venv\Scripts\python.exe -- run this from the project checkout.
    exit /b 1
)

echo Starting on http://localhost:8080   (Ctrl+C to stop)
".venv\Scripts\python.exe" "webapp\cut_app.py" %*
