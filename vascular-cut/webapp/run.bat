@echo off
REM The plain Streamlit version. run_cut.bat is the one you want.
REM
REM This predates webapp\cut_app.py and has neither the live cut planes nor the 3-D
REM patch view. It is kept because it needs no VTK render window, so it still runs
REM where the other one will not.
REM
REM Settings live in .streamlit\config.toml -- notably address=127.0.0.1, which is
REM what keeps Windows Firewall from prompting. If you ever see a firewall dialog
REM when starting this, something has bound a public interface; do not click Allow,
REM check the config first.

cd /d "%~dp0.."

if not exist ".venv\Scripts\python.exe" (
    echo Could not find .venv\Scripts\python.exe -- run this from the project checkout.
    exit /b 1
)

echo Starting on http://localhost:8501  (Ctrl+C to stop)
".venv\Scripts\python.exe" -m streamlit run "webapp\app.py" %*
