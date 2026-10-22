# Browser checks for the web app

These drive a real browser. Everything here exists because reading the source was
not enough, each of these caught a failure that source inspection had already
"confirmed" was fine:

| script | what it caught |
|---|---|
| `check_ui.py` | the stylesheet never reaching the browser at all, and opening rows rendering with no class of their own |
| `check_drag.py` | the slider thumb frozen through an entire drag, jumping only on release |
| `check_errs.py` | server-side exceptions on every drag event, invisible from the page as a bare "Exception raised" |
| `check_measure.py` | the measure tools shipping unreachable, and em dashes in the panel text |
| `check_shape.py` | `class_=` reaching no element, so half the stylesheet was dead; and CSS uppercasing δK into ΔK, which is a different symbol |
| `check_upload.py` | starts its own public-mode server: uploads, downloads, storage, and no path or filename leaking to the page |
| `check_privacy.py` | starts its own public-mode server: a second visitor refused before seeing anything, the next visitor's page clean, a stale token and link dead, uploads and loads that finish after their visitor left dropped, reload |

The pattern worth keeping: **measure the rendered result, not the intent.** Two of
my own measurements were wrong before the checks were, hashing the VTK canvas via
`toDataURL` reports "unchanged" even when the scene has moved, because the remote
view does not draw into a canvas JavaScript can read. `check_errs.py` screenshots
the viewport region instead.

## Running them

```powershell
# once
.venv\Scripts\pip install playwright
.venv\Scripts\python -m playwright install chromium

# start the app, then in another shell
.venv\Scripts\python tests\web\check_ui.py
.venv\Scripts\python tests\web\check_drag.py
.venv\Scripts\python tests\web\check_errs.py
.venv\Scripts\python tests\web\check_measure.py
.venv\Scripts\python tests\web\check_shape.py
```

Each prints `ALL PASS` or a list of failures, and drops a screenshot beside itself.

`check_shape.py` runs delta K for real, in a worker process, and takes about a
minute. It deliberately runs at 25,000 faces rather than the reference 100,000, so
the numbers it prints are not comparable to anything; `tests/test_shape_job.py`
checks the values.

Every check honours `VC_URL`, so a build can be checked on another port without
evicting whatever is on 8080:

```powershell
$env:VC_CUT_PORT=8087; .venv\Scripts\python webapp\cut_app.py
# then
$env:VC_URL="http://127.0.0.1:8087/"; .venv\Scripts\python tests\web\check_shape.py
```

## Three things that will waste your time if you do not know them

**Only one browser at a time.** trame's built-in server holds one state shared by
every tab on it, so a tab you left open drives the same case as the check and will
starve it. In public mode (`VC_PUBLIC=1`) the app goes further and admits one browser
per process: a second one gets "This tool is in use" instead of the app.

**Kill the server by its listening PID.** `ps -W` reports two ids for one process,
and killing the wrong one leaves the old build serving. Every "the fix did not take"
in this section so far has been that, not the fix. Use
`netstat -ano | findstr 8080` and `Stop-Process -Id <pid> -Force`.

**Check the build stamp.** The toolbar carries one (`ui-N`). If a check reports
stale-looking results, confirm the stamp on screen matches `cut_app.BUILD` before
diagnosing anything else, a cached page cost several rounds of chasing the wrong
bug.
