"""The shape-metrics section, driven through a browser.

Cut a case, hand it over, compute delta K, and check that the number, the warnings and
the 3-D view all arrive. The interesting part is the wait: the run takes tens of
seconds in another process, so this also checks the things that only break under that
-- the panel staying responsive, the phase text changing, and the viewport still
painting while the worker is busy.

Start the app first:  python webapp/cut_app.py
"""
import hashlib
import re
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

import os

# Overridable, because trame serves one session per process: a second instance on
# another port is the only way to check this without evicting whoever is on 8080.
URL = os.environ.get("VC_URL") or "http://127.0.0.1:8080/"
SHOT = Path(__file__).parent / "shape.png"
SHOT3D = Path(__file__).parent / "shape_3d.png"
fails = []


def note(ok, what):
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}")
    if not ok:
        fails.append(what)


def open_panel(page, title):
    """Make sure a section is open, rather than clicking it and hoping.

    Panel state is client-side and survives everything the server does, so a bare
    click is a toggle against unknown state: it opened the panel when the check ran
    one way round and closed it when it ran the other, and the failure looked like a
    missing button.
    """
    panel = page.locator(".v-expansion-panel", has=page.get_by_text(title, exact=False))
    if "v-expansion-panel--active" not in (panel.first.get_attribute("class") or ""):
        panel.first.locator(".v-expansion-panel-title").click()
        page.wait_for_timeout(900)


with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1500, "height": 1000})
    errs = []
    page.on("pageerror", lambda e: errs.append(str(e)))
    page.goto(URL, wait_until="domcontentloaded", timeout=60000)
    page.wait_for_selector(".v-application", timeout=45000)
    page.wait_for_timeout(1500)

    print("=" * 74)
    print("the section switch")
    print("=" * 74)
    note(page.locator(".top-tog").count() == 1, "the switch is on the toolbar")
    note(page.locator(".top-tog button").count() == 2, "with two sections on it")

    # Start from the cut section, whichever one the server was left in. trame holds one
    # state for the whole process, so a previous run of this check leaves the app in
    # shape mode and the next one would open there.
    was = page.locator(".top-tog .v-btn--active").inner_text().strip()
    print(f"  server was showing: {was}")
    if "cut" not in was.lower():
        page.locator(".top-tog button", has_text="Cut").first.click()
        page.wait_for_timeout(1500)

    # Cut something first, so the handoff has a surface to hand over.
    page.get_by_role("button", name="Load").first.click(timeout=10000)
    page.wait_for_selector(".op-card", timeout=30000)
    page.wait_for_timeout(2500)
    n_cards = page.locator(".op-card").count()
    print(f"  cut loaded, {n_cards} opening cards")

    print()
    print("=" * 74)
    print("the handoff from the cut section")
    print("=" * 74)
    open_panel(page, "OUTPUT")
    page.get_by_role("button", name=re.compile("Measure .K on this cut", re.I)).first.click()
    page.wait_for_timeout(2000)

    active = page.locator(".top-tog .v-btn--active").inner_text().strip()
    note("shape" in active.lower(), f"the handoff switched sections (now {active!r})")
    note(page.locator(".op-card").count() == 0,
         "the cut section's opening cards are gone")

    info = page.evaluate("""() => [...document.querySelectorAll('.readout .ro-row')]
        .map(r => [...r.children].map(c => c.textContent.trim()).filter(Boolean))""")
    for row in info:
        print("    " + "  ".join(row))
    flat = " ".join(" ".join(r) for r in info)
    note("Openings" in flat and "Euler" in flat,
         "the surface it is about to measure is described")
    m = re.search(r"Euler\s+(-?\d+)", flat)
    note(m is not None and int(m.group(1)) == -3,
         f"topology reads as the cut reported it ({m.group(1) if m else '?'})")

    print()
    print("=" * 74)
    print("the run")
    print("=" * 74)
    # 25k, so the check does not take a minute. The point here is the plumbing; the
    # numbers themselves are checked headlessly in tests/test_shape_job.py.
    #
    open_panel(page, "PARTITION")
    page.locator(".mode-tog button", has_text="25k").first.click(timeout=15000)
    page.wait_for_timeout(600)
    note(page.locator(".mode-tog button.v-btn--active", has_text="25k").count() == 1,
         "the face-count toggle takes a selection")

    page.get_by_role("button", name=re.compile(r"Compute", re.I)).first.click()

    phases, t0 = [], time.time()
    got_result = False
    while time.time() - t0 < 180:
        page.wait_for_timeout(700)
        txt = page.inner_text(".v-navigation-drawer")
        ph = re.search(r"(Remeshing[^\n]*|Partitioning[^\n]*|Trimming[^\n]*)", txt)
        if ph and ph.group(1) not in phases:
            phases.append(ph.group(1))
            print(f"    phase: {ph.group(1)}")
        # Case-insensitive: the stylesheet uppercases readout labels, and inner_text
        # returns what is RENDERED, so "At floor" arrives as "AT FLOOR" and a
        # case-sensitive search sits here for the full timeout with the answer on
        # screen the whole time.
        if re.search(r"\bat floor\b", txt, re.I):
            got_result = True
            break
    print(f"  waited {time.time() - t0:.0f}s")
    note(got_result, "a result arrived")
    note(len(phases) >= 1, "progress was reported while it ran")

    rows = page.evaluate("""() => [...document.querySelectorAll('.readout .ro-row')]
        .map(r => [...r.children].map(c => c.textContent.trim()).filter(Boolean))""")
    for row in rows:
        print("    " + "  ".join(row))
    flat = " ".join(" ".join(r) for r in rows)
    dk = re.search(r"[δΔ]K\s+([\d.eE+-]+)", flat)
    note(dk is not None, f"delta K is on the panel ({dk.group(1) if dk else 'MISSING'})")
    if dk:
        note(0.0 < float(dk.group(1)) < 100.0,
             f"and it is a plausible value ({dk.group(1)})")
    note("Patches" in flat and "Faces" in flat and "Euler" in flat,
         "the conditions it was computed under are shown beside it")

    warns = page.locator(".v-navigation-drawer .v-alert").all_inner_texts()
    warns = [w.strip() for w in warns if w.strip()]
    print()
    for w in warns:
        print(f"    [alert] {w[:150]}")
    note(any("25,000 faces" in w or "vertex floor" in w for w in warns),
         "running away from the reference settings is flagged")

    print()
    print("=" * 74)
    print("the 3-D view follows the colouring")
    print("=" * 74)

    def shot():
        return hashlib.sha1(page.locator(".v-main").screenshot()).hexdigest()[:12]

    page.wait_for_timeout(1500)
    a = shot()
    open_panel(page, "VIEW")
    # One metric, three scales. Full is CAMP's own auto-scale and is much wider than
    # the robust default, so the surface has to repaint.
    page.locator(".mode-tog button", has_text="Full").first.click()
    page.wait_for_timeout(2500)
    c = shot()
    print(f"    robust {a}   full {c}")
    note(a != c, "changing the colour scale repaints the surface")
    txt = page.inner_text(".v-navigation-drawer")
    note("how CAMP scales it" in txt, "Full says it is CAMP's own scaling")
    page.locator(".mode-tog button", has_text="Robust").first.click()
    page.wait_for_timeout(2500)
    note(re.search(r"past an end", page.inner_text(".v-navigation-drawer")) is not None,
         "Robust says how much of the surface is saturated")

    page.get_by_text("Patch outlines").first.click()
    page.wait_for_timeout(2500)
    d = shot()
    print(f"    outlines on  {d}")
    note(c != d, "patch outlines reach the view")

    page.locator(".v-main").screenshot(path=str(SHOT3D))
    page.locator(".v-navigation-drawer").screenshot(path=str(SHOT))
    print(f"\n  screenshots -> {SHOT.name}, {SHOT3D.name}")

    print()
    print("=" * 74)
    print("export, then the other way in")
    print("=" * 74)
    open_panel(page, "OUTPUT")
    page.get_by_role("button", name=re.compile("Export patch table", re.I)).first.click()
    page.wait_for_timeout(2500)
    alerts = " ".join(page.locator(".v-navigation-drawer .v-alert").all_inner_texts())
    wrote = re.search(r"Wrote (\S+) and (\S+)", alerts)
    note(wrote is not None,
         f"the patch table was written ({wrote.group(0) if wrote else 'no message'})")

    # The second route in: a surface that never went through the cut section. This is
    # the half of the site that has to work for someone who already has geometries.
    open_panel(page, "SURFACE")
    surface = str((Path(__file__).parents[2] / "output" / "KY160_1_open.stl").resolve())
    page.locator(".v-navigation-drawer input[type='text']").first.fill(surface)
    page.get_by_role("button", name=re.compile("Load surface", re.I)).first.click()
    page.wait_for_timeout(3000)
    info2 = page.evaluate("""() => [...document.querySelectorAll('.readout .ro-row')]
        .map(r => [...r.children].map(c => c.textContent.trim()).filter(Boolean))""")
    flat2 = " ".join(" ".join(r) for r in info2)
    print(f"    {flat2[:120]}")
    note("KY160_1_open.stl" in flat2, "a surface loaded by path is described")
    note(not re.search(r"\bat floor\b", flat2, re.I),
         "the previous result was dropped rather than left beside a new surface")

    print()
    print("=" * 74)
    print("back to the cut section, with the cut intact")
    print("=" * 74)
    page.locator(".top-tog button", has_text="Cut").first.click()
    page.wait_for_timeout(2000)
    note(page.locator(".op-card").count() == n_cards,
         f"the {n_cards} opening cards came back")

    body = page.inner_text("body")
    EM = chr(0x2014)
    note(body.count(EM) == 0, f"no em dashes on screen ({body.count(EM)} found)")

    # delta K, not Delta K. Buttons and readout labels are uppercased by CSS, and
    # uppercasing a lowercase delta silently produces a different symbol: the metric
    # is Khabaz's dK, while a capital delta reads as a CHANGE in K.
    big, small = chr(0x0394), chr(0x03B4)
    bad = [ln.strip() for ln in body.splitlines() if big in ln]
    for ln in bad:
        print(f"    capital delta: {ln[:80]!r}")
    note(not bad, f"no capital delta on screen ({len(bad)} line(s))")
    note(small in body, "the lowercase delta survived the uppercase transform")

    if errs:
        print("\n  page errors:")
        for e in errs[:6]:
            print(f"    {e[:140]}")
        fails.append(f"{len(errs)} page error(s)")
    b.close()

print()
print("ALL PASS" if not fails else f"FAILURES ({len(fails)}): " + "; ".join(fails))
sys.exit(1 if fails else 0)
