"""Look at the page the way the browser does, instead of guessing.

Three claims I have made and not verified: the stylesheet reaches the DOM, the
rules actually apply, and the build now running is the one on disk. Check all
three against computed styles.
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "webapp"))
sys.path.insert(0, str(ROOT))
# Read the expected build stamp from the app itself, so this check cannot go stale
# and report a false failure the way it just did.
import cut_app as _APP  # noqa: E402
BUILD = _APP.BUILD

import os

# Overridable: trame serves one session per process, so checking a build without
# evicting whoever is on 8080 means starting a second instance on another port.
URL = os.environ.get("VC_URL") or "http://127.0.0.1:8080/"
SHOT = Path(__file__).parent / "ui.png"
fails = []

with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1500, "height": 950})
    errors = []
    page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.goto(URL, wait_until="domcontentloaded", timeout=60000)
    # The app mounts over a websocket, so DOM-ready is far too early. Wait for a
    # real Vuetify root rather than a fixed sleep.
    try:
        page.wait_for_selector(".v-application", timeout=45000)
        page.wait_for_timeout(2500)
    except Exception:
        print("!! .v-application never appeared, the Vue app did not mount")
        print("   body text:", repr(page.inner_text("body"))[:300])
        print("   html head:", page.evaluate("() => document.body.innerHTML.length"),
              "bytes in body")

    print("=" * 74)
    print("is the running build the one on disk?")
    print("=" * 74)
    body = page.inner_text("body")
    print(f"  build stamp {BUILD!r} on screen : {BUILD in body}")
    if BUILD not in body:
        fails.append("build stamp missing, browser is running older code")

    print()
    print("=" * 74)
    print("did the stylesheet reach the browser?")
    print("=" * 74)
    # Served as a linked stylesheet by enable_module, so look at the sheet list
    # rather than for an injected <style> element.
    sheets = page.evaluate(
        "() => [...document.styleSheets].map(s => s.href || 'inline')")
    ours = [s for s in sheets if "app.css" in (s or "")]
    print(f"  stylesheets loaded           : {len(sheets)}")
    print(f"  ours                         : {ours}")
    if not ours:
        fails.append("app.css is not among the loaded stylesheets")

    print()
    print("=" * 74)
    print("do the rules actually apply?")
    print("=" * 74)
    tb = page.evaluate(
        "() => { const el = document.querySelector('.v-toolbar');"
        " return el ? getComputedStyle(el).backgroundColor : null; }")
    # --bar in app.css. Assert it is the styled dark, not Vuetify's default white,
    # rather than pinning an exact shade the palette is allowed to move.
    print(f"  toolbar background           : {tb}   (want a dark --bar, not white)")
    dark = tb and sum(int(n) for n in tb.strip("rgb() ").split(",")[:3]) < 250
    if not dark:
        fails.append(f"toolbar not styled (got {tb})")

    dr = page.evaluate(
        "() => { const el = document.querySelector('.v-navigation-drawer');"
        " return el ? getComputedStyle(el).backgroundColor : null; }")
    print(f"  drawer background            : {dr}   (transparent: the page shows through)")

    fs = page.evaluate(
        "() => { const el = document.querySelector('.v-application');"
        " return el ? getComputedStyle(el).fontSize : null; }")
    print(f"  base font size               : {fs}   (want 14px)")
    if fs != "14px":
        fails.append(f"base font not applied (got {fs})")

    print()
    print("=" * 74)
    print("load a case and check an opening card")
    print("=" * 74)
    try:
        page.get_by_role("button", name="Load").first.click(timeout=8000)
        page.wait_for_timeout(9000)
        cards = page.locator(".op-card")
        n = cards.count()
        print(f"  .op-card elements            : {n}")
        # What classes DO the opening rows carry? Find them by their slider instead.
        actual = page.evaluate(
            "() => [...document.querySelectorAll('.v-navigation-drawer .v-card')]"
            ".slice(0, 3).map(e => e.className)")
        print("  classes on the first cards   :")
        for c in actual:
            print(f"    {c}")
        if n:
            radius = page.evaluate(
                "() => getComputedStyle(document.querySelector('.op-card')).borderRadius")
            border = page.evaluate(
                "() => getComputedStyle(document.querySelector('.op-card')).borderColor")
            print(f"  card border radius           : {radius}   (want 2px)")
            print(f"  card border colour           : {border}")
            mono = page.evaluate(
                "() => { const e = document.querySelector('.op-dia');"
                " return e ? getComputedStyle(e).fontFamily : 'none'; }")
            print(f"  diameter font family         : {mono[:52]}")
            if radius == "0px":
                fails.append(f"card styling not applied (radius {radius})")
        else:
            fails.append("no opening cards rendered after Load")
    except Exception as e:
        print(f"  could not load a case: {type(e).__name__}: {e}")
        fails.append("Load did not work in the browser")

    page.screenshot(path=str(SHOT), full_page=False)
    print(f"\n  screenshot -> {SHOT}")

    if errors:
        print("\n  console errors:")
        for e in errors[:6]:
            print(f"    {e[:120]}")

    b.close()

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
