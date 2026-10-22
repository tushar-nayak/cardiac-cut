"""Do the measure tools actually work in the browser?

They were never ported, so this is the first check of them. Open View & measure,
press Measure, and read the numbers off the panel.
"""
import re
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

import os

# Overridable: trame serves one session per process, so checking a build without
# evicting whoever is on 8080 means starting a second instance on another port.
URL = os.environ.get("VC_URL") or "http://127.0.0.1:8080/"
SHOT = Path(__file__).parent / "measure.png"
fails = []

with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1500, "height": 1000})
    errs = []
    page.on("pageerror", lambda e: errs.append(str(e)))
    page.goto(URL, wait_until="domcontentloaded", timeout=60000)
    page.wait_for_selector(".v-application", timeout=45000)
    page.wait_for_timeout(1500)
    page.get_by_role("button", name="Load").first.click(timeout=10000)
    page.wait_for_selector(".op-card", timeout=30000)
    page.wait_for_timeout(3000)

    print("=" * 70)
    print("no em dashes anywhere on screen")
    print("=" * 70)
    body = page.inner_text("body")
    # chr(), not the character itself. A repo-wide em-dash strip rewrote this very
    # literal into a hyphen, and the check then counted every ordinary hyphen on
    # the page and reported nine failures that were not there.
    EM = chr(0x2014)
    dashes = body.count(EM)
    print(f"  em dashes in visible text : {dashes}")
    if dashes:
        for line in body.splitlines():
            if EM in line:
                print(f"    {line.strip()[:90]}")
        fails.append(f"{dashes} em dash(es) still on screen")

    print()
    print("=" * 70)
    print("measure")
    print("=" * 70)
    page.get_by_text("VIEW & MEASURE", exact=False).first.click()
    page.wait_for_timeout(1200)
    btn = page.get_by_role("button", name="Measure", exact=True).first
    btn.click(timeout=10000)
    try:
        page.wait_for_selector(".readout", timeout=60000)
    except Exception:
        fails.append("the readout never appeared")
    page.wait_for_timeout(1500)

    rows = page.evaluate("""() => [...document.querySelectorAll('.readout .ro-row')]
        .map(r => [...r.children].map(c => c.textContent.trim()))""")
    for r in rows:
        print("   ", "  ".join(x for x in r if x))
    if not rows:
        fails.append("no readout rows")
    else:
        flat = " ".join(" ".join(r) for r in rows)
        # Case-insensitive: CSS uppercases the label, textContent does not.
        mx = re.search(r"\bmax\s+([\d.]+)", flat, re.I)
        print()
        print(f"  max diameter parsed       : {mx.group(1) if mx else 'NOT FOUND'}")
        if not mx:
            fails.append("max diameter not shown")
        elif not (10.0 < float(mx.group(1)) < 200.0):
            fails.append(f"max diameter implausible: {mx.group(1)}")

    page.locator(".v-navigation-drawer").screenshot(path=str(SHOT))
    print(f"\n  screenshot -> {SHOT}")
    if errs:
        print("  page errors:")
        for e in errs[:5]:
            print(f"    {e[:120]}")
    b.close()

print()
print("ALL PASS" if not fails else f"FAILURES: {fails}")
