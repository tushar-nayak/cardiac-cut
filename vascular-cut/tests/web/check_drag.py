"""Does the slider thumb follow the mouse, or only catch up on release?

Drag the first opening's slider in steps, sampling the thumb position and the mm
box after each move but BEFORE releasing. If the thumb only moves on mouseup, the
local binding is not reactive and the number box will not track either.
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

import os

# Overridable: trame serves one session per process, so checking a build without
# evicting whoever is on 8080 means starting a second instance on another port.
URL = os.environ.get("VC_URL") or "http://127.0.0.1:8080/"
SHOT = Path(__file__).parent / "drag.png"

with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1500, "height": 950})
    page.goto(URL, wait_until="domcontentloaded", timeout=60000)
    page.wait_for_selector(".v-application", timeout=45000)
    page.wait_for_timeout(2000)
    page.get_by_role("button", name="Load").first.click(timeout=10000)
    page.wait_for_selector(".op-card", timeout=30000)
    page.wait_for_timeout(3000)

    card = page.locator(".op-card").first
    slider = card.locator(".v-slider-track__background").first
    thumb = card.locator(".v-slider-thumb").first
    box = card.locator("input[type=number]").first

    def sample(tag):
        tx = page.evaluate(
            "() => { const t = document.querySelector('.op-card .v-slider-thumb');"
            " return t ? Math.round(t.getBoundingClientRect().left) : null; }")
        val = box.input_value()
        print(f"    {tag:<22} thumb x={tx}   box={val!r}")
        return tx, val

    print("=" * 70)
    print("drag the first opening's depth slider")
    print("=" * 70)
    tb = slider.bounding_box()
    thb = thumb.bounding_box()
    print(f"  track: x={tb['x']:.0f}..{tb['x']+tb['width']:.0f}  thumb at {thb['x']:.0f}")

    start = sample("before")
    page.mouse.move(thb["x"] + thb["width"] / 2, thb["y"] + thb["height"] / 2)
    page.mouse.down()
    page.wait_for_timeout(200)

    during = []
    for frac in (0.25, 0.45, 0.65):
        page.mouse.move(tb["x"] + tb["width"] * frac, thb["y"] + thb["height"] / 2,
                        steps=8)
        page.wait_for_timeout(450)
        during.append(sample(f"dragging {frac:.0%}"))

    page.mouse.up()
    page.wait_for_timeout(2500)
    after = sample("after release")

    page.screenshot(path=str(SHOT))
    print(f"\n  screenshot -> {SHOT}")

    moved_during = any(d[0] is not None and start[0] is not None
                       and abs(d[0] - start[0]) > 4 for d in during)
    box_during = any(d[1] != start[1] for d in during)
    moved_after = (after[0] is not None and start[0] is not None
                   and abs(after[0] - start[0]) > 4)

    print()
    print(f"  thumb moved DURING the drag  : {moved_during}")
    print(f"  box  updated DURING the drag : {box_during}")
    print(f"  thumb moved after release    : {moved_after}")

    b.close()

print()
if moved_during and box_during:
    print("ALL PASS")
else:
    print("FAILURES: the control does not track the drag "
          f"(thumb={moved_during}, box={box_during}, after_release={moved_after})")
