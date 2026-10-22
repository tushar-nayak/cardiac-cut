"""What exactly is throwing during a drag?

page.on('pageerror') stringifies to 'Object' for non-Error throws, which hides the
message. Grab name/message/stack off the exception instead, and pull the real pixel
content off the canvas so 'did the 3-D change' is measured, not inferred.
"""
import hashlib

from playwright.sync_api import sync_playwright

import os

# Overridable: trame serves one session per process, so checking a build without
# evicting whoever is on 8080 means starting a second instance on another port.
URL = os.environ.get("VC_URL") or "http://127.0.0.1:8080/"
errs = []

with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1500, "height": 950})

    def on_err(e):
        try:
            errs.append(f"{getattr(e, 'name', '?')}: {getattr(e, 'message', str(e))}")
        except Exception:
            errs.append(repr(e))

    page.on("pageerror", on_err)
    page.on("console", lambda m: errs.append("CONSOLE " + m.text)
            if m.type == "error" else None)

    page.goto(URL, wait_until="domcontentloaded", timeout=60000)
    page.wait_for_selector(".v-application", timeout=45000)
    # capture unhandled rejections too -- trame triggers return promises
    page.evaluate("""() => {
      window.__errs = [];
      window.addEventListener('unhandledrejection', e => {
        window.__errs.push('REJECTION ' + (e.reason && (e.reason.message || e.reason)));
      });
      window.addEventListener('error', e => {
        window.__errs.push('ERROR ' + e.message);
      });
    }""")
    page.wait_for_timeout(1200)
    page.get_by_role("button", name="Load").first.click(timeout=10000)
    page.wait_for_selector(".op-card", timeout=30000)
    page.wait_for_timeout(3500)

    # Hashing the canvas via toDataURL kept returning the same digest even when the
    # scene had changed -- the remote view does not necessarily draw into the canvas
    # we can read. Screenshot the viewport region instead: that is the actual pixels.
    vp = page.evaluate("""() => {
      const el = document.querySelector('.v-main');
      if (!el) return null;
      const r = el.getBoundingClientRect();
      return {x: Math.round(r.x), y: Math.round(r.y),
              width: Math.round(r.width), height: Math.round(r.height)};
    }""")

    def canvas_hash():
        png = page.screenshot(clip=vp) if vp else page.screenshot()
        return hashlib.sha1(png).hexdigest()[:12], len(png)

    card = page.locator(".op-card").first
    track = card.locator(".v-slider-track__background").first
    thumb = card.locator(".v-slider-thumb").first
    tb, thb = track.bounding_box(), thumb.bounding_box()

    print("=" * 70)
    print("canvas content across a drag")
    print("=" * 70)
    h0 = canvas_hash()
    print(f"    before         {h0[0]}  ({h0[1]} chars)")
    page.mouse.move(thb["x"] + thb["width"] / 2, thb["y"] + thb["height"] / 2)
    page.mouse.down()
    page.wait_for_timeout(250)
    during = []
    for frac in (0.3, 0.5, 0.7):
        page.mouse.move(tb["x"] + tb["width"] * frac,
                        thb["y"] + thb["height"] / 2, steps=10)
        page.wait_for_timeout(800)
        h = canvas_hash()
        during.append(h[0])
        print(f"    dragging {frac:.0%}   {h[0]}  ({h[1]} chars)")
    page.mouse.up()
    page.wait_for_timeout(3000)
    hE = canvas_hash()
    print(f"    after release  {hE[0]}  ({hE[1]} chars)")

    page_errs = page.evaluate("() => window.__errs || []")
    b.close()

print()
print("=" * 70)
print("errors")
print("=" * 70)
allerrs = errs + page_errs
if not allerrs:
    print("  none")
for e in dict.fromkeys(allerrs):
    print(f"  {e[:200]}")

print()
changed = len(set(during) - {h0[0]}) > 0
print(f"  canvas changed during drag : {changed}")
print(f"  canvas changed after       : {hE[0] != h0[0]}")
