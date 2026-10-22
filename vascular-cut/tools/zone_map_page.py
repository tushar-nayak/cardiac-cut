"""Turn the zone maps into one page to scroll.

    .venv\\Scripts\\python tools\\zone_map_page.py --tag GT

Self-contained: every image is inlined, because the page is published as an artifact
and a strict CSP blocks any request off-host.

The page leads with the two properties that decide whether the zone set is right --
nothing claimed twice, nothing left over -- because those are checkable at a glance
across 30 cases in a way that per-zone areas are not.
"""
import argparse
import base64
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

CSS = """
.pairs{margin:28px 0 8px}
.pairs h2{font-size:1.15rem;margin:0 0 6px}
.pair{margin:0 0 18px;background:var(--surface);border:1px solid var(--line);
  border-radius:10px;padding:10px;box-shadow:var(--shadow)}
.pair img{width:100%;height:auto;display:block;border-radius:6px}
.pair figcaption{color:var(--muted);font-size:.82rem;margin-top:8px}
:root{
  --ground:#eceff3; --surface:#ffffff; --sunk:#e3e8ee;
  --ink:#141b23; --muted:#5d6b7a; --line:#d5dce4;
  --accent:#c2410c; --ok:#2f7d55; --warn:#b45309; --bad:#e11d48;
  --shadow:0 1px 2px rgba(20,27,35,.06),0 8px 24px rgba(20,27,35,.05);
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --ground:#11151a; --surface:#19202a; --sunk:#141b23;
    --ink:#e7ecf1; --muted:#93a2b3; --line:#2b3642;
    --accent:#e2673a; --ok:#4fa87a; --warn:#d08a2c; --bad:#f2506e;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.3);
  }
}
:root[data-theme="dark"]{
  --ground:#11151a; --surface:#19202a; --sunk:#141b23;
  --ink:#e7ecf1; --muted:#93a2b3; --line:#2b3642;
  --accent:#e2673a; --ok:#4fa87a; --warn:#d08a2c; --bad:#f2506e;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.3);
}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);
     font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
     -webkit-font-smoothing:antialiased}
.serif{font-family:"Iowan Old Style",Georgia,"Times New Roman",serif}
.mono{font-family:ui-monospace,SFMono-Regular,"SF Mono",Consolas,monospace;
      font-variant-numeric:tabular-nums}
.wrap{max-width:1120px;margin:0 auto;padding:0 20px}
header{position:sticky;top:0;z-index:10;background:var(--ground);
       border-bottom:1px solid var(--line)}
.head-top{display:flex;flex-wrap:wrap;gap:18px;align-items:baseline;
          padding:18px 0 8px}
h1{font-size:25px;margin:0;font-weight:600;letter-spacing:-.01em;
   text-wrap:balance}
h1 .sub{color:var(--muted);font-weight:400}
.tally{display:flex;gap:16px;margin-left:auto;flex-wrap:wrap}
.tally b{font-size:19px;font-weight:600;display:block;line-height:1.15}
.tally span{font-size:11px;text-transform:uppercase;letter-spacing:.07em;
            color:var(--muted)}
.key{display:flex;gap:14px;flex-wrap:wrap;padding:4px 0 12px;font-size:12.5px;
     color:var(--muted);align-items:center}
.k{display:inline-flex;align-items:center;gap:6px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block;
    box-shadow:0 0 0 1px rgba(0,0,0,.12) inset}
.jump{display:flex;gap:4px;overflow-x:auto;padding:0 0 12px;scrollbar-width:thin}
.jump a{font-size:11.5px;color:var(--muted);text-decoration:none;white-space:nowrap;
        padding:3px 7px;border-radius:5px;border:1px solid transparent}
.jump a:hover{color:var(--ink);border-color:var(--line);background:var(--surface)}
.jump a.flag{color:var(--bad)}
main{padding:20px 0 80px;display:grid;gap:16px;
     grid-template-columns:repeat(auto-fill,minmax(330px,1fr))}
.case{background:var(--surface);border:1px solid var(--line);border-radius:10px;
      box-shadow:var(--shadow);overflow:hidden;scroll-margin-top:170px}
.case-head{display:flex;gap:10px;align-items:baseline;flex-wrap:wrap;
           padding:11px 13px 9px;border-bottom:1px solid var(--line)}
.case-id{font-size:15.5px;font-weight:600}
.route{font-size:11.5px;color:var(--muted)}
.chip{font-size:10.5px;letter-spacing:.04em;text-transform:uppercase;
      padding:2px 7px;border-radius:20px;border:1px solid currentColor;
      margin-left:auto}
.chip.ok{color:var(--ok)} .chip.warn{color:var(--warn)} .chip.bad{color:var(--bad)}
img{display:block;width:100%;height:auto;background:var(--sunk)}
.bar{display:flex;height:9px;border-top:1px solid var(--line)}
.bar i{display:block;height:100%}
.zlist{padding:9px 13px 11px;font-size:11.5px;color:var(--muted);
       display:grid;gap:3px}
.zrow{display:flex;gap:8px;align-items:center}
.zrow .sw{width:9px;height:9px;flex:none}
.zrow .nm{color:var(--ink)}
.zrow .st{margin-left:auto}
.note{padding:9px 13px;font-size:11.5px;color:var(--bad);
      border-top:1px solid var(--line);background:var(--sunk)}
.legend{grid-column:1/-1;font-size:12.5px;color:var(--muted);max-width:64ch}
.legend b{color:var(--ink);font-weight:600}
footer{color:var(--muted);font-size:12px;padding:0 0 40px}
"""


def fmt_end(z):
    return "end" if z["s1"] is None else f'{z["s1"]:.0f}'


def b64(p):
    return base64.b64encode(p.read_bytes()).decode("ascii")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="GT")
    ap.add_argument("--label", default="thoracic")
    args = ap.parse_args()

    src = ROOT / "tools" / f"zone_map_{args.tag}"
    data = json.loads((src / "index.json").read_text(encoding="utf-8"))
    cols, cases = data["colours"], data["cases"]
    out = ROOT / "tools" / f"zone_map_{args.tag}.html"

    n_over = sum(1 for c in cases if c.get("overlap_faces"))
    n_gap = sum(1 for c in cases if c.get("gap_faces"))
    n_err = sum(1 for c in cases if "error" in c)
    n_cases = len(cases)
    n_forced = sum(c.get("stub_forced", 0) for c in cases)
    n_before = 25          # stubs carrying two zones, tools/zone_stub_bench.py
    n_wallchg = 0          # of them aortic by the radius test

    short = {"Ishimaru zone 0": "zone 0", "Ishimaru zone 1": "zone 1",
             "Ishimaru zone 2": "zone 2", "Ishimaru zone 3": "zone 3",
             "Ishimaru zone 4": "zone 4"}
    cards, jump = [], []
    for c in cases:
        cid = c["case"]
        if "error" in c:
            jump.append(f'<a class="flag" href="#{cid}">{cid}</a>')
            cards.append(f'<section class="case" id="{cid}"><div class="case-head">'
                         f'<span class="case-id">{html.escape(cid)}</span>'
                         f'<span class="chip bad">failed</span></div>'
                         f'<div class="note">{html.escape(c["error"])}</div></section>')
            continue
        bad = c["overlap_faces"] > 0
        warn = c["gap_faces"] > 0.01 * c["faces"] or not c["base_ok"]
        jump.append(f'<a class="{"flag" if (bad or warn) else ""}" '
                    f'href="#{cid}">{cid}</a>')
        # A ZONE IS MISSING BECAUSE TWO VESSELS SHARE AN ORIGIN, not because anything
        # failed. On a bovine arch the left carotid comes off the brachiocephalic, so
        # zone 1 has zero length and is never offered; 12 of these 30 are like that.
        # Naming it here matters because the zones are numbered, not counted: a case
        # showing 0, 2, 3, 4 is right, and drawing its zone 2 in zone 1's colour --
        # which is what colouring by position did -- is not.
        have = [z.get("zone") for z in c["zones"]]
        missing = [n for n in range(5) if n not in have and n < max(have or [0])]
        if missing:
            chip = ('<span class="chip warn">shared origin: no zone '
                    + ", ".join(str(n) for n in missing) + '</span>')
        elif bad:
            chip = f'<span class="chip bad">{c["overlap_faces"]} overlap</span>'
        elif warn:
            chip = (f'<span class="chip warn">'
                    f'{100 * c["gap_faces"] / c["faces"]:.0f}% unclaimed</span>')
        else:
            chip = f'<span class="chip ok">{len(c["zones"])} zones tile</span>'

        tot = max(sum(z["verts"] for z in c["zones"]), 1)
        # By the zone's own number: a bovine arch has no zone 1 and colouring by
        # position would paint zone 2 in zone 1's colour.
        def _zn(z, i):
            return int(z.get("zone", i))
        bar = "".join(
            f'<i style="background:{cols[_zn(z, i) % len(cols)]};'
            f'width:{100 * z["verts"] / tot:.2f}%"></i>'
            for i, z in enumerate(c["zones"]))
        zl = "".join(
            f'<div class="zrow"><span class="sw" '
            f'style="background:{cols[_zn(z, i) % len(cols)]}"></span>'
            f'<span class="nm">{html.escape(short.get(z["name"], z["name"]))}</span>'
            f'<span class="st mono">{z["s0"]:.0f}–{fmt_end(z)} mm</span></div>'
            for i, z in enumerate(c["zones"]))
        note = ("" if c["base_ok"] else
                f'<div class="note">Branch cut refused: '
                f'{html.escape(c["base_msg"])}</div>')
        cards.append(
            f'<section class="case" id="{cid}"><div class="case-head">'
            f'<span class="case-id">{html.escape(cid)}</span>'
            f'<span class="route mono">{c["route_mm"]:.0f} mm</span>{chip}</div>'
            f'<img alt="{html.escape(cid)}, both sides, coloured by zone" '
            f'src="data:image/png;base64,{b64(src / c["png"])}">'
            f'<div class="bar">{bar}</div><div class="zlist">{zl}</div>{note}'
            f'</section>')

    pairs_dir = src / "pairs"
    pairs = sorted(pairs_dir.glob("*.png")) if pairs_dir.is_dir() else []
    pairs_html = ""
    if pairs:
        # Close-ups, framed on the left subclavian, because that is where the question
        # is: does zone 3 begin immediately after it.
        strip = "".join(
            f'<figure class="pair"><img alt="{html.escape(q.stem)}, branches before '
            f'and after the stub rule" '
            f'src="data:image/png;base64,{b64(q)}">'
            f'<figcaption><b class="mono">{html.escape(q.stem)}</b> '
            f'— distal lip, then proximal lip + 8 mm, then the clearing plane'
            f'</figcaption></figure>'
            for q in pairs)
        pairs_html = (f'<div class="pairs">'
                      f'<h2 class="serif">Zone 3 begins where the plane first clears '
                      f'the subclavian’s mouth</h2>'
                      f'<p class="legend">Every attempt before the plane moved a '
                      f'<b>level set</b> of the harmonic coordinate. A level set '
                      f'follows the wall, so on a bend it runs oblique — no amount of '
                      f'shifting it makes it the transverse cut the clinical figure '
                      f'draws, which is why the boundary looked wrong wherever it was '
                      f'placed. Three separate placements were tried and rejected on '
                      f'that same misconception.</p>'
                      f'<p class="legend">It is now a plane, and the question was '
                      f'always WHICH plane. Two anchors were tried and both were a lip '
                      f'of the mouth plus a constant: the distal lip, and the proximal '
                      f'lip plus 8 mm. Each was fitted to the five hand-marked cases '
                      f'and each then failed on eight cases judged by eye — in both '
                      f'directions at once, overshooting the branch where the mouth is '
                      f'short and cutting through it where the mouth is long. The '
                      f'mouth’s own span runs from 1 to 20 mm of centerline '
                      f'across this cohort, so no constant can track it. The 8 mm was '
                      f'not a discovery; it was the span of the compact mouths it was '
                      f'fitted to.</p>'
                      f'<p class="legend">The quantity that does track the judgements '
                      f'is the plane’s <b>signed clearance</b> over the ostium: '
                      f'+0.05 mm on the case called correct, −5.2 and −7.3 on the two '
                      f'showing zone-3 faces on a branch, −17.6 on the one called '
                      f'“encroaching on the LSA”, +9.2 on “way too far”, '
                      f'+13.7 and +18.6 on the two called short. One number, monotone, '
                      f'the verdicts not overlapping. So the rule asks for clearance '
                      f'zero directly: sweep the centerline and keep the FIRST station '
                      f'whose normal plane leaves every rim vertex behind it. No '
                      f'percentile, no offset, no fitted constant — and no harmonic '
                      f'coordinate, which is why KK5_1 and KY183_1 now run the rule at '
                      f'all instead of silently keeping the old level set when the '
                      f'capped solve returned nothing.</p>'
                      f'<p class="legend">Measured over the cohort: clearance 0.00 to '
                      f'2.79 mm and <b>never negative</b>, so the plane no longer cuts '
                      f'any ostium; faces with one corner in a mouth and one corner in '
                      f'zone 3 — what renders as zone 3 sitting on a branch — fall from '
                      f'<b>677 to 29</b>, on one case instead of sixteen; zone 3 is a '
                      f'single connected piece on every case; zone islands <b>0</b>.</p>'
                      f'<p class="legend"><b>What it costs, stated plainly.</b> Against '
                      f'the hand marks it is worse: pooled mean 5.2 mm against 3.0, '
                      f'median 2.8 against 1.5. The three compact-mouthed marked cases '
                      f'do not move (KK32_0 2.9→2.8, KK24_1 0.7→0.8, DJ8_1 3.0→3.1); '
                      f'the two elongated ones do (DJ18_2 2.5→6.6, KK16_1 4.9→9.1). '
                      f'That is structural, not tuning: on those two the marks sit '
                      f'INSIDE the mouth’s span — the marked line cuts the '
                      f'subclavian — so what was drawn is not a plane. It hugs the '
                      f'ostium’s distal margin on one flank, and any plane doing '
                      f'that must slice the other. Satisfying both would need a '
                      f'boundary that follows the mouth and closes transversely '
                      f'elsewhere, which is a change of shape, not of anchor, and has '
                      f'not been built.</p>'
                      f'<p class="legend">Three anchors per case below, same framing: '
                      f'distal lip, then proximal lip + 8 mm, then the clearing plane. '
                      f'Red spheres are the marks where there are any; the black '
                      f'outline is the subclavian’s mouth.</p>'
                      f'{strip}</div>')

    keys = "".join(
        f'<span class="k"><span class="sw" style="background:{cols[i]}"></span>'
        f'zone {i}</span>' for i in range(5))
    doc = f"""<title>Aortic Zone Map</title>
<style>{CSS}</style>
<header><div class="wrap">
  <div class="head-top">
    <h1 class="serif">Zone map <span class="sub">— {args.label} cohort</span></h1>
    <div class="tally">
      <div><b class="mono">{len(cases)}</b><span>cases</span></div>
      <div><b class="mono">{n_over}</b><span>with overlap</span></div>
      <div><b class="mono">{n_gap}</b><span>with unclaimed wall</span></div>
    </div>
  </div>
  <div class="key">{keys}
    <span class="k"><span class="sw" style="background:{data['unclaimed']}"></span>
      unclaimed</span>
    <span class="k"><span class="sw" style="background:{data['overlap']}"></span>
      claimed twice</span>
  </div>
  <nav class="jump">{"".join(jump)}</nav>
</div></header>
<main class="wrap">
  <p class="legend"><b>Each branch now carries one zone.</b> A great vessel that
  curls back over the arch used to come out one colour for most of its length and
  then carry a streak of the next zone down the flank that hovers over it. The stub
  is now found from the crease at its own mouth — the join of two convex tubes is
  creased inward all the way round, wherever the vessel afterwards leans — and every
  vertex behind that crease takes the zone its ostium already chose.
  Across {n_cases} cases that repaints {n_forced:,} vertices and takes stubs carrying
  two zones from {n_before} to <b>0</b> of 89. The wall is untouched: a radius test the
  stub finder never sees — read off the vessel's own cross sections, so a branch
  cannot inflate it — calls {n_wallchg} of those vertices aortic, and boundary
  flatness is 2.44 mm median over 107 boundaries against 2.46 before, with overlap and
  unclaimed wall still nil. Taking a stub back can strand a few vertices of
  the zone it was borrowed from, so the ownership pass repairs connectivity with the
  stubs frozen: <b>no zone sits as an island inside another on any case</b>, and the
  zone-piece count is 137 before and after. The three cases marked <b>two stubs
  merged</b> are bovine arches where neighbouring vessels sit close enough that one
  region swallows two.</p>
  <p class="legend"><b>Zones are numbered, not counted.</b> On a bovine arch the left
  carotid arises from the brachiocephalic, so those two share an ostium, zone 1 has zero
  length and is never offered, and the case runs 0, 2, 3, 4. Twelve of these thirty are
  like that and one more skips zone 2. Each zone is drawn in its own number's colour and
  each branch is placed by its own ostium read off the centerline — so a bovine
  trunk's daughters both come out zone 0, where their shared origin is.</p>
  <p class="legend">The stub is grown a step at a time, and a step is kept only if what
  it ADDS is vessel rather than wall — judged by the median distance of its new
  vertices from the aortic centerline, in local aortic radii. Without that, a vessel
  with a wide opening could leak a splash of arch under the rim-length bound, which is
  what put an arrow-shaped tongue on KR47_1's ascending aorta. No single rim bound
  fixed it: loose enough for one vessel on that case leaked, tight enough to stop the
  leak reduced another to nothing but its rim.</p>
  <p class="legend">Two more things the eye caught that no check had. A zone boundary's
  collar reaches around the wall, and at a flared aortic root that is close to a
  boundary 30 mm up the vessel — on KK40_3 it handed 182 root vertices to zone 1 and
  then to zone 2, a band of the wrong zone around the inlet. The vessel's own two ends
  are now off limits to a boundary a zone away. And the aortic radius profile believed
  any ring it found: where a normal plane grazes the wall the ring runs ALONG the vessel
  and reads several calibres wide, 103 mm at KK4_2's root against a true 20. Rings are
  now rejected on shape, which is what starved that case's stubs.</p>
  <p class="legend">One change was made, measured, and <b>reverted on sight</b>.
  Placing each boundary on its landmark's own cross section takes the median
  displacement from the vessel it names to 0.00 mm and would have fixed KK40_3, whose
  zone 1 / zone 2 boundary sits 15 mm past the left carotid. It also grew zone 3 into a
  lobe over the arch apex and squeezed the bands either side of it, on every case. The
  metric improved and the map got worse — landing a ring on its landmark's plane says
  nothing about whether the band it bounds is the right shape. KK40_3's boundary stays
  open rather than trading every case's appearance for it.</p>
  {pairs_html}
  <p class="legend">Each case from both sides, every Ishimaru zone in its own colour.
  The zones <b>tile</b> the aorta — zone 1 ends where zone 2 begins — so a piece of
  wall claimed by two zones is a defect and is painted crimson. Membership is
  half-open, <span class="mono">s0 ≤ station &lt; s1</span>: a branch belongs to the
  zone it is the landmark <b>for</b>, and to that one only. The bar under each case is
  how the wall divides between them.</p>
  {"".join(cards)}
</main>
<footer class="wrap">Rendered by <span class="mono">tools/zone_map.py</span>.
Branches cut to nubs first, so a zone carries its vessels' origins rather than however
much of them the scan happened to trace.</footer>
"""
    out.write_text(doc, encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size / 1e6:.1f} MB)")
    print(f"  {len(cases)} cases, {n_over} with overlap, {n_gap} with unclaimed, "
          f"{n_err} failed")


if __name__ == "__main__":
    main()
