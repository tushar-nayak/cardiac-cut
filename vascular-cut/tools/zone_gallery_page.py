"""Turn the rendered strips into one page to scroll and check.

    .venv\\Scripts\\python tools\\zone_gallery_page.py                  # thoracic
    .venv\\Scripts\\python tools\\zone_gallery_page.py --tag JK

Reads tools/zone_gallery_<tag>/ and writes a single self-contained HTML file with
every image inlined, because the page is published as an artifact and a strict CSP
blocks any request off-host. Nothing is fetched at view time.

The page is for ONE job: looking at 30 cases quickly and spotting the ones that are
wrong. So it leads with what needs attention, keeps every case one scroll-stop tall,
and puts the numbers that decide "is this right" -- area, loops, the span asked for --
directly under the panel they belong to.
"""
import argparse
import base64
import html
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Warm renders on cool ground: the vessels are the only warm thing on the page, so
# they read as the subject and the chrome stays out of the way. The accent is not a
# taste choice -- it is the exact orange the zones are painted in.
CSS = """
:root{
  --ground:#eceff3; --surface:#ffffff; --sunk:#e3e8ee;
  --ink:#141b23; --muted:#5d6b7a; --line:#d5dce4;
  --accent:#c2410c; --ok:#2f7d55; --warn:#b45309; --bad:#b42318;
  --shadow:0 1px 2px rgba(20,27,35,.06),0 8px 24px rgba(20,27,35,.05);
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --ground:#11151a; --surface:#19202a; --sunk:#141b23;
    --ink:#e7ecf1; --muted:#93a2b3; --line:#2b3642;
    --accent:#e2673a; --ok:#4fa87a; --warn:#d08a2c; --bad:#e0685c;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.3);
  }
}
:root[data-theme="dark"]{
  --ground:#11151a; --surface:#19202a; --sunk:#141b23;
  --ink:#e7ecf1; --muted:#93a2b3; --line:#2b3642;
  --accent:#e2673a; --ok:#4fa87a; --warn:#d08a2c; --bad:#e0685c;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.3);
}
*{box-sizing:border-box}
body{
  margin:0; background:var(--ground); color:var(--ink);
  font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  -webkit-font-smoothing:antialiased;
}
.serif{font-family:"Iowan Old Style",Georgia,"Times New Roman",serif}
.mono{font-family:ui-monospace,SFMono-Regular,"SF Mono",Consolas,monospace;
      font-variant-numeric:tabular-nums}
.wrap{max-width:1240px;margin:0 auto;padding:0 20px}

header{
  position:sticky; top:0; z-index:10; background:var(--ground);
  border-bottom:1px solid var(--line);
}
.head-top{display:flex;flex-wrap:wrap;gap:18px;align-items:baseline;
          padding:18px 0 10px}
h1{font-size:25px;margin:0;font-weight:600;letter-spacing:-.01em;
   text-wrap:balance}
h1 .sub{color:var(--muted);font-weight:400}
.tally{display:flex;gap:16px;margin-left:auto;flex-wrap:wrap}
.tally b{font-size:19px;font-weight:600;display:block;line-height:1.15}
.tally span{font-size:11px;text-transform:uppercase;letter-spacing:.07em;
            color:var(--muted)}
.controls{display:flex;gap:8px;align-items:center;padding-bottom:10px;
          flex-wrap:wrap}
button{
  font:inherit;font-size:13px;color:var(--ink);background:var(--surface);
  border:1px solid var(--line);border-radius:6px;padding:5px 11px;cursor:pointer;
}
button:hover{border-color:var(--muted)}
button[aria-pressed="true"]{background:var(--ink);color:var(--ground);
                            border-color:var(--ink)}
button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.jump{display:flex;gap:4px;overflow-x:auto;padding:0 0 12px;
      scrollbar-width:thin}
.jump a{
  font-size:11.5px;color:var(--muted);text-decoration:none;white-space:nowrap;
  padding:3px 7px;border-radius:5px;border:1px solid transparent;
}
.jump a:hover{color:var(--ink);border-color:var(--line);background:var(--surface)}
.jump a.flag{color:var(--accent)}

main{padding:22px 0 80px;display:flex;flex-direction:column;gap:18px}
.case{
  background:var(--surface);border:1px solid var(--line);border-radius:10px;
  box-shadow:var(--shadow);overflow:hidden;scroll-margin-top:150px;
}
.case-head{display:flex;gap:14px;align-items:baseline;flex-wrap:wrap;
           padding:13px 16px 11px;border-bottom:1px solid var(--line)}
.case-id{font-size:17px;font-weight:600;letter-spacing:-.01em}
.route{font-size:12.5px;color:var(--muted)}
.chips{display:flex;gap:6px;margin-left:auto;flex-wrap:wrap}
.chip{font-size:11px;letter-spacing:.04em;text-transform:uppercase;
      padding:2px 8px;border-radius:20px;border:1px solid currentColor}
.chip.ok{color:var(--ok)} .chip.warn{color:var(--warn)} .chip.bad{color:var(--bad)}
.strip{display:block;width:100%;height:auto;cursor:zoom-in;
       background:var(--sunk)}
.case.zoom .strip{cursor:zoom-out}
.labels{display:grid;gap:0;border-top:1px solid var(--line)}
.lab{padding:9px 8px 11px;border-left:1px solid var(--line);text-align:center}
.lab:first-child{border-left:0}
.lab .n{font-size:12.5px;font-weight:600;display:block;margin-bottom:3px}
.lab .n.cut{color:var(--muted);font-weight:500}
.lab .d{font-size:11px;color:var(--muted);display:block;line-height:1.45}
.lab .isl{color:var(--warn)}
.note{padding:10px 16px;font-size:12.5px;color:var(--bad);
      border-top:1px solid var(--line);background:var(--sunk)}
.legend{font-size:12.5px;color:var(--muted);padding:2px 0 16px;max-width:62ch}
.legend b{color:var(--ink);font-weight:600}
.sw{display:inline-block;width:9px;height:9px;border-radius:2px;
    background:var(--accent);vertical-align:baseline}
footer{color:var(--muted);font-size:12px;padding:0 0 40px}
@media (max-width:720px){
  .lab .d{display:none}
  .lab .n{font-size:10.5px}
}
"""

JS = """
const cases=[...document.querySelectorAll('.case')];
const btns=[...document.querySelectorAll('[data-filter]')];
function apply(mode){
  btns.forEach(b=>b.setAttribute('aria-pressed',b.dataset.filter===mode));
  cases.forEach(c=>{
    c.hidden = (mode==='flag' && c.dataset.flag!=='1');
  });
}
btns.forEach(b=>b.addEventListener('click',()=>apply(b.dataset.filter)));
apply('all');
cases.forEach(c=>{
  const img=c.querySelector('.strip');
  img.addEventListener('click',()=>c.classList.toggle('zoom'));
});
"""


def b64(p):
    return base64.b64encode(p.read_bytes()).decode("ascii")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="GT")
    ap.add_argument("--label", default="thoracic")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    src = ROOT / "tools" / f"zone_gallery_{args.tag}"
    index = json.loads((src / "index.json").read_text(encoding="utf-8"))
    out = Path(args.out or (ROOT / "tools" / f"zone_gallery_{args.tag}.html"))

    n_case = len(index)
    n_zone = sum(max(len(r.get("panels", [])) - 1, 0) for r in index)
    n_flag = 0
    rows = []
    jump = []

    for r in index:
        cid = r["case"]
        if "error" in r:
            n_flag += 1
            jump.append(f'<a class="flag" href="#{cid}">{cid}</a>')
            rows.append(
                f'<section class="case" id="{cid}" data-flag="1">'
                f'<div class="case-head"><span class="case-id">{html.escape(cid)}</span>'
                f'<span class="chips"><span class="chip bad">failed</span></span></div>'
                f'<div class="note">{html.escape(r["error"])}</div></section>')
            continue

        panels = r["panels"]
        zones = panels[1:]
        isl = [p for p in zones
               if p.get("islands") and p["faces"] and
               p["islands"] > 0.02 * p["faces"]]
        flag = (not r["base_ok"]) or bool(isl)
        n_flag += int(flag)
        jump.append(f'<a class="{"flag" if flag else ""}" href="#{cid}">{cid}</a>')

        chips = [f'<span class="chip ok">{len(zones)} zones</span>']
        if isl:
            chips.append(f'<span class="chip warn">{len(isl)} with islands</span>')
        if not r["base_ok"]:
            chips.append('<span class="chip bad">branch cut refused</span>')

        labs = []
        for p in panels:
            nm = p["name"]
            if p.get("area_mm2") is None:
                labs.append('<div class="lab"><span class="n cut">as cut</span>'
                            '<span class="d">every branch cut,<br>whole vessel</span>'
                            '</div>')
                continue
            span = ("open distal" if p["s1"] is None
                    else f'{p["s1"] - p["s0"]:.0f} mm span')
            extra = ""
            if p.get("islands") and p["faces"] and p["islands"] > 0.02 * p["faces"]:
                extra = (f'<br><span class="isl">{p["islands"]} faces dropped '
                         f'({100 * p["islands"] / p["faces"]:.0f}%)</span>')
            labs.append(
                f'<div class="lab"><span class="n">'
                f'{html.escape(nm.replace("Ishimaru ", ""))}</span>'
                f'<span class="d mono">{p["area_mm2"]:,.0f} mm² · '
                f'{p["loops"]} loops<br>{span}{extra}</span></div>')

        note = ("" if r["base_ok"] else
                f'<div class="note">Branch cut refused: '
                f'{html.escape(r["base_msg"])}</div>')
        rows.append(
            f'<section class="case" id="{cid}" data-flag="{int(flag)}">'
            f'<div class="case-head"><span class="case-id">{html.escape(cid)}</span>'
            f'<span class="route mono">route {r["route_mm"]:.0f} mm</span>'
            f'<span class="chips">{"".join(chips)}</span></div>'
            f'<img class="strip" alt="{html.escape(cid)}: whole vessel then each zone"'
            f' src="data:image/png;base64,{b64(src / r["png"])}">'
            f'<div class="labels" style="grid-template-columns:repeat('
            f'{len(panels)},1fr)">{"".join(labs)}</div>{note}</section>')

    doc = f"""<title>Zone Isolation Review</title>
<style>{CSS}</style>
<header><div class="wrap">
  <div class="head-top">
    <h1 class="serif">Zone isolation <span class="sub">— {args.label} cohort</span></h1>
    <div class="tally">
      <div><b class="mono">{n_case}</b><span>cases</span></div>
      <div><b class="mono">{n_zone}</b><span>zones cut</span></div>
      <div><b class="mono">{n_flag}</b><span>need a look</span></div>
    </div>
  </div>
  <div class="controls">
    <button data-filter="all" aria-pressed="true">All cases</button>
    <button data-filter="flag" aria-pressed="false">Needs a look</button>
    <span class="legend" style="padding:0">Click any strip to enlarge it.</span>
  </div>
  <nav class="jump">{"".join(jump)}</nav>
</div></header>
<main class="wrap">
  <p class="legend">Each strip is one case under one camera: the whole vessel as the
  tool cuts it, then every zone it offers, <b><span class="sw"></span> the zone</b>
  against the rest of the aorta ghosted. Boundaries follow the surface at the
  station the anatomy gives — they are <b>not planar</b>, so these are surfaces to
  measure, not inlets to mesh for flow.</p>
  {"".join(rows)}
</main>
<footer class="wrap">Rendered by <span class="mono">tools/zone_gallery.py</span>
from <span class="mono">region_station.region_by_station</span>. Loops counts every
boundary of the isolated surface: two ends plus any branch kept inside the
zone.</footer>
<script>{JS}</script>
"""
    out.write_text(doc, encoding="utf-8")
    mb = out.stat().st_size / 1e6
    print(f"wrote {out}  ({mb:.1f} MB)")
    print(f"  {n_case} cases, {n_zone} zones, {n_flag} flagged")


if __name__ == "__main__":
    main()
