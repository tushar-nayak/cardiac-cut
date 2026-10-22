"""How many region presets does the TOOL actually accept, end to end?

    python tools/region_bench.py                       # abdominal
    python tools/region_bench.py Z:/path/to/cohort     # any cohort

This is the region number that means something to a user. anat_bench --cut counts
whether the two BOUNDARY PLANES leave two openings, which is a useful diagnostic and
a flattering proxy: it says 177 of 182 on the abdominal cohort where the tool accepts
146. The difference is every branch cut inside the region, plus validate_cut over the
result -- neither of which the proxy touches.

Grouped by failure so the next thing to fix is visible rather than buried.
"""
import sys
from pathlib import Path
ROOT = Path(r"c:\Users\cmdavis2\Projects\branches")
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/"tools"))
from vascular_api import CutSession
from anat_bench import find_cases

n = ok = 0
bad = []
for cid, cl, mp in find_cases(sys.argv[1] if len(sys.argv) > 1
                              else r"Y:/aorta/abdominal/MK/JK"):
    try:
        s = CutSession(str(cl), str(mp))
        names = [p["name"] for p in s.presets()]
    except Exception as e:
        print(f"  {cid:10s} SETUP {type(e).__name__}: {str(e)[:50]}")
        continue
    for nm in names:
        n += 1
        try:
            s.apply_preset(nm)
            st = s.recut()
        except Exception as e:
            bad.append((cid, nm, f"{type(e).__name__}: {str(e)[:44]}")); continue
        finally:
            s.clear_region()
        if st["ok"]:
            ok += 1
        else:
            bad.append((cid, nm, st["message"][:52]))
print(f"\n  {ok} of {n} region presets cut end to end ({100*ok/max(n,1):.0f}%)")
print(f"  {len(bad)} refused:")
for c, nm, why in bad:
    print(f"    {c:10s} {nm[:26]:26s} {why}")
