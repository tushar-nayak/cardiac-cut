"""The aortic arch, isolated, with its great vessels left as nubs. Thoracic cohort.

    python tools/arch_cohort.py run                 # cut every case, score, write
    python tools/arch_cohort.py run --k 0.75        # a different nub length
    python tools/arch_cohort.py run --only KY160_1,KK20_2
    python tools/arch_cohort.py report              # read the last run

WHAT IT MAKES, and why it is one geometry and not two

Each case comes out as ONE mesh: the arch between its two region stations, with the
brachiocephalic, carotid and subclavian cut back to stubs standing proud of the
aortic wall. That is the nub-in preset and the "Aortic arch" region applied together,
which is exactly what a curvature comparison of the arch wants -- the aorta trimmed to
the arch rather than to its own ends, and the branches present but not dominating the
surface.

HOW IT CUTS, and why it goes through CutSession

Every cut is made by the same session object the desktop window and the website use,
rather than by a batch path of its own. That matters more here than anywhere else in
the repository: a nub on a vessel that SPLITS -- a brachiocephalic feeding two great
vessels -- is a TRIM, and trims inside a region were silently dropped by the API until
they were fixed to match the window (tests/test_region_trim_parity.py). A batch that
cut its own way would have produced arch meshes with no nubs on them and reported
nothing wrong.

So the flow is: plan the nubs against the session's OWN tree, hand the plan back to
the session as trims and offsets, then isolate the arch and re-cut. What the batch
scores is what the window will open.

WHAT IS SCORED

    openings   how many the cut left, against how many were asked for
    flat       every rim within RIM_FLATNESS_MM of the plane that made it
    validate   what the tool itself says -- the bar an export has to clear
    proud      how far each nub stands past the aortic wall, in mm

A case is FLAGGED when any nub had to be rescued or clamped, when a bound had to be
moved to cut, or when validate refuses. Flagged does not mean wrong: it means the run
made a choice a person should look at, which is what the workbench is for.

Only the thoracic cohort. Meshes and centerlines are read; nothing is written into it.
"""
import argparse
import json
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np                                                 # noqa: E402

import nub_cut as NC                                               # noqa: E402
from anat_bench import find_cases                                  # noqa: E402
from branch_cut import RIM_FLATNESS_MM, boundary_loops             # noqa: E402
from vascular_api import CutSession, CutError                      # noqa: E402

COHORT = r"Z:/Martino/branches_CAMP/GT"
DEFAULT_OUT = ROOT / "output" / "arch_nubs"
DEFAULT_RUN = HERE / "arch_run.jsonl"
ARCH_PRESET = "Aortic arch"


def _arch_span(s):
    """The arch preset's two stations, or None if this case does not offer it."""
    for p in s.presets():
        if p["name"].strip().lower() == ARCH_PRESET.lower():
            return p["name"], p.get("s0"), p.get("s1")
    return None, None, None


def _apply_nubs(s, plan, used, offsets, terms):
    """Put a nub plan into the session as trims and offsets.

    Matched by SEGMENT, not by label. nub_terminals builds its own Terminal objects
    and a trim's label is minted by whoever made it -- the session's is
    "<segment>.trim" -- so the labels do not line up even though both were built from
    the same tree. Segment names do.
    """
    notes = []
    for t in terms:
        seg = t.segment.name
        off = float(offsets.get(t.label, 0.0) or 0.0)
        if getattr(t, "is_interior", False):
            try:
                label = s.add_trim(seg)
            except CutError as e:
                notes.append(f"{t.name}: {seg} is not offered as a trim here "
                             f"({str(e)[:50]}), its nub stays as separate openings")
                continue
            s.set_offset(label, off)
        else:
            hit = next((x for x in s.tree.terminals if x.segment.name == seg
                        and bool(x.at_end) == bool(t.at_end)), None)
            if hit is None:
                notes.append(f"{t.name}: no opening on {seg} in the session")
                continue
            hit.mode = t.mode
            if getattr(t, "band_mm", None):
                hit.band_mm = float(t.band_mm)
            s.set_offset(hit.label, off)
    return notes


def _proud(plan):
    """Per-vessel nub length past the aortic wall, from the plan."""
    out = {}
    for v in plan:
        nm = v.get("name") or v.get("segment")
        p = v.get("proud_mm", v.get("nub_mm"))
        if p is not None:
            out[str(nm)] = round(float(p), 2)
    return out


def run_case(name, csv, mat, k, outdir, write_mesh=True, export_dir=None):
    t0 = time.time()
    s = CutSession(str(csv), str(mat), family="thoracic")
    region_name, s0, s1 = _arch_span(s)
    if region_name is None:
        raise ValueError("this case offers no aortic arch region")

    # Planned against the SESSION's own tree, so every segment name and terminal the
    # plan refers to is one the session already has.
    _out, plan, nub_notes, info = NC.cut_nubs(s.base, s.tree, k=k,
                                              anat=getattr(s, "anat_vessels", None))
    notes = list(nub_notes)
    notes += _apply_nubs(s, plan, info.get("used"), info["offsets"], info["terminals"])

    s.apply_preset(region_name)
    st = s.recut()

    loops = boundary_loops(s.mesh) if s.mesh is not None else []
    active, offs = s._cut_list()
    flat = None
    if st["ok"] and s.mesh is not None:
        worst = 0.0
        verts = np.asarray(s.mesh.vertices)
        for t in active:
            pt, nrm = t.cut_plane(float(offs[t.label]))
            if not loops:
                break
            rim = min(loops, key=lambda L: float(
                np.linalg.norm(verts[np.asarray(L)].mean(axis=0) - pt)))
            worst = max(worst, float(np.abs(
                (verts[np.asarray(rim)] - pt) @ nrm).max()))
        flat = round(worst, 3)

    cuts = []
    for t in active:
        cuts.append({
            "label": t.label, "segment": t.segment.name, "name": str(t.name),
            "offset_mm": float(offs[t.label]), "mode": getattr(t, "mode", None),
            "interior": bool(getattr(t, "is_interior", False)),
            "bound": getattr(t, "_region_offset", None) is not None,
            "trunk": bool(getattr(t, "is_trunk", False)),
            "band_mm": (None if getattr(t, "band_mm", None) is None
                        else float(t.band_mm)),
            "tip": [float(x) for x in t.tip],
            "tilt_a_deg": float(getattr(t, "tilt_a_deg", 0.0) or 0.0),
            "tilt_b_deg": float(getattr(t, "tilt_b_deg", 0.0) or 0.0),
            "subtree": sorted(getattr(t, "subtree_segments", None) or []),
        })

    moved = []
    if getattr(s, "_trunk_moved", None):
        moved.append(s._trunk_moved)
    for b in (getattr(s, "_bound_moved", None) or []):
        moved.append(b)

    rec = {"case": name, "k": float(k), "secs": round(time.time() - t0, 1),
           "csv": str(csv), "mat": str(mat),
           "smooth_mm": float(getattr(s, "smooth_mm", 0.0)),
           "n_trunk_ends": 2,
           "region": {"name": region_name, "s0": None if s0 is None else float(s0),
                      "s1": None if s1 is None else float(s1)},
           "ok": bool(st["ok"]), "message": str(st.get("message", "")),
           "openings": len(loops), "expected": len(active),
           "euler": (None if s.mesh is None else int(s.mesh.euler_number)),
           "area_mm2": (None if s.mesh is None else round(float(s.mesh.area), 1)),
           "faces": (None if s.mesh is None else int(len(s.mesh.faces))),
           "flat_mm": flat, "moved": moved, "notes": notes,
           "proud_mm": _proud(plan), "cuts": cuts}

    if write_mesh and st["ok"]:
        outdir.mkdir(parents=True, exist_ok=True)
        stem = f"{name.split('/')[-1]}_arch_k{k:g}"
        s.mesh.export(str(outdir / f"{stem}.stl"))
        (outdir / f"{stem}.json").write_text(json.dumps(rec, indent=2))
        rec["mesh"] = f"{stem}.stl"

    # The same export the workbench window writes, in the same place and under the
    # same name, so a case re-cut in batch REPLACES the one a person exported by hand
    # rather than sitting beside it in a second format. Everything downstream --
    # arch_viewer, arch_bovine_plot, the CAMP run -- reads <case>_arch_cut_planes.json
    # and would otherwise keep reading the superseded cut.
    if export_dir is not None and st["ok"]:
        stl, side = s.export(export_dir, f"{name.split('/')[-1]}_arch")
        rec["export"] = str(stl)
        rec["export_sidecar"] = str(side)
    return rec


def flagged(rec):
    """Did the run make a choice a person should look at?"""
    if not rec.get("ok"):
        return True
    if rec.get("moved"):
        return True
    if rec.get("flat_mm") is not None and rec["flat_mm"] > RIM_FLATNESS_MM:
        return True
    return any("rescue" in n or "grown" in n or "clamped" in n
               for n in (rec.get("notes") or []))


def cmd_run(args):
    cases = find_cases(args.cohort)
    only = {x.strip() for x in (args.only or "").split(",") if x.strip()}
    if only:
        cases = [c for c in cases if c[0] in only]
    out = Path(args.run)
    outdir = Path(args.out)
    print(f"{len(cases)} cases, k={args.k}, meshes -> {outdir}\n", flush=True)
    recs = []
    with out.open("w", encoding="utf-8") as fh:
        for i, (cid, csv, mat) in enumerate(cases, 1):
            try:
                rec = run_case(cid, csv, mat, args.k, outdir,
                               write_mesh=not args.no_mesh,
                               export_dir=args.export)
            except Exception as e:                              # noqa: BLE001
                rec = {"case": cid, "csv": str(csv), "mat": str(mat),
                       "ok": False, "error": f"{type(e).__name__}: {e}",
                       "trace": traceback.format_exc()[-400:]}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            recs.append(rec)
            mark = "ok " if rec.get("ok") else "REFUSED"
            extra = rec.get("error") or rec.get("message", "")
            print(f"  [{i}/{len(cases)}] {cid:10s} {mark} "
                  f"{rec.get('openings', '-')}/{rec.get('expected', '-')} openings "
                  f"{str(extra)[:46]}", flush=True)
    _report(recs)
    return 0


def _report(recs):
    ok = [r for r in recs if r.get("ok")]
    flag = [r for r in ok if flagged(r)]
    print(f"\n{'=' * 70}\n  {len(ok)} of {len(recs)} cases cut the arch with nubs\n"
          f"{'=' * 70}")
    print(f"  clean and unflagged : {len(ok) - len(flag)}")
    print(f"  cut but flagged     : {len(flag)}")
    bad = [r for r in recs if not r.get("ok")]
    if bad:
        print(f"  refused / raised    : {len(bad)}")
        for r in bad:
            print(f"      {r['case']:10s} {str(r.get('error') or r.get('message'))[:60]}")
    if flag:
        print("\n  flagged, and why:")
        for r in flag:
            why = []
            if r.get("moved"):
                why.append(f"{len(r['moved'])} plane(s) moved by the tool")
            if r.get("flat_mm") is not None and r["flat_mm"] > RIM_FLATNESS_MM:
                why.append(f"worst rim {r['flat_mm']:.2f} mm")
            grown = [n for n in (r.get("notes") or [])
                     if "rescue" in n or "grown" in n or "clamped" in n]
            if grown:
                why.append(f"{len(grown)} nub(s) adjusted")
            print(f"      {r['case']:10s} {'; '.join(why)}")
    areas = [r["area_mm2"] for r in ok if r.get("area_mm2")]
    if areas:
        print(f"\n  surface area {min(areas):,.0f} to {max(areas):,.0f} mm2 "
              f"(median {sorted(areas)[len(areas)//2]:,.0f})")
    return 0


def cmd_report(args):
    recs = [json.loads(l) for l in Path(args.run).read_text().splitlines() if l.strip()]
    return _report(recs)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("run", "report"):
        p = sub.add_parser(nm)
        p.add_argument("--run", default=str(DEFAULT_RUN))
        if nm == "run":
            p.add_argument("--cohort", default=COHORT)
            p.add_argument("--out", default=str(DEFAULT_OUT))
            p.add_argument("--k", type=float, default=NC.NUB_K)
            p.add_argument("--only")
            p.add_argument("--no-mesh", action="store_true")
            p.add_argument("--export", default=None, help=(
                "also write <case>_arch.stl and <case>_arch_cut_planes.json here, "
                "the workbench's own export, replacing what is already there"))
    args = ap.parse_args()
    return cmd_run(args) if args.cmd == "run" else cmd_report(args)


if __name__ == "__main__":
    raise SystemExit(main())
