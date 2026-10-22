"""Open the desktop tool on a case with its nub-in cut already made.

The batch preset gets 85 of 85 nubs cut cleanly, and "cleanly" is not the same as
"right" -- a nub that had to grow is a nub nobody chose, and some vessels want a
human. This is the bridge: `tools/nub_workbench.py` walks the cohort, and this
module puts each case's batch cut into the window so the work starts from what the
preset produced rather than from a default 6 mm on everything.

Nothing here changes how the cutter or the window behave. It sets the same things a
person would set -- add these trims, put these sliders here, cut these ends
horizontally -- and then asks for a re-cut.

Three things about the window that this has to respect, all of which bite quietly:

  _build_rows writes the CURRENT row values back into _offset_memory before it
  rebuilds. So seeding the memory and then rebuilding hands your seeds straight back
  to whatever the sliders already showed. Set the rows, not the memory.

  set_value CLAMPS to the row's own limit, which for a trim is two radii short of
  the junction feeding it. Every cohort nub fits (the deepest is 0.884 of its
  vessel against a 0.95 ceiling), but a clamp would be silent, so it is checked and
  reported rather than assumed.

  Adding a trim rebuilds every row, so all the trims have to go in BEFORE any
  slider is touched, or the values are thrown away by the next rebuild.
"""
import json
from pathlib import Path

import numpy as np

from branch_cut import pre_bifurcation_candidates
from nub_cut import heal_split_segments

# How far apart two loads of the same case may place the same opening and still be the
# same opening. A weld or a split moves a tip by a fraction of a millimetre at most;
# the nearest DIFFERENT opening on these models is centimetres away.
TIP_MATCH_MM = 5.0


def load_record(path):
    return json.loads(Path(path).read_text())


def apply(win, rec, out_dir=None, name=None):
    """Put the cut described by *rec* into the open window. Returns a note list."""
    notes = []
    # Where exports go is set first and unconditionally. It is the one part that
    # still matters when there is no cut to restore, and a pass over the cohort that
    # quietly wrote its results to the previous default would be found out late.
    if out_dir is not None:
        win.out_dir = Path(out_dir)
        win.out_dir.mkdir(parents=True, exist_ok=True)
        if hasattr(win, "out_edit"):
            win.out_edit.setText(str(win.out_dir))
    if name:
        win.name_edit.setText(str(name))

    if win.tree is None:
        return notes + ["no centerline loaded, cannot restore a cut"]

    # The centerline must be the one the batch cut, not the one the window would have
    # chosen. The window PROBES the line as supplied and, if its turning exceeds the
    # jagged threshold, ticks its own smoothing box and reloads at whatever is in its
    # settings. Six of the thirty thoracic cases are jagged, so on those six the window
    # silently cut a different centerline from the one the batch scored: same offsets,
    # different planes. Three of them then refused to export a cut the batch had called
    # clean, and the other three looked fine while being different geometry -- the
    # worse outcome of the two. Set explicitly from the record, never inferred, and
    # never left to a per-user setting.
    reload_needed = False

    want = rec.get("smooth_mm")
    if want is not None and hasattr(win, "smooth_cb"):
        on = float(want) > 0.0
        now = (win.smooth_mm.value() if win.smooth_cb.isChecked() else 0.0)
        if abs(float(now) - float(want)) > 1e-6:
            win.smooth_cb.blockSignals(True)
            win.smooth_cb.setChecked(on)
            win.smooth_cb.blockSignals(False)
            if on:
                win.smooth_mm.blockSignals(True)
                win.smooth_mm.setValue(float(want))
                win.smooth_mm.blockSignals(False)
            reload_needed = True
            notes.append(
                f"centerline set to {'smoothed %.1f mm' % want if on else 'as supplied'}"
                f", matching the batch (the window would have used "
                f"{'smoothed %.1f mm' % now if now else 'as supplied'})")

    # Which anatomy this is. The window infers it from the relative size of the
    # largest openings and reads two of the thirty as abdominal, which does more than
    # misname them: the count decides which terminals are TRUNKS, so an aortic end
    # gets classified as a branch and nubbed instead of squared off. Declared from the
    # record, so the window classifies the case the same way the batch did.
    n_trunk = rec.get("n_trunk_ends")
    if n_trunk and hasattr(win, "family"):
        idx = win.family.findData(int(n_trunk))
        if idx >= 0 and win.family.currentIndex() != idx:
            was = win.family.currentText()
            win.family.blockSignals(True)
            win.family.setCurrentIndex(idx)
            win.family.blockSignals(False)
            reload_needed = True
            notes.append(f"anatomy declared \"{win.family.currentText()}\" as the batch "
                         f"had it (the window had chosen \"{was}\")")

    # One reload for both, not one each: each costs a full re-read and re-cut.
    if reload_needed:
        win.load_model(keep_smooth_choice=True)
    if win.tree is None:
        return notes + ["the model failed to reload"]

    # Applied whether or not there is a cut to restore, because it decides whether the
    # case is cuttable AT ALL. Seven of the thirty arrive with the centerline broken at
    # a non-junction, and the stub left at the end of the model measures its radius
    # against the mesh's end cap; DJ18_2's descending aorta then refuses at every
    # offset in both plane modes, and cuts with a rim 0.000 mm off plane once repaired.
    for nm, was, now, src in heal_split_segments(win.tree, win.base):
        notes.append(f"split-line repair: {nm} radius {was:.2f} -> {now:.2f} mm, "
                     f"taken from {src} which it continues into")

    cuts = rec.get("cuts") or []
    if not cuts:
        win._on_change()
        return notes + ["no cut list to restore, the case is opened as it loads"]

    # 1. trims first: each one rebuilds every row, which would discard slider values
    want_trim = [c for c in cuts if c["interior"]]
    if want_trim:
        cands = pre_bifurcation_candidates(win.tree, win.base)
        for c in want_trim:
            if any(tr["cand"]["segment"] == c["segment"] for tr in win.trims):
                continue
            cand = next((x for x in cands if x["segment"] == c["segment"]), None)
            if cand is None:
                notes.append(f"{c['name']}: {c['segment']} is not offered as a trim "
                             f"here, its nub is left as separate openings")
                continue
            win._add_trim(cand)

    # 1b. the REGION, if this record isolates one. After the trims and before any
    #     slider is touched, for the same reason the trims come first: _set_region
    #     rebuilds every row, so a value set before it is thrown away.
    #
    #     A region record turns a nub cut into an arch-and-nubs cut -- the great
    #     vessels left as stubs, the aorta squared off at the arch's two stations
    #     rather than at its own ends. Restored by STATION, not by preset name: the
    #     batch may have moved a bound to clear an ostium it was undercutting, and
    #     re-deriving the preset here would put it back where it did not cut.
    reg = rec.get("region")
    if reg and reg.get("s0") is not None:
        s0, s1 = reg["s0"], reg.get("s1")
        win._set_region((float(s0), None if s1 is None else float(s1)),
                        name=reg.get("name") or "")
        span = ("open at the far end" if s1 is None
                else f"{float(s1) - float(s0):.0f} mm long")
        notes.append(f"region \"{reg.get('name') or 'custom'}\" restored at "
                     f"{float(s0):.1f}"
                     f"{'' if s1 is None else ' to %.1f' % float(s1)} mm ({span})")

    # 2. plane modes. The preset squares the aorta's ends off horizontally, and drops
    #    back to perpendicular where the vessel is too oblique for that; both are
    #    recorded per terminal, so neither has to be worked out again.
    def _match(t):
        """The cut this row is, by label, then segment, then by where its tip is.

        Segment names are not stable between two loads of the same file: the weld
        step renames a host to <name>_A/<name>_B when it splits it, and did so in the
        window and not the batch on KP35_1. Position is stable, so it is the last and
        most trustworthy word.
        """
        interior = bool(getattr(t, "is_interior", False))
        for c in cuts:
            if c.get("label") == t.label and c["interior"] == interior:
                return c
        for c in cuts:
            if c["segment"] == t.segment.name and c["interior"] == interior:
                return c
        best, best_d = None, TIP_MATCH_MM
        for c in cuts:
            tip = c.get("tip")
            if tip is None or c["interior"] != interior:
                continue
            d = float(np.linalg.norm(np.asarray(tip, float) - np.asarray(t.tip, float)))
            if d < best_d:
                best, best_d = c, d
        return best
    for r in win.rows:
        t = r.terminal
        c = _match(t)
        if c is not None:
            if t.mode != c["mode"]:
                t.mode = c["mode"]
            # How far the cut reaches around its plane -- the Plane size control.
            # Left automatic it is scaled from the segment's measured radius, and
            # where that radius is wrong the band cannot span the vessel and the
            # opening reads as torn at every offset. Whatever the batch had to widen
            # it to has to be set here too, or the window cuts a vessel the batch
            # never cut.
            if c.get("band_mm") is not None:
                t.band_mm = float(c["band_mm"])
                notes.append(f"{t.name}: plane size set to {t.band_mm:.1f} mm, "
                             f"as the batch cut it")
            # Tilt off centerline-perpendicular, where nothing flat would cut. Without
            # it the window re-cuts that vessel flat and gets the torn opening the
            # batch avoided.
            ta = float(c.get("tilt_a_deg") or 0.0)
            tb = float(c.get("tilt_b_deg") or 0.0)
            if abs(ta) > 1e-9 or abs(tb) > 1e-9:
                t.tilt_a_deg, t.tilt_b_deg = ta, tb
                notes.append(f"{t.name}: plane tilted {ta:+.0f}/{tb:+.0f} deg, "
                             f"as the batch cut it")

    # 3. offsets. Matched by SEGMENT and kind, not by label: the batch calls an
    #    interior nub "<seg>.nub" and the window calls the same cut "<seg>.trim".
    placed, missed = 0, []
    for r in win.rows:
        t = r.terminal
        want = _match(t)
        if want is None:
            missed.append(t.name)
            continue
        r.blockSignals(True)
        r.set_value(float(want["offset_mm"]))
        r.blockSignals(False)
        got = float(r.value())
        if abs(got - float(want["offset_mm"])) > 0.05:
            notes.append(f"{t.name}: asked for {want['offset_mm']:.1f} mm, the slider "
                         f"stops at {got:.1f} mm")
        else:
            placed += 1
    if missed:
        notes.append(f"not in the batch cut, left at the default: {', '.join(missed)}")

    win._offset_memory.update({r.terminal.label: r.value() for r in win.rows})
    win._on_change()                      # one re-cut, after everything is in place
    notes.insert(0, f"nub-in cut restored: {placed} opening(s) placed"
                    f"{f', {len(win.trims)} trimmed group(s)' if win.trims else ''}")
    return notes


def apply_from_file(win, path, out_dir=None, name=None):
    return apply(win, load_record(path), out_dir=out_dir, name=name)


def plan_here(win, k=None, anat=None):
    """Work out the nub-in cut for the model ALREADY OPEN in *win*.

    Returns a record in exactly the shape the batch writes, so it goes straight to
    apply() and gets the same trim handling, the same tip-based matching, and the
    same notes. Until this existed the preset was reachable only by running the
    cohort tool and then reopening the window on its output -- the rule was validated
    on 30 cases and then not offered to anyone cutting a single one.

    The cut is performed on a COPY here and only its offsets are kept, rather than
    the rule being re-derived from k. That is not wasted work: a rescue rewrites an
    offset when a nub will not cut cleanly, and an obliquity fallback rewrites a
    trunk end's plane mode, and neither is recoverable from k alone. Re-deriving
    would hand the window a plan the batch would not have produced.
    """
    import nub_cut as NC

    if win.tree is None or win.base is None:
        raise ValueError("load a mesh and a centerline first")
    k = NC.NUB_K if k is None else float(k)
    if anat is None:
        anat = dict(getattr(win, "anat_vessels", None) or {})

    # Seven of the thirty thoracic cases arrive with the centerline broken at a
    # non-junction, and the stub left there measures its radius against the mesh's end
    # cap. Every nub length is a multiple of a measured radius, so this comes first.
    heal = list(heal_split_segments(win.tree, win.base))

    out, _plan, notes, info = NC.cut_nubs(win.base, win.tree, k=k, anat=anat)
    cuts = []
    for t in info["terminals"]:
        if t.label not in info["offsets"]:
            continue
        cuts.append({
            "label": t.label, "segment": t.segment.name, "name": t.name,
            "offset_mm": float(info["offsets"][t.label]), "mode": t.mode,
            "interior": bool(getattr(t, "is_interior", False)),
            "trunk": bool(getattr(t, "is_trunk", False)),
            "band_mm": (None if getattr(t, "band_mm", None) is None
                        else float(t.band_mm)),
            "tip": [float(x) for x in t.tip],
            "tilt_a_deg": float(getattr(t, "tilt_a_deg", 0.0) or 0.0),
            "tilt_b_deg": float(getattr(t, "tilt_b_deg", 0.0) or 0.0),
            "subtree": sorted(getattr(t, "subtree_segments", None) or []),
        })
    # smooth_mm and n_trunk_ends are deliberately absent: apply() reloads the model
    # when either disagrees with the window, and here the window IS the source. Its
    # own settings are already what this plan was computed against.
    return {"case": "this window", "k": float(k), "cuts": cuts,
            "notes": [f"split-line repair: {nm} radius {was:.2f} -> {now:.2f} mm"
                      for nm, was, now, _src in heal] + list(notes),
            "openings": len(out.faces) and int(out.euler_number),
            "validate": info.get("validate")}


def apply_here(win, k=None):
    """Compute the nub-in cut for the open model and put it in the window."""
    rec = plan_here(win, k=k)
    return list(rec.get("notes") or []) + apply(win, rec)
