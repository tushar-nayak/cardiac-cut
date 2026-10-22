# delta K for the web app: run the job, hold the answer, colour the patches.
#
# The counterpart of vascular_api.CutSession -- everything the shape-metrics section
# needs that is not trame. Keeping it here means the section can be driven and checked
# without a browser, which is the only practical way to test a 30-second computation.
#
# The metric itself is in shape_metrics.py. This module is about running it somewhere
# safe and turning the result into something a 3-D view can draw.

import itertools
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import shape_metrics as SM                                  # noqa: E402

WORKER = ROOT / "shape_worker.py"
_RUNS = itertools.count(1)          # numbers every run this process starts

# The reference cohort ran at 100,000. The others are here because a coarser run is
# useful while you are still deciding where to cut, and because seeing delta K move
# with face count is the fastest way to understand that it is not a shape constant.
FACE_TARGETS = (25_000, 50_000, 100_000, 200_000)

# 0.0025 is the reference, and the rest of the ladder is here so that the sensitivity
# can be seen rather than asserted. 0 is offered because it is the honest untouched
# number and because someone will ask for it; it is also the setting at which two runs
# of the SAME surface differ by 2.6x on average. 0.05 was dropped: it discards a tenth
# of all patches, which stops being a trim and starts being a different metric.
TRIM_OPTIONS = (0.0, 0.001, 0.0025, 0.005, 0.01)

# How many independent partitions to average. Error falls as 1/sqrt(n) and runtime
# rises as n, so the useful range is short: 4 halves it, 8 halves it again, and past
# that the compute is better spent elsewhere.
SEED_OPTIONS = (1, 2, 4, 8)


class ShapeError(RuntimeError):
    pass


# ── how the patches are coloured ─────────────────────────────────────────────
#
# Integrated K, and only that. It is what the lab's own CAMP figures show, so the
# picture on screen can be laid beside them. Derived views -- deviation from the
# area-weighted mean, Gaussian curvature, fit quality, patch identity -- were briefly
# offered here and were a distraction: they invite reading a shape off a colouring
# that no published figure uses, and the thing that actually needed fixing was the
# SCALE, not the choice of metric.
COLOUR_TITLE = "Integrated K"
COLOUR_MAP = "plasma"


class ShapeJob:
    """One delta K run.

    The computation happens in a separate process. Not for isolation's sake: ACVD and
    CAMP are VTK underneath, and running them on a thread beside this app's live
    render window segfaults the interpreter within seconds. See shape_worker.py.
    """

    def __init__(self, workdir=None):
        self.dir = Path(workdir or tempfile.mkdtemp(prefix="vcshape_"))
        self.dir.mkdir(parents=True, exist_ok=True)
        self._blank()

    def reset(self):
        """Forget the surface and the run, keeping the working folder.

        For the web app's wipe between visitors: a fresh ShapeJob would make a fresh
        temp folder every visit and leave the old one behind, so the object is emptied
        in place instead. Stops a running worker first; its files stay in the folder
        (nothing here deletes), and nothing reads them again once the fields are blank.

        Raises if the worker will not stop. The caller -- the wipe between visitors --
        then fails closed: nobody new is let in while a run of the last person's is
        still writing.
        """
        if self.running:
            self.proc.kill()
            self.proc.wait(timeout=5)       # TimeoutExpired propagates, on purpose
        self._close_log()
        self._blank()

    def _blank(self):
        self.source = None            # Path to the surface being measured
        self.source_name = ""
        self.source_kind = ""         # "cut" (from this session) or "file"
        self.source_faces = 0
        self.source_euler = None
        self.source_area = None
        self.source_mesh = None
        self.source_defects = None

        self.proc = None
        self.log_path = None
        self.npz_path = None
        self.started_at = 0.0
        self.elapsed = 0.0

        self.summary = None           # scalars, once it finishes
        self.arrays = None            # per-patch and per-face arrays
        self.error = ""
        self._edges = None            # patch outlines, built once on demand

    # ── where the surface comes from ─────────────────────────────────────────
    def _inspect(self, path):
        """Face count, area and topology, before spending 30 seconds on the metric.

        Worth the second it costs: a closed surface or a mesh that is far too coarse
        is knowable up front, and saying so before the run is the difference between
        a warning and a wasted wait.
        """
        import trimesh
        m = trimesh.load(str(path), process=False)
        w = SM.weld(m)
        self.source_faces = int(len(w.faces))
        self.source_euler = int(w.euler_number)
        self.source_area = float(w.area)
        # Kept for the viewport. Whichever surface is about to be measured should be
        # the one on screen while you choose the settings, not the previous one.
        self.source_mesh = w
        # A second, against a run of half a minute, and it is the only chance to say
        # "this surface has geometry the vessel does not have" BEFORE someone reads a
        # number off it.
        try:
            self.source_defects = SM.surface_defects(w)
        except Exception:                                   # noqa: BLE001
            self.source_defects = None

    def adopt_mesh(self, mesh, name, expect_euler=None, kind="cut"):
        """Take the surface straight from a cut in this session.

        Written out as STL and read back, rather than passed in memory, because the
        worker is another process and because this is the exact file a person would
        have exported and re-uploaded. Measuring anything else would make the two
        routes able to disagree.
        """
        if mesh is None:
            raise ShapeError("There is no cut surface to measure yet.")
        stem = "".join(c for c in str(name or "cut")
                       if c.isalnum() or c in "-_.").strip(".") or "cut"
        path = self.dir / f"{stem}.stl"
        mesh.export(str(path))
        self.source = path
        self.source_name = path.name
        # "cut" from the Cut section, "heart" from the Cardiac section. A ventricle
        # wall is not compared with the aortic reference cohort.
        self.source_kind = kind
        self._inspect(path)
        if expect_euler is not None and self.source_euler != int(expect_euler):
            raise ShapeError(f"The exported surface has euler number "
                             f"{self.source_euler}, but the cut reported "
                             f"{int(expect_euler)}. Export changed the topology.")
        self.clear_result()
        return self

    def adopt_path(self, path):
        """Take a surface that already exists on this machine."""
        p = Path(str(path).strip().strip('"'))
        if not p.exists():
            raise ShapeError(f"{p} does not exist.")
        if p.is_dir():
            raise ShapeError(f"{p} is a folder, not a surface.")
        try:
            self._inspect(p)
        except ShapeError:
            raise
        except Exception as e:                              # noqa: BLE001
            raise ShapeError(f"Could not read {p.name}: {type(e).__name__}: {e}")
        self.source = p
        self.source_name = p.name
        self.source_kind = "file"
        self.clear_result()
        return self

    def source_info(self):
        if self.source is None:
            return None
        d = self.source_defects or {}
        return {
            "name": self.source_name,
            "kind": self.source_kind,
            "faces": self.source_faces,
            "euler": self.source_euler,
            "area_mm2": self.source_area,
            "openings": (2 - self.source_euler) if self.source_euler is not None
                        else None,
            "closed": self.source_euler == 2,
            "folds": d.get("folds"),
            "max_dihedral": d.get("max_dihedral"),
            "fold_span_mm": d.get("fold_span_mm"),
            "nonmanifold_edges": d.get("nonmanifold_edges"),
        }

    # ── running ──────────────────────────────────────────────────────────────
    @property
    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def clear_result(self):
        self.summary = None
        self.arrays = None
        self.error = ""
        self._edges = None

    def start(self, target_faces, trim_frac, n_seeds=1):
        if self.source is None:
            raise ShapeError("Choose a surface first.")
        if self.running:
            raise ShapeError("A run is already in progress.")
        self.clear_result()

        n_seeds = max(int(n_seeds), 1)
        seeds = list(range(42, 42 + n_seeds))
        # Unique per run, not only per setting. The folder outlives each web visitor,
        # and with settings-only names a crashed run could pick up the result file
        # or error text an EARLIER run -- another person's -- left under the same
        # name. Nothing here deletes, so the names have to keep runs apart.
        stamp = (f"{int(target_faces)}_{trim_frac:g}_s{n_seeds}_"
                 f"{time.strftime('%Y%m%d%H%M%S')}_{os.getpid()}_{next(_RUNS)}")
        self.npz_path = self.dir / f"result_{stamp}.npz"
        self.log_path = self.dir / f"run_{stamp}.log"
        job = self.dir / f"job_{stamp}.json"
        job.write_text(json.dumps({
            "mesh": str(self.source),
            "out": str(self.npz_path),
            "target_faces": int(target_faces),
            "trim_frac": float(trim_frac),
            # One partition is one draw: re-seeding the same mesh moves delta K by
            # ~16%, because the number is carried by a handful of patches wherever
            # curvature changes fastest and k-means places those boundaries
            # arbitrarily. seeds=[42] alone is the historic behaviour exactly.
            "seeds": seeds,
            # ACVD may move every vertex, but it must not open or close a hole. Handing
            # the worker the topology we already measured turns that into an assertion.
            "expect_euler": self.source_euler,
        }, indent=2), encoding="utf-8")

        # Same interpreter, so the worker gets the same venv without anyone having to
        # arrange it. stderr folded into the log: CAMP's traceback and its progress
        # belong in one place, in order.
        self._log_fh = open(self.log_path, "w", encoding="utf-8", errors="replace")
        self.proc = subprocess.Popen(
            [sys.executable, str(WORKER), str(job)],
            stdout=self._log_fh, stderr=subprocess.STDOUT, cwd=str(ROOT))
        self.started_at = time.monotonic()
        self.elapsed = 0.0

    def cancel(self):
        if self.running:
            self.proc.kill()
            self.proc.wait(timeout=5)
        self._close_log()
        self.proc = None

    def _close_log(self):
        fh = getattr(self, "_log_fh", None)
        if fh is not None and not fh.closed:
            fh.close()

    def _log_text(self):
        if self.log_path is None or not self.log_path.exists():
            return ""
        try:
            return self.log_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return ""

    def phase(self):
        """What the run is doing now.

        CAMP prints once when it starts and again when it finishes, and nothing at all
        during the split loop -- so a log tail alone would sit on "Loaded mesh" for
        most of the run and look stuck. Read the markers as phase boundaries instead,
        and say what happens between them.
        """
        log = self._log_text()
        if "[done]" in log:
            return "Trimming and integrating curvature"
        if "[input]" in log:
            return "Partitioning: splitting patches until every Monge fit converges"
        return "Remeshing to a uniform triangulation"

    def poll(self):
        """Where the run is. Returns one of idle / running / done / failed."""
        if self.proc is None:
            if self.error:
                return {"state": "failed", "message": self.error}
            if self.summary is not None:
                return {"state": "done", "message": "Done."}
            return {"state": "idle", "message": ""}

        self.elapsed = time.monotonic() - self.started_at
        code = self.proc.poll()
        if code is None:
            return {"state": "running", "message": self.phase(),
                    "elapsed": self.elapsed}

        self._close_log()
        self.proc = None
        side = Path(str(self.npz_path) + ".json")
        payload = {}
        if side.exists():
            try:
                payload = json.loads(side.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                payload = {}

        if code != 0 or "error" in payload or not self.npz_path.exists():
            self.error = payload.get("error") or self._failure_text(code)
            return {"state": "failed", "message": self.error}

        with np.load(self.npz_path) as z:
            self.arrays = {k: z[k] for k in z.files}
        self.summary = payload
        self.summary["seconds"] = round(self.elapsed, 1)
        return {"state": "done", "message": "Done.",
                "elapsed": self.elapsed}

    def _failure_text(self, code):
        """A worker that dies without writing its side file leaves only the log.

        The last line is nearly always the exception, and showing it beats "exit 1",
        which tells a person nothing they can act on.
        """
        lines = [ln.strip() for ln in self._log_text().splitlines() if ln.strip()]
        tail = lines[-1] if lines else ""
        if tail:
            return f"The run failed (exit {code}). {tail}"
        return (f"The run failed (exit {code}) without writing anything. It may have "
                f"run out of memory; try a smaller target face count.")

    # ── the answer ───────────────────────────────────────────────────────────
    @property
    def done(self):
        return self.summary is not None and self.arrays is not None

    def mean_K(self):
        """The area-weighted mean that delta K is the variance about."""
        a = self.arrays
        keep = a["keep"]
        q, w = a["integrated_K"][keep], a["patch_area"][keep]
        if w.sum() <= 0:
            return 0.0
        w = w / w.sum()
        return float(np.sum(w * q))

    # What the colour bar spans, and why the default is not the full range.
    #
    # Integrated K is avg_kg x patch_area, and its distribution is long-tailed: on
    # KY160_1 it runs -11.94 to 3.367 while the great majority of patches sit within a
    # few tenths of zero. Handed the full range -- which is what CAMP does, since
    # campsrc passes cmin=None, cmax=None and lets plotly auto-scale -- every one of
    # those patches lands in the top few per cent of the map and the whole vessel comes
    # out one flat colour. The picture is faithful and it shows nothing.
    #
    # So the default clips a few per cent off each end and lets the extremes saturate.
    # CAMP has the same escape hatch: campsrc.camp takes color_bounds and clips to it,
    # so setting limits is its own mechanism rather than a departure from it.
    #
    # The clip is ADJUSTABLE because no constant carries. The range that first showed
    # structure across the aorta body on KY160_1 was about -0.2 to 0.2, saturating
    # 6.2%; a median +/- 3 MAD rule reproduces that almost exactly there and saturates
    # 21% of JK5_1, whose range is six times wider and whose maximum is 270. Measured
    # across percentile rules on both, the middle half of the surface never uses more
    # than a quarter of the map. Too variable to decide for someone.
    SCALES = ("robust", "full", "manual")
    DEFAULT_CLIP_PCT = 5.0

    def _clim(self, vals, scale, lo, hi, clip_pct=None):
        full = (float(np.min(vals)), float(np.max(vals)))
        wide = full if full[1] > full[0] else (full[0], full[0] + 1.0)
        if scale == "manual" and lo is not None and hi is not None:
            a, b = float(lo), float(hi)
            if b > a:
                return (a, b), full
        if scale == "full":
            return wide, full
        pct = float(np.clip(
            self.DEFAULT_CLIP_PCT if clip_pct is None else clip_pct, 0.0, 45.0))
        a, b = (float(x) for x in np.percentile(vals, [pct, 100.0 - pct]))
        return ((a, b) if b > a else wide), full

    def face_values(self, scale="robust", lo=None, hi=None, clip_pct=None):
        """Per-face integrated K, its colour map, and the range the bar should span.

        *scale* is "robust" (clip *clip_pct* off each end), "full" (CAMP's own
        auto-scale) or "manual" with explicit *lo* and *hi*.
        """
        if not self.done:
            return None
        a = self.arrays
        per = a["integrated_K"]
        vals = per[a["face_labels"]]
        clim, full = self._clim(vals, scale, lo, hi, clip_pct=clip_pct)
        return {"values": vals, "clim": clim, "cmap": COLOUR_MAP, "n_colors": 256,
                "title": COLOUR_TITLE, "bar": True, "range": full,
                "clipped": (clim[0] > full[0] + 1e-12) or (clim[1] < full[1] - 1e-12)}

    def kept_face_mask(self):
        """True for faces whose patch survived the trim."""
        if not self.done:
            return None
        return self.arrays["keep"][self.arrays["face_labels"]]

    def patch_edges(self):
        """Vertex pairs along patch boundaries, plus the mesh rim.

        Drawing these is what makes the partition legible. Colour alone does not: two
        adjacent patches with similar integrated K are the same shade, so a surface
        can look like one smooth field when it is in fact a hundred separate fits.

        An interior edge is shared by two faces and sorts next to its twin, so the
        boundary edges are the pairs whose two owners carry different labels.
        """
        if self._edges is not None or not self.done:
            return self._edges
        f = np.asarray(self.arrays["faces"], dtype=np.int64)
        lab = np.asarray(self.arrays["face_labels"])

        e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]], axis=0)
        owner = np.tile(np.arange(len(f), dtype=np.int64), 3)
        key = np.sort(e, axis=1)
        order = np.lexsort((key[:, 1], key[:, 0]))
        key, owner, e = key[order], owner[order], e[order]

        paired = np.flatnonzero(np.all(key[1:] == key[:-1], axis=1))
        between = e[paired][lab[owner[paired]] != lab[owner[paired + 1]]]

        # Edges with no twin are the open rim, which is worth outlining too: it is
        # where the cut was, and it is the boundary CAMP's patches had to respect.
        lone = np.ones(len(key), dtype=bool)
        lone[paired] = False
        lone[paired + 1] = False

        self._edges = np.vstack([between, e[lone]]) if lone.any() else between
        return self._edges

    def polydata_arrays(self):
        """The partitioned surface, or the source while there is nothing to colour."""
        if self.done:
            return self.arrays["vertices"], self.arrays["faces"]
        if self.source_mesh is not None:
            return (np.asarray(self.source_mesh.vertices, dtype=float),
                    np.asarray(self.source_mesh.faces, dtype=np.int64))
        return None, None

    # ── how to read the number ───────────────────────────────────────────────
    def notes(self):
        """Everything that would make this delta K incomparable to the cohort.

        delta K is extensive and partition dependent. The value is always computable
        and is only meaningful next to another value produced the same way, so the
        conditions that differ are the most important thing on the panel.
        """
        out = []
        s = self.summary or {}
        if not s:
            return out

        # First, and always, because it is the one thing that changes how every other
        # number on the panel should be read. A delta K is not a measurement of the
        # anatomy alone -- it is a measurement of the anatomy AND of a partition that
        # gets rebuilt from scratch each run.
        n_seeds = int(s.get("n_seeds", 1) or 1)
        band = SM.partition_band(n_seeds)
        if n_seeds > 1 and s.get("k_fluct_se"):
            out.append(("info",
                        f"Averaged over {n_seeds} independent partitions: delta K = "
                        f"{s['k_fluct']:.5g} +/- {s['k_fluct_se']:.2g} "
                        f"({100 * s['k_fluct_se'] / s['k_fluct']:.0f}%). The patches "
                        f"drawn are the FIRST partition ({s['k_fluct_seed']:.5g}); the "
                        f"number is the mean of all {n_seeds}."))
        else:
            out.append(("warn",
                        f"ONE partition. delta K is carried by ~22 patches out of "
                        f"~850, sited where curvature changes fastest, and where "
                        f"k-means puts those boundaries is arbitrary -- re-running "
                        f"with a different seed moves the answer by about "
                        f"{100 * SM.PARTITION_CV:.0f}%. Raise the partition count to "
                        f"average that down: it falls as 1/sqrt(n), so 4 gives "
                        f"+/-{100 * SM.partition_band(4):.0f}% and 8 gives "
                        f"+/-{100 * SM.partition_band(8):.0f}%."))
        out.append(("info",
                    f"Comparable only at a fixed mesh density. delta K drifts with "
                    f"target_faces (about {100 * SM.DENSITY_CV:.0f}% across a 4x "
                    f"change) because min_patch_size is a vertex count, so the "
                    f"smallest patch covers a different physical area at each "
                    f"density. Averaging partitions does not remove that -- it is a "
                    f"change of scale, not noise."))

        f = self.source_defects or {}
        if f.get("folds"):
            span = f.get("fold_span_mm")
            where = (f", clustered within {span:.0f} mm of each other"
                     if span is not None and span < 25.0 else "")
            out.append(("warn",
                        f"{f['folds']} folded edge(s) in the source surface, worst "
                        f"dihedral {f['max_dihedral']:.0f}deg{where}. That is the "
                        f"surface doubling back through itself -- geometry the vessel "
                        f"does not have, inherited from the segmentation rather than "
                        f"introduced here. Curvature fitted across a fold is not "
                        f"anatomy. It is NOT a correction to apply: across the "
                        f"reference cohort fold count does not predict how much delta "
                        f"K moves. Fix it upstream if the value matters."))
        if f.get("nonmanifold_edges"):
            out.append(("warn",
                        f"{f['nonmanifold_edges']} non-manifold edge(s): more than two "
                        f"faces meet along them, so the surface has no consistent "
                        f"normal there."))

        if self.source_kind == "heart":
            # The reference cohort is 30 thoracic AORTAS. Its area range, floor
            # fraction, face count and trim say nothing about a ventricle wall, so none
            # of those comparisons is made; what is left applies to any surface.
            out.append(("info",
                        "A ventricle wall, from the Cardiac section. There is no "
                        "reference cohort for ventricles, so compare this value only "
                        "with other ventricle walls run at the same face count and "
                        "trim. The partition-noise figures above were measured on "
                        "aortas; ventricles have not been measured yet."))
            return self._generic_notes(out, s)
        return self._reference_notes(out, s)

    def _reference_notes(self, out, s):
        lo, hi = SM.REFERENCE["area_mm2"]
        if not s.get("in_reference_range", True):
            out.append(("warn",
                        f"This surface is {s['A_mm2']:,.0f} mm2, outside the "
                        f"reference cohort's {lo:,.0f} to {hi:,.0f} mm2. delta K is "
                        f"still computed, but there is no reference curve under it "
                        f"to place the value against."))

        f_lo, f_hi = SM.REFERENCE["pct_at_min_size"]
        pct = s.get("pct_at_min_size")
        if pct is not None and not (f_lo <= pct <= f_hi):
            side = "coarser" if pct > f_hi else "finer"
            out.append(("warn",
                        f"{pct:.1f}% of patches sit at the vertex floor, outside the "
                        f"reference's {f_lo} to {f_hi}%. The partition is {side} than "
                        f"the cohort's, so this value is not on their scale. Change "
                        f"the target face count and run it again."))

        if int(s.get("target_faces", 0)) != SM.REFERENCE["target_faces"]:
            out.append(("warn",
                        f"Run at {int(s['target_faces']):,} faces. The reference "
                        f"cohort ran at {SM.REFERENCE['target_faces']:,}, and delta K "
                        f"moves with face count."))

        tf = float(s.get("trim_frac", 0))
        if abs(tf - SM.REFERENCE["trim_frac"]) > 1e-12:
            out.append(("warn",
                        f"Trim fraction {tf:g}, not the reference's "
                        f"{SM.REFERENCE['trim_frac']:g}. The trim is a quantile on "
                        f"integrated K and it takes BOTH tails, so {tf:g} discards "
                        f"{200 * tf:.2g}% of patches, not {100 * tf:.2g}%."))

        return self._generic_notes(out, s)

    def _generic_notes(self, out, s):
        """The notes that hold for any surface, aorta or not."""
        chi = s.get("euler_number")
        if chi == 2:
            out.append(("warn",
                        "This surface is closed. CAMP was built for open geometries, "
                        "and a closed one has no rim, so the partition near where the "
                        "openings would be is not comparable. Cut it first."))

        pf = s.get("pct_failed_quality")
        if pf is not None and pf > 20.0:
            out.append(("warn",
                        f"{pf:.1f}% of patches never met the fit thresholds "
                        f"(R2 >= {SM.PARAMS['R2_thresh']}, eta <= "
                        f"{SM.PARAMS['eta_thresh']}). They are still counted, and they "
                        f"are deliberately never dropped: a patch that will not take a "
                        f"quadratic is where the surface is least like a quadratic, "
                        f"which is the thing delta K exists to measure. Excluding them "
                        f"does not clean the number up, it removes most of it. A high "
                        f"figure usually means the surface is noisy or has slivers; "
                        f"try the sliver-free rim when cutting."))
        return out

    def repeatability(self):
        """(median, worst) ratio, and the trim it was measured at, for this run."""
        return SM.repeatability_band((self.summary or {}).get("trim_frac"))

    def headline(self):
        """The conditions the number is only valid under, in one line."""
        s = self.summary or {}
        if not s:
            return ""
        return (f"{int(s.get('faces', 0)):,} faces, {s.get('n_patches', 0):,} patches, "
                f"trim {s.get('trim_frac', 0):g}, CAMP {s.get('campsrc', '?')}")

    def csv(self):
        """The per-patch table, as the download in the Streamlit app produces it."""
        if not self.done:
            return ""
        a = self.arrays
        return SM.patch_table({
            "integrated_K": a["integrated_K"], "patch_area": a["patch_area"],
            "avg_kg": a["avg_kg"], "eta_patch": a["eta_patch"],
            "R2_patch": a["R2_patch"], "keep": a["keep"],
        }).to_csv(index=False)

    def export(self, out_dir, stem):
        """The patch table, and the numbers that say how it was produced.

        The sidecar is not optional decoration. A delta K with no record of face
        count, trim and CAMP build is a number nobody can reproduce or compare.
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        table = out_dir / f"{stem}_patches.csv"
        side = out_dir / f"{stem}_delta_k.json"
        # newline="" or every row ends \r\r\n: pandas already terminates with \r\n and
        # the default text mode translates the \n again. Excel copes; anything that
        # splits on \r sees a blank row between every patch.
        table.write_text(self.csv(), encoding="utf-8", newline="")
        side.write_text(json.dumps({
            **{k: v for k, v in (self.summary or {}).items()},
            "params": SM.PARAMS,
            "reference": SM.REFERENCE,
            "source_kind": self.source_kind,
        }, indent=2, default=str), encoding="utf-8")
        return table, side
