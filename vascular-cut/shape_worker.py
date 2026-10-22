# Run one delta K job in a process of its own, and write the answer to disk.
#
#   python shape_worker.py <job.json>
#
# This exists because of a crash, not because of tidiness. The web app holds a live
# VTK render window; ACVD and CAMP are VTK and numpy underneath. Running the metric
# on a worker THREAD beside that render window segfaults the interpreter -- measured,
# repeatably, exit 139 within a few seconds of the first concurrent render. A thread
# cannot be made safe here, so the work goes somewhere with its own address space.
#
# It also buys two things worth having on their own: a CAMP failure becomes a non-zero
# exit code instead of a dead server, and progress is readable, because CAMP's own
# log lines land in a file the caller can tail while the job runs.
#
# The job file:
#   {"mesh": "<path>", "out": "<path.npz>", "target_faces": 100000,
#    "trim_frac": 0.01, "expect_euler": -3, "seeds": [42, 43, 44, 45]}
#
# "seeds" is optional and defaults to [42], which is the historic single-partition
# behaviour byte for byte. With more than one, CAMP runs once per seed and the reported
# delta K is their mean, with a standard error -- one partition is a draw, not the
# answer. The ARRAYS written are the first seed's, because patch labels cannot be
# averaged across partitions that carve the surface differently.
#
# On success: <out> holds the arrays, <out>.json holds the scalars, exit 0.
# On failure: <out>.json holds {"error": "..."}, exit 1.

import json
import sys
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def run(job):
    import trimesh
    import shape_metrics as SM

    mesh = trimesh.load(str(job["mesh"]), process=False)
    seeds = [int(x) for x in job.get("seeds") or [42]]

    def _progress(i, n):
        # Straight to the log the caller is already tailing, so a multi-seed run does
        # not look like a hung single-seed one.
        print(f"[seeds] partition {i + 1} of {n} (seed {seeds[i]})", flush=True)

    res = SM.delta_k_averaged(
        mesh,
        seeds=seeds,
        progress=_progress if len(seeds) > 1 else None,
        target_faces=int(job["target_faces"]),
        trim_frac=float(job["trim_frac"]),
        # The plotly figure is CAMP's own view of the partition. This app draws the
        # patches itself in the 3-D scene, so building it would cost seconds and be
        # thrown away.
        make_figure=False,
        expect_euler=job.get("expect_euler"),
    )

    # The partitioned mesh, NOT the remeshed one. CAMP repairs connectivity before it
    # partitions, so face_labels is indexed against the mesh it returns. Colouring the
    # pre-CAMP faces with these labels would look right and be wrong.
    cmesh = res["mesh"]
    verts = np.asarray(cmesh.vertices, dtype=np.float32)
    faces = np.asarray(cmesh.faces, dtype=np.int32)
    labels = np.asarray(res["face_labels"], dtype=np.int32).ravel()
    if len(labels) != len(faces):
        raise ValueError(f"CAMP returned {len(labels)} face labels for "
                         f"{len(faces)} faces; they must be one per face")

    out = Path(job["out"])
    np.savez_compressed(
        out,
        vertices=verts, faces=faces, face_labels=labels,
        integrated_K=np.asarray(res["integrated_K"], dtype=np.float64),
        patch_area=np.asarray(res["patch_area"], dtype=np.float64),
        avg_kg=np.asarray(res["avg_kg"], dtype=np.float64),
        eta_patch=np.asarray(res["eta_patch"], dtype=np.float64),
        R2_patch=np.asarray(res["R2_patch"], dtype=np.float64),
        keep=np.asarray(res["keep"], dtype=bool),
    )

    # Everything the caller needs to render a readout without opening the npz.
    summary = {k: v for k, v in res.items()
               if isinstance(v, (int, float, bool, str, type(None)))}
    summary["seeds"] = res.get("seeds", [42])
    summary["k_fluct_all"] = res.get("k_fluct_all", [res["k_fluct"]])
    summary["campsrc"] = SM.campsrc_origin()
    summary["source"] = Path(job["mesh"]).name
    return summary


def main(argv):
    if len(argv) != 2:
        print(__doc__ or "usage: shape_worker.py <job.json>", file=sys.stderr)
        return 2
    job = json.loads(Path(argv[1]).read_text(encoding="utf-8"))
    side = Path(str(job["out"]) + ".json")
    try:
        summary = run(job)
    except Exception as e:                                  # noqa: BLE001
        # The traceback goes to the log the caller is already tailing; the message
        # goes somewhere the caller can read without parsing it back out of the log.
        traceback.print_exc()
        side.write_text(json.dumps({"error": f"{type(e).__name__}: {e}"}, indent=2),
                        encoding="utf-8")
        return 1
    side.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
