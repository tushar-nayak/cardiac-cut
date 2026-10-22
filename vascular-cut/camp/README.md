# CAMP, vendored

`campsrc.py` here is the Pocivavsek Lab's CAMP implementation, the convergent
adaptive patch partition behind the δK (`k_fluct`) shape metric. It was copied in so
the web app can compute δK on a machine that has only this repository.

**It is now also the CANONICAL copy** — decided by CAMP's author on 31 August 2026.
It began as a vendored snapshot, but the four partition changes recorded in
`CHANGES.md` were made here and this is the lineage that carries forward. Where this
disagrees with another copy, including `Y:\aorta\python_code\camp`, this one wins.
Mirrored to `https://github.com/SurgBioMech/CAMP` and to a sister folder on the lab
share.

## Provenance

| | |
|---|---|
| copied from | `c:\Users\cmdavis2\Projects\camp\UpdatedCamp\campsrc.py` |
| source modified | 22 July 2026, 15:43 |
| lines | 812, identical to upstream |
| git blob SHA-1 | `893d0d2d028ee9cd9274cd7f4b0a9ffc75376bd2` |
| edited here? | **yes, since 19 Aug 2026 — see [CHANGES.md](CHANGES.md)** |

Keeping it unmodified *was* deliberate: a diff against upstream stays meaningful, and
nobody has to wonder whether a local edit changed a number. That is no longer true, so
the next best thing applies — every change is enumerated in `CHANGES.md` with its
measured effect, and **each one defaults to the previous behaviour or is verified
byte-identical to it**.

## What is different from upstream

| change | default | reproduces upstream when |
|---|---|---|
| connectivity split | on | `split_disconnected=False` |
| fan gate | `fan_thresh=90.0` | `fan_thresh=None` |
| orphan re-homing | off | always (off is upstream) |
| `seed` parameter | 42 | `seed=42`, verified to 1e-12 |

Running with `split_disconnected=False, fan_thresh=None, orphan_adjacency=False,
seed=42` gives upstream's numbers exactly. That is how the golden value in
`tests/test_shape_job.py` was re-derived when it moved.

**These changes should go past their author before anything built on them is
published** — they alter δK for every branched geometry, by up to 85% on cases with
disconnected patches.

**On hashes.** `.gitattributes` stores `*.py` with LF endings, and the upstream file
is CRLF, so a plain SHA-256 of the two files does *not* match even though the code
is identical, 31,979 bytes there against 31,168 here, purely line endings. Compare
the **git blob SHA-1** above instead, which is computed after normalisation and is
therefore the same on every platform:

```powershell
git rev-parse HEAD:camp/campsrc.py      # must equal the value in the table
```

For the record, the upstream file as it sits on disk is SHA-256
`6C236F288C413AEB5EBCE37686ACC697064479B90629E9B2D24067795FF2D959`. Verified at
copy time that the two differ only in line endings, 812 lines each, byte-identical
once CRLF is normalised to LF.

## Why a copy at all

There were three files named `campsrc.py` in the source project, and they are not
the same method:

| path | size | |
|---|---|---|
| `camp/campsrc.py` | 34,833 B | older, `seed_patches=250`, different trimming |
| `camp/UpdatedCamp/campsrc.py` | 31,979 B | **this one** |
| `camp/UpdatedCamp/camp/campsrc.py` | 36,452 B | older still (June) |

A plain `import campsrc` resolves to whichever lands on `sys.path` first, and the
wrong one returns a different δK without complaining. One copy in one place removes
that failure mode.

## Refreshing it

When the lab's version moves on:

```powershell
Copy-Item <upstream>\campsrc.py camp\campsrc.py -Force
git add camp\campsrc.py
git rev-parse :camp/campsrc.py          # the new blob SHA-1, for the table above
```

Then re-run the δK checks, the partition is sensitive to its parameters, and a
changed patch count or `k_fluct` should be a decision, not a surprise.

## Dependencies

Self-contained apart from third-party packages: numpy, pandas, trimesh, scipy,
plotly, scikit-learn. All are pinned in the project's `requirements.txt`.

## Overriding it

`shape_metrics.py` uses this vendored copy by default. To run against a different
build without touching the repository, point `CAMP_SRC` at it:

```powershell
$env:CAMP_SRC = "C:\path\to\other\campsrc.py"
```
