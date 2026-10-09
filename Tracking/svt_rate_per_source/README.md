# SVT reconstructed-hit rates by source

Read EICrecon `events` trees, assign each SVT reconstructed hit to its MC source,
and produce surface maps and rate tables. Copy this entire folder anywhere;
it needs no files or imports from the parent analysis project, ROOT bindings,
PODIO bindings, or detector installation. Input data and Python packages are
the only external requirements.

| File | Purpose |
| --- | --- |
| `plot_hit_rate_per_source.py` | Command-line workflow: chunked input → counts → rates → PDFs/CSV. |
| `svt_sources.py` | Explicit SVT windows, relation names, source mapping, and ancestry tracing. |
| `requirements.txt` | Four required Python packages; Python 3.9 or newer. |
| `test_svt_rate.py` | Offline numerical, relation, chunking, and output checks. |

## Start here

Use an existing analysis environment, or create one from this folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
# Only for root:// inputs; omit when using local files:
python -m pip install XRootD
python plot_hit_rate_per_source.py --help
```

A short check on a local file, followed by a full run:

```bash
python -u plot_hit_rate_per_source.py \
  --input-files /path/to/sample.eicrecon.edm4eic.root \
  --label check --entry-stop 10 --chunk-size 5

python -u plot_hit_rate_per_source.py \
  --input-files /path/to/sample.eicrecon.edm4eic.root \
  --label full --chunk-size 50
```

Outputs default to `outputs/LABEL/` beside the script. Use a new label for
each run; nonempty output folders are protected from overwriting. Save your
production command in a shell script before running it.

## Several files or remote input

Pass several paths after `--input-files`, or use `--file-list inputs.txt`.
The text list contains one local path or full URL per line; blank lines and
lines beginning with `#` are ignored. Local paths in a list are relative to
the list's directory; command-line paths are relative to your working directory.
Both forms accept an optional `:events` suffix. Duplicate inputs are rejected.

JLab's current XRootD layout is:

```text
root://dtn2304.jlab.org:8443//jlab-osdf-ro/eic/EPIC/volatile/RECO/<campaign>/<dataset>/<file>.root
```

Use an actual file URL from the [EIC production file catalog/access guide](https://eic.github.io/software/production_file.html);
availability depends on the campaign. BNL and other full URLs are accepted
unchanged. This release opens exactly the supplied locations; it does not
switch servers automatically.

```bash
python -u plot_hit_rate_per_source.py \
  --file-list inputs.txt --label production --chunk-size 50

# Same inputs and normalization, finer spatial bins:
python -u plot_hit_rate_per_source.py \
  --file-list inputs.txt --label production_tile --chunk-size 50 --tile
```

| Option | Interpretation |
| --- | --- |
| `--entry-start A --entry-stop B` | Process `[A, B)` separately in every file; stop is clipped to its length. |
| `--max-files N` | First N listed files; convenient for a quick check. |
| `--chunk-size N` | Memory/IO control only; default 10 events. Reduce for busy background frames. |
| `--frame-us T` | Duration of each equal-weight frame; default 2 µs. Controls rate normalization. |
| `--tile` | 3.5 × 9.8 mm bins; default is 20 × 20 mm, called RSU here. |
| `--min-surface-hits N` | Omit sparse PDF pages; preserves all counts and CSV rows. Default 1. |
| `--output-dir DIR` | Write under `DIR/LABEL/`. |

## Physics conventions — read before interpreting rates

Each selected event must represent the same physical duration and have equal
statistical weight. Zero-hit frames count in the denominator. Files are pooled
by event count, not averaged with equal weight per file.

```text
N = total number of selected event frames across all files
T = duration per frame in µs
rate_scale = 1000 / (N × T)                  [kHz per counted hit]
bin_rate = accumulated hits in bin × rate_scale
surface_rate = sum of in-map bin rates
peak_bin_rate = maximum of the accumulated bin rates
```

At T = 2 µs, this reproduces the original `500/N` scale. It estimates hit
arrival rates, not the probability of occupancy, event rate, or a rate per area.
Single-particle or triggered/weighted DIS samples need an appropriate physical
normalization; the default 2 µs assumption does not establish one.

The source path is `RecHit → RawHit → RawHitAssociation → SimHit → MCParticle`.
For multiple associations, choose the largest weight; equal weights choose
the first row. For `generatorStatus == 0`, repeatedly follow the first MC
parent to a nonzero status. Cycles, missing parents, invalid parent ranges,
and chains beyond 100 links remain unresolved and enter `other`. ObjectID
indices are event-local within the explicitly named collections; collection
IDs are not remapped to alternative detector collections.

| MC generatorStatus | Source label |
| --- | --- |
| 1, 2 | DIS |
| 2001, 2002 | SR |
| 3001, 3002 | Bremstrahlung (spelling retained from the parent workflow) |
| 4001, 4002 | Coulomb |
| 5001, 5002 | Touschek |
| 6001, 6002 | Proton beam gas |
| Other values, including unresolved 0 | other |

`unmatched` means a valid raw-hit index has no association. `invalid_relation`
covers null rec-hit relations, invalid selected sim/particle indices, and
nonfinite or negative weights. Diagnostics are always saved in `*_run.json`.
`ALL` includes every in-map rec hit, including these categories.
`other background` is the inclusive sum of the four non-SR beam backgrounds;
do not add it or `ALL` to the exclusive sources when computing a total.
The program checks exclusive-source conservation against `ALL` bin by bin.

Barrels use `(z, rφ)`, with `φ = atan2(y, x)` in radians and r in mm; disks
use `(x, y)`. The five barrel and ten disk windows are listed in
`svt_sources.py`; window endpoints are excluded. These fixed projected grids
approximate RSU/tile dimensions, not the physical sensor segmentation or
cellID layout. Bins are anchored at the negative map edge and have full nominal
widths; the final bin may extend beyond the nominal map extent. Selected hits
outside the grids appear in `outside_map_hits` and are excluded from rates.
Update the windows when using a different geometry.

## Outputs and reuse

| `LABEL_…` output | Contents |
| --- | --- |
| `SOURCE_svt_hit_rate.pdf` | One page per populated surface; color is bin rate in kHz. |
| `summary.csv` | Long table by `lname`, `source`, and `rec_branch`; rates, raw counts, occupied-frame counts, and out-of-map counts. |
| `input_files.csv` | Exact URLs/paths, available events, actual `[start, stop)` ranges, and processed counts, including skipped files. |
| `maps.npz` | Integer raw histograms and bin edges in mm; axes are `(x, y)` as defined above. |
| `run.json` | Normalization, arguments, versions, code hashes, windows, relations, diagnostics, and archive-key index. |

`rsu_or_tile_max_rate_khz` is the peak of the time-averaged map;
`max_raw_hits_per_event_tile` is the largest bin count in any single frame
(the historical column name is retained for RSU runs too).
Absent source/surface rows mean no hits were accumulated for that source.

CSV headers beginning with `#` contain normalization metadata. Example:

```python
import pandas as pd  # Optional; not required to run this release.
summary = pd.read_csv("outputs/full/full_summary.csv", comment="#")
```

Recover a rate map without rereading ROOT:

```python
import json
import numpy as np
run = json.load(open("outputs/full/full_run.json"))
with np.load("outputs/full/full_maps.npz") as maps:
    entry = next(m for m in run["maps"] if m["lname"] == "L0" and m["source"] == "ALL")
    rate_khz = maps[entry["key"]] * run["rate_scale_khz_per_count"]
    xedges_mm = maps[entry["xedges_key"]]
    yedges_mm = maps[entry["yedges_key"]]
```

