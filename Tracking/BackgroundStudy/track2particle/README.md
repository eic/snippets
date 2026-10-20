# Hit-based track-to-particle analysis for background study
Shujie Li, 08.2026, code generated with assistance from OpenAI Codex.


## Goal

This package traces reconstructed central-tracker trajectories to Monte Carlo
particles, then measures:

- signal-particle reconstruction efficiency;
- candidate-track purity and background/fake-track fractions;
- hits per trajectory;
- momentum and polar-angle resolutions.

The expensive ROOT/PODIO event loop is run once. Later changes to eta or
momentum plotting bins reuse the compressed per-file CSV chunks.

## Code structure

- `track2particle.py`: event processing, per-file chunks, aggregation, and
  efficiency/purity/track-quality plots.
- `track2particle_resol.py`: resolution fits and the fixed seven-eta-panel
  plots, using existing chunks only.
- `epic_analysis_base.py`: uproot readers, shared cuts, coordinate utilities,
  and Gaussian fitting.
- `epic_analysis_podio.py`: PODIO relation tracing used during event
  processing. This file is a required runtime dependency and should be
  distributed with the other three files.

`track2particle_resol.py` depends only on `track2particle.py` and
`epic_analysis_base.py`; it does not depend on the separate single-particle
performance package.

## Environment

Run within eic-shell to get PODIO suport. The
Python environment needs PODIO/EDM4hep plus `numpy`, `pandas`, `awkward`,
`uproot`, `matplotlib`, `seaborn`, `particle`, and `lmfit`.

The default input and output roots are relative paths named `rootfiles/` and
`track2particle_output/`. For reproducible work, pass both paths explicitly.

## Input

The input is an EICrecon ROOT file containing at least:

- `CentralCKFTrajectories` and `CentralCKFTrackParameters`;
- `MCParticles`;
- the central tracking hit-to-sim-hit-to-MC-particle relation chain expected
  by `epic_analysis_podio.py`.

For a numbered campaign, files are expected under
`INPUT_ROOT/CONFIG/` with names of the form:

```text
rec_pythia8NCDIS_..._1.0001_CONFIG_n99_skip0.root
```

Use repeated `--input-file /full/path/file.root` arguments for files that do
not follow this convention.

## Minimal use

Process one file first:

```bash
python track2particle.py \
    --input-file /path/to/rec_file.root \
    --output-root /path/to/track2particle_output \
    --analysis-label my_configuration \
    --events 10 --workers 1 \
    --process-only --require-all-events
```

Process and aggregate a 100-file campaign:

```bash
python track2particle.py \
    --config epic \
    --input-root /path/to/rootfiles \
    --output-root /path/to/track2particle_output \
    --file-start 1 --file-stop 100 \
    --events 99 --workers 8 \
    --eta-range -1 1 0.1 \
    --require-all-events
```

Rebuild performance plots without reopening ROOT files:

```bash
python track2particle.py \
    --config epic \
    --input-root /path/to/rootfiles \
    --output-root /path/to/track2particle_output \
    --file-start 1 --file-stop 100 \
    --eta-range -1 1 0.1 \
    --aggregate-only --require-all-events
```

Build the resolution tables and fixed seven-eta-panel plots:

```bash
python track2particle_resol.py \
    --config epic \
    --input-root /path/to/rootfiles \
    --chunk-root /path/to/track2particle_output \
    --file-start 1 --file-stop 100
```

Run either script with `--help` for all options.

## Physics definitions

The default cuts are imported from `epic_analysis_base.py`. At the time of
this benchmark they are:

```text
generated momentum > 0.2 GeV
|vertex r| < 1 mm
|vertex z| < 100 mm
at least 4 hits per candidate track
strict majority: dominant hit fraction > 0.5
```

Every trajectory hit enters the hit-fraction denominator. Invalid or noise
relations use particle ID `-1`.

- `good_signal`: strict-majority match to generator status 1 or 2.
- `background_track`: strict-majority match to another real MC particle.
- `fake_or_ghost`: no strict-majority real-particle match.

These three classes are mutually exclusive and exhaustive for candidate
tracks. A selected signal particle is reconstructed when at least one
candidate track has a strict-majority match to that particle.

Efficiency uses truth momentum and eta:

```text
selected signal particles with a matched track / selected signal particles
```

Purity uses reconstructed momentum and eta:

```text
good signal tracks / all candidate tracks
```

Resolution bins use truth momentum and eta. The stored residuals are:

```text
dp/p [%]      = 100 * (p_reco - p_truth) / p_truth
delta theta   = theta_reco - theta_truth
theta_reco    = 2 * atan(exp(-eta_reco))
```

Gaussian widths are fitted with the robust clipping procedure in
`epic_analysis_base.hist_gaus`. Empty or under-populated bins are written as
`NaN`, not zero.

## Output layout

```text
track2particle_output/
├── per_run_analysis/
│   └── rec_<input-file-tag>/
│       ├── analysis_settings.json
│       ├── analysis.log
│       ├── manifest.csv
│       ├── summary.csv
│       └── chunks/
│           ├── particles_<event-range>.csv.gz
│           └── trajectories_<event-range>.csv.gz
├── <analysis-label>/performance/
│   └── eta_<range>_step_<width>/
│       ├── file_summary.csv
│       ├── efficiency_bins.csv
│       ├── purity_bins.csv
│       ├── performance_summary.csv
│       ├── track_class_bins.csv
│       ├── track_class_summary.csv
│       ├── hits_per_trajectory.csv
│       └── *.png
└── resolution/<configuration>/
    ├── input_file_summary.csv
    ├── resolution_tracks.csv.gz
    ├── resolution_summary.csv
    └── *.png, *.pdf
```

Manifest chunk paths are relative to the per-file directory, so the complete
output tree can be moved. Older manifests with absolute paths are still
readable.

Important chunk columns are:

- particle chunks: `event`, `particle_id`, `mom`, `theta`, `eta`, `PDG`, and
  `generatorStatus`;
- trajectory chunks: `event`, `traj_id`, `reco_mom`, `reco_eta`, `reco_pt`,
  `total_count`, `max_count`, `max_fraction`, `most_common_source`, and
  `part_status`.

The binned efficiency/purity tables contain bin edges, numerator,
denominator, value, and a 68% Wilson interval. `resolution_summary.csv`
contains the fitted mean, sigma, and sigma uncertainty for every one of the
seven eta regions and every truth-momentum bin.

## Resume and changing cuts

`--resume` reuses an existing per-file result only when its recorded input,
event range, and event-level cuts match and every named chunk exists.

Changing `--momentum-min` changes which MC particles are saved during event
processing; use a new output root and rerun the event loop. Eta ranges,
momentum plot bins, and track hit-count/fraction cuts are applied during
aggregation and can reuse existing chunks, provided the needed particles were
saved.

## Performance benchmark

Measured on NERSC Perlmutter CPU nodes with files on local NERSC storage:

| Stage | Data | Wall time | Peak memory |
|---|---:|---:|---:|
| Event processing, one worker | 99 events in one file | median 57.4 s; 95th percentile 58.6 s | median 1.69 GiB; 95th percentile 1.70 GiB |
| Aggregation | 100 files, 9,900 events | 12.2 s | 226 MiB |
| Resolution join, fits, and plots | 100 files, 9,900 events | 9.0 s | 242 MiB |

The event-processing statistics use 200 independently measured noise-sample
files. Memory scales approximately with the number of simultaneous file
workers, so choose `--workers` from both available CPU cores and roughly
1.8 GiB of memory per worker. Aggregation and resolution are single-process
and much lighter.
