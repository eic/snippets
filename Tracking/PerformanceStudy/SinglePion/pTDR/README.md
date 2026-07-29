# Single Particle Performance Workflow

`ShujieLi lbl gov. March 2026`

This module provides a workflow to scan single pion simulation files to produce the __pTDR style__ tracking resolution v.s. momentum summary plot (7 eta bins). 
The tracking efficiency, resolution, and pull summaries are extracted using the same logic as `plot_single_resol.py` and the deprecated __7-eta-bin__ folder.


## Requirements

- Python environment with `uproot`, `awkward`, `pandas`, `matplotlib`, `numpy`, `lmfit`.

## Inputs and naming

Each input simulation file is expected to cover a combination of: 
* Default momentum: 0.5, 1, 2, 5, 10, 15 GeV
* Default eta range: ("-3.5 -3.0" "-3.0 -2.5" "-2.5 -1" "-1 1"  "1 2.5" "2.5 3.0" "3.0 3.5")

Example rootfile generation scripts: `submit_batch.sh`, `run_sim_nohup.sh`

The analysis scripts read eta and momentum from the filename tag. Files should start with `rec_` and include:
- `eta_<eta_lo>_<eta_hi>`
- `mom_<mom>GeV` (or `<mom>GeV` / `<mom>MeV`)

Example:

```
rec_GoldCoating5um_eta_-1_1_mom_5GeV_n10000.root
```



## Run the study

`run_study.py` extracts tracking efficiency, resolution, and pull distribution for each rootfile using the same logic as `plot_single_resol.py`. Results are appended as one row to the CSV table `performance_table.csv`. Corresponding plots are saved under `./plots`:

```bash

python run_study.py \
  /path/to/rec_*.root \
  --plots-dir plots \
```

Batch run with skip-existing:

```bash
python run_study.py \
  --pattern "*GoldCoating*.root" \
  --skip-existing \
```

## Make summary plots

`run_plot.py` reads the `performance_table.csv` table and produces a single multi-page PDF with the summary efficiency + dp/theta/phi/dca plots in pTDR format:

```bash
python run_plot.py \
  --setting GoldCoating5um
```

Compare two settings:

```bash
python run_plot.py \
  --setting GoldCoating5um \
  --setting2 GoldCoating10um
```

## Outputs

- `performance_table.csv`: append-only summary table with per-file stats and binned efficiency arrays.
- `plots/`: per-file PDFs from `run_study.py` and combined multi-page PDF from `run_plot.py`.

## Notes

- `pwg_requirements.txt` is loaded from the current working directory, the module folder, or the legacy snippets path.
- Theta/phi resolutions are converted from mrad to rad in plotting.
- Files not starting with `rec_` are skipped.
