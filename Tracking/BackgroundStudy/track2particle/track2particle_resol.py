#!/usr/bin/env python3
"""Plot DIS momentum and polar-angle resolutions from saved CSV chunks.

This script never reopens reconstructed ROOT files. It joins each good signal
track to its dominant MC particle using ``(event, particle_id)``, calculates
residuals in truth momentum and eta bins, fits Gaussian widths, and writes the
standard seven-eta-panel plots.

Example:

    python track2particle_resol.py \
        --config epic \
        --input-root /path/to/rootfiles \
        --chunk-root /path/to/track2particle_output
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Keep local analysis modules importable when this file is launched directly.
MODULE_DIRECTORY = Path(__file__).resolve().parent
if str(MODULE_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(MODULE_DIRECTORY))

import epic_analysis_base as ana
import track2particle as track_analysis


DEFAULT_INPUT_ROOT = Path("rootfiles")
DEFAULT_CHUNK_ROOT = Path("track2particle_output")
DEFAULT_MOMENTUM_BINS = np.array([0.0, 0.5, 1.0, 5.0, 1000.0])

# These are the standard tracking-performance eta regions. They are fixed so
# results from different configurations always have directly comparable rows.
ETA_BIN_RANGES = [
    (-3.5, -3.0),
    (-3.0, -2.5),
    (-2.5, -1.0),
    (-1.0, 1.0),
    (1.0, 2.5),
    (2.5, 3.0),
    (3.0, 3.5),
]

RESOLUTION_SPECS = {
    "dp": {
        "column": "resol_dp",
        "ylabel": r"$\delta p/p$ [%]",
        "scale": 1.0,
        "default_y_hi": 16.0,
        "filename": "tracking_dp_over_p_resolution.png",
    },
    "theta": {
        "column": "resol_theta",
        "ylabel": r"$\theta$ [rad]",
        "scale": 1.0 / 1000.0,
        "default_y_hi": 0.01,
        "filename": "tracking_theta_resolution.png",
    },
}


def load_matched_tracks(input_files, chunk_root, cuts):
    """Load good signal tracks and join each to its dominant MC particle."""
    matched_frames = []
    status_rows = []

    for input_file in input_files:
        output_directory = track_analysis.get_output_directory(
            input_file, chunk_root
        )
        manifest_file = output_directory / "manifest.csv"
        if not manifest_file.is_file():
            raise FileNotFoundError(f"Missing manifest: {manifest_file}")

        manifest = pd.read_csv(manifest_file)
        trajectories, particles = track_analysis.load_active_chunks(
            manifest, manifest_file
        )
        _, good_signal_tracks = track_analysis.select_valid_tracks(
            trajectories, cuts
        )

        particle_key = ["event", "particle_id"]
        if particles.duplicated(particle_key).any():
            raise ValueError(f"Particle keys are not unique in {manifest_file}")

        truth_columns = [
            "event",
            "particle_id",
            "mom",
            "theta",
            "eta",
            "PDG",
            "generatorStatus",
        ]
        missing_columns = [
            column for column in truth_columns if column not in particles
        ]
        if missing_columns:
            raise KeyError(
                f"{manifest_file} is missing particle columns "
                f"{missing_columns}"
            )

        # A strict-majority source is unique, so this is a many-tracks-to-one-
        # particle join. An inner join excludes tracks whose saved truth row is
        # unavailable and records that loss in the per-file summary.
        matched = good_signal_tracks.merge(
            particles[truth_columns],
            left_on=["event", "most_common_source"],
            right_on=particle_key,
            how="inner",
            validate="many_to_one",
        )
        matched["source_file"] = str(input_file)
        matched_frames.append(matched)
        status_rows.append(
            {
                "input_file": str(input_file),
                "events": len(manifest),
                "events_ok": int(manifest["status"].eq("ok").sum()),
                "good_signal_tracks": len(good_signal_tracks),
                "truth_matched_tracks": len(matched),
            }
        )

    matched_tracks = (
        pd.concat(matched_frames, ignore_index=True)
        if matched_frames
        else pd.DataFrame()
    )
    return matched_tracks, pd.DataFrame(status_rows)


def calculate_residuals(matched_tracks):
    """Calculate reconstructed-minus-truth momentum and theta residuals."""
    tracks = matched_tracks.copy()

    # eta = -log(tan(theta/2)); reconstructed eta is stored in the trajectory
    # chunks, so this recovers the fitted polar angle without reopening ROOT.
    tracks["reco_theta"] = 2.0 * np.arctan(np.exp(-tracks["reco_eta"]))
    tracks["resol_dp"] = (
        (tracks["reco_mom"] - tracks["mom"]) / tracks["mom"] * 100.0
    )
    tracks["resol_theta"] = (
        (tracks["reco_theta"] - tracks["theta"]) * 1000.0
    )

    required = [
        "mom",
        "eta",
        "reco_mom",
        "reco_eta",
        "resol_dp",
        "resol_theta",
    ]
    finite = np.isfinite(tracks[required]).all(axis=1)
    finite &= tracks["mom"] > 0
    return tracks[finite].copy()


def summarize_resolution_bins(tracks, momentum_bins, setting):
    """Fit residual widths in each truth eta and momentum bin."""
    rows = []
    for eta_lo, eta_hi in ETA_BIN_RANGES:
        eta_selected = tracks[
            tracks["eta"].between(eta_lo, eta_hi, inclusive="left")
        ]
        for momentum_index, (momentum_lo, momentum_hi) in enumerate(
            zip(momentum_bins[:-1], momentum_bins[1:])
        ):
            selected = eta_selected[
                (eta_selected["mom"] >= momentum_lo)
                & (eta_selected["mom"] < momentum_hi)
            ]
            row = {
                "setting": setting,
                "eta_lo": eta_lo,
                "eta_hi": eta_hi,
                "momentum_bin": momentum_index,
                "momentum_lo": momentum_lo,
                "momentum_hi": momentum_hi,
                # Truth momentum is used on the x axis to avoid a resolution-
                # dependent bin-position bias.
                "mom_gev": (
                    selected["mom"].mean() if not selected.empty else np.nan
                ),
                "n_tracks": len(selected),
            }
            for column in ("resol_dp", "resol_theta"):
                mean, sigma, sigma_error = ana.hist_gaus(
                    selected[column], ax=None, bins=101
                )
                row[f"{column}_mean"] = mean
                row[f"{column}_sigma"] = sigma
                row[f"{column}_sigma_err"] = sigma_error
            rows.append(row)
    return pd.DataFrame(rows)


def load_pwg_requirements(pwg_file):
    """Read an optional PWG requirement table for the momentum plot."""
    if pwg_file is None:
        return None
    if not Path(pwg_file).is_file():
        raise FileNotFoundError(f"Missing PWG requirement table: {pwg_file}")
    return pd.read_csv(pwg_file, sep=r"\s+", skiprows=1)


def pwg_dp_requirement(pwg_table, eta, momentum):
    """Evaluate the PWG dp/p requirement for one eta region."""
    if pwg_table is None:
        return None
    selected = pwg_table[
        (pwg_table["eta_lo"] <= eta) & (pwg_table["eta_hi"] > eta)
    ]
    if selected.empty:
        return None
    coefficient_a = selected["dp_par1"].iloc[0]
    coefficient_b = selected["dp_par2"].iloc[0]
    return np.sqrt((coefficient_a * momentum) ** 2 + coefficient_b**2)


def adaptive_y_limits(default_y_hi, plotted_values):
    """Choose a power-of-two y range that contains and resolves all points."""
    values = np.asarray(plotted_values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return -0.05 * default_y_hi, default_y_hi

    data_max = float(values.max())
    y_hi = float(default_y_hi)
    while data_max > y_hi:
        y_hi *= 2.0
    while data_max <= 0.5 * y_hi:
        y_hi /= 2.0
    return -0.05 * y_hi, y_hi


def resolution_x_limits(summary):
    """Return one x limit per eta panel, including every truth-momentum point."""
    limits = []
    for eta_lo, eta_hi in ETA_BIN_RANGES:
        panel = summary[
            np.isclose(summary["eta_lo"], eta_lo)
            & np.isclose(summary["eta_hi"], eta_hi)
        ]
        finite_momentum = panel["mom_gev"].replace(
            [np.inf, -np.inf], np.nan
        ).dropna()
        data_limit = 1.10 * finite_momentum.max() if len(finite_momentum) else 0
        limits.append(max(15.0, data_limit))
    return limits


def plot_resolution(summary, variable, setting, pwg_table=None):
    """Create one standard seven-panel resolution figure."""
    spec = RESOLUTION_SPECS[variable]
    sigma_column = f"{spec['column']}_sigma"
    error_column = f"{spec['column']}_sigma_err"
    x_limits = resolution_x_limits(summary)

    figure, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.ravel()
    legend_handle = None
    requirement_handle = None

    for panel_index, (eta_lo, eta_hi) in enumerate(ETA_BIN_RANGES):
        axis = axes[panel_index]
        panel = summary[
            np.isclose(summary["eta_lo"], eta_lo)
            & np.isclose(summary["eta_hi"], eta_hi)
        ].dropna(subset=["mom_gev", sigma_column, error_column])
        panel = panel[panel[sigma_column] > 0]

        y_values = panel[sigma_column] * spec["scale"]
        if not panel.empty:
            legend_handle = axis.errorbar(
                panel["mom_gev"],
                y_values,
                yerr=panel[error_column] * spec["scale"],
                color="tab:blue",
                linestyle="none",
                marker="o",
                label=setting,
            )

        if variable == "dp" and pwg_table is not None:
            x_line = np.linspace(0.001, x_limits[panel_index], 1000)
            y_line = pwg_dp_requirement(pwg_table, eta_lo, x_line)
            if y_line is not None:
                requirement_handle, = axis.plot(
                    x_line, y_line, "k--", label="PWG requirement"
                )

        y_lo, y_hi = adaptive_y_limits(spec["default_y_hi"], y_values)
        axis.set_ylim(y_lo, y_hi)
        axis.set_xlim(0, 1.05 * x_limits[panel_index])
        axis.text(
            0.08,
            0.9,
            f"{eta_lo}<$\\eta$<{eta_hi}",
            fontsize=14,
            transform=axis.transAxes,
        )

    # The eighth panel is reserved for a shared legend.
    axes[-1].axis("off")
    handles = [handle for handle in (legend_handle, requirement_handle) if handle]
    if handles:
        axes[-1].legend(handles=handles, frameon=False, loc="upper left")
    for axis in axes[4:7]:
        axis.set_xlabel("momentum [GeV/c]")
    axes[0].set_ylabel(spec["ylabel"])
    axes[4].set_ylabel(spec["ylabel"])
    figure.tight_layout()
    return figure


def make_resolution_plots(summary, output_directory, setting, pwg_file):
    """Write individual PNG plots and a two-page combined PDF."""
    output_directory.mkdir(parents=True, exist_ok=True)
    pwg_table = load_pwg_requirements(pwg_file)
    combined_pdf = output_directory / f"tracking_resolutions_{setting}.pdf"

    with PdfPages(combined_pdf) as pdf:
        for variable, spec in RESOLUTION_SPECS.items():
            figure = plot_resolution(summary, variable, setting, pwg_table)
            figure.savefig(output_directory / spec["filename"], dpi=160)
            pdf.savefig(figure)
            plt.close(figure)
    return combined_pdf


def parse_arguments():
    """Parse and validate the command-line interface."""
    parser = argparse.ArgumentParser(
        description="Plot dp/p and theta resolutions from track2particle chunks."
    )
    parser.add_argument("--config", default="epic")
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--chunk-root", type=Path, default=DEFAULT_CHUNK_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--file-start", type=int, default=1)
    parser.add_argument("--file-stop", type=int, default=100)
    parser.add_argument(
        "--momentum-bins",
        nargs="+",
        type=float,
        default=DEFAULT_MOMENTUM_BINS.tolist(),
        metavar="EDGE",
    )
    parser.add_argument("--pwg-file", type=Path, default=None)
    args = parser.parse_args()

    momentum_bins = np.asarray(args.momentum_bins, dtype=float)
    if (
        len(momentum_bins) < 2
        or not np.all(np.isfinite(momentum_bins))
        or not np.all(np.diff(momentum_bins) > 0)
    ):
        parser.error("--momentum-bins must be finite increasing edges")
    if args.file_start < 1 or args.file_stop < args.file_start:
        parser.error("--file-start/--file-stop define an invalid range")
    args.momentum_bins = momentum_bins
    return args


def main():
    """Load existing chunks, fit resolution widths, and write tables/plots."""
    args = parse_arguments()
    output_directory = args.output_dir or (
        args.chunk_root / "resolution" / args.config
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    input_files = track_analysis.build_input_files(
        args.config, args.input_root, args.file_start, args.file_stop
    )
    matched_tracks, file_summary = load_matched_tracks(
        input_files, args.chunk_root, track_analysis.DEFAULT_CUTS
    )
    tracks = calculate_residuals(matched_tracks)
    summary = summarize_resolution_bins(
        tracks, args.momentum_bins, args.config
    )

    # Preserve both the fitted products and enough joined rows to audit them.
    file_summary.to_csv(output_directory / "input_file_summary.csv", index=False)
    tracks.to_csv(output_directory / "resolution_tracks.csv.gz", index=False)
    summary.to_csv(output_directory / "resolution_summary.csv", index=False)
    combined_pdf = make_resolution_plots(
        summary, output_directory, args.config, args.pwg_file
    )

    print(f"Input files:            {len(input_files)}")
    print(
        "Events OK/requested:    "
        f"{file_summary['events_ok'].sum()}/{file_summary['events'].sum()}"
    )
    print(f"Good signal tracks:     {file_summary['good_signal_tracks'].sum()}")
    print(f"Joined to saved truth:  {file_summary['truth_matched_tracks'].sum()}")
    print(f"Finite residual rows:   {len(tracks)}")
    print(f"Eta bins:               {ETA_BIN_RANGES}")
    print(f"Momentum bins [GeV]:    {args.momentum_bins.tolist()}")
    print(
        "Successful bin fits:    "
        f"dp/p={summary['resol_dp_sigma'].notna().sum()}, "
        f"theta={summary['resol_theta_sigma'].notna().sum()} "
        f"(of {len(summary)})"
    )
    print(f"Summary table:          {output_directory / 'resolution_summary.csv'}")
    print(f"Plots:                  {combined_pdf}")


if __name__ == "__main__":
    main()
