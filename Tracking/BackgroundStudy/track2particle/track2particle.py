#!/usr/bin/env python3
"""Trace CKF trajectories to MC particles for a reconstructed DIS file set.

The expensive PODIO reader is opened once per input file.  Parallelism is
therefore across independent files, rather than across event ranges within one
file.  This is substantially more efficient for the 100-file DIS samples.

Run this script inside a matching EIC environment. This example processes a
numbered reconstructed-file campaign:

    python track2particle.py \
        --config epic_noise \
        --input-root /path/to/rootfiles \
        --file-start 1 --file-stop 100 \
        --events 99 --workers 8 \
        --output-root /path/to/track2particle_output \
        --eta-range -1 1 0.1 \
        --require-all-events

Omit the momentum options to use the shared defaults from
``epic_analysis_base.py``. Use a new output root when changing an event-level
cut, because that cut affects the saved per-file particle chunks.
"""

from __future__ import annotations

import argparse
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import io
import json
from multiprocessing import get_context
from pathlib import Path
import resource
import signal
import sys
import time

import awkward as ak
import numpy as np
import pandas as pd

import epic_analysis_base as ana


# Relative defaults make a cloned copy usable without editing personal paths.
DEFAULT_INPUT_ROOT = Path("rootfiles")
DEFAULT_OUTPUT_ROOT = Path("track2particle_output")
DEFAULT_MODULE_DIRECTORY = Path(__file__).resolve().parent
PER_RUN_DIRECTORY = "per_run_analysis"

DEFAULT_CUTS = {
    "momentum_min": ana.TRACK_MOM_MIN,  # GeV
    "vertex_r_max": ana.VERTEX_CUT_R_MAX,  # mm
    "vertex_z_max": ana.VERTEX_CUT_Z_MAX,  # mm
    "track_hit_count_min": ana.TRACK_HIT_COUNT_MIN,
    "track_hit_fraction_min": ana.TRACK_HIT_FRACTION_MIN,
}

DEFAULT_MOMENTUM_BINS = np.array([0.0, 0.5, 1.0, 5.0, 1000.0])
DEFAULT_ETA_BINS = np.round(np.arange(-4.0, 4.01, 0.1), 1)
TRACK_CLASS_ORDER = ("good_signal", "background_track", "fake_or_ghost")
TRACK_CLASS_LABELS = {
    "good_signal": "good signal",
    "background_track": "background track",
    "fake_or_ghost": "fake/ghost",
}


def normalize_eta_bins(eta_bins):
    """Return increasing eta-bin edges.

    A three-value, non-increasing input such as ``[-1, 1, 0.1]`` is
    interpreted as ``[eta_min, eta_max, step]``. An increasing sequence is
    treated as explicit bin edges.
    """
    eta_bins = np.asarray(eta_bins, dtype=float)
    if (
        len(eta_bins) == 3
        and eta_bins[0] < eta_bins[1]
        and eta_bins[2] > 0
        and not np.all(np.diff(eta_bins) > 0)
    ):
        eta_min, eta_max, eta_step = eta_bins
        eta_bins = np.arange(eta_min, eta_max, eta_step)
        eta_bins = np.append(eta_bins, eta_max)

    if (
        eta_bins.ndim != 1
        or len(eta_bins) < 2
        or not np.all(np.isfinite(eta_bins))
        or not np.all(np.diff(eta_bins) > 0)
    ):
        raise ValueError(
            "eta_bins must be increasing bin edges or "
            "[eta_min, eta_max, step]"
        )
    return eta_bins


class EventTimeout(Exception):
    """Raised when one event exceeds the configured analysis time limit."""


def stop_slow_event(signum, frame):
    """Convert a worker alarm into an exception handled per event."""
    del signum, frame
    raise EventTimeout("event exceeded the configured time limit")


def get_output_directory(input_file, output_root):
    """Return the per-file output directory for one reconstructed input.

    New analyses live under ``per_run_analysis``. The legacy flat location is
    returned only when it already exists, which keeps older output roots
    readable without affecting the new layout.
    """
    file_tag = Path(input_file).name
    for suffix in (".root", ".edm4eic", ".edm4hep"):
        if file_tag.endswith(suffix):
            file_tag = file_tag[: -len(suffix)]
    output_root = Path(output_root)
    output_directory = output_root / PER_RUN_DIRECTORY / file_tag
    legacy_directory = output_root / file_tag
    if not output_directory.exists() and legacy_directory.exists():
        return legacy_directory
    return output_directory


def manifest_chunk_paths(manifest, manifest_file, column):
    """Resolve unique chunk paths from one manifest column.

    Current manifests store paths relative to their per-file directory, so an
    output tree can be moved intact. Absolute paths from older manifests remain
    supported.
    """
    manifest_directory = Path(manifest_file).parent
    paths = []
    for value in manifest[column].dropna().unique():
        if not str(value):
            continue
        path = Path(value)
        paths.append(path if path.is_absolute() else manifest_directory / path)
    return paths


def build_input_files(config, input_root, file_start, file_stop):
    """Build and validate the expected numbered reconstructed-file paths."""
    input_directory = Path(input_root) / config
    input_files = [
        input_directory
        / (
            "rec_pythia8NCDIS_10x275_minQ2=1_beamEffects_"
            f"xAngle=-0.025_hiDiv_1.{file_index:04d}_{config}_n99_skip0.root"
        )
        for file_index in range(file_start, file_stop + 1)
    ]
    missing_files = [path for path in input_files if not path.is_file()]
    if missing_files:
        preview = "\n".join(f"  {path}" for path in missing_files[:10])
        if len(missing_files) > 10:
            preview += f"\n  ... and {len(missing_files) - 10} more"
        raise FileNotFoundError(
            f"Missing {len(missing_files)} expected input files:\n{preview}"
        )
    return input_files


def manifest_is_reusable(manifest_file, first_event, number_of_events):
    """Check whether a prior manifest covers exactly the requested event range."""
    manifest_file = Path(manifest_file)
    if not manifest_file.is_file():
        return False
    try:
        manifest = pd.read_csv(manifest_file)
    except Exception:
        return False
    required_columns = {
        "event",
        "status",
        "trajectory_file",
        "particle_file",
    }
    if not required_columns.issubset(manifest.columns):
        return False
    expected_events = np.arange(
        first_event,
        first_event + number_of_events,
        dtype=int,
    )
    actual_events = np.sort(manifest["event"].to_numpy(dtype=int))
    if not np.array_equal(actual_events, expected_events):
        return False

    # A manifest is authoritative only if every chunk path it names still
    # exists. Empty paths are allowed for events that produced no output rows.
    for column in ("trajectory_file", "particle_file"):
        paths = manifest_chunk_paths(manifest, manifest_file, column)
        if any(not path.is_file() for path in paths):
            return False
    return True


def settings_are_reusable(settings_file, requested_settings):
    """Return whether saved processing settings match the requested values."""
    settings_file = Path(settings_file)
    if not settings_file.is_file():
        return False
    try:
        with settings_file.open() as stream:
            saved_settings = json.load(stream)
    except (OSError, ValueError):
        return False
    return saved_settings == requested_settings


def analyze_input_file(task):
    """Analyze one complete input file inside one clean worker process."""
    input_file = Path(task["input_file"])
    output_root = Path(task["output_root"])
    module_directory = Path(task["module_directory"])
    first_event = int(task["first_event"])
    number_of_events = int(task["number_of_events"])
    timeout_seconds = int(task["timeout_seconds"])
    cuts = dict(task["cuts"])
    resume = bool(task["resume"])

    output_directory = get_output_directory(input_file, output_root)
    chunk_directory = output_directory / "chunks"
    manifest_file = output_directory / "manifest.csv"
    settings_file = output_directory / "analysis_settings.json"
    worker_log = output_directory / "analysis.log"
    output_directory.mkdir(parents=True, exist_ok=True)
    chunk_directory.mkdir(parents=True, exist_ok=True)

    requested_settings = {
        "schema_version": 1,
        "input_file": str(input_file.resolve()),
        "first_event": first_event,
        "number_of_events": number_of_events,
        "cuts": cuts,
    }

    if resume and manifest_file.exists() and not settings_are_reusable(
        settings_file,
        requested_settings,
    ):
        raise ValueError(
            "Existing manifest was produced with different or unrecorded "
            f"analysis settings: {manifest_file}. Use a different "
            "--output-root for the new cuts."
        )

    if (
        resume
        and settings_are_reusable(settings_file, requested_settings)
        and manifest_is_reusable(
            manifest_file,
            first_event,
            number_of_events,
        )
    ):
        manifest = pd.read_csv(manifest_file)
        return {
            "input_file": str(input_file),
            "output_directory": str(output_directory),
            "manifest_file": str(manifest_file),
            "status": "reused",
            "events_ok": int((manifest["status"] == "ok").sum()),
            "events_problem": int((manifest["status"] != "ok").sum()),
            "elapsed_seconds": 0.0,
            "max_rss_mb": 0.0,
        }

    if manifest_file.exists() and not resume:
        raise FileExistsError(
            f"Output already exists; use --resume or a different output root: "
            f"{manifest_file}"
        )

    start_time = time.monotonic()
    trajectory_tables = []
    particle_tables = []
    status_rows = []

    with settings_file.open("w") as stream:
        json.dump(requested_settings, stream, indent=2, sort_keys=True)
        stream.write("\n")

    # Keep verbose PODIO/helper output in a file beside this input's manifest.
    # The parent process prints only one concise progress line per input file.
    with worker_log.open("w") as log_stream:
        with contextlib.redirect_stdout(log_stream), contextlib.redirect_stderr(
            log_stream
        ):
            print(f"Input: {input_file}")
            print(f"UTC start: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

            if str(module_directory) not in sys.path:
                sys.path.insert(0, str(module_directory))
            import epic_analysis_base as ana_worker
            import epic_analysis_podio as pod_worker

            uproot_tree = ana_worker.read_ur(str(input_file), "events")
            available_events = int(uproot_tree.num_entries)
            stop_event = min(
                first_event + number_of_events,
                available_events,
            )
            event_numbers = np.arange(first_event, stop_event, dtype=int)
            if len(event_numbers) != number_of_events:
                raise ValueError(
                    f"{input_file}: requested {number_of_events} events from "
                    f"{first_event}, but only {len(event_numbers)} are available"
                )

            podio_events = pod_worker.read_podio(str(input_file))
            signal.signal(signal.SIGALRM, stop_slow_event)

            for event_number in event_numbers:
                try:
                    if timeout_seconds > 0:
                        signal.alarm(timeout_seconds)

                    event_number = int(event_number)
                    event = podio_events[event_number]

                    # Suppress repetitive per-event helper messages while
                    # retaining exceptions and the event status in manifest.csv.
                    quiet_output = io.StringIO()
                    with contextlib.redirect_stdout(quiet_output):
                        mc_particles = ak.to_dataframe(
                            ana_worker.get_part(
                                uproot_tree,
                                entry_start=event_number,
                                entry_stop=event_number + 1,
                                kprimary=1,
                            )
                        ).reset_index()
                        track_parameters = ak.to_dataframe(
                            ana_worker.get_params(
                                uproot_tree,
                                "CentralCKFTrackParameters",
                                entry_start=event_number,
                                entry_stop=event_number + 1,
                            )
                        ).reset_index()
                        (
                            selected_particles,
                            trajectories,
                        ) = pod_worker.get_part_traj_counts(
                            event,
                            mc_particles,
                            ksignal=1,
                        )

                    track_parameters = track_parameters.rename(
                        columns={
                            "subentry": "traj_id",
                            "mom": "reco_mom",
                            "eta": "reco_eta",
                            "pt": "reco_pt",
                        }
                    )[["traj_id", "reco_mom", "reco_eta", "reco_pt"]]

                    particle_cut = (
                        (
                            selected_particles["vertex_r"].abs()
                            < cuts["vertex_r_max"]
                        )
                        & (
                            selected_particles["vertex.z"].abs()
                            < cuts["vertex_z_max"]
                        )
                        & (
                            selected_particles["mom"]
                            > cuts["momentum_min"]
                        )
                    )
                    selected_particles = selected_particles[
                        particle_cut
                    ].copy()
                    particle_id_column = (
                        "orig_subentry"
                        if "orig_subentry" in selected_particles.columns
                        else "subentry"
                    )
                    selected_particles["particle_id"] = selected_particles[
                        particle_id_column
                    ].astype(int)

                    trajectories = trajectories.merge(
                        track_parameters,
                        on="traj_id",
                        how="left",
                        validate="one_to_one",
                    )
                    trajectories["event"] = event_number
                    selected_particles["event"] = event_number

                    if not trajectories.empty:
                        trajectory_tables.append(trajectories)
                    if not selected_particles.empty:
                        particle_tables.append(selected_particles)

                    status_rows.append(
                        {
                            "event": event_number,
                            "status": "ok",
                            "n_trajectories": int(len(trajectories)),
                            "n_selected_particles": int(
                                len(selected_particles)
                            ),
                            "error": "",
                        }
                    )
                except EventTimeout:
                    status_rows.append(
                        {
                            "event": int(event_number),
                            "status": "timeout",
                            "n_trajectories": 0,
                            "n_selected_particles": 0,
                            "error": (
                                f"event exceeded the {timeout_seconds}-second "
                                "time limit"
                            ),
                        }
                    )
                except Exception as error:
                    status_rows.append(
                        {
                            "event": int(event_number),
                            "status": "error",
                            "n_trajectories": 0,
                            "n_selected_particles": 0,
                            "error": f"{type(error).__name__}: {error}",
                        }
                    )
                finally:
                    signal.alarm(0)

            first_processed_event = int(event_numbers[0])
            last_processed_event = int(event_numbers[-1])
            chunk_label = (
                f"{first_processed_event:06d}_{last_processed_event:06d}"
            )

            trajectory_file = ""
            if trajectory_tables:
                trajectory_file = (
                    chunk_directory / f"trajectories_{chunk_label}.csv.gz"
                )
                pd.concat(
                    trajectory_tables,
                    ignore_index=True,
                ).to_csv(trajectory_file, index=False)

            particle_file = ""
            if particle_tables:
                particle_file = (
                    chunk_directory / f"particles_{chunk_label}.csv.gz"
                )
                pd.concat(
                    particle_tables,
                    ignore_index=True,
                ).to_csv(particle_file, index=False)

            # Save portable paths relative to this per-file directory. The
            # manifest may then be moved together with its chunks.
            trajectory_manifest_path = (
                str(trajectory_file.relative_to(output_directory))
                if trajectory_file
                else ""
            )
            particle_manifest_path = (
                str(particle_file.relative_to(output_directory))
                if particle_file
                else ""
            )
            for row in status_rows:
                row["trajectory_file"] = (
                    trajectory_manifest_path
                )
                row["particle_file"] = particle_manifest_path

            manifest = pd.DataFrame(status_rows).sort_values("event")
            manifest.to_csv(manifest_file, index=False)
            print(f"UTC end: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    elapsed_seconds = time.monotonic() - start_time
    max_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    return {
        "input_file": str(input_file),
        "output_directory": str(output_directory),
        "manifest_file": str(manifest_file),
        "status": "processed",
        "events_ok": int((manifest["status"] == "ok").sum()),
        "events_problem": int((manifest["status"] != "ok").sum()),
        "elapsed_seconds": elapsed_seconds,
        "max_rss_mb": max_rss_mb,
    }


def process_input_files(
    input_files,
    output_root,
    module_directory,
    first_event,
    number_of_events,
    timeout_seconds,
    cuts,
    number_of_workers,
    resume,
):
    """Process independent input files concurrently."""
    tasks = [
        {
            "input_file": str(input_file),
            "output_root": str(output_root),
            "module_directory": str(module_directory),
            "first_event": first_event,
            "number_of_events": number_of_events,
            "timeout_seconds": timeout_seconds,
            "cuts": cuts,
            "resume": resume,
        }
        for input_file in input_files
    ]
    worker_count = min(max(1, int(number_of_workers)), len(tasks))
    print(
        f"Processing {len(tasks)} files with {worker_count} file workers; "
        f"{number_of_events} events per file"
    )

    results = []
    if worker_count == 1:
        for task_index, task in enumerate(tasks, start=1):
            result = analyze_input_file(task)
            results.append(result)
            print_file_progress(task_index, len(tasks), result)
        return results

    # The standalone parent never imports PODIO or ROOT. Forked workers
    # therefore start without inherited readers and import them independently.
    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=get_context("fork"),
    ) as executor:
        futures = {
            executor.submit(analyze_input_file, task): task["input_file"]
            for task in tasks
        }
        for completed_count, future in enumerate(
            as_completed(futures),
            start=1,
        ):
            input_file = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = {
                    "input_file": input_file,
                    "output_directory": str(
                        get_output_directory(input_file, output_root)
                    ),
                    "manifest_file": "",
                    "status": "worker_error",
                    "events_ok": 0,
                    "events_problem": number_of_events,
                    "elapsed_seconds": np.nan,
                    "max_rss_mb": np.nan,
                    "error": f"{type(error).__name__}: {error}",
                }
            results.append(result)
            print_file_progress(completed_count, len(tasks), result)
    return results


def print_file_progress(completed_count, total_count, result):
    """Print one concise progress line for a completed input file."""
    print(
        f"[{completed_count:03d}/{total_count:03d}] "
        f"{Path(result['input_file']).name}: {result['status']}, "
        f"ok={result['events_ok']}, problem={result['events_problem']}, "
        f"time={result['elapsed_seconds']:.1f}s, "
        f"maxRSS={result['max_rss_mb']:.0f} MB",
        flush=True,
    )
    if result.get("error"):
        print(f"    {result['error']}", flush=True)


def load_active_chunks(manifest, manifest_file=None):
    """Load chunk files named by successful rows in one manifest.

    ``manifest_file`` is required for portable relative chunk paths. It is
    optional only for compatibility with callers holding old absolute paths.
    """
    good_events = manifest[manifest["status"] == "ok"]
    manifest_file = Path(manifest_file or "manifest.csv")
    trajectory_files = manifest_chunk_paths(
        good_events, manifest_file, "trajectory_file"
    )
    particle_files = manifest_chunk_paths(
        good_events, manifest_file, "particle_file"
    )
    trajectories = (
        pd.concat(
            [pd.read_csv(path) for path in trajectory_files],
            ignore_index=True,
        )
        if trajectory_files
        else pd.DataFrame()
    )
    particles = (
        pd.concat(
            [pd.read_csv(path) for path in particle_files],
            ignore_index=True,
        )
        if particle_files
        else pd.DataFrame()
    )
    return trajectories, particles


def select_valid_tracks(trajectories, cuts):
    """Classify reconstructed candidate tracks by their dominant hit source.

    Every trajectory hit contributes to ``total_count``. Noise or invalid
    hits have particle ID -1, while detector-background hits retain their
    real particle ID and generator status. A track is matched only when a
    real particle supplies strictly more than half of all trajectory hits.
    """
    if trajectories.empty:
        return pd.DataFrame(), pd.DataFrame()

    valid_tracks = trajectories[
        trajectories["total_count"].abs() >= cuts["track_hit_count_min"]
    ].copy()
    has_real_source = valid_tracks["most_common_source"] >= 0
    has_majority = (
        valid_tracks["max_fraction"]
        > cuts["track_hit_fraction_min"]
    )
    from_signal = valid_tracks["part_status"].isin([1, 2])
    valid_tracks["is_majority_matched"] = has_real_source & has_majority
    valid_tracks["is_good_signal"] = (
        valid_tracks["is_majority_matched"] & from_signal
    )
    valid_tracks["track_class"] = np.select(
        [
            valid_tracks["is_good_signal"],
            valid_tracks["is_majority_matched"],
        ],
        TRACK_CLASS_ORDER[:2],
        default=TRACK_CLASS_ORDER[2],
    )
    good_signal_tracks = valid_tracks[
        valid_tracks["is_good_signal"]
    ].copy()
    return valid_tracks, good_signal_tracks


def mark_particles_with_tracks(particles, valid_tracks):
    """Mark particles supplying a strict majority of a candidate track's hits."""
    particles = particles.copy()
    if particles.empty:
        particles["has_valid_track"] = pd.Series(dtype=bool)
        return particles

    majority_tracks = valid_tracks[
        valid_tracks["is_majority_matched"]
    ]
    valid_particle_keys = (
        majority_tracks[["event", "most_common_source"]]
        .rename(columns={"most_common_source": "particle_id"})
        .drop_duplicates()
        if not majority_tracks.empty
        else pd.DataFrame(columns=["event", "particle_id"])
    )
    particles = particles.merge(
        valid_particle_keys.assign(has_valid_track=True),
        on=["event", "particle_id"],
        how="left",
    )
    particles["has_valid_track"] = (
        particles["has_valid_track"].fillna(False).astype(bool)
    )
    return particles


def build_binned_performance(
    data,
    passed_column,
    momentum_column,
    eta_column,
    metric_name,
    momentum_bins,
    eta_bins,
):
    """Return one row per momentum/eta bin with a 68% Wilson interval."""
    eta_bins = normalize_eta_bins(eta_bins)
    required_columns = {passed_column, momentum_column, eta_column}
    missing_columns = required_columns - set(data.columns)
    if missing_columns:
        raise KeyError(
            f"{metric_name}: missing columns {sorted(missing_columns)}"
        )

    finite_data = data.dropna(
        subset=[momentum_column, eta_column]
    ).copy()
    finite_data[passed_column] = finite_data[passed_column].astype(bool)
    denominator, _, _ = np.histogram2d(
        finite_data[momentum_column],
        finite_data[eta_column],
        bins=[momentum_bins, eta_bins],
    )
    numerator, _, _ = np.histogram2d(
        finite_data.loc[finite_data[passed_column], momentum_column],
        finite_data.loc[finite_data[passed_column], eta_column],
        bins=[momentum_bins, eta_bins],
    )

    rows = []
    for momentum_index in range(len(momentum_bins) - 1):
        for eta_index in range(len(eta_bins) - 1):
            total = int(denominator[momentum_index, eta_index])
            passed = int(numerator[momentum_index, eta_index])
            value = np.nan
            interval_low = np.nan
            interval_high = np.nan
            if total > 0:
                value = passed / total
                scale = 1.0 + 1.0 / total
                center = (value + 0.5 / total) / scale
                half_width = (
                    np.sqrt(
                        value * (1.0 - value) / total
                        + 0.25 / total**2
                    )
                    / scale
                )
                interval_low = max(0.0, center - half_width)
                interval_high = min(1.0, center + half_width)

            rows.append(
                {
                    "metric": metric_name,
                    "momentum_bin": momentum_index,
                    "momentum_low_GeV": momentum_bins[momentum_index],
                    "momentum_high_GeV": momentum_bins[
                        momentum_index + 1
                    ],
                    "eta_bin": eta_index,
                    "eta_low": eta_bins[eta_index],
                    "eta_high": eta_bins[eta_index + 1],
                    "eta_center": 0.5
                    * (eta_bins[eta_index] + eta_bins[eta_index + 1]),
                    "numerator": passed,
                    "denominator": total,
                    "value": value,
                    "interval_low": interval_low,
                    "interval_high": interval_high,
                }
            )
    return pd.DataFrame(rows)


def _plot_fraction_distribution(
    denominator_data,
    numerator_data,
    metric_label,
    beam,
    setting,
    nhit_cut,
    bins,
    momentum_bins,
    output_directory,
):
    """Plot numerator and denominator eta spectra in each momentum bin."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    bins = normalize_eta_bins(bins)
    eta_min, eta_max = bins[0], bins[-1]
    n_panels = max(1, len(momentum_bins) - 1)
    ncols = 2 if n_panels > 1 else 1
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        sharex=True,
        sharey=True,
        figsize=(3 * ncols, 3 * nrows),
    )
    axes = np.atleast_1d(axes).ravel()
    total_denominator = 0
    total_numerator = 0

    for panel, (mom_min, mom_max) in enumerate(
        zip(momentum_bins[:-1], momentum_bins[1:])
    ):
        ax = axes[panel]
        denominator = denominator_data[
            (denominator_data.mom >= mom_min)
            & (denominator_data.mom < mom_max)
            & (denominator_data.eta >= eta_min)
            & (denominator_data.eta <= eta_max)
        ]
        numerator = numerator_data[
            (numerator_data.mom >= mom_min)
            & (numerator_data.mom < mom_max)
            & (numerator_data.eta >= eta_min)
            & (numerator_data.eta <= eta_max)
        ]
        sns.histplot(denominator, x="eta", bins=bins, ax=ax)
        sns.histplot(numerator, x="eta", bins=bins, ax=ax)

        denominator_count = len(denominator)
        numerator_count = len(numerator)
        total_denominator += denominator_count
        total_numerator += numerator_count
        fraction = (
            numerator_count / denominator_count
            if denominator_count
            else np.nan
        )
        panel_label = "Eff" if metric_label == "Efficiency" else metric_label
        ax.text(
            0.15,
            0.8,
            f"{mom_min}<p<{mom_max} GeV\n{panel_label}: "
            f"{fraction:.3f}={numerator_count}/{denominator_count}",
            fontsize=10,
            transform=ax.transAxes,
            ha="left",
        )
        ax.set_xlim(eta_min, eta_max)
        ax.set_xlabel(r"$\eta$")

    for ax in axes[len(momentum_bins) - 1 :]:
        ax.set_visible(False)

    total_fraction = (
        total_numerator / total_denominator
        if total_denominator
        else np.nan
    )
    fig.suptitle(
        f"{metric_label} ({nhit_cut} hits) | {beam}, {setting} | "
        f"total={total_fraction:.3f} "
        f"({total_numerator}/{total_denominator})"
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    metric_tag = "eff" if metric_label == "Efficiency" else "purity"
    fig.savefig(
        output_directory / f"bg_{beam}_{setting}_{nhit_cut}_{metric_tag}.png"
    )
    plt.close(fig)
    return total_fraction, total_numerator, total_denominator


def plot_purity_distribution(
    valid_track_params,
    good_signal_track_params,
    beam,
    setting,
    nhit_cut,
    bins=None,
    mom_bin=None,
    outdir=None,
):
    """Plot good-signal tracks over all candidate tracks versus eta."""
    bins = np.arange(-4, 4.01, 0.1) if bins is None else bins
    momentum_bins = (
        DEFAULT_MOMENTUM_BINS if mom_bin is None else np.asarray(mom_bin)
    )
    return _plot_fraction_distribution(
        valid_track_params,
        good_signal_track_params,
        "Purity",
        beam,
        setting,
        nhit_cut,
        bins,
        momentum_bins,
        outdir,
    )


def plot_efficiency_distribution(
    particles,
    beam,
    setting,
    nhit_cut,
    bins=None,
    mom_bin=None,
    outdir=None,
):
    """Plot reconstructed signal particles over selected particles versus eta."""
    bins = np.arange(-4, 4.01, 0.1) if bins is None else bins
    momentum_bins = (
        DEFAULT_MOMENTUM_BINS if mom_bin is None else np.asarray(mom_bin)
    )
    particles_with_track = particles[particles["has_valid_track"]]
    return _plot_fraction_distribution(
        particles,
        particles_with_track,
        "Efficiency",
        beam,
        setting,
        nhit_cut,
        bins,
        momentum_bins,
        outdir,
    )


def plot_track_quality_metrics(
    track_class_bins,
    tracks,
    momentum_bins,
    analysis_label,
    output_directory,
):
    """Plot track-class fractions and the physical hit-count distribution."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import seaborn as sns

    class_colors = sns.color_palette(n_colors=len(TRACK_CLASS_ORDER))
    panel_count = len(momentum_bins) - 1
    ncols = 2 if panel_count > 1 else 1
    nrows = int(np.ceil(panel_count / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        sharex=True,
        sharey=True,
        figsize=(3.5 * ncols, 3 * nrows),
    )
    axes = np.atleast_1d(axes).ravel()
    for momentum_index, ax in enumerate(axes[:panel_count]):
        panel_data = track_class_bins[
            track_class_bins["momentum_bin"] == momentum_index
        ]
        for class_name in TRACK_CLASS_ORDER:
            class_data = panel_data[panel_data["metric"] == class_name]
            ax.plot(
                class_data["eta_center"],
                class_data["value"],
                marker=".",
                linewidth=1,
                label=TRACK_CLASS_LABELS[class_name],
            )
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(r"$\eta$")
        ax.set_ylabel("track fraction")
        ax.set_title(
            f"{momentum_bins[momentum_index]} < p < "
            f"{momentum_bins[momentum_index + 1]} GeV"
        )
    for ax in axes[panel_count:]:
        ax.set_visible(False)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Track classification | {analysis_label}")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(Path(output_directory) / "track_class_fractions.png")
    plt.close(fig)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        sharex=True,
        sharey=True,
        figsize=(3.5 * ncols, 3 * nrows),
    )
    axes = np.atleast_1d(axes).ravel()
    for momentum_index, ax in enumerate(axes[:panel_count]):
        momentum_low = momentum_bins[momentum_index]
        momentum_high = momentum_bins[momentum_index + 1]
        panel_tracks = tracks[
            (tracks["reco_mom"] >= momentum_low)
            & (tracks["reco_mom"] < momentum_high)
        ].copy()
        if panel_tracks.empty:
            ax.text(
                0.5,
                0.5,
                "no tracks",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_xlim(3, 15)
            ax.set_xticks(np.arange(3, 16, 2))
            ax.set_xlabel("hits per trajectory")
            ax.set_ylabel("tracks")
            ax.set_title(f"{momentum_low} < p < {momentum_high} GeV")
            continue
        panel_tracks["track_class_label"] = panel_tracks[
            "track_class"
        ].map(TRACK_CLASS_LABELS)
        sns.histplot(
            panel_tracks,
            x="total_count",
            hue="track_class_label",
            hue_order=[TRACK_CLASS_LABELS[name] for name in TRACK_CLASS_ORDER],
            palette=class_colors,
            discrete=True,
            element="step",
            fill=False,
            common_norm=False,
            legend=False,
            ax=ax,
        )
        ax.set_xlim(3, 15)
        ax.set_xticks(np.arange(3, 16, 2))
        ax.set_xlabel("hits per trajectory")
        ax.set_ylabel("tracks")
        ax.set_title(f"{momentum_low} < p < {momentum_high} GeV")

    for ax in axes[panel_count:]:
        ax.set_visible(False)
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=class_colors[index],
        label=TRACK_CLASS_LABELS[class_name],
        )
        for index, class_name in enumerate(TRACK_CLASS_ORDER)
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=len(TRACK_CLASS_ORDER),
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(
        f"Candidate-track hit counts | {analysis_label}",
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(Path(output_directory) / "hits_per_trajectory.png")
    plt.close(fig)


def aggregate_results(
    input_files,
    output_root,
    analysis_label,
    beam_label,
    cuts,
    momentum_bins,
    eta_bins,
    performance_tag=None,
):
    """Summarize manifests/chunks and write combined performance products."""
    eta_bins = normalize_eta_bins(eta_bins)
    summary_rows = []
    particle_frames = []
    track_frames = []

    for input_file in input_files:
        output_directory = get_output_directory(input_file, output_root)
        manifest_file = output_directory / "manifest.csv"
        if not manifest_file.is_file():
            raise FileNotFoundError(
                f"Missing manifest needed for aggregation: {manifest_file}"
            )
        manifest = pd.read_csv(manifest_file)
        trajectories, particles = load_active_chunks(manifest, manifest_file)
        valid_tracks, good_signal_tracks = select_valid_tracks(
            trajectories,
            cuts,
        )
        particles = mark_particles_with_tracks(particles, valid_tracks)

        n_particles = len(particles)
        n_particles_with_track = int(particles["has_valid_track"].sum())
        n_valid_tracks = len(valid_tracks)
        n_good_signal_tracks = len(good_signal_tracks)
        n_background_tracks = int(
            valid_tracks["track_class"].eq("background_track").sum()
        )
        n_fake_or_ghost_tracks = int(
            valid_tracks["track_class"].eq("fake_or_ghost").sum()
        )
        summary = {
            "input_file": str(input_file),
            "events_requested": int(len(manifest)),
            "events_ok": int((manifest["status"] == "ok").sum()),
            "events_timeout": int(
                (manifest["status"] == "timeout").sum()
            ),
            "events_error": int(
                manifest["status"].isin(["error", "worker_error"]).sum()
            ),
            "selected_particles": int(n_particles),
            "particles_with_valid_track": int(n_particles_with_track),
            "valid_tracks": int(n_valid_tracks),
            "good_signal_tracks": int(n_good_signal_tracks),
            "background_tracks": n_background_tracks,
            "fake_or_ghost_tracks": n_fake_or_ghost_tracks,
            "efficiency": (
                n_particles_with_track / n_particles
                if n_particles
                else np.nan
            ),
            "purity": (
                n_good_signal_tracks / n_valid_tracks
                if n_valid_tracks
                else np.nan
            ),
            "background_track_fraction": (
                n_background_tracks / n_valid_tracks
                if n_valid_tracks
                else np.nan
            ),
            "fake_or_ghost_fraction": (
                n_fake_or_ghost_tracks / n_valid_tracks
                if n_valid_tracks
                else np.nan
            ),
            "mean_hits_per_trajectory": (
                valid_tracks["total_count"].mean()
                if n_valid_tracks
                else np.nan
            ),
            "median_hits_per_trajectory": (
                valid_tracks["total_count"].median()
                if n_valid_tracks
                else np.nan
            ),
            "output_directory": str(output_directory),
        }
        summary_rows.append(summary)
        pd.DataFrame([summary]).to_csv(
            output_directory / "summary.csv",
            index=False,
        )

        if not particles.empty:
            particle_frames.append(
                particles.assign(source_file=str(input_file))
            )
        if not valid_tracks.empty:
            track_frames.append(
                valid_tracks.assign(source_file=str(input_file))
            )

    combined_particles = (
        pd.concat(particle_frames, ignore_index=True)
        if particle_frames
        else pd.DataFrame()
    )
    combined_valid_tracks = (
        pd.concat(track_frames, ignore_index=True)
        if track_frames
        else pd.DataFrame()
    )
    if combined_particles.empty or combined_valid_tracks.empty:
        raise RuntimeError(
            "Aggregation found no selected particles or valid tracks"
        )

    efficiency_bins = build_binned_performance(
        combined_particles,
        passed_column="has_valid_track",
        momentum_column="mom",
        eta_column="eta",
        metric_name="efficiency",
        momentum_bins=momentum_bins,
        eta_bins=eta_bins,
    )
    purity_bins = build_binned_performance(
        combined_valid_tracks,
        passed_column="is_good_signal",
        momentum_column="reco_mom",
        eta_column="reco_eta",
        metric_name="purity",
        momentum_bins=momentum_bins,
        eta_bins=eta_bins,
    )

    performance_directory = Path(output_root) / analysis_label / "performance"
    if performance_tag:
        performance_directory = performance_directory / performance_tag
    performance_directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(
        performance_directory / "file_summary.csv",
        index=False,
    )
    efficiency_bins.to_csv(
        performance_directory / "efficiency_bins.csv",
        index=False,
    )
    purity_bins.to_csv(
        performance_directory / "purity_bins.csv",
        index=False,
    )

    eta_min, eta_max = eta_bins[0], eta_bins[-1]
    momentum_min, momentum_max = momentum_bins[0], momentum_bins[-1]
    tracks_in_analysis_range = combined_valid_tracks[
        (combined_valid_tracks["reco_eta"] >= eta_min)
        & (combined_valid_tracks["reco_eta"] <= eta_max)
        & (combined_valid_tracks["reco_mom"] >= momentum_min)
        & (combined_valid_tracks["reco_mom"] < momentum_max)
    ].copy()
    track_columns = [
        "source_file",
        "event",
        "traj_id",
        "reco_mom",
        "reco_eta",
        "reco_pt",
        "total_count",
        "max_count",
        "max_fraction",
        "most_common_source",
        "part_status",
        "track_class",
    ]
    tracks_in_analysis_range[track_columns].to_csv(
        performance_directory / "hits_per_trajectory.csv",
        index=False,
    )

    track_class_bin_frames = []
    for class_name in TRACK_CLASS_ORDER:
        class_data = combined_valid_tracks.assign(
            in_track_class=combined_valid_tracks["track_class"].eq(class_name)
        )
        track_class_bin_frames.append(
            build_binned_performance(
                class_data,
                passed_column="in_track_class",
                momentum_column="reco_mom",
                eta_column="reco_eta",
                metric_name=class_name,
                momentum_bins=momentum_bins,
                eta_bins=eta_bins,
            )
        )
    track_class_bins = pd.concat(track_class_bin_frames, ignore_index=True)
    track_class_bins.to_csv(
        performance_directory / "track_class_bins.csv",
        index=False,
    )

    candidate_tracks = len(tracks_in_analysis_range)
    track_class_counts = tracks_in_analysis_range[
        "track_class"
    ].value_counts()
    track_class_summary = {
        "eta_min": eta_min,
        "eta_max": eta_max,
        "candidate_tracks": candidate_tracks,
        "good_signal_tracks": int(track_class_counts.get("good_signal", 0)),
        "background_tracks": int(track_class_counts.get("background_track", 0)),
        "fake_or_ghost_tracks": int(track_class_counts.get("fake_or_ghost", 0)),
        "mean_hits_per_trajectory": tracks_in_analysis_range[
            "total_count"
        ].mean(),
        "median_hits_per_trajectory": tracks_in_analysis_range[
            "total_count"
        ].median(),
        "hits_per_trajectory_q10": tracks_in_analysis_range[
            "total_count"
        ].quantile(0.10),
        "hits_per_trajectory_q90": tracks_in_analysis_range[
            "total_count"
        ].quantile(0.90),
    }
    for class_name, count_name in (
        ("good_signal", "purity"),
        ("background_track", "background_track_fraction"),
        ("fake_or_ghost", "fake_or_ghost_fraction"),
    ):
        count = int(track_class_counts.get(class_name, 0))
        track_class_summary[count_name] = (
            count / candidate_tracks if candidate_tracks else np.nan
        )
    pd.DataFrame([track_class_summary]).to_csv(
        performance_directory / "track_class_summary.csv",
        index=False,
    )
    plot_track_quality_metrics(
        track_class_bins,
        tracks_in_analysis_range,
        momentum_bins,
        analysis_label,
        performance_directory,
    )

    valid_track_params = combined_valid_tracks.rename(
        columns={
            "reco_mom": "mom",
            "reco_eta": "eta",
            "reco_pt": "pt",
        }
    )
    good_signal_track_params = valid_track_params[
        valid_track_params["is_good_signal"]
    ].copy()
    nhit_cut = cuts["track_hit_count_min"]
    (
        efficiency_total,
        efficiency_numerator,
        efficiency_denominator,
    ) = plot_efficiency_distribution(
        combined_particles,
        beam_label,
        analysis_label,
        nhit_cut,
        bins=eta_bins,
        mom_bin=momentum_bins,
        outdir=performance_directory,
    )
    (
        purity_total,
        purity_numerator,
        purity_denominator,
    ) = plot_purity_distribution(
        valid_track_params,
        good_signal_track_params,
        beam_label,
        analysis_label,
        nhit_cut,
        bins=eta_bins,
        mom_bin=momentum_bins,
        outdir=performance_directory,
    )

    performance_summary = pd.DataFrame(
        [
            {
                "metric": "efficiency",
                "numerator": efficiency_numerator,
                "denominator": efficiency_denominator,
                "value": efficiency_total,
            },
            {
                "metric": "purity",
                "numerator": purity_numerator,
                "denominator": purity_denominator,
                "value": purity_total,
            },
            {
                "metric": "background_track_fraction",
                "numerator": track_class_summary["background_tracks"],
                "denominator": candidate_tracks,
                "value": track_class_summary["background_track_fraction"],
            },
            {
                "metric": "fake_or_ghost_fraction",
                "numerator": track_class_summary["fake_or_ghost_tracks"],
                "denominator": candidate_tracks,
                "value": track_class_summary["fake_or_ghost_fraction"],
            },
        ]
    )
    performance_summary["analysis_label"] = analysis_label
    performance_summary["input_files"] = len(input_files)
    performance_summary.to_csv(
        performance_directory / "performance_summary.csv",
        index=False,
    )
    return performance_summary, performance_directory


def parse_arguments():
    """Parse and validate the command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Analyze reconstructed DIS files with one PODIO reader per file."
        )
    )
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Detector/campaign tag in numbered input filenames, e.g. "
            "no_L4_HD4_noise. Required unless --input-file is used."
        ),
    )
    parser.add_argument(
        "--input-file",
        action="append",
        type=Path,
        default=None,
        help=(
            "Exact reconstructed ROOT file to process or aggregate. Repeat "
            "for multiple files; bypasses --config and the numbered pattern."
        ),
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help="Directory containing one subdirectory per configuration.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Root directory for per-file chunks and combined performance.",
    )
    parser.add_argument("--file-start", type=int, default=1)
    parser.add_argument("--file-stop", type=int, default=100)
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--events", type=int, default=99)
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of input files processed concurrently.",
    )
    parser.add_argument("--timeout-seconds", type=int, default=120)
    parser.add_argument("--beam-label", default="10x275")
    parser.add_argument(
        "--analysis-label",
        default=None,
        help="Combined-output label; defaults to --config.",
    )
    parser.add_argument(
        "--eta-range",
        nargs=3,
        type=float,
        metavar=("MIN", "MAX", "STEP"),
        default=None,
        help="Eta range for combined tables and plots, e.g. -1 1 0.1.",
    )
    parser.add_argument(
        "--momentum-min",
        type=float,
        default=None,
        help=(
            "Minimum generated-particle momentum in GeV. The default comes "
            "from epic_analysis_base.TRACK_MOM_MIN."
        ),
    )
    parser.add_argument(
        "--momentum-bins",
        nargs="+",
        type=float,
        default=None,
        metavar="EDGE",
        help=(
            "Increasing momentum-bin edges in GeV. Defaults to the script's "
            "standard bins."
        ),
    )
    parser.add_argument(
        "--performance-tag",
        default=None,
        help=(
            "Optional subdirectory under performance/. With --eta-range, "
            "a non-overwriting eta tag is generated automatically."
        ),
    )
    parser.add_argument(
        "--module-directory",
        type=Path,
        default=DEFAULT_MODULE_DIRECTORY,
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse complete existing manifests instead of overwriting them.",
    )
    parser.add_argument(
        "--require-all-events",
        action="store_true",
        help="Exit with an error if any requested event times out or fails.",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--process-only",
        action="store_true",
        help="Generate per-file manifests/chunks without combined plots.",
    )
    mode.add_argument(
        "--aggregate-only",
        action="store_true",
        help="Reuse existing manifests/chunks and only rebuild summaries.",
    )
    args = parser.parse_args()

    if args.file_start < 1 or args.file_stop < args.file_start:
        parser.error("--file-start/--file-stop define an invalid range")
    if args.first_event < 0 or args.events < 1:
        parser.error("--first-event must be >=0 and --events must be >=1")
    if args.workers < 1:
        parser.error("--workers must be >=1")
    if args.momentum_min is not None and (
        not np.isfinite(args.momentum_min) or args.momentum_min < 0
    ):
        parser.error("--momentum-min must be finite and >=0")
    if args.momentum_bins is not None:
        momentum_bins = np.asarray(args.momentum_bins, dtype=float)
        if (
            len(momentum_bins) < 2
            or not np.all(np.isfinite(momentum_bins))
            or not np.all(np.diff(momentum_bins) > 0)
        ):
            parser.error("--momentum-bins must be finite increasing edges")
    if args.input_file is None and args.config is None:
        parser.error("provide --config or at least one --input-file")
    return args


def value_tag(value):
    """Return a filename-safe compact representation of one bin value."""
    return f"{value:g}".replace("-", "m").replace(".", "p")


def main():
    """Process reconstructed files and/or aggregate their saved chunks."""
    args = parse_arguments()
    cuts = dict(DEFAULT_CUTS)
    if args.momentum_min is not None:
        cuts["momentum_min"] = args.momentum_min
    momentum_bins = (
        np.asarray(args.momentum_bins, dtype=float)
        if args.momentum_bins is not None
        else DEFAULT_MOMENTUM_BINS
    )
    if args.input_file:
        input_files = [path.resolve() for path in args.input_file]
        missing_files = [path for path in input_files if not path.is_file()]
        if missing_files:
            raise FileNotFoundError(
                "Missing explicit input file(s):\n"
                + "\n".join(f"  {path}" for path in missing_files)
            )
        inferred_label = input_files[0].parent.name
    else:
        input_files = build_input_files(
            args.config,
            args.input_root,
            args.file_start,
            args.file_stop,
        )
        inferred_label = args.config

    analysis_label = args.analysis_label or args.config or inferred_label
    eta_bins = (
        normalize_eta_bins(args.eta_range)
        if args.eta_range is not None
        else DEFAULT_ETA_BINS
    )
    performance_tag = args.performance_tag
    if performance_tag is None and args.eta_range is not None:
        performance_tag = (
            f"eta_{value_tag(eta_bins[0])}_{value_tag(eta_bins[-1])}_"
            f"step_{value_tag(args.eta_range[2])}"
        )
    args.output_root.mkdir(parents=True, exist_ok=True)

    print(f"Configuration: {args.config or inferred_label}")
    print(f"Input files:   {len(input_files)}")
    if not args.input_file:
        print(f"File range:    {args.file_start} through {args.file_stop}")
    print(f"Event range:   {args.first_event} through "
          f"{args.first_event + args.events - 1}")
    print(f"Output root:   {args.output_root}")
    print(f"Cuts:          {cuts}")
    print(f"Momentum bins: {momentum_bins.tolist()} GeV")

    if not args.aggregate_only:
        run_results = process_input_files(
            input_files=input_files,
            output_root=args.output_root,
            module_directory=args.module_directory,
            first_event=args.first_event,
            number_of_events=args.events,
            timeout_seconds=args.timeout_seconds,
            cuts=cuts,
            number_of_workers=args.workers,
            resume=args.resume,
        )
        run_summary = pd.DataFrame(run_results)
        if len(input_files) == 1:
            run_summary_file = (
                get_output_directory(input_files[0], args.output_root)
                / "run_summary.csv"
            )
        else:
            run_summary_file = (
                args.output_root / f"{analysis_label}_run_summary.csv"
            )
        run_summary.to_csv(run_summary_file, index=False)
        print(f"Saved run summary: {run_summary_file}")

        worker_failures = run_summary["status"].eq("worker_error")
        if worker_failures.any():
            failed_files = run_summary.loc[
                worker_failures,
                "input_file",
            ].tolist()
            raise RuntimeError(
                f"{len(failed_files)} file workers failed; rerun with "
                "--resume after inspecting the per-file logs"
            )
        if args.require_all_events and run_summary["events_problem"].sum() > 0:
            problem_count = int(run_summary["events_problem"].sum())
            raise RuntimeError(
                f"{problem_count} requested event(s) failed or timed out; "
                "inspect the per-file manifests before aggregation"
            )

    if not args.process_only:
        performance_summary, performance_directory = aggregate_results(
            input_files=input_files,
            output_root=args.output_root,
            analysis_label=analysis_label,
            beam_label=args.beam_label,
            cuts=cuts,
            momentum_bins=momentum_bins,
            eta_bins=eta_bins,
            performance_tag=performance_tag,
        )
        print(f"Saved performance products: {performance_directory}")
        print(performance_summary.to_string(index=False))


if __name__ == "__main__":
    main()
