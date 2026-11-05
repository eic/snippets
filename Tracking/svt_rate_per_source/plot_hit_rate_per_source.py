#!/usr/bin/env python3
"""Plot SVT reconstructed-hit rates by MC source from chunked EDM4eic events.

Run with --help. All project imports are local to this folder; ROOT and PODIO
Python bindings are unnecessary. See README.md for normalization assumptions.
"""

import argparse
import csv
import hashlib
import json
import math
import re
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import awkward as ak
import matplotlib

matplotlib.use("Agg")  # Batch jobs and remote terminals do not need a display.
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.backends.backend_pdf import PdfPages

from svt_sources import (BACKGROUND_COMPONENTS, RELATIONS, SOURCE_ORDER,
                         STATUS_TO_SOURCE, SURFACES, assign_sources, read_branch)


HERE = Path(__file__).resolve().parent
SUMMARY_COLUMNS = ["lname", "source", "rec_branch", "raw_hits",
                   "scaled_total_rate_khz", "rsu_or_tile_max_rate_khz",
                   "raw_rsu_or_tile_max", "events_with_hits",
                   "max_raw_hits_per_event_tile", "outside_map_hits"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input-files", nargs="+", help="Local ROOT files or full root:// URLs; optional :events suffix.")
    inputs.add_argument("--file-list", type=Path, help="One input per line; blank lines and # comments are ignored.")
    parser.add_argument("--label", default="svt_rate", help="Name of a fresh output subdirectory (default: svt_rate).")
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs", help="Parent output directory (default: outputs beside this script).")
    parser.add_argument("--entry-start", type=int, default=0, help="First event per file, inclusive (default: 0).")
    parser.add_argument("--entry-stop", type=int, help="Last event per file, exclusive (default: end of file).")
    parser.add_argument("--max-files", type=int, help="Use only the first N files for a short check.")
    parser.add_argument("--chunk-size", type=int, default=10, help="Events read at once (default: 10); does not change normalization.")
    parser.add_argument("--frame-us", type=float, default=2.0, help="Physical duration of each equal-weight event frame in us (default: 2).")
    parser.add_argument("--tile", action="store_true", help="Use 3.5 x 9.8 mm bins; otherwise use 20 x 20 mm RSU bins.")
    parser.add_argument("--min-surface-hits", type=int, default=1, help="Minimum in-map source hits for a PDF page; CSV counts are unaffected.")
    args = parser.parse_args()
    if args.entry_start < 0 or (args.entry_stop is not None and args.entry_stop <= args.entry_start):
        parser.error("Require 0 <= entry-start < entry-stop (when supplied)")
    if args.chunk_size <= 0 or (args.max_files is not None and args.max_files <= 0):
        parser.error("chunk-size and max-files must be positive")
    if not math.isfinite(args.frame_us) or args.frame_us <= 0 or args.min_surface_hits < 1:
        parser.error("frame-us must be finite and positive; min-surface-hits must be >= 1")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.label):
        parser.error("label must start with a letter/digit and contain only letters, digits, _, -, or .")
    return args


def build_inputs(args):
    """Resolve local paths and preserve full remote URLs, in the supplied order."""
    if args.file_list:
        base = args.file_list.resolve().parent
        items = [line.strip() for line in args.file_list.read_text().splitlines()
                 if line.strip() and not line.lstrip().startswith("#")]
    else:
        base = Path.cwd()
        items = args.input_files
    if args.max_files:
        items = items[:args.max_files]
    paths = []
    for item in items:
        path = item[:-len(":events")] if item.endswith(":events") else item
        if not path.endswith(".root"):
            raise ValueError(f"Expected a ROOT file with optional :events suffix: {item}")
        if "://" not in path:
            path = str((base / Path(path).expanduser()).resolve())
        paths.append(path)
    if not paths:
        raise ValueError("The input list is empty")
    if len(set(paths)) != len(paths):
        raise ValueError("Repeated input files would double-count events; remove duplicate entries")
    return paths


def make_surface_maps(tile=False):
    """Build fixed grids; keep complete nominal-size bins at the positive edge."""
    dx, dy = (3.5, 9.8) if tile else (20.0, 20.0)
    maps = {}
    for name, branch, barrel, window in SURFACES:
        extent = (int(window[1] * np.pi // 100) + 1) * 100 if barrel else 500
        # The original workflow stopped one bin short; include the full extent.
        xedges = -extent + dx * np.arange(math.ceil(2 * extent / dx) + 1)
        yedges = -extent + dy * np.arange(math.ceil(2 * extent / dy) + 1)
        shape = (len(xedges) - 1, len(yedges) - 1)
        maps[name] = {"branch": branch, "barrel": barrel, "window": window,
                      "xedges": xedges, "yedges": yedges,
                      "counts": {"ALL": np.zeros(shape, dtype=np.int64)},
                      "counters": defaultdict(Counter)}
    return maps


def fill_event(surface_maps, branch, hits, sources):
    """Fill ALL, exclusive sources, and the inclusive other-background group."""
    if len(hits) != len(sources):
        raise ValueError(f"{branch}: rec-hit and raw-hit relation lengths disagree")
    x, y, z = (ak.to_numpy(hits[field]) for field in ("position.x", "position.y", "position.z"))
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y)) and np.all(np.isfinite(z))):
        raise ValueError(f"{branch}: nonfinite reconstructed-hit coordinates")
    radius = np.hypot(x, y)
    for surface in surface_maps.values():
        if surface["branch"] != branch:
            continue
        low, high = surface["window"]
        selector = radius if surface["barrel"] else z
        selected = (selector > low) & (selector < high)
        if not np.any(selected):
            continue
        px = z[selected] if surface["barrel"] else x[selected]
        py = (radius * np.arctan2(y, x))[selected] if surface["barrel"] else y[selected]
        labels = sources[selected]
        masks = {"ALL": np.ones(len(labels), dtype=bool)}
        masks.update({source: labels == source for source in np.unique(labels)})
        background = np.isin(labels, list(BACKGROUND_COMPONENTS))
        if np.any(background):
            masks["other background"] = background
        for source, mask in masks.items():
            counts = np.histogram2d(px[mask], py[mask], bins=(surface["xedges"], surface["yedges"]))[0].astype(np.int64)
            if source not in surface["counts"]:
                surface["counts"][source] = np.zeros_like(counts)
            surface["counts"][source] += counts
            counter = surface["counters"][source]
            in_map = int(counts.sum())
            counter["events_with_hits"] += int(in_map > 0)
            counter["max_raw_hits_per_event_tile"] = max(counter["max_raw_hits_per_event_tile"], int(counts.max()))
            counter["outside_map_hits"] += int(mask.sum()) - in_map


def process_files(paths, args):
    """Keep one event chunk and accumulated maps in memory; share MC truth."""
    maps = make_surface_maps(args.tile)
    diagnostics = {config["rec"]: Counter() for config in RELATIONS}
    inputs = []
    total_events = 0
    started = time.perf_counter()
    for file_id, path in enumerate(paths):
        print(f"[{file_id + 1}/{len(paths)}] {path}", flush=True)
        # Disable array caching: completed chunks need not remain in memory.
        with uproot.open(path, array_cache=None) as root_file:
            tree = root_file["events"]
            stop = tree.num_entries if args.entry_stop is None else min(args.entry_stop, tree.num_entries)
            start = min(args.entry_start, tree.num_entries)
            inputs.append({"file_id": file_id, "input_url": path + ":events",
                           "tree_events": tree.num_entries, "entry_start": start,
                           "entry_stop": max(start, stop), "processed_events": max(0, stop - start)})
            for chunk_start in range(start, stop, args.chunk_size):
                chunk_stop = min(chunk_start + args.chunk_size, stop)
                print(f"  events [{chunk_start}, {chunk_stop})", flush=True)
                mc = read_branch(tree, "MCParticles", chunk_start, chunk_stop,
                                 ["generatorStatus", "parents_begin", "parents_end"])
                parents = read_branch(tree, "_MCParticles_parents", chunk_start, chunk_stop, ["index"])
                for config in RELATIONS:
                    branch = config["rec"]
                    hits = read_branch(tree, branch, chunk_start, chunk_stop,
                                       ["position.x", "position.y", "position.z"])
                    relations = {key: read_branch(tree, config[key], chunk_start, chunk_stop,
                                                  ["weight"] if key == "association" else ["index"])
                                 for key in ("rec_raw", "association", "association_raw", "association_sim", "sim_particle")}
                    for event in range(chunk_stop - chunk_start):
                        labels, counts = assign_sources(
                            relations["rec_raw"][event], relations["association"][event]["weight"],
                            relations["association_raw"][event], relations["association_sim"][event],
                            relations["sim_particle"][event], mc[event], parents[event])
                        diagnostics[branch].update(counts)
                        fill_event(maps, branch, hits[event], labels)
                total_events += chunk_stop - chunk_start
    if total_events == 0:
        raise ValueError("No events were selected; check input sizes and the requested event range")
    # ALL must conserve hit counts bin-by-bin, including unresolved categories.
    for name, surface in maps.items():
        exclusive = np.zeros_like(surface["counts"]["ALL"])
        for source, counts in surface["counts"].items():
            if source not in ("ALL", "other background"):
                exclusive += counts
        if not np.array_equal(exclusive, surface["counts"]["ALL"]):
            raise RuntimeError(f"Source-count conservation failed on {name}")
    return maps, inputs, diagnostics, total_events, time.perf_counter() - started


def summarize(maps, scale):
    """Return per-surface/source counts and rates, without integer truncation."""
    rows = []
    for name, surface in maps.items():
        for source in SOURCE_ORDER:
            if source not in surface["counts"]:
                continue
            counts = surface["counts"][source]
            raw_hits, raw_peak = int(counts.sum()), int(counts.max())
            row = {"lname": name, "source": source, "rec_branch": surface["branch"],
                   "raw_hits": raw_hits, "scaled_total_rate_khz": raw_hits * scale,
                   "rsu_or_tile_max_rate_khz": raw_peak * scale, "raw_rsu_or_tile_max": raw_peak}
            for key in SUMMARY_COLUMNS[-3:]:
                row[key] = int(surface["counters"][source][key])
            rows.append(row)
    return rows


def write_csv(path, rows, columns, metadata):
    with path.open("w", newline="") as output:
        for key, value in metadata.items():
            output.write(f"# {key}: {value}\n")
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def plot_maps(maps, scale, args, output):
    """Plot already accumulated counts; no rebinning or integer rate rounding."""
    pdfs = []
    for source in SOURCE_ORDER:
        populated = [s for s in maps.items() if source in s[1]["counts"]
                     and s[1]["counts"][source].sum() >= args.min_surface_hits]
        if not populated:
            continue
        filename = f"{args.label}_{source.replace(' ', '_')}_svt_hit_rate.pdf"
        with PdfPages(output / filename) as pdf:
            for name, surface in populated:
                rates = surface["counts"][source] * scale
                fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
                image = ax.pcolormesh(surface["xedges"], surface["yedges"],
                                     np.ma.masked_where(rates.T == 0, rates.T), shading="flat")
                fig.colorbar(image, ax=ax, label="Hit rate per bin [kHz]")
                ax.set_aspect("equal")
                ax.set_xlabel("z [mm]" if surface["barrel"] else "x [mm]")
                ax.set_ylabel("rφ [mm]" if surface["barrel"] else "y [mm]")
                ax.set_title(f"{name} | {source} | {args.label}")
                bin_label = "tile" if args.tile else "RSU"
                ax.text(0.02, 0.98, f"Total: {rates.sum():.4g} kHz\nMaximum per {bin_label}: {rates.max():.4g} kHz",
                        transform=ax.transAxes, va="top", bbox={"facecolor": "white", "alpha": 0.85})
                pdf.savefig(fig)
                plt.close(fig)
        pdfs.append(filename)
    return pdfs


def save_maps(maps, output, label):
    """Archive raw histograms and mm bin edges with an explicit JSON key index."""
    arrays, index = {}, []
    for surface_id, (name, surface) in enumerate(maps.items()):
        prefix = f"surface_{surface_id:02d}"
        arrays[prefix + "_xedges_mm"] = surface["xedges"]
        arrays[prefix + "_yedges_mm"] = surface["yedges"]
        for source in SOURCE_ORDER:
            if source not in surface["counts"]:
                continue
            key = prefix + "_" + source.replace(" ", "_")
            arrays[key] = surface["counts"][source]
            index.append({"key": key, "lname": name, "source": source,
                          "xedges_key": prefix + "_xedges_mm", "yedges_key": prefix + "_yedges_mm"})
    np.savez_compressed(output / f"{label}_maps.npz", **arrays)
    return index


def main():
    args = parse_args()
    paths = build_inputs(args)
    output = args.output_dir / args.label
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output}. Use a new --label to retain earlier results.")
    # Create output only after processing succeeds; failed reads leave no partial tables.
    maps, inputs, diagnostics, events, elapsed = process_files(paths, args)
    output.mkdir(parents=True, exist_ok=True)
    scale = 1000.0 / (args.frame_us * events)
    metadata = {"label": args.label, "processed_events": events, "frame_us": args.frame_us,
                "bin_mode": "tile" if args.tile else "RSU", "rate_scale_khz_per_count": scale}
    write_csv(output / f"{args.label}_summary.csv", summarize(maps, scale), SUMMARY_COLUMNS, metadata)
    write_csv(output / f"{args.label}_input_files.csv", inputs, list(inputs[0]), metadata)
    map_index = save_maps(maps, output, args.label)
    pdfs = plot_maps(maps, scale, args, output)
    record = {
        **metadata, "created_utc": datetime.now(timezone.utc).isoformat(),
        "processing_seconds": elapsed, "input_files": inputs,
        "arguments": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "versions": {"numpy": np.__version__, "awkward": ak.__version__, "uproot": uproot.__version__, "matplotlib": matplotlib.__version__},
        "code_sha256": {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                        for name in ("plot_hit_rate_per_source.py", "svt_sources.py")},
        "status_to_source": STATUS_TO_SOURCE, "relations": RELATIONS,
        "surface_windows_mm": [{"lname": n, "rec_branch": b, "is_barrel": barrel, "range": window}
                               for n, b, barrel, window in SURFACES],
        "diagnostics": {branch: dict(counts) for branch, counts in diagnostics.items()},
        "maps": map_index, "pdfs": pdfs,
    }
    (output / f"{args.label}_run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Wrote {output}: {events} events, {scale:.6g} kHz/count, processing {elapsed:.2f} s", flush=True)


if __name__ == "__main__":
    main()
