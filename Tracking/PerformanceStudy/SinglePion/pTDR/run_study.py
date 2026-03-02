import argparse
import csv
import os
import sys
from datetime import datetime
from pathlib import Path
import glob

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

import epic_analysis_base as ana
from single_particle_performance.io import build_track_mc_df, load_tree, parse_filename_metadata


def _ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def _plot_efficiency(pdf, mc_primary, df_matched, tag, eta_range=None):
    # Diagnostic plot for a single file (not used by run_plot.py).
    bins = np.arange(-4, 4, 0.1)
    centers, eff, err = compute_binned_efficiency(mc_primary, df_matched, bins)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(7, 6), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    ax_top.hist(mc_primary["eta_mc"].to_numpy(), bins=bins, histtype="step", color="black", label="Gen")
    ax_top.hist(df_matched["eta_mc"].to_numpy(), bins=bins, histtype="step", color="tab:blue", label="Reco matched")
    ax_top.set_ylabel("Entries")
    ax_top.legend(frameon=False)
    ax_top.grid(True, alpha=0.3)

    ax_bot.errorbar(centers, eff, yerr=err, fmt="o", ms=3, lw=1, color="tab:blue", ecolor="gray", capsize=2)
    ax_bot.axhline(1.0, color="gray", lw=1, ls="--")
    ax_bot.set_xlabel(r"$\eta$")
    ax_bot.set_ylabel("Efficiency")
    ax_bot.set_xlim(bins[0], bins[-1])
    ax_bot.set_ylim(0, 1.05)
    ax_bot.grid(True, alpha=0.3)
    if eta_range:
        ax_bot.axvspan(eta_range[0], eta_range[1], color="tab:blue", alpha=0.1)

    fig.suptitle(tag)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def compute_binned_efficiency(mc_primary, df_matched, bins):
    # Return bin centers + efficiency/error arrays for plotting later.
    eta_gen = mc_primary["eta_mc"].to_numpy()
    eta_rec = df_matched["eta_mc"].to_numpy()

    gen_counts, bin_edges = np.histogram(eta_gen, bins=bins)
    rec_counts, _ = np.histogram(eta_rec, bins=bins)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    eff = np.divide(rec_counts, gen_counts, out=np.zeros_like(rec_counts, dtype=float), where=gen_counts > 0)
    err = np.zeros_like(eff)
    mask = gen_counts > 0
    err[mask] = np.sqrt(eff[mask] * (1.0 - eff[mask]) / gen_counts[mask])
    return centers, eff, err


def _plot_distributions(pdf, df, tag):
    plots = [
        ("pull_theta", "Pull distribution(theta)"),
        ("pull_phi", "Pull distribution(phi)"),
        ("pull_qoverp", "Pull distribution(q/p)"),
        ("resol_theta", "Resolution (theta [mrad])"),
        ("resol_phi", "Resolution (phi [mrad])"),
        ("resol_dp", "Resolution (dp/p [%])"),
    ]
    if "resol_dca" in df.columns:
        plots.append(("resol_dca", "Resolution (DCA$_r$ [mm])"))

    cols = 3
    rows = int(np.ceil(len(plots) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows))
    axes_flat = axes.flatten()

    for ax, (col, title) in zip(axes_flat, plots):
        mean, sigma, sigma_err = ana.hist_gaus(df[col], ax=ax, bins=101)
        label = title
        if np.isfinite(sigma) and np.isfinite(sigma_err):
            label = f"{title}\n$\\sigma$={sigma:.3g} ± {sigma_err:.3g}"
        ax.set_title(label)
        ax.set_ylabel("Entries")
        ax.grid(True, alpha=0.3)

    for ax in axes_flat[len(plots):]:
        ax.axis("off")

    fig.suptitle(tag)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _write_row(table_path, row):
    # Append to a single CSV to keep a growing summary across runs.
    columns = _ordered_columns(row)
    exists = os.path.exists(table_path)
    with open(table_path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        if not exists:
            writer.writeheader()
        writer.writerow(_format_row(row))


def _format_row(row, sig_digits=5):
    formatted = {}
    for key, value in row.items():
        if isinstance(value, float):
            formatted[key] = f"{value:.{sig_digits}g}"
        elif isinstance(value, (list, tuple)):
            formatted[key] = "[" + ", ".join(f"{vv:.{sig_digits}g}" for vv in value) + "]"
        else:
            formatted[key] = value
    return formatted


def _ordered_columns(row):
    columns = list(row.keys())
    if "timestamp" in columns:
        columns.remove("timestamp")
        columns.insert(1, "timestamp")
    if "source_file" in columns:
        columns.remove("source_file")
        columns.append("source_file")
    return columns


def select_by_eta_mom(df, eta_range=None, mom_gev=None, mom_tol=0.02):
    # Apply eta/momentum selection in MC truth space.
    mask = np.ones(len(df), dtype=bool)
    if eta_range is not None:
        eta_lo, eta_hi = eta_range
        mask &= df["eta_mc"].between(eta_lo, eta_hi, inclusive="left")
    if mom_gev is not None:
        lo = mom_gev * (1.0 - mom_tol)
        hi = mom_gev * (1.0 + mom_tol)
        mask &= (df["p_mc"] >= lo) & (df["p_mc"] <= hi)
    return df[mask]


def compute_metrics(df):
    # Compute pull/resolution quantities in the same units as plot_single_resol.py.
    valid = (df["cov_5"] > 0) & (df["cov_9"] > 0) & (df["cov_14"] > 0)
    valid &= np.isfinite(df["theta_rec"]) & np.isfinite(df["phi_rec"]) & np.isfinite(df["qOverP_rec"])
    valid &= np.isfinite(df["theta_mc"]) & np.isfinite(df["phi_mc"]) & np.isfinite(df["qOverP_true"])
    df = df[valid].copy()

    df["pull_theta"] = (df["theta_rec"] - df["theta_mc"]) / np.sqrt(df["cov_9"])
    df["pull_phi"] = (df["phi_rec"] - df["phi_mc"]) / np.sqrt(df["cov_5"])
    df["pull_qoverp"] = (np.abs(df["qOverP_rec"]) - np.abs(df["qOverP_true"])) / np.sqrt(df["cov_14"])

    df["resol_theta"] = (df["theta_rec"] - df["theta_mc"]) * 1000.0
    df["resol_phi"] = (df["phi_rec"] - df["phi_mc"]) * 1000.0
    df["p_rec"] = 1.0 / np.abs(df["qOverP_rec"])
    df["resol_dp"] = (df["p_rec"] - df["p_mc"]) / df["p_mc"] * 100.0

    if "loc.a" in df.columns:
        df["resol_dca"] = df["loc.a"]

    return df


def compute_efficiency(mc, df_matched, eta_range=None, mom_gev=None, mom_tol=0.02, pid=211):
    # Compute inclusive efficiency for the selected eta/momentum region.
    mc_primary = mc[mc["generatorStatus"] == 1].copy()
    if pid is not None and "PDG" in mc_primary.columns:
        mc_primary = mc_primary[mc_primary["PDG"] == pid]

    if eta_range is not None:
        eta_lo, eta_hi = eta_range
        mc_primary = mc_primary[mc_primary["eta_mc"].between(eta_lo, eta_hi, inclusive="left")]

    if mom_gev is not None:
        lo = mom_gev * (1.0 - mom_tol)
        hi = mom_gev * (1.0 + mom_tol)
        mc_primary = mc_primary[(mc_primary["p_mc"] >= lo) & (mc_primary["p_mc"] <= hi)]

    n_gen = len(mc_primary)
    matched = df_matched.copy()
    if eta_range is not None:
        eta_lo, eta_hi = eta_range
        matched = matched[matched["eta_mc"].between(eta_lo, eta_hi, inclusive="left")]
    if mom_gev is not None:
        lo = mom_gev * (1.0 - mom_tol)
        hi = mom_gev * (1.0 + mom_tol)
        matched = matched[(matched["p_mc"] >= lo) & (matched["p_mc"] <= hi)]

    matched = matched.drop_duplicates(subset=["entry", "mc_index"])
    n_rec = len(matched)

    if n_gen <= 0:
        return 0.0, 0.0, n_gen, n_rec

    eff = n_rec / n_gen
    eff_err = np.sqrt(eff * (1.0 - eff) / n_gen)
    return eff, eff_err, n_gen, n_rec


def summarize_metrics(df, columns, bins=101):
    # Fit gaussians for summary stats stored in the CSV.
    summary = {}
    for col in columns:
        if col not in df.columns:
            summary[col] = {"mean": np.nan, "sigma": np.nan, "sigma_err": np.nan}
            continue
        mean, sigma, sigma_err = ana.hist_gaus(df[col], ax=None, bins=bins)
        summary[col] = {"mean": mean, "sigma": sigma, "sigma_err": sigma_err}
    return summary


def run_single_file(
    fname,
    s3_dir="",
    eta_range=None,
    mom_gev=None,
    mom_tol=0.02,
    pid=211,
    track_params="CentralCKFTrackParameters",
    assoc_name="CentralCKFTrackAssociations",
    track_collection="CentralCKFTracks",
    out_dir=".",
    plots_dir="./plots",
    table_path="./performance_table.csv",
    entry_stop=None,
):
    # Enforce the rec_ naming convention so metadata parsing is consistent.
    if not os.path.basename(fname).startswith("rec_"):
        print(f"Skipping non-rec file: {fname}")
        return None
    meta = parse_filename_metadata(fname)
    tag = meta["tag"]

    if eta_range is None and meta["eta_lo"] is not None and meta["eta_hi"] is not None:
        eta_range = (meta["eta_lo"], meta["eta_hi"])
    if mom_gev is None and meta["mom_gev"] is not None:
        mom_gev = meta["mom_gev"]

    if eta_range is None or mom_gev is None:
        raise ValueError("eta_range and mom_gev must be provided or found in the filename tag")

    tree = load_tree(fname, s3_dir=s3_dir, entry_stop=entry_stop)
    df, mc = build_track_mc_df(
        tree,
        track_params=track_params,
        assoc_name=assoc_name,
        track_collection=track_collection,
    )

    df = compute_metrics(df)
    df_sel = select_by_eta_mom(df, eta_range=eta_range, mom_gev=mom_gev, mom_tol=mom_tol)

    mc_primary = mc[mc["generatorStatus"] == 1].copy()
    if pid is not None and "PDG" in mc_primary.columns:
        mc_primary = mc_primary[mc_primary["PDG"] == pid]
    if mom_gev is not None:
        lo = mom_gev * (1.0 - mom_tol)
        hi = mom_gev * (1.0 + mom_tol)
        mc_primary = mc_primary[(mc_primary["p_mc"] >= lo) & (mc_primary["p_mc"] <= hi)]

    eff, eff_err, n_gen, n_rec = compute_efficiency(
        mc,
        df,
        eta_range=eta_range,
        mom_gev=mom_gev,
        mom_tol=mom_tol,
        pid=pid,
    )

    # Bin-by-bin efficiency is stored for run_plot.py aggregation.
    eff_bins, eff_vals, eff_errs = compute_binned_efficiency(mc_primary, df_sel, np.arange(-4, 4, 0.1))

    metric_columns = [
        "pull_theta",
        "pull_phi",
        "pull_qoverp",
        "resol_theta",
        "resol_phi",
        "resol_dp",
        "resol_dca",
    ]
    summary = summarize_metrics(df_sel, metric_columns)

    row = {
        "tag": tag,
        "source_file": os.path.abspath(fname),
        "setting": meta.get("setting"),
        "eta_lo": eta_range[0],
        "eta_hi": eta_range[1],
        "mom_gev": mom_gev,
        "pid": pid,
        "mom_tol": mom_tol,
        "n_gen": n_gen,
        "n_rec": n_rec,
        "eff": eff,
        "eff_err": eff_err,
        "eff_bins": eff_bins.tolist(),
        "eff_values": eff_vals.tolist(),
        "eff_errors": eff_errs.tolist(),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }

    for col, stats in summary.items():
        row[f"{col}_mean"] = stats["mean"]
        row[f"{col}_sigma"] = stats["sigma"]
        row[f"{col}_sigma_err"] = stats["sigma_err"]

    _ensure_dir(out_dir)
    _ensure_dir(plots_dir)
    _write_row(table_path, row)

    pdf_path = os.path.join(plots_dir, f"{tag}_performance.pdf")
    with PdfPages(pdf_path) as pdf:
        _plot_efficiency(pdf, mc_primary, df_sel, tag, eta_range=eta_range)
        _plot_distributions(pdf, df_sel, tag)

    return row


def main():
    parser = argparse.ArgumentParser(description="Run single-particle performance study on ROOT files.")
    parser.add_argument("inputs", nargs="*", help="Input ROOT files")
    parser.add_argument("--pattern", default=None, help="Glob pattern for input ROOT files")
    parser.add_argument("--s3-dir", default="", help="Remote s3/xrootd directory")
    parser.add_argument("--eta", nargs=2, type=float, metavar=("ETA_LO", "ETA_HI"), help="Override eta range")
    parser.add_argument("--mom", type=float, help="Override momentum in GeV")
    parser.add_argument("--mom-tol", type=float, default=0.02, help="Momentum tolerance as fraction")
    parser.add_argument("--pid", type=int, default=211, help="PDG id (default: pi+ 211)")
    parser.add_argument("--track-params", default="CentralCKFTrackParameters")
    parser.add_argument("--assoc-name", default="CentralCKFTrackAssociations")
    parser.add_argument("--track-collection", default="CentralCKFTracks")
    parser.add_argument("--out-dir", default=".")
    parser.add_argument("--plots-dir", default="./plots")
    parser.add_argument("--table-path", default="./performance_table.csv")
    parser.add_argument("--entry-stop", type=int, default=None)
    parser.add_argument("--skip-existing", action="store_true", help="Skip files already in the CSV table")

    args = parser.parse_args()

    ana.configure_analysis_environment()

    inputs = list(args.inputs)
    if args.pattern:
        inputs.extend(sorted(glob.glob(args.pattern)))

    if not inputs:
        raise SystemExit("No input files provided. Use positional files or --pattern.")

    if args.skip_existing and os.path.exists(args.table_path):
        existing_tags = set()
        with open(args.table_path, newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                tag = row.get("tag")
                if tag:
                    existing_tags.add(tag)
        if existing_tags:
            inputs = [f for f in inputs if parse_filename_metadata(f)["tag"] not in existing_tags]

    if not inputs:
        print("No new files to process after applying --skip-existing.")
        return

    for fname in inputs:
        run_single_file(
            fname,
            s3_dir=args.s3_dir,
            eta_range=tuple(args.eta) if args.eta else None,
            mom_gev=args.mom,
            mom_tol=args.mom_tol,
            pid=args.pid,
            track_params=args.track_params,
            assoc_name=args.assoc_name,
            track_collection=args.track_collection,
            out_dir=args.out_dir,
            plots_dir=args.plots_dir,
            table_path=args.table_path,
            entry_stop=args.entry_stop,
        )


if __name__ == "__main__":
    main()
