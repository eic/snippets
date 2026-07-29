import os
import re
from pathlib import Path

import numpy as np

import epic_analysis_base as ana


def _col_id(name):
    # Look up collection IDs in the podio metadata table.
    return next(k for k, v in ana.COL_TABLE.items() if v == name)


def parse_filename_metadata(path):
    """Parse filename for eta range, momentum, and tag/setting.

    Returns dict with keys: tag, setting, eta_lo, eta_hi, mom_gev.
    """
    fname = Path(path).name
    tag = fname
    if tag.endswith(".root"):
        tag = tag[:-5]

    setting = None
    if "rec_" in tag and "_eta" in tag:
        setting = tag.split("rec_", 1)[1].split("_eta", 1)[0]
    elif tag.startswith("rec_"):
        setting = tag.split("rec_", 1)[1].split("_", 1)[0]

    eta_lo = eta_hi = None
    eta_match = re.search(r"eta_(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)", tag)
    if eta_match:
        eta_lo = float(eta_match.group(1))
        eta_hi = float(eta_match.group(2))

    mom_gev = None
    mom_match = re.search(r"mom_([0-9]+(?:\.[0-9]+)?)GeV", tag)
    if mom_match:
        mom_gev = float(mom_match.group(1))
    else:
        gev_match = re.search(r"([0-9]+(?:\.[0-9]+)?)GeV", tag)
        if gev_match:
            mom_gev = float(gev_match.group(1))
        else:
            mev_match = re.search(r"([0-9]+(?:\.[0-9]+)?)MeV", tag)
            if mev_match:
                mom_gev = float(mev_match.group(1)) / 1000.0

    return {
        "tag": tag,
        "setting": setting,
        "eta_lo": eta_lo,
        "eta_hi": eta_hi,
        "mom_gev": mom_gev,
    }


def load_tree(fname, s3_dir="", entry_stop=None):
    # Update the global collection table for downstream lookups.
    ana.COL_TABLE = ana.get_col_table(fname, s3_dir=s3_dir)
    tree = ana.read_ur(fname, "events", s3_dir=s3_dir, entry_stop=entry_stop)
    return tree


def build_track_mc_df(
    tree,
    track_params="CentralCKFTrackParameters",
    assoc_name="CentralCKFTrackAssociations",
    track_collection="CentralCKFTracks",
):
    # Mirror the matching logic used in plot_single_resol.py.
    cov_branch = f"{track_params}/{track_params}.covariance.covariance[21]"

    params = ana.get_branch_df(tree, track_params)
    cov_df = ana.get_branch_df(tree, cov_branch)
    cov_wide = cov_df.pivot_table(index=["entry", "subentry"], columns="subsubentry", values="values")
    cov_wide.columns = [f"cov_{c}" for c in cov_wide.columns]
    cov_wide = cov_wide.reset_index()

    params = params.merge(cov_wide, on=["entry", "subentry"], how="left")
    params = params.rename(
        columns={
            "subentry": "track_index",
            "theta": "theta_rec",
            "phi": "phi_rec",
            "qOverP": "qOverP_rec",
        }
    )

    rec_assoc = ana.get_branch_df(tree, f"_{assoc_name}_rec")
    sim_assoc = ana.get_branch_df(tree, f"_{assoc_name}_sim")
    assoc = rec_assoc[["entry", "subentry", "index", "collectionID"]].rename(
        columns={"subentry": "assoc_index", "index": "track_index", "collectionID": "track_colID"}
    )
    assoc["mc_index"] = sim_assoc["index"]
    assoc["mc_colID"] = sim_assoc["collectionID"]

    track_col_id = _col_id(track_collection)
    mc_col_id = _col_id("MCParticles")
    assoc = assoc[(assoc["track_colID"] == track_col_id) & (assoc["mc_colID"] == mc_col_id)]

    df = assoc.merge(params, on=["entry", "track_index"], how="inner")

    mc = ana.get_branch_df(tree, "MCParticles")
    mc = mc.rename(columns={"subentry": "mc_index"})
    mc["p_mc"] = np.sqrt(mc["momentum.x"] ** 2 + mc["momentum.y"] ** 2 + mc["momentum.z"] ** 2)
    mc["theta_mc"] = np.arccos(mc["momentum.z"] / mc["p_mc"])
    mc["phi_mc"] = np.arctan2(mc["momentum.y"], mc["momentum.x"])
    mc["eta_mc"] = ana.theta2eta(mc["theta_mc"])
    mc["qOverP_true"] = mc["charge"] / mc["p_mc"]

    df = df.merge(
        mc[["entry", "mc_index", "theta_mc", "phi_mc", "eta_mc", "qOverP_true", "p_mc", "PDG", "generatorStatus"]],
        on=["entry", "mc_index"],
        how="inner",
    )

    return df, mc
