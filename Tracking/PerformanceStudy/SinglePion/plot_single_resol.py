import os
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

_PDF_BASENAME = 'plot_single_resol'
_PDF_COVER_ADDED = False

def _cover_title_from_globals():
    for key in ('tag', 'fname'):
        if key in globals() and globals()[key] is not None:
            val = str(globals()[key])
            return os.path.basename(val)
    return _PDF_BASENAME

def _add_pdf_cover(fname, pdf_pages):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.6, 'source file: '+fname+'.root', ha='center', va='center', fontsize=16)
    fig.text(0.5, 0.5, datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
             ha='center', va='center', fontsize=12)
    pdf_pages.savefig(fig)
    plt.close(fig)

_PDF_PATH = f'{_PDF_BASENAME}_figures.pdf'
_pdf_pages = PdfPages(_PDF_PATH)

def _save_new_figs(figs_before):
    global _PDF_COVER_ADDED
    if not _PDF_COVER_ADDED:
        _add_pdf_cover(_cover_title_from_globals(), _pdf_pages)
        _PDF_COVER_ADDED = True
    for num in plt.get_fignums():
        if num not in figs_before:
            _pdf_pages.savefig(plt.figure(num))

# Minimal example for single particle analysis:
# * pull distributions for theta, phi, and q/p (from qOverP). Definitions follow Tracking_Performances.C (cov[5], cov[9], cov[14]).
# * resolutions (theta, phi, dp/p)


import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

repo_root = Path('/global/homes/s/shujie/eic_dir/worksim/scripts/modules')
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

import epic_analysis_base as ana
ana.configure_analysis_environment()

# ## load local file
# fname = "/global/cfs/cdirs/m3763/shujie/worksim/matmap/rec__mom_0.5GeV_Central_2601_full_map.root"
# ana.COL_TABLE = ana.get_col_table(fname)
# tree = ana.read_ur(fname, 'events')

## load simulation campaign file
deg_list=['3to50deg']#, '45to135deg','130to177deg']
mom_list=['500MeV', '1GeV', '2GeV', '5GeV', '10GeV', '20GeV'] ## omit 100 and 200 MeV
deg_name = deg_list[0]
mom_name = mom_list[0]
s3_path = f"EPIC/RECO/25.12.0/epic_craterlake/SINGLE/pi+/{mom_name}/{deg_name}/"
fname    = f"pi+_{mom_name}_{deg_name}.0001.eicrecon.edm4eic.root"

ana.COL_TABLE = ana.get_col_table(fname,s3_dir=s3_path)
tree = ana.read_ur(fname, 'events',s3_dir=s3_path)

tag=fname.split('/')[-1]
tag=tag.split('.root')[0]

track_params = 'CentralCKFTrackParameters'
assoc_name = 'CentralCKFTrackAssociations'
cov_branch = f'{track_params}/{track_params}.covariance.covariance[21]'

def _col_id(name):
    return next(k for k, v in ana.COL_TABLE.items() if v == name)

track_col_id = _col_id('CentralCKFTracks')
mc_col_id = _col_id('MCParticles')

params = ana.get_branch_df(tree, track_params)
cov_df = ana.get_branch_df(tree, cov_branch)
cov_wide = cov_df.pivot_table(index=['entry', 'subentry'], columns='subsubentry', values='values')
cov_wide.columns = [f'cov_{c}' for c in cov_wide.columns]
cov_wide = cov_wide.reset_index()

params = params.merge(cov_wide, on=['entry', 'subentry'], how='left')
params = params.rename(columns={'subentry': 'track_index', 'theta': 'theta_rec', 'phi': 'phi_rec', 'qOverP': 'qOverP_rec'})

rec_assoc = ana.get_branch_df(tree, f'_{assoc_name}_rec')
sim_assoc = ana.get_branch_df(tree, f'_{assoc_name}_sim')
assoc = rec_assoc[['entry', 'subentry', 'index', 'collectionID']].rename(
    columns={'subentry': 'assoc_index', 'index': 'track_index', 'collectionID': 'track_colID'}
)
assoc['mc_index'] = sim_assoc['index']
assoc['mc_colID'] = sim_assoc['collectionID']
assoc = assoc[(assoc['track_colID'] == track_col_id) & (assoc['mc_colID'] == mc_col_id)]

df = assoc.merge(params, on=['entry', 'track_index'], how='inner')

mc = ana.get_branch_df(tree, 'MCParticles')
mc = mc.rename(columns={'subentry': 'mc_index'})
mc['p_mc'] = np.sqrt(mc['momentum.x']**2 + mc['momentum.y']**2 + mc['momentum.z']**2)
mc['theta_mc'] = np.arccos(mc['momentum.z'] / mc['p_mc'])
mc['phi_mc'] = np.arctan2(mc['momentum.y'], mc['momentum.x'])
mc['qOverP_true'] = mc['charge'] / mc['p_mc']

df = df.merge(mc[['entry', 'mc_index', 'theta_mc', 'phi_mc', 'qOverP_true', 'p_mc']], on=['entry', 'mc_index'], how='inner')

valid = (df['cov_5'] > 0) & (df['cov_9'] > 0) & (df['cov_14'] > 0)
valid &= np.isfinite(df['theta_rec']) & np.isfinite(df['phi_rec']) & np.isfinite(df['qOverP_rec'])
valid &= np.isfinite(df['theta_mc']) & np.isfinite(df['phi_mc']) & np.isfinite(df['qOverP_true'])
df = df[valid]

df['pull_theta'] = (df['theta_rec'] - df['theta_mc']) / np.sqrt(df['cov_9'])
df['pull_phi'] = (df['phi_rec'] - df['phi_mc']) / np.sqrt(df['cov_5'])
df['pull_qoverp'] = (np.abs(df['qOverP_rec']) - np.abs(df['qOverP_true'])) / np.sqrt(df['cov_14'])



df['resol_theta'] = (df['theta_rec'] - df['theta_mc']) *1000
df['resol_phi']   = (df['phi_rec'] - df['phi_mc']) *1000
df['p_rec']     = (1./np.abs(df['qOverP_rec']))
df['resol_dp']   = (df['p_rec'] - df['p_mc'])/df['p_mc']*100

_fig_nums_before = set(plt.get_fignums())
bins = np.arange(-4, 4, 0.1)
eta_mc_in_rec = ana.theta2eta(df['theta_mc'])
eta_mc = ana.theta2eta(mc[mc.generatorStatus == 1]['theta_mc'])

ha_counts, bin_edges = np.histogram(eta_mc_in_rec, bins=bins)
hb_counts, _ = np.histogram(eta_mc, bins=bins)

bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
eff = np.divide(ha_counts, hb_counts, out=np.zeros_like(ha_counts, dtype=float), where=hb_counts > 0)
eff_err = np.zeros_like(eff)
mask = hb_counts > 0
eff_err[mask] = np.sqrt(eff[mask] * (1.0 - eff[mask]) / hb_counts[mask])

fig, (ax_top, ax_bot) = plt.subplots(2,1, figsize=(7, 6), sharex=True, gridspec_kw={'height_ratios': [3, 2]})
ax_top.hist(eta_mc, bins=bins, histtype='step', color='black', label='Gen')
ax_top.hist(eta_mc_in_rec, bins=bins, histtype='step', color='tab:blue', label='Reco matched')
ax_top.set_ylabel('Entries')
ax_top.legend(frameon=False)
ax_top.grid(True, alpha=0.3)

ax_bot.errorbar(bin_centers, eff, yerr=eff_err, fmt='o', ms=3, lw=1, color='tab:blue', ecolor='gray', capsize=2)
ax_bot.axhline(1.0, color='gray', lw=1, ls='--')
ax_bot.set_xlabel(r'$\eta$')
ax_bot.set_ylabel('Efficiency')
ax_bot.set_xlim(bin_edges[0], bin_edges[-1])
ax_bot.set_ylim(0, 1.05)
ax_bot.grid(True, alpha=0.3)
plt.tight_layout()
## save to pdf
_save_new_figs(_fig_nums_before)

_fig_nums_before = set(plt.get_fignums())
import math
plots = [
    ('pull_theta' , 'Pull distribution(theta)'),
    ('pull_phi'   , 'Pull distribution(phi)'),
    ('pull_qoverp', 'Pull distribution(q/p)'),
    ('resol_theta', 'Resolution (theta[mrad])'),
    ('resol_phi'  , 'Resolution (phi[mrad])'),
    ('resol_dp'   , 'Resolution (dp/p[%])'),
    ('loc.a'   , 'Resolution (DCA$_r$[mm])'),
]

cols = 3
rows = math.ceil(len(plots) / cols)
fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows))

def _label_from_title(title):
    ss=title.split("(")
    label1=ss[0]+":"
    label2 = ss[1].replace(')', '')
    return label1, label2

axes_flat = axes.flatten()
for ax, (col, title) in zip(axes_flat, plots):
    bins = 101
    mean, sigma, sigma_err = ana.hist_gaus(df[col], ax, bins=bins)
    label1,label2 = _label_from_title(title)
    if np.isfinite(sigma) and np.isfinite(sigma_err):
        label2 = f'{label2} = {sigma:.3g} +/- {sigma_err:.3g}'
    ax.text(
        0.05,
        0.9,
        label1+"\n"+label2,
        transform=ax.transAxes,
        ha='left',
        va='center',
        fontsize=12,
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none'),
    )
    ax.set_ylabel('Entries')
    ax.set_xlabel('')
    ax.grid(True, alpha=0.3)

for ax in axes_flat[len(plots):]:
    ax.axis('off')

plt.tight_layout()

## save to pdf
_save_new_figs(_fig_nums_before)


_pdf_pages.close()
