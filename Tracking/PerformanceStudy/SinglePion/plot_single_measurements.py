import os
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

_PDF_BASENAME = 'plot_single_measurements'
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

# inspect single particle hit quality: measurement, outlier, chi2, residual

## load modules

import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import awkward as ak
# Add repo root to sys.path
path_to_epic_analysis_base_py="/global/homes/s/shujie/eic_dir/worksim/scripts/modules/"
repo_root = Path(path_to_epic_analysis_base_py)
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

import epic_analysis_base as ana

## local file
# fname     = "/global/cfs/cdirs/m3763/shujie/worksim/matmap/"+"rec__mom_0.5GeV_Central_2601_athavan_map.root"
# s3_dir    = ''

## simulation campaign
mom_name="500MeV"
deg_name="3to50deg"
s3_dir = f"EPIC/RECO/25.12.0/epic_craterlake/SINGLE/pi+/{mom_name}/{deg_name}/"
fname    = f"pi+_{mom_name}_{deg_name}.0001.eicrecon.edm4eic.root"
tag=fname.split('/')[-1]
tag=tag.split('.root')[0]

ana.COL_TABLE = ana.get_col_table(fname,s3_dir,0)
tree          = ana.read_ur(fname,"events",s3_dir,entry_stop=100)

_fig_nums_before = set(plt.get_fignums())

parts = ana.get_part(tree)
fig,axs=plt.subplots(1,3, figsize=(12,3))
_=axs[0].hist(parts["eta"],bins=100)
axs[0].set_title("eta")
_=axs[1].hist(parts["phi"],bins=100)
axs[1].set_title("phi")
_=axs[2].hist(parts["mom"],bins=np.arange(0,20,0.1))
axs[2].set_title("momentum [GeV]")
## save to pdf
_save_new_figs(_fig_nums_before)

traj_meas_hits=ana.get_traj_hits_particle(tree)
traj_out_hits=ana.get_traj_hits_particle(tree,measurement_name="outliers_deprecated")

## append chi2 values
vname = "measurementChi2"
chi2=ana.get_relation(tree,"CentralCKFTrajectories",vname)
traj_meas_hits[vname]=np.hstack(chi2["values"])

vname = "outlierChi2"
chi2=ana.get_relation(tree,"CentralCKFTrajectories",vname)
traj_out_hits=traj_out_hits[traj_out_hits.nOutliers>0].reset_index()
traj_out_hits[vname]=np.hstack(chi2["values"])

# traj=ana.get_traj_relations(tree,l_var=["measurementChi2","outlierChi2"])

# chi2 (for barrel region)

_fig_nums_before = set(plt.get_fignums())
for nn in ["barrel"]:#, "endcap"]:
    plt.figure()
    if nn=="barrel":
        ax='x'
        ay='y'
    elif nn=="endcap":
        ax='z'
        ay='r'
    plt.xlabel(ax+' [mm]')
    plt.ylabel(ay+' [mm]')
    ax = 'position.'+ax
    ay = 'position.'+ay
    cond = traj_meas_hits["hits_colName"].str.contains(nn, case=False, na=False)
    df = traj_meas_hits[cond]
    plt.scatter(df[ax],df[ay],s=0.1,c='b',label="measurement hits")

    cond = traj_out_hits["hits_colName"].str.contains(nn, case=False, na=False)
    df = traj_out_hits[cond]
    plt.scatter(df[ax],df[ay],s=1,c='r',label="outlier hits")
    plt.legend()
    plt.title(nn+" region")
    ## save to pdf
    
_save_new_figs(_fig_nums_before)

_fig_nums_before = set(plt.get_fignums())
## Rec hits v.s. R
## ---------Barrel------------##
## number of svt hits per layer
### events have sim/rec hits from given layer 
plt.figure()
scounts=np.zeros(len(ana.barrel_range))
rcounts=np.zeros(len(ana.barrel_range))
# for ss,rr in zip(name_sim_barrel,name_rec_barrel):
sim_col = pd.Series(ana.name_sim_barrel).drop_duplicates().tolist()
rec_col = pd.Series(ana.name_rec_barrel).drop_duplicates().tolist()
for ss,rr in zip(sim_col, rec_col):
    sim_hits= ana.get_branch_df(tree,ss)
    sim_hits['position.r']=np.sqrt(sim_hits["position.x"]**2+sim_hits["position.y"]**2)
    r_sim = sim_hits["position.r"]
    rec_hits= ana.get_branch_df(tree,rr)
    rec_hits['position.r']=np.sqrt(rec_hits["position.x"]**2+rec_hits["position.y"]**2)
    r_rec = rec_hits["position.r"]
    
    for ii,(r1,r2) in enumerate(ana.barrel_range):
        n_sim=len(sim_hits[(r_sim>r1)&(r_sim<r2)]["entry"].unique())
        n_rec=len(rec_hits[(r_rec>r1)&(r_rec<r2)]["entry"].unique())
        if n_sim>0:
            print(ss, r1,r2, n_sim, n_rec)
            scounts[ii]+=n_sim            
            rcounts[ii]+=n_rec
# plt.scatter((r1+r2)/2, n_meas)

## plot rec hits v.s. measurement
bin_centers=[(a+b)/2 for (a,b) in ana.barrel_range]
_=plt.hist(bin_centers, bins=np.hstack(ana.barrel_range), weights=scounts, histtype="step", label="sim hits")
_=plt.hist(bin_centers, bins=np.hstack(ana.barrel_range), weights=rcounts, histtype="step",label="digi hits")


### events have measurements from given layer 
def check_ranges(arr, ranges):
    return [int(((arr >= lower) & (arr < upper))) for (lower, upper) in ranges]

# per row: list of 0/1 for each range
def get_range_counts(traj_hits,axis='r'):
    if axis=='r':
        my_range=ana.barrel_range
        my_var  ='position.r'
    elif axis=='z':
        my_range=ana.disk_range
        my_var  ='position.z'
    rbins = traj_hits[my_var].apply(lambda arr: check_ranges(arr, my_range))
    # build a DataFrame where each column is a range bin
    rbins_df = pd.DataFrame(rbins.tolist(), index=traj_hits.index)
    # group by entry/subentry and take max (any hit in range for that pair)
    rbins_unique = (
        rbins_df
        .groupby([traj_hits['entry'], traj_hits['subentry']])
        .max()
    )
    # total counts per range (one per unique entry/subentry)
    return rbins_unique.sum()

# sum the counts per range
# rbins = np.array(rbins.to_list()).sum(axis=0)# counts = range_check_matrix.sum(axis=0)
rbins = get_range_counts(traj_meas_hits, 'r')
_=plt.hist(bin_centers, bins=np.hstack(ana.barrel_range), weights=rbins,histtype="step",label="measurements")

rbins = get_range_counts(traj_out_hits, 'r')
_=plt.hist(bin_centers, bins=np.hstack(ana.barrel_range), weights=rbins,label="outliers")
plt.legend(frameon=0)

# _=plt.hist(np.hstack(traj_hits.measurements_r),bins=np.hstack(barrel_range))
plt.xlabel("R[mm]")
plt.xlim(0,900)
plt.title("Barrel hits in tracking")
plt.ylabel("# of events")
## save to pdf
_save_new_figs(_fig_nums_before)

_fig_nums_before = set(plt.get_fignums())
plt.figure()
bins=np.linspace(0,20,41)
df=traj_meas_hits
for nn in set(ana.name_rec_barrel):
    print(nn)
    cond=df["hits_colName"]==nn
    dfc = df[cond]
    if len(dfc)>0:
        _=plt.hist(np.hstack(dfc.measurementChi2),bins=bins,histtype="step",label=nn)
plt.legend()
plt.yscale('log')
plt.xlabel("chi2")
plt.ylabel("entries")
plt.title("measurement chi2")
## save to pdf
_save_new_figs(_fig_nums_before)

# residuals

seg = ana.get_branch_df(tree, '_CentralTrackSegments_points')
seg['position.r'] = np.sqrt(seg['position.x']**2 + seg['position.y']**2)
seg['position.phi'] = np.arctan2(seg['position.y'], seg['position.x'])

for lower, upper in ana.barrel_range:
    in_range = (seg['position.r'] >= lower) & (seg['position.r'] < upper)
    if not in_range.any():
        continue
    top_surfaces = seg.loc[in_range, 'surface'].value_counts().head(2).index ## drop apprach 1 and 2 surfaces
    seg = seg[~(in_range & seg['surface'].isin(top_surfaces))]

from scipy.spatial import cKDTree

def _match_closest_seg(traj_df, seg_df):
    if traj_df.empty or seg_df.empty:
        out = traj_df.copy()
        out['closest_seg_x'] = np.nan
        out['closest_seg_y'] = np.nan
        out['closest_seg_z'] = np.nan
        out['closest_seg_dist'] = np.nan
        return out

    seg_xyz = np.column_stack([seg_df['position.x'], seg_df['position.y'], seg_df['position.z']])
    traj_xyz = np.column_stack([traj_df['position.x'], traj_df['position.y'], traj_df['position.z']])
    ckdtree = cKDTree(seg_xyz)
    dist, idx = ckdtree.query(traj_xyz, k=1)

    matched = traj_df.copy()
    matched['closest_seg_x'] = seg_df.iloc[idx]['position.x'].to_numpy()
    matched['closest_seg_y'] = seg_df.iloc[idx]['position.y'].to_numpy()
    matched['closest_seg_z'] = seg_df.iloc[idx]['position.z'].to_numpy()
    matched['closest_seg_dist'] = dist
    return matched

traj_meas_hits = (
    traj_meas_hits
    .groupby('entry', group_keys=False)
    .apply(lambda g: _match_closest_seg(g, seg[seg['entry'] == g.name]))
)

traj_out_hits = (
    traj_out_hits
    .groupby('entry', group_keys=False)
    .apply(lambda g: _match_closest_seg(g, seg[seg['entry'] == g.name]))
)

# ientry = 1
# iseg = seg[seg['entry'] == ientry]
# itraj = traj_meas_hits[traj_meas_hits['entry'] == ientry]

# plt.scatter(itraj['position.x'], itraj['position.y'], s=10, label='traj')
# plt.scatter(iseg['position.x'], iseg['position.y'], s=3, label='seg')
# plt.scatter(itraj['closest_seg_x'], itraj['closest_seg_y'], s=6, label='closest seg')
# plt.legend(frameon=False)
# itraj

_fig_nums_before = set(plt.get_fignums())
bins = np.linspace(0, 1, 41)
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(7, 8), gridspec_kw={'hspace': 0.08})
for ii, rr in enumerate(ana.barrel_range):
    df = traj_meas_hits
    cond = (df['position.r'] > rr[0]) & (df['position.r'] < rr[1])
    counts_meas, _ = np.histogram(df[cond].closest_seg_dist, bins=bins)
    ax[0].hist(df[cond].closest_seg_dist, bins=bins, histtype='step', label=ana.barrel_name[ii])

    df = traj_out_hits
    cond = (df['position.r'] > rr[0]) & (df['position.r'] < rr[1])
    counts_out, _ = np.histogram(df[cond].closest_seg_dist, bins=bins)
    ax[1].hist(df[cond].closest_seg_dist, bins=bins, histtype='step', label=ana.barrel_name[ii])

    sum_counts = counts_meas + counts_out
    ax[2].hist(bins[:-1], bins=bins, weights=sum_counts, histtype='step', label=ana.barrel_name[ii])

ax[0].legend()

label_loc = (0.45, 0.92)
ax[0].text(*label_loc, 'measurements', transform=ax[0].transAxes, va='top')
ax[1].text(*label_loc, 'outliers', transform=ax[1].transAxes, va='top')
ax[2].text(*label_loc, 'all', transform=ax[2].transAxes, va='top')

ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[2].set_yscale('log')
ax[0].title('Hit residuals by category')
plt.tight_layout()
## save to pdf

# cond = (df['position.r'] > 400) & (df['position.r'] < 500)
# _=plt.hist(df[cond].closest_seg_dist,bins=100)
_save_new_figs(_fig_nums_before)


_pdf_pages.close()
