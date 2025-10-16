'''
    check tracking eff, resol for EIC-ePIC single particle simulation. 
    input: recon rootfiles with default name format: rec_{setting}_eta_{eta_list[ii]:g}_{eta_list[ii+1]:g}_{mom}GeV_{nev}.root
    Shujie Li, Sept 2024
'''

## this block of functions are copied from epic_analysis.ipynb 
from matplotlib.backends.backend_pdf import PdfPages
from lmfit.models import GaussianModel

from importlib import reload
import os
import sys
import types
import pandas as pd
import scipy
from scipy.signal import find_peaks
pd.set_option('display.max_rows', 500)
pd.options.display.max_rows = 40
pd.options.display.min_rows = 20
pd.options.display.max_columns = 100

import awkward as ak
# import ROOT
import uproot as ur
print('Uproot version: ' + ur.__version__) ## this script assumed version 4 
ur.default_library="pd" ## does not work???
from scipy import stats
import numpy as np
import argparse
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pylab as plt
plt.rcParams['figure.figsize'] = [8.0, 6.0]
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title]

deg2rad = np.pi/180.0
## convert theta to eta
def theta2eta(xx, inverse=0):
    xx = np.array(xx)
    if inverse==1:
        return np.arctan((np.e)**(-xx))*2
    else:
        return -np.log(np.tan(xx/2.))
        
## read a root tree with uproot. provide s3_dir to read from a remote server, xrootd required.
# dir = 'EPIC/RECO/23.11.0/epic_craterlake/DIS/NC/18x275/minQ2=10/'
# file = 'pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.tree.edm4eic.root'
def read_ur(fname, tname, s3_dir=""):
    if len(s3_dir)>1: # read from JLab server
        server = 'root://dtn-eic.jlab.org//work/eic2/'
        fname = server+s3_dir+fname
    tree     = ur.open(fname)[tname]
    print(f"read_ur: read {fname}:{tname}. {tree.num_entries} events in total")
    return tree

## read a branch from uproot tree with an option to flatten, return pandas dataframe
#  use kflatten=0 if want to get the nested dataframe
#  iev: event index, default=-1: get all events
def get_branch(tree,bname="",kflatten=1):
    if bname not in tree.keys():
        sys.exit("ERROR(get_branch): can't find branch "+bname)
    df    = tree[bname].array(library="ak")
    df    = ak.to_dataframe(df)
    if isinstance(df,pd.Series):
        return df #df.to_frame(name=bname.split("_")[-1])
    
    colls = df.columns
    if len(colls)<1:
        print('cannot find any leaf under branch {}'.format(bname))
        return pd.DataFrame()

    #remove prefix in name
    cols = [str(c) for c in colls if str(c).startswith('{}.'.format(bname))]
    #if it's an array member, the only column is "value"
    if not cols:
        return df
    # drop nested entry
    for cname in list(cols):
        if "[" in cname:
            cols.remove(cname) ## can not convert array to python. drop for now
        elif "covariance.covariance" in cname: ## for TrackParameters
            cols.remove(cname) 
    # rename and flat
    df = df[cols]
    df.rename(columns={c: c.replace(bname + '.', '') for c in df.columns}, inplace=True)
    if kflatten:
        df   = df.reset_index()
    # df.rename(columns={c: c[0].upper() + c[1:] for c in df.columns}, inplace=True)
    return  df

def hist_gaus(dataset, ax, bins=100, klog=0, header=None):
    ## select data in range if bins is provided as an array
    if not np.isscalar(bins):
        c1 = dataset <= bins[-1]
        c2 = dataset >= bins[0]
        dataset = dataset[c1 & c2]
    ## ----1st fit------
    n, bins, _ = ax.hist(dataset, bins, density=False, facecolor='b', alpha=0.3)
    xx    = bins[0:-1] + (bins[1] - bins[0]) / 2.0
    ymax  = np.max(n)
    std   = np.std(dataset)
    mean  = np.mean(dataset)

    c1 = xx <= (mean + 2 * std)
    c2 = xx >= (mean - 2 * std)
    cond  = c1 & c2

    ii = 0
    while len(n[cond]) < len(bins) / 2.0:
        # ax.cla()
        ax.clear()
        diff = (bins[-1] - bins[0]) / 2.0 / 2.0
        n, bins, _ = ax.hist(
            dataset,
            np.linspace(
                bins[0] + diff,
                bins[-1] - diff,
                len(bins)
            ),
            density=False,
            facecolor='b',
            alpha=0.3
        )
        xx = bins[0:-1] + (bins[1]-bins[0]) / 2.0
        c1 = xx <= (mean + 2 * std)
        c2 = xx >= (mean - 2 * std)
        cond = c1 & c2
        ii += 1
        if ii > 10:
            print("ERROR(hist_gaus): can not adjust the range.")
            return -1, -1, -1

    model = GaussianModel()
    # create parameters with initial guesses:
    params = model.make_params(center=np.median(xx[cond]), amplitude=np.max(n), sigma=np.std(xx[cond]))  
    result = model.fit(n, params, x=xx)
    
    # -----2nd fit--------
    std = result.params['sigma']
    # std = result.params['sigma'].get_value
    # print(mean,std)
    c1 = xx <= (mean + 2 * std)
    c2 = xx >= (mean - 2 * std)
    cond = c1 & c2
    if len(xx[cond]) < 10:
        print("Fit failed")
        return -1, -1, -1        
    model = GaussianModel()
    params = model.make_params(center=np.median(xx[cond]), amplitude=np.max(n[cond]), sigma=np.std(xx[cond]))  
    try: 
        result = model.fit(n[cond], params, x=xx[cond])
    except TypeError:
        print("Fit failed")
        return -1, -1,-1
    if result.params['sigma'].stderr ==  None:
        print("Fit failed")
        return -1, -1, -1
    #     print(result.fit_report())
        
    #     print (result.best_values)
    # plt.plot(xx, result.best_fit, 'r-', label='best fit')
    rv_mean = float(result.params['center'])
    rv_std = float(result.params['sigma'])
    rv_std_err = float(result.params['sigma'].stderr)
    rv_ampl = float(result.params['amplitude'])

    if len(result.best_fit) > 0:
        ax.plot(xx[cond], result.best_fit, 'r-', label='sigma=%g, err=%g' %(result.params['sigma'], result.params['sigma'].stderr))

    ax.legend(title=header, frameon=False, loc='upper left')

    ymax  = np.max(n)
    if klog:
        ax.set_yscale('log')
        ax.set_ylim(1, ymax * 10)
    else:
        ax.set_ylim(0, ymax * 1.3)
    return float(result.params['sigma']), float(result.params['sigma'].stderr)
    # return float(result.params['center']), float(result.params['sigma']), float(result.params['sigma'].stderr)


def hist_gaus(dataset,ax, bins=100,klog=0,header=None):
    ## select data in range if bins is provided as an array
    if not np.isscalar(bins):
        c1 = dataset<=bins[-1]
        c2 = dataset>=bins[0]
        dataset=dataset[c1&c2]
    ## ----1st fit------
    n, bins, patches = ax.hist(dataset, bins,density=False, facecolor='b', alpha=0.3)
    xx    = bins[0:-1]+(bins[1]-bins[0])/2.0
    ymax  = np.max(n)
    std   = np.std(dataset)
    mean  = np.mean(dataset)
    c1 = xx<=(mean+2*std)
    c2 = xx>=(mean-2*std)
    cond  = c1&c2

    ii=0
    while len(n[cond])<(len(bins)/2.0):
        # ax.cla()
        ax.clear()
        diff = (bins[-1]-bins[0])/2.0/2.0
        n, bins, patches = ax.hist(dataset, np.linspace(bins[0]+diff,bins[-1]-diff,len(bins)),density=False, facecolor='b', alpha=0.3)
        xx    = bins[0:-1]+(bins[1]-bins[0])/2.0
        c1 = xx<=(mean+2*std)
        c2 = xx>=(mean-2*std)
        cond  = c1&c2
        ii+=1
        if ii>5:
            print("ERROR(hist_gaus): can not adjust the range.")
            return -1,-1

    model = GaussianModel()
    # create parameters with initial guesses:
    params = model.make_params(center=np.median(xx[cond]), amplitude=np.max(n), sigma=np.std(xx[cond]))  
    result = model.fit(n, params, x=xx)
    
    # -----2nd fit--------
    std = result.params['sigma'].value
    # print(mean,std)
    c1 = xx<=(mean+2*std)
    c2 = xx>=(mean-2*std)
    cond = c1&c2
    if len(xx[cond])<10:
        print("Fit failed")
        return -1, -1        
    model = GaussianModel()
    params = model.make_params(center=np.median(xx[cond]), amplitude=np.max(n[cond]), sigma=np.std(xx[cond]))  
    try: 
        result = model.fit(n[cond], params, x=xx[cond])
    except TypeError:
        print("Fit failed")
        return -1,-1
    if result.params['sigma'].stderr ==  None:
        print("Fit failed")
        return -1, -1
    #     print(result.fit_report())
        
    #     print (result.best_values)
    # plt.plot(xx, result.best_fit, 'r-', label='best fit')
    if len(result.best_fit)>0:
        ax.plot(xx[cond], result.best_fit, 'r-', label='sigma=%g,err=%g' %(result.params['sigma'].value,result.params['sigma'].stderr))
    ax.legend(title=header, frameon=False,loc='upper left')

    ymax  = np.max(n)
    if klog:
        ax.set_yscale('log')
        ax.set_ylim(1,ymax*10)
    else:
        ax.set_ylim(0,ymax*1.3)
    return float(result.params['sigma'].value),float(result.params['sigma'].stderr)


def plot_eff(pion_o, pion,eta_bins=np.linspace(-4, 4, 21)):
    fig, ax = plt.subplots(1,1,figsize=[6,6])
    plt.title("")
    ## eff
    # original eta of all particle
    sim_eta, _ = np.histogram(pion_o['eta'].values, bins=eta_bins)
    # original eta of particles get reconstruted
    rec_eta, _ = np.histogram(pion['eta'], bins=eta_bins)
    track_eff_total = np.sum(rec_eta)/np.sum(sim_eta)

    eta_centers = (eta_bins[1:] + eta_bins[:-1])/2.
    eta_binsize = np.mean(np.diff(eta_centers))
    track_eff = np.nan_to_num(np.array(rec_eta)/np.array(sim_eta))
    
    # binary distribution, pq*sqrt(N)
    # TODO check the errors
    # eff = np.mean(track_eff)
    track_err = np.nan_to_num(track_eff*(1. - track_eff)*np.reciprocal(np.sqrt(sim_eta)))
    # rec_err = eff*(1. - eff)*np.sqrt(rec_eta)
    # track_eff_lower = track_eff - np.maximum(np.zeros(shape=rec_eta.shape), (rec_eta - rec_err)/sim_eta)
    # track_eff_upper = np.minimum(np.ones(shape=rec_eta.shape), (rec_eta + rec_err)/sim_eta) - track_eff
    track_eff_lower = track_eff - np.maximum(np.zeros(shape=rec_eta.shape), track_eff - track_err)
    track_eff_upper = np.minimum(np.ones(shape=rec_eta.shape), track_eff + track_err) - track_eff

    ax.errorbar(eta_centers, track_eff, xerr=eta_binsize/2., yerr=[track_eff_lower, track_eff_upper],
                fmt='o', capsize=3)
    ax.set_ylim(0., 1.1)
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylabel('Tracking Efficiency')#, fontsize=20)
    ax.set_xlabel(r'$\eta$')#, fontsize=20)
    ax.text(-4,1.04,"recon/generated events= %d / %d =%.2f" %(len(pion),len(pion_o),len(pion)/len(pion_o)))
    ax.axhline(1,ls='--',color='grey')
    return track_eff,track_err, eta_centers, fig


def plot_resol(pion,params):
    fig, axs = plt.subplots(2,2, figsize=(10,6),dpi=300)
    plt.title("")

    ## calculate resolutions
    dp_lim=5*2 #%
    th_lim=0.005*2 #rad
    ph_lim=0.03*2
    dca_lim = 3*2
    nbins  = 200

    # ----------------momentum resolution----------------------
    i  = 0
    ax = axs.flat[i]
    # obs    = 'theta'
    rec   = 1./np.array(params['qOverP'])
    sim   = np.array(pion['mom'])
    delta = (rec - sim)/sim*100 # in %
    sig_mom,err_mom=hist_gaus(delta,ax,np.linspace(-dp_lim, dp_lim, nbins+1),klog=0,header=None)
    ax.set_xlabel(r'$\delta p/p$ [%]')#, fontsize=20)

    # ----------------theta resolution----------------------
    i+=1
    ax = axs.flat[i]
    obs    = 'theta'
    rec    = np.array(params[obs])
    sim    = np.array(pion[obs])
    delta    = (rec - sim)
    sig_th,err_th=hist_gaus(delta,ax,np.linspace(-th_lim, th_lim, nbins+1),klog=0,header=None)
    ax.set_xlabel(r'$d\theta$ [rad]')#, fontsize=20)

    # ----------------phi resolution----------------------
    i+=1
    ax = axs.flat[i]
    obs   = 'phi'
    rec   = np.array(params[obs])
    sim   = np.array(pion[obs])
    delta = (rec - sim)
    sig_ph,err_ph=hist_gaus(delta,ax,np.linspace(-ph_lim, ph_lim, nbins+1),klog=0,header=None)
    ax.set_xlabel(r'$d\phi$ [rad]')#, fontsize=20)

    # ----------------theta resolution----------------------
    i+=1
    ax = axs.flat[i]
    obs    = 'loc.a'
    rec    = np.array(params[obs])
    sim    = 0#np.array(pion[obs])
    delta    = (rec - sim)
    sig_dca,err_dca=hist_gaus(delta,ax,np.linspace(-dca_lim, dca_lim, nbins+1),klog=0,header=None)
    ax.set_xlabel(r'DCA$_r$ [mm]')#, fontsize=20)
    return sig_mom, err_mom, sig_th, err_th, sig_ph, err_ph, sig_dca, err_dca, fig


## this only works for single particle simulation, no hit-based track-particle matching
## set eff_eta_bins to [] to disable eff plots. same for resol
## degrees of 3,10,40,140,170,177 correspond to eta bin of -3.5,-2.5,-1,1,2.5,3.5
def performance_plot(fname,rootfile_path,eff_eta_bins=np.linspace(-4,4,21),resol_deg_bins=np.array([3,10,40,140,170,177]),pid=211, out_dir="."):
    out_dir=out_dir+"/"
    ## read events tree
    tag     = fname.split('/')[-1][:-5]
    tree    = read_ur(fname,"events",rootfile_path)

    bname   = "CentralCKFTrackParameters"
    params  = get_branch(tree,bname)
    ## to deal with the nested subsubentry created by the covariance column
    ## also keep only the first track if more than one are reconstructed
    params  = params[(params.subentry==0)&(params.subsubentry==0)].reset_index() 
    params["eta"] = theta2eta(params.theta)
    params  = params.drop(['index','subsubentry'],axis=1)

    bname   = "MCParticles"
    part    = get_branch(tree,bname)
    # primary particle
    cond1   = part.generatorStatus==1
    # pid (defaut=pi+)
    cond2   = part.PDG==pid

    pion    = part[cond1&cond2].reset_index()
    x,y,z   = pion[["momentum.x","momentum.y","momentum.z"]].to_numpy().T
    r       = np.sqrt(x**2 + y**2 + z**2)  # Magnitude of the vector (distance to origin)
    pion["theta"]= np.arccos(z / r) 
    pion["phi"]  = np.arctan2(y, x)
    pion["eta"]  = theta2eta(pion.theta)
    pion["mom"]  = r
    # select particles that get reconstructed
    pion_o = pion #save all generaged pions
    cond = pion.entry.isin(params.entry)
    pion = pion.drop('index',axis=1)
    pion = pion[cond].reset_index()
    # now pion and params should be one-to-one

    ## eff plot
    if len(eff_eta_bins)>0:
        dump=plot_eff(pion_o,pion,eff_eta_bins)
        dump[-1].axes[0].set_title(fname)
        dump[-1].savefig(f'{out_dir}/eff_{tag}.png')
        plt.close()
        formatted_string = f"{tag};{dump[0].tolist()};{dump[1].tolist()};{dump[2].tolist()}"
        with open(f'{out_dir}eff_out.txt', 'a') as eff_file:
            eff_file.write(formatted_string + '\n')

    ## resolutions    
    if len(resol_deg_bins)>0:
        ## use this for Joe's rootfiles which is generated with eta_bin=0.5
        if len(resol_deg_bins)==1:
            dump=plot_resol(pion, params)
            dump[-1].axes[0].set_title(f"{tag}")
            dump[-1].savefig(f'../output/resol_{tag}.png')
            plt.close()
            temp = list(dump[0:-1])
            temp = ' '.join(map(str,temp))
            formatted_string = f"{tag} {temp}"
            with open(f'{out_dir}resol_out_whole.txt', 'a') as resol_file:
                resol_file.write(formatted_string + '\n')


        ## need to make slices of eta/theta for simulation campaign data
        else:
            for dd in np.arange(len(resol_deg_bins)-1): 
                deg_lo = resol_deg_bins[dd]
                deg_hi = resol_deg_bins[dd+1]
                cond1  = (pion.theta/deg2rad)>deg_lo
                cond2  = (pion.theta/deg2rad)<=deg_hi
                cond   = cond1&cond2
                pion_slice   = pion[cond].reset_index()
                params_slice = params[cond].reset_index()
                ## only proceed with good stats
                if len(pion_slice)>100:
                    dump=plot_resol(pion_slice, params_slice)
                    dump[-1].axes[0].set_title(f"{deg_lo} to {deg_hi} in "+tag)
                    dump[-1].savefig(f'{out_dir}resol_{tag}_theta_{deg_lo}_{deg_hi}.png')
                    plt.close()
                    temp = list(dump[0:-1])
                    temp = ' '.join(map(str,temp))
                    formatted_string = f"{tag} {deg_lo} {deg_hi} {temp}"
                    with open(f'{out_dir}resol_out_slices.txt', 'a') as resol_file:
                        resol_file.write(formatted_string + '\n')


if __name__ == "__main__":
    eff_eta_bins   = np.linspace(-4,4,41)
    # resol_deg_bins =[0] 
    resol_deg_bins = np.array([3,10,40,140,170,177])
    mom_list = [0.5, 1, 2, 5, 10, 15]
    eta_list = [-3.5,-2.5,-1,1,2.5,3.5]
    setting="default"
    nev=5000
    dirname = ""
    for mom in mom_list:
        for ii in np.arange(len(eta_list)-1):
            filename = dirname+f"rec_{setting}_eta_{eta_list[ii]:g}_{eta_list[ii+1]:g}_{mom}GeV_{nev}.root"
            print(filename)
            performance_plot(filename,"",eff_eta_bins=eff_eta_bins, resol_deg_bins=resol_deg_bins)
            break 
