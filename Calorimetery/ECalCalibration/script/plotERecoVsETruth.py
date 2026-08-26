"""Reproduce the single uncalibrated CALOROC energy-response plot."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from script.common import (
    LOAD_CHUNK_SIZE,
    SliceResult,
    applyPlotStyle,
    boundsCrystalBallPositive,
    buildSliceHistogramFromChunks,
    crystal_ball,
    curveFitReducer,
    loadData,
    maskPositiveAboveThreshold,
    nERecoBins,
    nETruthBins,
    p0CrystalBall,
    plotThetaDistribution,
    reduceSlices,
)


applyPlotStyle()


def linear(x, a, b):
    return a * x + b


def main(chunk_iter, output, title, e_truth_min, e_truth_max, e_reco_min,
         e_reco_max, txt=None, ytitle=''):
    """Accumulate response histograms, fit each truth-energy slice, and plot."""
    e_truth_bins = np.linspace(e_truth_min, e_truth_max, nETruthBins)
    e_reco_bins = np.linspace(e_reco_min, e_reco_max, nERecoBins)
    theta_bins = np.linspace(0, 180, 180)

    hist = np.zeros((len(e_truth_bins) - 1, len(e_reco_bins) - 1))
    theta_counts = np.zeros(len(theta_bins) - 1)
    for chunk in chunk_iter:
        h, _, _ = np.histogram2d(
            chunk['ETruth'], chunk['EReco'], bins=(e_truth_bins, e_reco_bins)
        )
        hist += h
        tc, _ = np.histogram(chunk['theta'], bins=theta_bins)
        theta_counts += tc

    err = hist / np.sqrt(hist)
    xc = 0.5 * (e_truth_bins[1:] + e_truth_bins[:-1])
    yc = 0.5 * (e_reco_bins[1:] + e_reco_bins[:-1])
    mean, _, std, _ = reduceSlices(
        xc, yc, hist, err,
        reducer=curveFitReducer(
            crystal_ball, p0CrystalBall, boundsCrystalBallPositive
        ),
        sliceMask=maskPositiveAboveThreshold(0.05, 0.4),
    )

    x, y = np.meshgrid(e_truth_bins, e_reco_bins)
    fig, (ax, ax2) = plt.subplots(ncols=2, figsize=(14, 7))
    pc = ax.pcolormesh(x, y, hist.T, cmap='magma_r')
    ax.errorbar(xc, mean, yerr=np.abs(std), marker='o')
    try:
        finite = np.isfinite(mean) & np.isfinite(std)
        popt, _ = curve_fit(
            linear, xc[finite], mean[finite], sigma=np.abs(std[finite]),
            absolute_sigma=True, p0=[1, 0]
        )
        fine_bins = np.linspace(e_truth_bins[0], e_truth_bins[-1], 100)
        ax.plot(
            fine_bins, linear(fine_bins, *popt),
            label='y = %.3e x + %.3e' % (popt[0], popt[1])
        )
    except Exception:
        pass
    ax.legend(frameon=False)

    if txt is not None:
        np.savetxt(txt, np.column_stack((xc, mean, np.abs(std))), fmt='%.6f')
    ax.set_ylim(e_reco_min, e_reco_max)
    try:
        fig.colorbar(pc, ax=ax)
    except Exception:
        pass
    ax.set_xlabel(r'$E_{truth}$ (GeV)')
    ax.set_ylabel(ytitle)
    plotThetaDistribution(ax2, theta_counts, theta_bins)
    fig.suptitle(title)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output')
    parser.add_argument('data')
    parser.add_argument('pdg', type=int)
    parser.add_argument('title')
    parser.add_argument('--eTruthMin', type=float, required=True)
    parser.add_argument('--eTruthMax', type=float, required=True)
    parser.add_argument('--eRecoMin', type=float, required=True)
    parser.add_argument('--eRecoMax', type=float, required=True)
    parser.add_argument('--ytitle')
    parser.add_argument('--txt', default=None)
    args = parser.parse_args()
    columns = [
        'mcId', 'PDG', 'theta', 'EReco', 'ETruth', 'PhiReco', 'EtaReco',
        'PhiTruth', 'EtaTruth', 'status', 'nevent'
    ]
    main(
        loadData(args.data, args.pdg, columns), args.output, args.title,
        args.eTruthMin, args.eTruthMax, args.eRecoMin, args.eRecoMax,
        args.txt, args.ytitle
    )
