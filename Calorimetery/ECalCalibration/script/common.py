"""Shared helpers used by the standalone energy-response reproduction."""

import os
import sys
from dataclasses import dataclass
from typing import Callable, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

LOAD_CHUNK_SIZE = 1_000_000
nETruthBins = 30
nERecoBins = 150


def applyPlotStyle():
    """Apply 2× scaling to matplotlib font and line sizes for readability."""
    plt.rcParams.update({
        'axes.titlesize': 24, 'axes.labelsize': 22,
        'xtick.labelsize': 20, 'ytick.labelsize': 20,
        'legend.fontsize': 16, 'figure.titlesize': 18,
        'lines.linewidth': 1.3, 'lines.markersize': 7.8,
    })


def crystal_ball(x, A, mu, sigma, alpha, n):
    """Evaluate the un-normalized Crystal Ball model."""
    x = np.array(x)
    t = (x - mu) / sigma
    abs_alpha = np.abs(alpha)
    a = (n / abs_alpha) ** n * np.exp(-0.5 * abs_alpha ** 2)
    b = n / abs_alpha - abs_alpha
    result = np.zeros_like(t)
    gaussian = t < abs_alpha
    result[gaussian] = np.exp(-0.5 * t[gaussian] ** 2)
    result[~gaussian] = a * (b + t[~gaussian]) ** (-n)
    return A * result


def loadData(path, pdg, columns, chunksize=LOAD_CHUNK_SIZE):
    """Yield filtered DataFrame chunks from a space-separated data file.

    Each chunk is filtered to status >= 0 and the requested PDG value; empty
    chunks are skipped so downstream histogram accumulation can stream safely.
    """
    reader = pd.read_csv(path, sep=' ', names=columns, chunksize=chunksize)
    for raw_chunk in reader:
        chunk = raw_chunk[(raw_chunk['status'] >= 0) & (raw_chunk['PDG'] == pdg)]
        if len(chunk):
            yield chunk


def buildSliceHistogramFromChunks(chunkIter, xcol, ycol, xbins, ybins,
                                  xTransform=None, yTransform=None):
    """Stream-accumulate a 2D histogram from DataFrame chunks.

    Optional transforms replace the named columns and are useful for derived
    quantities. Returns centers, histogram, and Poisson-style errors.
    """
    hist = np.zeros((len(xbins) - 1, len(ybins) - 1))
    for chunk in chunkIter:
        x = xTransform(chunk) if xTransform is not None else chunk[xcol].to_numpy()
        y = yTransform(chunk) if yTransform is not None else chunk[ycol].to_numpy()
        h, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
        hist += h
    err = hist / np.sqrt(hist)
    xcenters = 0.5 * (xbins[1:] + xbins[:-1])
    ycenters = 0.5 * (ybins[1:] + ybins[:-1])
    return xcenters, ycenters, hist, err


def plotThetaDistribution(ax, thetaCounts, bins):
    """Draw the theta-distribution side panel from accumulated bin counts."""
    ax.stairs(thetaCounts, bins)
    ax.set_xlabel(r'$\theta$ (deg)')


def p0CrystalBall(ycenters, histSlice):
    """Build robust initial Crystal Ball parameters for one histogram slice."""
    cdf = np.cumsum(histSlice)
    cdf = cdf / cdf[-1]
    mu0 = float(np.interp(0.5, cdf, ycenters))
    y_low = float(np.interp(0.1587, cdf, ycenters))
    y_high = float(np.interp(0.8413, cdf, ycenters))
    sigma0 = 0.5 * (y_high - y_low)
    if sigma0 <= 0 or np.isnan(sigma0):
        yrange = ycenters[-1] - ycenters[0]
        sigma0 = yrange / 6.0 if yrange > 0 else 1.0
    return [np.max(histSlice), mu0, sigma0, 1.5, 2.0]


def boundsCrystalBallPositive(ycenters, histSlice, p0):
    """Return positive-peak Crystal Ball parameter bounds."""
    return (
        [0.1 * p0[0], 0.5 * p0[1], 0.1 * p0[2], 0.1, 1.01],
        [5 * p0[0], max(p0[1], ycenters[-1]), 5 * p0[2], 10, 100],
    )


def maskPositive(xcenter, ycenters, histSlice):
    """Keep positive-count bins."""
    return histSlice > 0


def maskPositiveAboveThreshold(thresholdFrac, lowXCutoff):
    """Make the historical positive-response slice mask."""
    def _mask(xcenter, ycenters, histSlice):
        if xcenter < lowXCutoff:
            return histSlice > 0
        threshold = thresholdFrac * max(ycenters)
        return (histSlice > 0) & (ycenters > threshold)
    return _mask


@dataclass
class SliceResult:
    """Per-slice center and width estimates returned by a reducer."""
    center: float
    centerErr: float
    width: float
    widthErr: float
    render: Optional[Callable] = None
    method: str = 'fit'


def reduceSlices(xcenters, ycenters, hist, err, reducer, sliceMask=None, outDir=None):
    """Reduce every x-slice and return centers, errors, widths, and errors.

    Slices with too few usable bins receive NaNs. If `outDir` is supplied, a
    reducer's optional render callback writes one diagnostic image per slice.
    """
    if sliceMask is None:
        sliceMask = maskPositive
    centers, center_errors, widths, width_errors = [], [], [], []
    for i, (hist_slice, err_slice) in enumerate(zip(hist, err)):
        selected = sliceMask(xcenters[i], ycenters, hist_slice)
        hm = hist_slice[selected]
        em = err_slice[selected]
        ym = ycenters[selected]
        result = reducer(ym, hm, em, ycenters) if len(hm) > 3 and np.sum(hm) > 0 else None
        if result is None:
            centers.append(np.nan); center_errors.append(np.nan)
            widths.append(np.nan); width_errors.append(np.nan)
            continue
        centers.append(result.center); center_errors.append(result.centerErr)
        widths.append(result.width); width_errors.append(result.widthErr)
        if outDir is not None and result.render is not None:
            fig, ax = plt.subplots(figsize=(7, 5))
            result.render(ax)
            fig.tight_layout()
            fig.savefig(os.path.join(outDir, 'img_%d.png' % i))
            plt.close(fig)
    return (np.array(centers), np.array(center_errors),
            np.array(widths), np.array(width_errors))


def _mode_and_crossing(ycenters, histSlice, heightFrac):
    """Estimate histogram mode and half-width at a chosen height fraction."""
    ycenters, histSlice = np.array(ycenters), np.array(histSlice)
    peak_index = np.argmax(histSlice)
    mode = float(ycenters[peak_index])
    height = heightFrac * float(histSlice[peak_index])
    left = np.nan
    for i in range(peak_index - 1, -1, -1):
        if histSlice[i] <= height:
            denom = histSlice[i + 1] - histSlice[i]
            frac = (height - histSlice[i]) / denom if denom else 0.5
            left = ycenters[i] + np.clip(frac, 0, 1) * (ycenters[i + 1] - ycenters[i])
            break
    right = np.nan
    for i in range(peak_index + 1, len(histSlice)):
        if histSlice[i] <= height:
            denom = histSlice[i] - histSlice[i - 1]
            frac = (height - histSlice[i - 1]) / denom if denom else 0.5
            right = ycenters[i - 1] + np.clip(frac, 0, 1) * (ycenters[i] - ycenters[i - 1])
            break
    width = (right - left) / 2 if np.isfinite(left) and np.isfinite(right) else np.nan
    return mode, width, left, right, height


def _bootstrap_mode_crossing(ycenters, histSlice, n=200):
    """Estimate mode/HWHM uncertainties with deterministic multinomial bootstrap."""
    total = int(round(np.sum(histSlice)))
    if total < 5:
        return np.nan, np.nan
    rng = np.random.default_rng(12345)
    probabilities = histSlice / histSlice.sum()
    centers, widths = [], []
    for _ in range(n):
        result = _mode_and_crossing(ycenters, rng.multinomial(total, probabilities), 0.5)
        centers.append(result[0])
        if np.isfinite(result[1]):
            widths.append(result[1])
    center_error = np.std(centers, ddof=1) if len(centers) >= 2 else np.nan
    width_error = np.std(widths, ddof=1) if len(widths) >= 10 else np.nan
    return center_error, width_error


def curveFitReducer(fitFunc, p0Func, boundsFunc=None):
    """Create a per-slice Crystal Ball fit reducer with historical fallback.

    Fits use scipy's `curve_fit`; low-statistics, failed, or physically
    implausible fits fall back to a mode plus HWHM estimate with bootstrap
    uncertainties, matching the parent pipeline's behavior.
    """
    def reduce_one(ycMasked, histMasked, errMasked, ycFull):
        total = np.sum(histMasked)
        if total < 20:
            return fallback(ycMasked, histMasked, errMasked, ycFull, 'low_stats')
        p0 = p0Func(ycMasked, histMasked)
        kwargs = dict(sigma=errMasked, absolute_sigma=True, p0=p0, maxfev=5000)
        if boundsFunc is not None:
            kwargs['bounds'] = boundsFunc(ycFull, histMasked, p0)
        try:
            popt, pcov = curve_fit(fitFunc, ycMasked, histMasked, **kwargs)
        except Exception as exc:
            return fallback(ycMasked, histMasked, errMasked, ycFull, f'fit exception {type(exc).__name__}')
        yrange = ycFull[-1] - ycFull[0]
        diag = np.diag(pcov)
        reason = None
        if popt[2] <= 0 or popt[2] > yrange:
            reason = 'invalid sigma'
        elif popt[1] < ycFull[0] or popt[1] > ycFull[-1]:
            reason = 'invalid mu'
        elif np.any(~np.isfinite(diag)) or np.any(diag < 0):
            reason = 'invalid covariance'
        elif pcov[1, 1] > yrange ** 2 or pcov[2, 2] > yrange ** 2:
            reason = 'uncertain fit'
        if reason:
            return fallback(ycMasked, histMasked, errMasked, ycFull, reason)
        return SliceResult(float(popt[1]), float(np.sqrt(pcov[1, 1])),
                           float(popt[2]), float(np.sqrt(pcov[2, 2])), method=fitFunc.__name__)

    def fallback(yc, histSlice, errSlice, ycFull, reason):
        """Return the historical mode/HWHM fallback result."""
        print(f'[curveFitReducer] slice fallback -> mode+HWHM: {reason}', file=sys.stderr)
        mode, width, left, right, height = _mode_and_crossing(yc, histSlice, 0.5)
        center_error, width_error = _bootstrap_mode_crossing(yc, histSlice)
        return SliceResult(float(mode), center_error, float(width) if np.isfinite(width) else np.nan,
                           width_error, method='mode_hwhm')

    return reduce_one
