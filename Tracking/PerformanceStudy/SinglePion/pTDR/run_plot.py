import argparse
import os
import sys
from pathlib import Path

import ast
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

plt.rcParams['figure.figsize'] = [8.0, 6.0]
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=MEDIUM_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=MEDIUM_SIZE)
plt.rc('ytick', labelsize=MEDIUM_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)

deg2rad = np.pi / 180.0


def theta2eta(xx, inverse=0):
    xx = np.array(xx)
    if inverse == 1:
        return np.arctan((np.e) ** (-xx)) * 2
    return -np.log(np.tan(xx / 2.0))


def _load_pwg_table(pwg_file):
    candidates = []
    if pwg_file:
        candidates.append(pwg_file)
    candidates.append(str(Path.cwd() / "pwg_requirements.txt"))
    candidates.append(str(Path(__file__).with_name("pwg_requirements.txt")))
    candidates.append(
        "/global/cfs/cdirs/m3763/shujie/worksim/snippets/Tracking/PerformanceStudy/SinglePion/7etabins/"
        "pwg_requirements.txt"
    )
    for path in candidates:
        if path and os.path.exists(path):
            return pd.read_csv(path, sep=r"\s+", skiprows=1)
    return None


def pwg_value(pwg, varname, eta, mom):
    if pwg is None:
        return None
    if varname not in ("dca", "dp"):
        return None
    cond = (pwg.eta_lo <= eta) & (pwg.eta_hi > eta)
    a = pwg[varname + "_par1"].values[cond][0]
    b = pwg[varname + "_par2"].values[cond][0]
    x = mom
    if varname == "dca":
        x *= np.sin(theta2eta(eta, 1))
        return np.sqrt((a / 1000.0 / x) ** 2 + (b / 1000.0) ** 2)
    if varname == "dp":
        return np.sqrt((a * x) ** 2 + b ** 2)
    return None


def _eta_center(df):
    return 0.5 * (df["eta_lo"] + df["eta_hi"])


def _parse_list(value):
    # Parse CSV list columns and coerce single NaN entries to 0.
    if value is None:
        return None
    parsed = value
    if not isinstance(value, (list, tuple, np.ndarray)):
        try:
            text = str(value)
            text = text.replace("nan", "None").replace("NaN", "None")
            parsed = ast.literal_eval(text)
        except Exception:
            return None
    if not isinstance(parsed, (list, tuple, np.ndarray)):
        return None
    cleaned = [0.0 if v is None or (isinstance(v, float) and np.isnan(v)) else v for v in parsed]
    arr = np.array(cleaned, dtype=float)
    return np.nan_to_num(arr, nan=0.0)


def plot_eff(df, out_path=None, setting="", title=None):
    if setting:
        df = df[df["setting"] == setting]

    df = df.dropna(subset=["eta_lo", "eta_hi", "mom_gev", "eff"]).copy()
    if df.empty:
        print("No data to plot for efficiency.")
        return

    mom_list = [0.5, 1, 2, 5, 10, 15]
    line_styles = [(0, (3, 3, 1, 2)), '--', '-.', ':', '-', (0, (3, 1, 1, 1))]

    eta_bins = None
    if "eff_bins" in df.columns:
        for val in df["eff_bins"]:
            parsed = _parse_list(val)
            if parsed is not None:
                eta_bins = parsed
                break

    plt.figure()
    for ii, mom in enumerate(mom_list):
        dft = df[np.isclose(df["mom_gev"], mom)].copy()
        if dft.empty:
            continue
        if "eff_values" in dft.columns and "eff_errors" in dft.columns and eta_bins is not None:
            ys = None
            cnt = None
            for _, row in dft.iterrows():
                yy = _parse_list(row["eff_values"])
                ee = _parse_list(row["eff_errors"])

                if yy is None or ee is None:
                    continue
                if ys is None:
                    ys = yy.astype(float)
                    cnt = (yy != 0).astype(int)
                else:
                    ys = ys + yy
                    cnt = cnt + (yy != 0).astype(int)
            if ys is None:
                continue
            cnt[cnt == 0] = 1
            ys = ys / cnt
            plt.plot(eta_bins, ys, ls=line_styles[ii], label=f"{mom} GeV")
        else:
            dft["eta_center"] = _eta_center(dft)
            dft = dft.groupby("eta_center")["eff"].mean().reset_index()
            dft = dft.sort_values("eta_center")
            plt.plot(dft["eta_center"], dft["eff"], ls=line_styles[ii], label=f"{mom} GeV")

    plt.legend(frameon=0, loc="upper left", ncol=2, fontsize=13)
    plt.ylim(0.0, 1.4)
    plt.xlim(-4, 4)
    plt.xlabel("$\\eta$")
    plt.ylabel("efficiency")
    plt.grid()
    if title:
        plt.title(title)
    if out_path:
        plt.savefig(out_path)
        plt.close()
        return None
    fig = plt.gcf()
    return fig


def plot_resol(df, varname, out_path, setting1="", setting2="", pwg=None, save=True):
    # Match the legacy 7-eta-bin layout and PWG overlays.
    df = df.dropna(subset=["eta_lo", "eta_hi", "mom_gev"]).copy()
    df["name"] = df["setting"].fillna("")
    df["mom"] = df["mom_gev"]

    eta_lo_pwg = [-3.5, -3.0, -2.5, -1.0, 1.0, 2.5, 3.0]
    eta_hi_pwg = [-3.0, -2.5, -1.0, 1.0, 2.5, 3.0, 3.5]

    if varname == "th":
        y_hi = 0.01
        yname = r"$\theta$ [rad]"
        xname = "momentum [GeV]"
        x_hi = [20, 20, 20, 20, 20, 20, 20]
        sig_col = "resol_theta_sigma"
        err_col = "resol_theta_sigma_err"
        scale = 1.0 / 1000.0
    elif varname == "ph":
        y_hi = 0.025
        yname = r"$\phi$ [rad]"
        xname = "momentum [GeV]"
        x_hi = [20, 20, 20, 20, 20, 20, 20]
        sig_col = "resol_phi_sigma"
        err_col = "resol_phi_sigma_err"
        scale = 1.0 / 1000.0
    elif varname == "dp":
        y_hi = 12
        yname = r"$\delta p/p$ [%]"
        xname = "momentum [GeV/c]"
        x_hi = [20, 20, 20, 20, 20, 20, 20]
        sig_col = "resol_dp_sigma"
        err_col = "resol_dp_sigma_err"
        scale = 1.0
    elif varname == "dca":
        y_hi = 1
        yname = "DCA$_r$ [mm]"
        xname = "pT [GeV]"
        x_hi = [1.5, 2.5, 5, 10, 5, 2.5, 1.5]
        sig_col = "resol_dca_sigma"
        err_col = "resol_dca_sigma_err"
        scale = 1.0
    else:
        print("ERROR(plot_resol): please use a valid varname: th, ph, dp, dca")
        return

    fig, axs = plt.subplots(2, 4, figsize=(16, 8))
    axs = axs.flat

    for ii, e_lo in enumerate(eta_lo_pwg):
        e_hi = eta_hi_pwg[ii]
        ax = axs[ii]

        c1 = df.eta_lo >= e_lo - 0.01
        c2 = df.eta_hi <= e_hi + 0.01
        dft = df[c1 & c2]

        if len(dft) == 0:
            continue

        dft = dft[["name", "mom", sig_col, err_col]].groupby(["name", "mom"]).mean().reset_index()

        if setting1:
            cond = (dft.name == setting1) & (dft[sig_col] > 0)
        else:
            cond = dft[sig_col] > 0

        if varname == "dca":
            xdata = dft.mom[cond] * np.sin(theta2eta((e_lo + e_hi) / 2, 1))
        else:
            xdata = dft.mom[cond]
        ax.errorbar(
            xdata,
            dft[sig_col][cond] * scale,
            yerr=dft[err_col][cond] * scale,
            color="b",
            ls="none",
            marker="o",
        )

        if setting2:
            cond = (dft.name == setting2) & (dft[sig_col] > 0)
            if varname == "dca":
                xdata = dft.mom[cond] * np.sin(theta2eta((e_lo + e_hi) / 2, 1))
            else:
                xdata = dft.mom[cond]
            ax.errorbar(
                xdata,
                dft[sig_col][cond] * scale,
                yerr=dft[err_col][cond] * scale,
                color="r",
                ls="none",
                marker="x",
            )

        xline = np.arange(0.001, 20, 0.001)
        y_pwg = pwg_value(pwg, varname, e_lo, xline) if pwg is not None else None
        if y_pwg is not None:
            ax.plot(xline, y_pwg, 'k--', zorder=10)

        ax.set_ylim(-y_hi * 0.05, y_hi)
        ax.set_xlim(0, x_hi[ii] * 1.05)
        ax.text(x_hi[ii] * 0.1, y_hi * 0.9, f"{e_lo}<$\\eta$<{e_hi}", fontsize=14)

    ax = axs[7]
    ax.axis('off')
    xline = np.arange(0.001, 20, 0.001)
    y_pwg = pwg_value(pwg, varname, eta_lo_pwg[0], xline) if pwg is not None else None
    sim_label = setting1 or "Simulation"
    if y_pwg is not None:
        ax.plot(xline, y_pwg - 100000, "k--", label="PWG Requirements")
        ax.errorbar(xline, y_pwg - 100000, ls="none", marker="o", color="blue", label=sim_label)
    else:
        ax.errorbar(xline, xline * 0 - 100000, ls="none", marker="o", color="blue", label=sim_label)
    if setting2:
        ax.errorbar(xline, xline * 0 - 100000, ls="none", marker="x", color="r", label=setting2)
    ax.set_ylim(0, 1)
    ax.legend(frameon=0, loc="upper left", fontsize=16)

    axs[4].set_xlabel(xname)
    axs[5].set_xlabel(xname)
    axs[6].set_xlabel(xname)

    axs[0].set_ylabel(yname)
    axs[4].set_ylabel(yname)

    plt.tight_layout()
    if save:
        plt.savefig(out_path)
        plt.close()
        return None
    return fig


def main():
    parser = argparse.ArgumentParser(description="Plot single-particle performance summaries.")
    parser.add_argument("--table", default="performance_table.csv")
    parser.add_argument("--out-dir", default="./plots")
    parser.add_argument("--setting", default="")
    parser.add_argument("--setting2", default="")
    # Always plot efficiency; no flag needed.
    parser.add_argument("--plot-resol", nargs="*", default=["dp", "th", "ph", "dca"])
    parser.add_argument("--pwg-file", default="")

    args = parser.parse_args()

    df = pd.read_csv(args.table)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pwg = _load_pwg_table(args.pwg_file)

    tag = args.setting or "all"
    if args.setting2:
        tag = f"{tag}_{args.setting2}"

    combined_path = str(out_dir / f"tracking_single_particle_perf_{tag}.pdf")
    print(f"Plotting efficiency and resolutions -> {combined_path}")
    with PdfPages(combined_path) as pdf:
        eff_fig = plot_eff(df, out_path=None, setting=args.setting, title=args.setting or "all")
        if eff_fig is not None:
            pdf.savefig(eff_fig)
            plt.close(eff_fig)

        if args.setting2:
            eff_fig = plot_eff(df, out_path=None, setting=args.setting2, title=args.setting2)
            if eff_fig is not None:
                pdf.savefig(eff_fig)
                plt.close(eff_fig)

        if args.plot_resol:
            for varname in args.plot_resol:
                fig = plot_resol(
                    df,
                    varname,
                    combined_path,
                    setting1=args.setting,
                    setting2=args.setting2,
                    pwg=pwg,
                    save=False,
                )
                if fig is not None:
                    pdf.savefig(fig)
                    plt.close(fig)


if __name__ == "__main__":
    main()
