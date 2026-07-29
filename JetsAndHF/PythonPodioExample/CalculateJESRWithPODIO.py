# =============================================================================
## @file   CalculateJESRWithPODIO.py
#  @author Derek Anderson, building on work by Brian Page
#  @date   07.28.2025
# -----------------------------------------------------------------------------
## @brief Example python script that calculates the Jet Energy Scale
#    (JES) and Resolution (JER) from EICrecon output using PODIO
#    interface. This builds on the calculation originally implemented
#    by Brian Page at
#
#      https://github.com/eic/physics_benchmarks/blob/master/benchmarks/Jets-HF/jets/analysis/jets.cxx
#
#  @usage Inside the eic-shell:
#
#      python CalculateJESRWithPODIO.py
#
#    Or, if you want to specify input/output files:
#
#      python CalculateJESRWithPODIO.py -i <my input file> -o <my output file>
#
# =============================================================================

import argparse as ap
import ROOT
import numpy as np

from podio.reading import get_reader

# default input/output files
IFileDefault = "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.4/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x130/q2_10to100/pythia6.428-1.0_NC_noRad_ep_10x130_q2_10to100_ab.0625.eicrecon.edm4eic.root"
OFileDefault = "test_podio.root"


def Calculate(
    ifile = IFileDefault,
    ofile = OFileDefault,
):
    """Calculate JES/R

    Calculate JES/R for a collection of jets. Excludes
    any generator-level jets with electrons to avoid
    DIS electron.

    Args:
      ifile: input file produced by EICrecon
      ofile: output file
    Returns:
      tuple of the JES and JER
    """

    # use podio reader to open example EICrecon output
    reader = get_reader(ifile)
    print("\n  Starting calculation:")
    print(f"      input  = {ifile}")
    print(f"      output = {ofile}")

    # create histogram for extracting JES/R
    hres = ROOT.TH1D("hEneResPODIO", ";(E_{jet}^{rec} - E_{jet}^{gen}) / E_{jet}^{gen}", 50, -2., 3.)
    hres.Sumw2()

    nframes = len(reader.get("events"))
    print(f"    Starting event loop, {nframes} events to process.")

    # now loop through all events in output
    for iframe, frame in enumerate(reader.get("events")):

        # print progress
        if (iframe % 100 == 0):
            print(f"      Processing event {iframe + 1}/{nframes}...")

        # grab collection of jets
        # available options:
        #   - ReconstructedChargedJets
        #   - GeneratedChargedJets
        #   - ReconstructedJets
        #   - GeneratedJets
        #   - ReconstructedCentauroJets
        #   - GeneratedCentauroJets
        gen_jets = frame.get("GeneratedChargedJets")
        rec_jets = frame.get("ReconstructedChargedJets")

        # loop through generated jets
        for gjet in gen_jets:

            # calculate pseudorapidity, azimuth
            gpt = np.sqrt((gjet.getMomentum().x**2) + (gjet.getMomentum().y**2))
            gphi = np.arctan2(gjet.getMomentum().y, gjet.getMomentum().x)
            gtheta = np.arctan2(gpt, gjet.getMomentum().z)
            geta = -np.log(np.tan(gtheta/2))

            # select only jets w/ pt > 5 (jet
            # reco below this gets tricky)
            if gjet.getEnergy() < 5:
                continue

            # select only jets in |eta_max| - R =
            # 3.5 - 1 = 2.5 (to avoid jets cut off
            # by acceptance)
            if np.abs(geta) > 2.5:
                continue 

            # filter out generated jets with electrons
            # in them (to avoid DIS electron)
            has_electron = False
            for gcst in gjet.getParticles():
                if gcst.getPDG() == 11:
                    has_electron = True
                    break

            if has_electron:
                continue

            # loop through reconstructed jets to find closest
            # one in eta-phi space
            match = None
            min_dist = 9999.
            for rjet in rec_jets:

                # calculate pseudorapidity, azimuth
                rpt = np.sqrt((rjet.getMomentum().x**2) + (rjet.getMomentum().y**2))
                rphi = np.arctan2(rjet.getMomentum().y, rjet.getMomentum().x)
                rtheta = np.arctan2(rpt, rjet.getMomentum().z)
                reta = -np.log(np.tan(rtheta/2))

                # calculate distance in eta-phi space between
                # gen & reco jets
                deta = reta - geta
                dphi = rphi - gphi
                dist = np.sqrt((deta**2) + (dphi**2))

                # if distance is outside tolerance, continue
                if dist > 0.05:
                    continue

                # if reco jet is closest so far, set as match
                if dist < min_dist:
                    match = rjet
                    min_dist = dist

            # if no match found, continue 
            if match is None:
                continue

            # otherwise calculate %-difference 
            eres = (match.getEnergy() - gjet.getEnergy())/gjet.getEnergy()
            hres.Fill(eres)

    # extract JES/R
    fres = ROOT.TF1("fEneRes", "gaus(0)", -0.5, 0.5)
    fres.SetParameter(0, hres.Integral())
    fres.SetParameter(1, hres.GetMean())
    fres.SetParameter(2, hres.GetRMS())
    hres.Fit(fres, "r")
    JES = fres.GetParameter(1) + 1
    JER = fres.GetParameter(2)
    print(f"    Finished event loop:")
    print(f"      JES = {JES}, JER = {JER}")

    # save histogram + fit for reference
    with ROOT.TFile(ofile, "recreate") as ofile:
        ofile.WriteObject(hres, "hEneResPODIO")
        ofile.Close()

    print("  Finished calculation!\n")
    return (JES, JER)


# main ========================================================================

if __name__ == "__main__":

    # set up argments
    parser = ap.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        help = "Input file",
        nargs = '?',
        const = IFileDefault,
        default = IFileDefault,
        type = str
    )
    parser.add_argument(
        "-o",
        "--output",
        help = "Output file",
        nargs = '?',
        const = OFileDefault,
        default = OFileDefault,
        type = str
    )
    args = parser.parse_args()

    # run calculation
    Calculate(args.input, args.output)

# end =========================================================================
