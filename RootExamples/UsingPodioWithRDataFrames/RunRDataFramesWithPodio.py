# =============================================================================
## @file   RunRDataFramesWithPodio.py
#  @author Derek Anderson
#  @date   01.30.2026
# -----------------------------------------------------------------------------
## @brief
#    Example python script illustrating how to use PODIO
#    interfaces (like relations) in an RDataFrame-based
#    analysis.
#
#  @usage
#    Has to be run inside eic-shell! Run with: 
#
#        $ python RunRDataFramesWithPodio.py
# ============================================================================

from dataclasses import dataclass
import argparser
import ROOT

@dataclass
class Options:
    """Dataclass to hold options

    Members:
      ofile: output file
      ifile: input file
      ijets: input jet collection
      ipars: input particle collection
    """
    ofile: str
    ifile: str
    ijets: str
    ipars: str

# Set default options
DefaultOpts = Options(
    ofile = "test_podio_frames.root",
    ifile = "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.4/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x130/q2_10to100/pythia6.428-1.0_NC_noRad_ep_10x130_q2_10to100_ab.0625.eicrecon.edm4eic.root",
    ijets = "ReconstructedCentauroJets",
    ipars = "ReconstructedBreitFrameParticles",
)


def RunRDFAnalysis(opts: Options):

# JIT a C++ function from Python
ROOT.gInterpreter.Declare("""
bool myFilter(float x) {
    return x > 10;
}
""")
 
    df = ROOT.RDataFrame("myTree", "myFile.root")
# Use the function in an RDF operation
sum = df.Filter("myFilter(x)").Sum("y")
print(sum.GetValue())

# main ========================================================================

if __name__ == "__main__":

# end =========================================================================
