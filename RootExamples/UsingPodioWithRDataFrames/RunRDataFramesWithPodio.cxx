// ============================================================================
//! \file   RunRDataFramesWithPodio.cxx
//! \author Derek Anderson
//! \date   01.30.2026
// ----------------------------------------------------------------------------
//! \brief
//!   Example ROOT macro illustrating how to use PODIO
//!   interfaces (like relations) in an RDataFrame-based
//!   analysis.
//!
//! \usage
//!   Has to be run inside eic-shell! Run with: 
//!
//!       $ root -b -q RunRDataFramesWithPodio.cxx
// ============================================================================

#ifndef RunRDataFramesWithPodio_cxx
#define RunRDataFramesWithPodio_cxx

#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4hep/utils/vector_utils.h>
#include <Math/Vector3D.h>
#include <podio/DataSource.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
#include <TFile.h>
#include <TH1.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// aliases for convenience
using RDF    = ROOT::RDataFrame;
using TH1Def = ROOT::RDF::TH1DModel;
using TH2Def = ROOT::RDF::TH2DModel;

R__LOAD_LIBRARY(/opt/local/lib/libedm4eic.so)
R__LOAD_LIBRARY(/opt/local/lib/libedm4hep.so)
R__LOAD_LIBRARY(/opt/local/lib/libedm4hepRDF.so)
R__LOAD_LIBRARY(/opt/local/lib/libedm4hepUtils.so)

// default options
static const std::string OFileDefault = "test_podio_frames.root";
static const std::string IFileDefault = "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.4/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x130/q2_10to100/pythia6.428-1.0_NC_noRad_ep_10x130_q2_10to100_ab.0625.eicrecon.edm4eic.root";
static const std::string IJetsDefault = "ReconstructedCentauroJets";
static const std::string IParsDefault = "ReconstructedBreitFrameParticles";



// ============================================================================
//! Struct to consolidate user options
// ============================================================================
struct Options {
  std::string ofile;  // output file
  std::string ifile;  // input file
  std::string ijets;  // input jet collection
  std::string ipars;  // input particle collection
} DefaultOptions = {
  OFileDefault,
  IFileDefault,
  IJetsDefault,
  IParsDefault
};



// ============================================================================
//! Run RDataFrame analysis using PODIO interfaces 
// ============================================================================
void RunRDataFramesWithPodio(const Options& opt = DefaultOptions) {

  // turn on histogram errors
  TH1::SetDefaultSumw2(true);
  std::cout << "\n  Starting PODIO + RDataFrames example!" << std::endl; 

  // open input/output --------------------------------------------------------

  // create output file
  TFile* output = new TFile(opt.ofile.data(), "recreate");
  if (!output) {
    std::cerr << "PANIC: couldn't open output file!\n" << std::endl;
    exit(-1);
  }
  std::cout << "    Opened output file:\n"
            << "      " << opt.ofile
            << std::endl;

  // initialize dataframe
//ROOT::RDataFrame rdf(std::make_unique<podio::DataSource>(filePath));
  ROOT::RDataFrame frame(std::make_unique<podio::DataSource>(opt.ifile));
  if (frame.Count() == 0) {
    std::cerr << "PANIC: no events found!" << std::endl;
    exit(-1);
  }
  std::cout << "    Opened RDataFrame input:\n"
            << "      " << opt.ifile
            << std::endl;

  // define histograms --------------------------------------------------------

  // helper struct for histogram binning
  struct Bin {
    uint32_t num;
    float    start;
    float    stop;
  };

  // define histogram binnings
  std::vector<Bin> bins = {
    {51, -0.5, 40.5},  // number
    {102, -1.0, 50.},  // energy
    {24, -0.1, 1.1}    // z
  };

  // create hist models for dataframe
  //   -- NOTE "cst" = "constituents"
  TH1Def hdEvtNPars = TH1Def("hEvtNPars", "No. of particles per event;N_{par}", bins[0].num, bins[0].start, bins[0].stop);
  TH1Def hdEvtNJets = TH1Def("hEvtNJets", "No. of jets per event;N_{jet}", bins[0].num, bins[0].start, bins[0].stop);
  TH1Def hdJetNCst  = TH1Def("hJetNCst", "No. of csts per jet;N_{jet}", bins[0].num, bins[0].start, bins[0].stop);
  TH1Def hdJetMom   = TH1Def("hJetMom", "Jet momentum;E_{jet} [GeV]", bins[1].num, bins[1].start, bins[1].stop);
  TH1Def hdCstZ     = TH1Def("hCstZ", "Cst z;z_{cst}", bins[2].num, bins[2].start, bins[2].stop);

  // lambdas for analysis -----------------------------------------------------

  // check if particle collection is empty
  auto hasPars = [](const edm4eic::ReconstructedParticleCollection& pars) {
    return pars.size() > 0;
  };

  // check if jet collection is empty
  //   -- TODO update to edm4eic::JetCollection when ready
  auto hasJets = [](const edm4eic::ReconstructedParticleCollection& jets) {
    return jets.size() > 0;
  };

  auto getNPars = [](const edm4eic::ReconstructedParticleCollection& pars) {
    return pars.size();
  };

  auto getNJets = [](const edm4eic::ReconstructedParticleCollection& jets) {
    return jets.size();
  };

  // get no. of cst.s per jet
  auto getNCsts = [](const edm4eic::ReconstructedParticleCollection& jets) {
    std::vector<uint32_t> ncsts;
    for (const auto& jet : jets) {
      ncsts.push_back(
        jet.getParticles().size()
      );
    }
    return ncsts; 
  };

  auto getJetP = [](const edm4eic::ReconstructedParticle& jet) {
    return edm4hep::utils::magnitude(jet.getMomentum());
  };

  auto getCstZ = [](const edm4eic::ReconstructedParticle& jet) {
    ROOT::Math::XYZVector pjet = ROOT::Math::XYZVector(
        jet.getMomentum().x,
        jet.getMomentum().y,
        jet.getMomentum().z
    );
    std::vector<double> zvals;
    for (const auto& cst : jet.getParticles()) {
      ROOT::Math::XYZVector pcst = ROOT::Math::XYZVector(
          cst.getMomentum().x,
          cst.getMomentum().y,
          cst.getMomentum().z
      );
      zvals.push_back(
        pjet.Dot(pcst) / pjet.Mag2()
      );
    }
    return zvals;
  };

  // run analysis -------------------------------------------------------------

  auto analysis = frame.Filter(hasPars, {opt.ipars})
                       .Filter(hasJets, {opt.ijets})
                       .Define("npars", getNPars, {opt.ipars})
                       .Define("njets", getNJets, {opt.ijets})
                       .Define("ncsts", getNCsts, {opt.ijets})
                       .Define("pjet", getJetP, {opt.ijets})
                       .Define("zcst", getCstZ, {opt.ijets});

  // extract histograms
  auto hEvtNPars = analysis.Histo1D(hdEvtNPars, "npars");
  auto hEvtNJets = analysis.Histo1D(hdEvtNJets, "njets");
  auto hJetNCst  = analysis.Histo1D(hdJetNCst, "ncsts");
  auto hJetMom   = analysis.Histo1D(hdJetMom, "pjet");
  auto hCstZ     = analysis.Histo1D(hdCstZ, "zcst");

  // save & close -------------------------------------------------------------

  // save histograms
  output    -> cd();
  hEvtNPars -> Write();
  hEvtNJets -> Write();
  hJetNCst  -> Write();
  hJetMom   -> Write();
  hCstZ     -> Write();

  // save dataframe for further inspection
  analysis.Snapshot("OutTree", opt.ofile, {"npars", "njets", "pjet", "zcst"});

  // close files
  output -> cd();
  output -> Close();
  std::cout << "    Closed output file\n"
            << "  PODIO + RDataFrame example finished!\n"
            << std::endl;

  // exit without error
  return;

}

#endif

// end ========================================================================
