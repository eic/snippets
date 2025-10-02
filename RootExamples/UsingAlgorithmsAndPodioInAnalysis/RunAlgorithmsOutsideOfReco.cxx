// ============================================================================
//! \file   RunAlgorithmsOutsideOfReco.cxx
//! \author Derek Anderson
//! \date   10.01.2025
// ----------------------------------------------------------------------------
//! \brief
//!   Example ROOT macro illustrating how to run an
//!   EICrecon algorithm outside of EICrecon, and how to
//!   write out your own PODIO collections.
//!
//! \usage
//!   Has to be run inside eic-shell! run the `setup.sh`
//!   script beforehand to make sure relevant include
//!   locations are in your paths.
//!
//!       $ source setup.sh
//!       $ root -b -q RunAlgorithmsOutsideOfReco.cxx
// ============================================================================

#ifndef RunAlgorithmsOutsideOfReco_cxx
#define RunAlgorithmsOutsideOfReco_cxx

#include <algorithms/logger.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/utils/vector_utils.h>
#include <edm4eic/InclusiveKinematicsCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <EICrecon/algorithms/interfaces/WithPodConfig.h>
#include <EICrecon/algorithms/reco/ElectronReconstruction.h>
#include <EICrecon/algorithms/reco/ElectronReconstructionConfig.h>
#include <EICrecon/algorithms/reco/InclusiveKinematicsElectron.h>
#include <EICrecon/algorithms/reco/JetReconstruction.h>
#include <EICrecon/algorithms/reco/JetReconstructionConfig.h>
#include <EICrecon/algorithms/reco/ScatteredElectronsEMinusPz.h>
#include <EICrecon/algorithms/reco/ScatteredElectronsEMinusPzConfig.h>
#include <podio/CollectionBase.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/ROOTWriter.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

R__LOAD_LIBRARY(EICrecon/plugins/reco.so)



// user options ---------------------------------------------------------------

struct Options {
  std::string in_file;   // input eicrecon file
  std::string out_hist;  // output file for histograms
  std::string out_tree;  // output file for podio tree
} DefaultOptions = {
  "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.06.1/epic_craterlake/DIS/NC/10x100/minQ2=10/pythia8NCDIS_10x100_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.1287.eicrecon.edm4eic.root",
  "test_algo.hists.root",
  "test_algo.tree.root"
};



// macro body -----------------------------------------------------------------

void RunAlgorithmsOutsideOfReco(const Options& opt = DefaultOptions) {

  // announce start
  std::cout << "\n  Starting algorithm example!" << std::endl;

  // open inputs/output files -------------------------------------------------

  // open input/output with podio reader/writer
  podio::ROOTWriter writer = podio::ROOTWriter(opt.out_tree);
  podio::ROOTReader reader = podio::ROOTReader();
  reader.openFile(opt.in_file);
  std::cout << "    Opened PODIO frame reader/writer." << std::endl;

  // open output histogram file
  TFile* outHist = new TFile(opt.out_hist.data(), "recreate");
  if (!outHist) {
    std::cerr << "PANIC: couldn't open output histogram file!" << std::endl;
    return;
  }
  std::cout << "    Opened output file." << std::endl;

  // define output histograms -------------------------------------------------

  // make sure hist errors on on
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  // helper struct for histogram binning
  struct Bin {
    uint32_t num;
    float    start;
    float    stop;
  };

  // define histogram binning
  std::vector<Bin> bins = {
    {101, -10., 1e3},  // Q^2
    {102, -1., 50.},   // energy
    {200 , -5., 5.}    // eta
  };

  // define histograms
  TH1D* hReconQ2  = new TH1D("hReconQ2", "EICrecon Q^{2} (e-method)", bins[0].num, bins[0].start, bins[0].stop);
  TH1D* hRerunQ2  = new TH1D("hRerunQ2", "Recalculated Q2 (e-method)", bins[0].num, bins[0].start, bins[0].stop);
  TH2D* hReconExH = new TH2D("hReconExH", "EICrecon jet E vs. #eta (anti-k_{T})", bins[2].num, bins[2].start, bins[2].stop, bins[1].num, bins[1].start, bins[1].stop);
  TH2D* hRerunExH = new TH2D("hRerunExH", "EICrecon jet E vs. #eta (ee-k_{T})", bins[2].num, bins[2].start, bins[2].stop, bins[1].num, bins[1].start, bins[1].stop);
  std::cout << "    Defined histograms." << std::endl;

  // initialize algorithms to rerun -------------------------------------------

  // modify configurations for algorithms we want to rerun
  eicrecon::JetReconstructionConfig cfgJetReco {
    .jetAlgo = "ee_kt_algorithm"
  };

  eicrecon::ElectronReconstructionConfig cfgElecReco {
    .min_energy_over_momentum = 0.5,
    .max_energy_over_momentum = 1.5
  };

  eicrecon::ScatteredElectronsEMinusPzConfig cfgDISElecSelect {
    .minEMinusPz = 1.0,
    .maxEMinusPz = 100.0
  };

  // set up new jet reconstruction algorithm
/* FIXME figure out correct templating for jet reco
  eicrecon::JetReconstruction<edm4eic::ReconstructedParticleCollection> algoJetReco("ReconstructedChargedEEKtJets");
  algoJetReco.level(algorithms::LogLevel::kInfo);
  algoJetReco.applyConfig(cfgJetReco);
  algoJetReco.init();
*/

  // set up new DIS/kinematic algorithms here
  eicrecon::ElectronReconstruction algoElecReco("LooselyReconstructedElectronsForDIS");
  algoElecReco.level(algorithms::LogLevel::kInfo);
  algoElecReco.applyConfig(cfgElecReco);
  algoElecReco.init();

  /* TODO set up new DIS electron selection & kinematic calculation here */

  // event loop ---------------------------------------------------------------

  // get no. of frames and announce start of event loop
  const uint64_t nEvts = reader.getEntries(podio::Category::Event);
  std::cout << "    Starting event loop: " << nEvts << " events to process." << std::endl;

  // iterate through frames (i.e. events in this case)
  for (uint64_t iEvt = 0; iEvt < nEvts; ++iEvt) {

    // announce progress
    uint64_t iProg = iEvt + 1;
    if (iProg % 100 == 0) {
      std::cout << "      Processing event " << iProg << "/" << nEvts << "..." << std::endl;
    }

    // open input frame
    auto iFrame = podio::Frame(reader.readNextEntry(podio::Category::Event));

    // grab & analyze default collections -------------------------------------

    // grab necessary collections from EICrecon output
    auto& mcPars       = iFrame.get<edm4hep::MCParticleCollection>("MCParticles");
    auto& recoPars     = iFrame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
    auto& recoChrgs    = iFrame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
    auto& recoChrgJets = iFrame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedJets");
    auto& recoElecs    = iFrame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedElectronsForDIS");  // maybe not...
    auto& disElecs     = iFrame.get<edm4eic::ReconstructedParticleCollection>("ScatteredElectronsEMinusPz");  // maybe not...
    auto& kineElecs    = iFrame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsElectron");
    auto& hadronFS     = iFrame.get<edm4eic::HadronicFinalStateCollection>("HadronicFinalState");

    // fill histograms from default collections
    for (const auto& kine : kineElecs) {
      hReconQ2 -> Fill( kine.getQ2() );
    }
    for (const auto& jet : recoChrgJets) {
      hReconExH -> Fill(
        edm4hep::utils::eta( jet.getMomentum() ),
        jet.getEnergy()
      );
    }

    // run modified algorithms & analyze/save output --------------------------

    // create output collections
    auto oHeaders      = std::make_unique<edm4hep::EventHeaderCollection>();
    auto oElectrons    = std::make_unique<edm4eic::ReconstructedParticleCollection>();
    auto oDISElectrons = std::make_unique<edm4eic::ReconstructedParticleCollection>();

    // set header info
    auto header = oHeaders -> create();
    header.setEventNumber(iEvt);
    header.setRunNumber(0);
    header.setWeight(1.0);

    /* TODO rerun jet reconstruction here */

    // rerun electron reconstruction 
    auto inElecReco  = std::make_tuple(&recoPars);
    auto outElecReco = std::make_tuple(oElectrons.get());
    algoElecReco.process(inElecReco, outElecReco);


    /* TODO rerun DIS electron selection & kinematic calculation here */

    // collect output collections into a frame
    // and write to TTree
    auto oFrame = podio::Frame();
    oFrame.put(std::move(oHeaders), "EventHeader");
    oFrame.put(std::move(oElectrons), "LooselyReconstructedElectronsForDIS");
    writer.writeFrame(oFrame, "events");

  }  // end event loop
  std::cout << "    Event loop finished." << std::endl;

  // save output & close ------------------------------------------------------

  // save histograms output
  outHist   -> cd();
  hReconQ2  -> Write();
  hRerunQ2  -> Write();
  hReconExH -> Write();
  hRerunExH -> Write();
  outHist   -> Close();
  std::cout << "    Saved and closed histogram output." << std::endl;

  // save output collections
  writer.finish();

  // announce end & exit
  std::cout << "  Finished algorithm example!\n" << std::endl;
  return;

}

#endif

// end ========================================================================
