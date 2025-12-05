// ============================================================================
//! \file   PlotDetectorTypes.cxx 
//! \author Derek Anderson
//! \date   11.20.2025
// ----------------------------------------------------------------------------
//! \brief
//!   Short macro to extract and plot detector 
//!   type for various sim hit collections.
//!
//! \usage
//!   Depends on PODIO and EDM4hep. Can be run
//!   out of the box inside eic-shell.
//!
//!       $ root -b -q PlotDetectorTypes.cxx
// ============================================================================

#define PlotDetectorTypes_cxx

#include <DD4hep/Detector.h>
#include <DD4hep/DetType.h>
#include <DDRec/CellIDPositionConverter.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>
#include <podio/CollectionBase.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH2.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>



// ============================================================================
//! Struct to consolidate user options
// ============================================================================
struct Options {
  std::string ifile;    // input file
  std::string ofile;    // output file
  bool        progress; // print progress through frame loop
} DefaultOptions = {
  "./forEMxHCalFlags.withTrackerTypesPropagated.py8ncdis18x275q1Knevt1K.edm4hep.root",
  "test.root",
  true
};



// ============================================================================
//! Extract and plot detector type for various hit collections
// ============================================================================
void PlotDetectorTypes(const Options& opt = DefaultOptions) {

  gErrorIgnoreLevel = kError;
  std::cout << "\n  Beginning detector types extraction/plotting." << std::endl;

  // --------------------------------------------------------------------------
  // open inputs/outputs
  // --------------------------------------------------------------------------

  // open file w/ frame reader
  podio::ROOTReader reader = podio::ROOTReader();
  reader.openFile( opt.ifile );
  std::cout << "    Opened PODIO reader:\n"
            << "      " << opt.ifile
            << std::endl;

  // open output file
  TFile* output = new TFile(opt.ofile.data(), "recreate");
  if (!output) {
    std::cerr << "PANIC: couldn't open output file!\n" << std::endl;
    exit(-1);
  } else {
    std::cout << "    Opened output file:\n"
              << "      " << opt.ofile
              << std::endl;
  }

  // --------------------------------------------------------------------------
  // create histograms
  // --------------------------------------------------------------------------

  // turn on errors
  TH2::SetDefaultSumw2(true);

  // helper struct for histogram binning
  struct Bin {
    uint32_t num;
    float    start;
    float    stop;
  };

  // histogram binnings to use
  std::vector<Bin> bins = {
      {20000, -10, 10},  // eta
      {550,  -50, 500}   // r [cm]
  };

  // axis titles
  std::string axes(";#eta; r [cm]");

  // define histograms
  TH2D* hPIDRxZ    = new TH2D("hPIDRxZ", axes.data(), bins[0].num, bins[0].start, bins[0].stop, bins [1].num, bins[1].start, bins[1].stop);
  TH2D* hTrkrRxZ   = new TH2D("hTrkrRxZ", axes.data(), bins[0].num, bins[0].start, bins[0].stop, bins[1].num, bins[1].start, bins[1].stop);
  TH2D* hECalRxZ   = new TH2D("hECalRxZ", axes.data(), bins[0].num, bins[0].start, bins[0].stop, bins[1].num, bins[1].start, bins[1].stop);
  TH2D* hHCalRxZ   = new TH2D("hHCalRxZ", axes.data(), bins[0].num, bins[0].start, bins[0].stop, bins[1].num, bins[1].start, bins[1].stop);
  TH2D* hBarrelRxZ = new TH2D("hBarrelRxZ", axes.data(), bins[0].num, bins[0].start, bins[0].stop, bins[1].num, bins[1].start, bins[1].stop);
  TH2D* hEndcapRxZ = new TH2D("hEndcapRxZ", axes.data(), bins[0].num, bins[0].start, bins[0].stop, bins[1].num, bins[1].start, bins[1].stop);
  std::cout << "    Created histograms." << std::endl;

  // --------------------------------------------------------------------------
  // initialize DD4hep interfaces
  // --------------------------------------------------------------------------

  // initialize detector description
  auto detector = dd4hep::Detector::make_unique("");
  try {
    auto* detConfig = std::getenv("DETECTOR_CONFIG");
    auto* detPath   = std::getenv("DETECTOR_PATH");
    if (detConfig && detPath) {
      detector -> fromCompact(
        std::string(detPath) + "/" + std::string(detConfig) + ".xml"
      );
    } else {
      std::cerr << "PANIC: check that DETECTOR_CONFIG and DETECTOR_PATH are set!\n" << std::endl;
      exit(-1);
    }
  } catch(const std::runtime_error& err) {
    std::cerr << "PANIC: error initializing detector!\n" << std::endl;
    exit(-1);
  }
  std::cout << "    Initialized detector." << std::endl;

  // set up id-position converter
  //   --> this will be needed to look up volume
  //       context for a hit
  dd4hep::rec::CellIDPositionConverter converter(*detector);

  // set which calo hit collections to plot
  std::vector<std::string> calCollects = {
    "B0ECalHits",
    "EcalBarrelImagingHits",
    "EcalBarrelScFiHits",
    "EcalEndcapNHits",
    "EcalEndcapPHits",
    "EcalEndcapPInsertHits",
    "EcalFarForwardZDCHits",
    "EcalLumiSpecHits",
    "HcalBarrelHits",
    "HcalEndcapNHits",
    "HcalEndcapPInsertHits",
    "HcalFarForwardZDCHits",
    "LFHCALHits",
    "LumiDirectPCALHits"
  };

  // set which tracker hit collections to plot
  //   --> (includes PID detectors)
  std::vector<std::string> trkCollects = {
    "B0TrackerHits",
    "BackwardMPGDEndcapHits",
    "ForwardMPGDEndcapHits",
    "ForwardOffMTrackerHits",
    "ForwardRomanPotHits",
    "LumiSpecTrackerHits",
    "MPGDBarrelHits",
    "OuterMPGDBarrelHits",
    "SiBarrelHits",
    "TaggerTrackerHits",
    "TOFBarrelHits",
    "TOFEndcapHits",
    "TrackerEndcapHits",
    "VertexBarrelHits"
  };

  // set which pid hit collections to plot
  std::vector<std::string> pidCollects = {
    "DIRCBarHits",
    "DRICHHits",
    "RICHEndcapNHits"
  };

  // announce which collections are going to be plotted
  std::cout << "    Plotting these collections:" << std::endl;
  for (const auto& collect : calCollects) {
    std::cout << "      " << collect << std::endl;
  }
  for (const auto& collect : trkCollects) {
    std::cout << "      " << collect << std::endl;
  }
  for (const auto& collect : pidCollects) {
    std::cout << "      " << collect << std::endl;
  }

  // --------------------------------------------------------------------------
  // frame loop 
  // --------------------------------------------------------------------------

  // announce start of frame loop
  const uint64_t nFrames = reader.getEntries(podio::Category::Event);
  std::cout << "    Starting frame loop: " << nFrames << " frames to process." << std::endl;

  // iterate through frames
  for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {

    // announce progress
    uint64_t iProg = iFrame + 1;
    if (opt.progress) {
      if (iProg % 100 == 0) {
        std::cout << "      Processing frame " << iProg << "/" << nFrames << "..." << std::endl;
      }
    }

    // grab frame for process 
    auto frame = podio::Frame( reader.readNextEntry(podio::Category::Event) );

    // calo hit collection loop -----------------------------------------------

    for (const auto& collect : calCollects) {

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collect);
      for (const auto& hit : hits) {

        // grab cell ID
        //   --> we'll need this to lookup the volume context
        const auto id = hit.getCellID();

        // calculate hit position in (eta, r) space
        const double eta = edm4hep::utils::eta(hit.getPosition()); 
        const double rad = std::hypot(hit.getPosition().x, hit.getPosition().y) / 10.;

        // grab grab volume context via CellID-position converter
        //   --> the volume context maps CellIDs onto the actual
        //       det elements
        const dd4hep::VolumeManagerContext* context = converter.findContext(id);
        if (context) {

          const dd4hep::DetElement det = context -> element;
          const dd4hep::DetType    type(det.typeFlag());
          if (type.is(dd4hep::DetType::CALORIMETER)) {
            if (type.is(dd4hep::DetType::ELECTROMAGNETIC)) {
              hECalRxZ -> Fill(eta, rad);
            }
            if (type.is(dd4hep::DetType::HADRONIC)) {
              hHCalRxZ -> Fill(eta, rad);
            }
            if (type.is(dd4hep::DetType::BARREL)) {
              hBarrelRxZ -> Fill(eta, rad);
            }
            if (type.is(dd4hep::DetType::ENDCAP)) {
              hEndcapRxZ -> Fill(eta, rad);
            }
          }
        }
      }  // end hit loop
    }  // end calo collection loop

    // tracker hit collection loop --------------------------------------------

    for (const auto& collect : trkCollects) {

      const auto& hits = frame.get<edm4hep::SimTrackerHitCollection>(collect);
      for (const auto& hit : hits) {

        // grab cell ID
        //   --> we'll need this to lookup the volume context
        const auto id = hit.getCellID();

        // calculate hit position in (eta, r) space
        const double eta = edm4hep::utils::eta(hit.getPosition()); 
        const double rad = std::hypot(hit.getPosition().x, hit.getPosition().y) / 10.;

        // grab grab volume context via CellID-position converter
        //   --> the volume context maps CellIDs onto the actual
        //       det elements
        const dd4hep::VolumeManagerContext* context = converter.findContext(id);
        if (context) {

          const dd4hep::DetElement det = context -> element;
          const dd4hep::DetType    type(det.typeFlag());
          if (type.is(dd4hep::DetType::TRACKER)) {
            if (type.is(dd4hep::DetType::BARREL)) {
              hBarrelRxZ -> Fill(eta, rad);
            }
            if (type.is(dd4hep::DetType::ENDCAP)) {
              hEndcapRxZ -> Fill(eta, rad);
            }
            hTrkrRxZ -> Fill(eta, rad);
          }
        }
      }  // end hit loop
    }  // end trk collection loop

    // pid hit collection loop ------------------------------------------------

    for (const auto& collect : pidCollects) {

      const auto& hits = frame.get<edm4hep::SimTrackerHitCollection>(collect);
      for (const auto& hit : hits) {

        // grab cell ID
        //   --> we'll need this to lookup the volume context
        const auto id = hit.getCellID();

        // calculate hit position in (eta, r) space
        const double eta = edm4hep::utils::eta(hit.getPosition()); 
        const double rad = std::hypot(hit.getPosition().x, hit.getPosition().y) / 10.;

        // grab grab volume context via CellID-position converter
        //   --> the volume context maps CellIDs onto the actual
        //       det elements
        const dd4hep::VolumeManagerContext* context = converter.findContext(id);
        if (context) {

          const dd4hep::DetElement det = context -> element;
          const dd4hep::DetType    type(det.typeFlag());
          if (type.is(dd4hep::DetType::CHERENKOV)) {
            if (type.is(dd4hep::DetType::BARREL)) {
              hBarrelRxZ -> Fill(eta, rad);
            }
            if (type.is(dd4hep::DetType::ENDCAP)) {
              hEndcapRxZ -> Fill(eta, rad);
            }
            hPIDRxZ -> Fill(eta, rad);
          }
        }
      }  // end hit loop
    }  // end pid collection loop

  }  // end frame loop
  std::cout << "    Frame loop complete!" << std::endl;

  // --------------------------------------------------------------------------
  // Draw hit distributions
  // --------------------------------------------------------------------------

  // set histogram colors
  const uint32_t kPID    = 859;
  const uint32_t kTrkr   = 819;
  const uint32_t kECal   = 899;
  const uint32_t kHCal   = 879;
  const uint32_t kBarrel = 809;
  const uint32_t kEndcap = 859;

  hPIDRxZ    -> SetMarkerColor(kPID);
  hPIDRxZ    -> SetLineColor(kPID);
  hPIDRxZ    -> SetFillColor(kPID);
  hTrkrRxZ   -> SetMarkerColor(kTrkr);
  hTrkrRxZ   -> SetLineColor(kTrkr);
  hTrkrRxZ   -> SetFillColor(kTrkr);
  hECalRxZ   -> SetMarkerColor(kECal);
  hECalRxZ   -> SetLineColor(kECal);
  hECalRxZ   -> SetFillColor(kECal);
  hHCalRxZ   -> SetMarkerColor(kHCal);
  hHCalRxZ   -> SetLineColor(kHCal);
  hHCalRxZ   -> SetFillColor(kHCal);
  hBarrelRxZ -> SetMarkerColor(kBarrel);
  hBarrelRxZ -> SetLineColor(kBarrel);
  hBarrelRxZ -> SetFillColor(kBarrel);
  hEndcapRxZ -> SetMarkerColor(kEndcap);
  hEndcapRxZ -> SetLineColor(kEndcap);
  hEndcapRxZ -> SetFillColor(kEndcap);

  // plot to compare
  //   --> tracker vs. ecal vs. hcal
  //   --> barrel vs. endcap
  TCanvas* cType = new TCanvas("cType", "", 950, 950);
  cType    -> SetTopMargin(0.02);
  cType    -> SetRightMargin(0.02);
  cType    -> cd(); 
  hTrkrRxZ -> Draw("BOX");
  hPIDRxZ  -> Draw("BOX SAME");
  hECalRxZ -> Draw("BOX SAME");
  hHCalRxZ -> Draw("BOX SAME");
  output   -> cd();
  cType    -> Write();
  cType    -> Close();

  TCanvas* cRegion = new TCanvas("cRegion", "", 950, 950);
  cRegion    -> SetTopMargin(0.02);
  cRegion    -> SetRightMargin(0.02);
  cRegion    -> cd();
  hBarrelRxZ -> Draw("BOX");
  hEndcapRxZ -> Draw("BOX SAME");
  output     -> cd();
  cRegion    -> Write();
  cRegion    -> Close();

  // --------------------------------------------------------------------------
  // save & exit
  // --------------------------------------------------------------------------

  // save output & close
  output     -> cd();
  hPIDRxZ    -> Write();
  hTrkrRxZ   -> Write();
  hECalRxZ   -> Write();
  hHCalRxZ   -> Write();
  hBarrelRxZ -> Write();
  hEndcapRxZ -> Write();
  output     -> Close();
  std::cout << "    Saved histograms and closed file." << std::endl;

  // exit
  std::cout << "  Macro complete!\n" << std::endl;
  return;

}

// end ========================================================================
