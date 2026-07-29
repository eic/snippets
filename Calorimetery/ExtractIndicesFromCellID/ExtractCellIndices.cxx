// ============================================================================
//! \file   ExtractCellIndices.cxx 
//! \author Derek Anderson
//! \date   01.29.2026
// ----------------------------------------------------------------------------
//! \brief
//!   Short macro to illustrate how to extract 
//!   indices from a CellID. Uses calorimeter
//!   hits as the example, but applies to any
//!   object that has a CellID.
//!
//! \usage
//!   Depends on PODIO, EDM4eic, and EDM4hep.
//!   Can be run out of the box inside eic-shell.
//!
//!       $ root -b -q ExtractCellIndices.cxx
// ============================================================================

#define ExtractCellIndices_cxx

#include <DD4hep/BitFieldCoder.h>
#include <DD4hep/Detector.h>
#include <DD4hep/IDDescriptor.h>
#include <edm4eic/CalorimeterHitCollection.h>
#include <podio/CollectionBase.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <iostream>
#include <string>

// default options
static const std::string IFileDefault = "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.4/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x130/q2_10to100/pythia6.428-1.0_NC_noRad_ep_10x130_q2_10to100_ab.0625.eicrecon.edm4eic.root";
static const std::string IReadDefault = "HcalBarrelHits";
static const std::string IHitsDefault = "HcalBarrelRecHits";



// ============================================================================
//! Struct to consolidate user options
// ============================================================================
struct Options {
  std::string ifile;  // input file
  std::string iread;  // hit collection defining readout (e.g. HcalBarrelHits)
  std::string ihits;  // calorimeter hit collection to process (e.g. HcalBarrelRecHits)
  uint32_t    nhits;  // max no. of hits to print out
} DefaultOptions = {
  IFileDefault,
  IReadDefault,
  IHitsDefault,
  5
};



// ============================================================================
//! Extract indices from cell ID's 
// ============================================================================
void ExtractCellIndices(const Options& opt = DefaultOptions) {

  gErrorIgnoreLevel = kError;
  std::cout << "\n  Beginning cell index extraction!" << std::endl;

  // --------------------------------------------------------------------------
  // open input
  // --------------------------------------------------------------------------

  // open file w/ frame reader
  podio::ROOTReader reader = podio::ROOTReader();
  reader.openFile( opt.ifile );
  std::cout << "    Opened PODIO reader:\n"
            << "      " << opt.ifile
            << std::endl;

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

  // set up the ID descriptor
  dd4hep::IDDescriptor descriptor;
  try {
    descriptor  = detector -> readout(opt.iread.data()).idSpec();
  } catch (const std::runtime_error& err) {
    std::cerr << "PANIC: readout class is not in output!\n" << std::endl;
    exit(-1);
  }
  std::cout << "    Initialized descriptor.\n"
            << "    Available fields:"
            << std::endl;

  // extract field map and print it out
  //   --> NOTE the CellID is a bit field that stores all of
  //       the various indices and IDs, i.e. fields, that
  //       define a unique cell.
  //   --> The no. of bits allocated to each field are set
  //       in the "readouts" element in each detector's
  //       compact (xml) file
  //   --> To extract the value of a field from the CellID:
  //         1) First DD4hep masks out all bits NOT
  //            allocated for the field, and
  //         2) Then it right shifts the remaining bits so
  //            that we can convert the value to decimal
  //       The "offset" here is how many bits we need to
  //       right shift.
  auto fields = descriptor.fields(); 
  for (const auto& field : fields) {
      std::cout << "      name  = "  << field.first << "\n"
                << "      offset = " << field.second
                << std::endl;
  }
  std::cout << "    Field indices:" << std::endl;

  // grab decoder and test it
  //   --> NOTE we'll need the decoder to extract
  //       individual fields from a CellID
  dd4hep::DDSegmentation::BitFieldCoder* decoder = NULL;
  try {
    decoder = descriptor.decoder();
    for (const auto& field : fields) {
      std::cout << "      " << field.first << " = " << decoder -> index(field.first) << std::endl;
    }
  } catch (const std::runtime_error &err) {
    std::cerr << "PANIC: something went wrong grabbing the decoder!\n" << std::endl;
    exit(-1);
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
    bool     print = iProg % 100 == 0;
    if (print) {
      std::cout << "      Processing frame " << iProg << "/" << nFrames << "..." << std::endl;
    }

    // grab frame and collection to process
    auto        frame = podio::Frame( reader.readNextEntry(podio::Category::Event) );
    const auto& hits  = frame.get<edm4eic::CalorimeterHitCollection>(opt.ihits);

    // print up to nhits out -------------------------------------------------- 

    if (print) {
      for (uint32_t iHit = 0; const auto& hit : hits) {

        if (iHit >= opt.nhits) {
          break;
        }

        const auto idObj  = hit.getObjectID();
        const auto idCell = hit.getCellID();
        std::cout << "        [Hit " << idObj.index << "] cell ID = " << idCell << std::endl;


        for (const auto& field : fields) {
            std::cout << "          " << field.first << " = " << decoder -> get(idCell, field.first) << std::endl;
        }
        ++iHit;

      }
    }
  }  // end frame loop
  std::cout << "    Frame loop complete!" << std::endl;

  // --------------------------------------------------------------------------
  // exit macro 
  // --------------------------------------------------------------------------

  std::cout << "  Extraction complete!\n" << std::endl;
  return;

}

// end ========================================================================
