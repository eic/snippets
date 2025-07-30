// ----------------------------------------------------------------------------
// 'MatchProjectionsAndClusters.cxx'
// Derek Anderson
// 05.20.2024
//
// Quick ROOT macro to match calorimeter clusters
// to track projections.
//
// Usage:
//   root -b -q MatchProjectionsAndClusters.cxx++
// ----------------------------------------------------------------------------

#define MATCHPROJECTIONSANDCLUSTERS_CXX

// c++ utilities
#include <cmath>
#include <string>
#include <limits>
#include <optional>
#include <iostream>
// root libraries
#include <TSystem.h>
// podio libraries
#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/CollectionBase.h>
// edm4eic types
#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/TrackPoint.h>
#include <edm4eic/TrackSegment.h>
#include <edm4eic/TrackSegmentCollection.h>
// edm4hep types
#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>



// user options ---------------------------------------------------------------

struct Options {
  std::string in_file;      // input file
  std::string out_file;     // output file
  std::string clusters;     // name of collection of clusters to match to
  std::string projections;  // name of collection of track projections
  uint64_t    system;       // system id of calo to match to
} DefaultOptions = {
  "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.06.1/epic_craterlake/DIS/NC/10x100/minQ2=10/pythia8NCDIS_10x100_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.1287.eicrecon.edm4eic.root",
  "test.root",
  "EcalBarrelClusters",
  "CalorimeterTrackProjections",
  101 //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
};



// macro body -----------------------------------------------------------------

void MatchProjectionsAndClusters(const Options& opt = DefaultOptions) {

  // announce start of macro
  std::cout << "\n  Beginning cluster-track projection matching macro!" << std::endl;

  // open file w/ frame reader
  podio::ROOTReader reader = podio::ROOTReader();
  reader.openFile( opt.in_file );
  std::cout << "    Opened ROOT-based frame reader." << std::endl;

  // get no. of frames and annoucne
  const uint64_t nFrames = reader.getEntries(podio::Category::Event);
  std::cout << "    Starting frame loop: " << reader.getEntries(podio::Category::Event) << " frames to process." << std::endl;

  // iterate through frames (i.e. events in this case)
  uint64_t nClustTotal   = 0;
  uint64_t nClustMatched = 0;
  for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {

    // announce progress
    std::cout << "      Processing frame " << iFrame + 1 << "/" << nFrames << "..." << std::endl;

    // grab frame
    auto frame = podio::Frame( reader.readNextEntry(podio::Category::Event) );

    // grab collections
    auto& clusters = frame.get<edm4eic::ClusterCollection>( opt.clusters );
    auto& segments = frame.get<edm4eic::TrackSegmentCollection>( opt.projections );

    // loop over clusters
    for (size_t iClust = 0; edm4eic::Cluster cluster : clusters) {

      // grab eta/phi of cluster
      const double etaClust  = edm4hep::utils::eta( cluster.getPosition() );
      const double phiClust  = edm4hep::utils::angleAzimuthal( cluster.getPosition() );

      // match based on eta/phi dstiance
      double distMatch = std::numeric_limits<double>::max(); 

      // loop over projections to find matching one
      std::optional<edm4eic::TrackPoint> match;
      for (edm4eic::TrackSegment segment : segments) {
        for (edm4eic::TrackPoint projection : segment.getPoints()) {

          // ignore if not pointing to calo or at face of calo
          const bool isInSystem = (projection.system == opt.system);
          const bool isAtFace   = (projection.surface == 1);
          if (!isInSystem || !isAtFace) continue;

          // grab eta/phi of projection
          const double etaProject = edm4hep::utils::eta( projection.position );
          const double phiProject = edm4hep::utils::angleAzimuthal( projection.position );

          // get distance to projection
          const double distProject = std::hypot(
            etaClust - etaProject,
            phiClust - phiProject
          );

          // if smallest distance found, update variables accordingly
          if (distProject < distMatch) {
            distMatch = distProject;
            match     = projection;
          }

        }  // end point loop
      }  // end segment loop 

      // do analysis if match was found...
      if ( match.has_value() ) {
        edm4eic::TrackPoint matchedProject = match.value();
        std::cout << "        Found match! match = " << matchedProject << std::endl;
        ++nClustMatched;
      }
      ++nClustTotal;

    }  // end cluster loop
  }  // end frame loop
  std::cout << "    Finished frame loop: " << nClustMatched << "/" << nClustTotal << " clusters matched." << std::endl;

  // announce end and exit
  std::cout << "  End of macro!\n" << std::endl;
  return;

}

// end ------------------------------------------------------------------------
