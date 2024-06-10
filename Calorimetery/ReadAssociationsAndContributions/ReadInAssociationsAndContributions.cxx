// ----------------------------------------------------------------------------
// 'ReadInAssociationsAndContributions.cxx'
// Derek Anderson
// 06.10.2024
//
// Quick ROOT macro to read in both the
// associations of and contributions
// to a reconstructed cluster.
//
// Usage:
//   root -b -q ReadInAssociationsAndContributions.cxx++
// ----------------------------------------------------------------------------

#define READINASSOCIATIONSANDCONTRIBUTIONS_CXX

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
#include <podio/CollectionBase.h>
#include <podio/ROOTFrameReader.h>
// edm4hep types
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/CaloHitContribution.h>
// edm4eic types
#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociation.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>
#include <edm4eic/CalorimeterHit.h>



// user options ---------------------------------------------------------------

struct Options {
  std::string in_file;       // input file
  std::string out_file;      // output file
  std::string clusters;      // name of collection of clusters
  std::string associations;  // name of collection of cluster-particle associations
  std::string sim_hits;      // name of collection of mc calorimeter hits
} DefaultOptions = {
  "../normal/forQuickAssociationAnswers_normalMode_nEvt500.py8s18x275minq100ang0025div1.run0file0.podio.root",
  "test.root",
  "LFHCALClusters",
  "LFHCALClusterAssociations",
  "LFHCALHits"
};



// macro body -----------------------------------------------------------------

void ReadInAssociationsAndContributions(const Options& opt = DefaultOptions) {

  // announce start of macro
  std::cout << "\n  Beginning cluster-track projection matching macro!" << std::endl;

  // open file w/ frame reader
  podio::ROOTFrameReader reader = podio::ROOTFrameReader();
  reader.openFile( opt.in_file );
  std::cout << "    Opened ROOT-based frame reader." << std::endl;

  // get no. of frames and annoucne
  const uint64_t nFrames = reader.getEntries(podio::Category::Event);
  std::cout << "    Starting frame loop: " << reader.getEntries(podio::Category::Event) << " frames to process." << std::endl;

  // iterate through frames (i.e. events in this case)
  for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {

    // announce progress
    std::cout << "      Processing frame " << iFrame + 1 << "/" << nFrames << "..." << std::endl;

    // grab frame
    auto frame = podio::Frame( reader.readNextEntry(podio::Category::Event) );

    // grab collections
    auto& clusters     = frame.get<edm4eic::ClusterCollection>( opt.clusters );
    auto& associations = frame.get<edm4eic::MCRecoClusterParticleAssociationCollection>( opt.associations );
    auto& sim_hits     = frame.get<edm4hep::SimCalorimeterHitCollection>( opt.sim_hits );

    // look at associations ---------------------------------------------------

    // loop over associations
    for (size_t iAssoc = 0; edm4eic::MCRecoClusterParticleAssociation assoc : associations) {

      // grab associated cluster, particle
      auto cluster  = assoc.getRec();
      auto particle = assoc.getSim();

      // print IDs of both objects
      std::cout << "        [Assoc. #" << iAssoc << "]: (cluster, particle) ID = (" << assoc.getRecID() << ", " << assoc.getSimID() << ")" << std::endl;

      // now print some additional info from each
      std::cout << "          Cluster:  energy = " << cluster.getEnergy() << ", nHits = " << cluster.getNhits() << "\n"
                << "          Particle: energy = " << particle.getEnergy() << ", pdg = " << particle.getPDG() << ", gen status = " << particle.getGeneratorStatus()
                << std::endl;
      ++iAssoc;

    }  // end association loop

    // look at contributions --------------------------------------------------

    // loop over clusters
    for (size_t iClust = 0; edm4eic::Cluster cluster : clusters) {

      // loop through reconstructed cluster hits (cells)
      for (size_t iRec = 0; edm4eic::CalorimeterHit rec : cluster.getHits()) {

        // grab cell ID
        const uint64_t recCellID = rec.getCellID();

        // loop over sim hits
        for (size_t iSim = 0; edm4hep::SimCalorimeterHit sim : sim_hits) {

          // match sim-to-reco hit based on cell ID
          const uint64_t simCellID  = sim.getCellID();
          const bool     isSameCell = (recCellID == simCellID);
          if (!isSameCell) continue;

          // now loop over contributions
          for (size_t iContrib = 0; edm4hep::CaloHitContribution contrib : sim.getContributions()) {

            // grab corresponding particle
            auto particle = contrib.getParticle();

            // print IDs of all 5 objects in hierarchy
            std::cout << "        [Contrib #" << iContrib << "] (cluster, reco hit, sim hit, contrib, particle) ID = ("
                      << cluster.getObjectID().index  << ", "
                      << rec.getObjectID().index      << ", "
                      << sim.getObjectID().index      << ", "
                      << contrib.getObjectID().index  << ", "
                      << particle.getObjectID().index << ")"
                      << std::endl;

            // and now print some info from each layer
            std::cout << "          Cluster:  energy = " << cluster.getEnergy() << ", nHits = << " << cluster.getNhits() << "\n"
                      << "          Rec Hit:  energy = " << rec.getEnergy() << ", time = " << rec.getTime() << "\n"
                      << "          Sim Hit:  energy = " << sim.getEnergy() << "\n"
                      << "          Contrib:  energy = " << contrib.getEnergy() << ", time = " << contrib.getTime() << "\n"
                      << "          Particle: energy = " << particle.getEnergy() << ", pdg = " << particle.getPDG() << ", gen status = " << particle.getGeneratorStatus()
                      << std::endl;
            ++iContrib;

          }  // end contrib loop
          ++iSim;

        }  // end sim hit loop
        ++iRec;

      }  // end reco hit (cell) loop
      ++iClust;

    }  // end cluster loop
  }  // end frame loop
  std::cout << "    Finished frame loop"  << std::endl;

  // announce end and exit
  std::cout << "  End of macro!\n" << std::endl;
  return;

}

// end ------------------------------------------------------------------------
