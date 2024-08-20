// ----------------------------------------------------------------------------
// 'CheckClustParEPProcessor.h'
// Derek Anderson
// 05.07.2024
//
// A JANA plugin to check the E/p of calorimeter
// clusters and their associated MCParticles
// ----------------------------------------------------------------------------

#ifndef CHECKCLUSTPAREPPROCESSOR_H
#define CHECKCLUSTPAREPPROCESSOR_H

// c++ utilities
#include <string>
#include <vector>
#include <utility>
// root libraries
#include <TH1.h>
#include <TFile.h>
// jana libraries
#include <JANA/JEventProcessor.h>
#include <JANA/JEventProcessorSequentialRoot.h>
// services
#include "services/rootfile/RootFile_service.h"
// edm types
#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociation.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>

using namespace std;



// CheckClustParEPProcessor definition ----------------------------------------

class CheckClustParEPProcessor: public JEventProcessorSequentialRoot {

  // struct to hold user options
  struct Config {
    string sPlugin;
    string sMCPars;
    vector<pair<string, string>> sClustAndAssocs;
  } m_config = {
    "CheckClustParEP",
    "MCParticles",
    {
      make_pair("HcalEndcapNClusters", "HcalEndcapNClusterAssociations"),
      make_pair("EcalEndcapNClusters", "EcalEndcapNClusterAssociations"),
      make_pair("HcalBarrelClusters",  "HcalBarrelClusterAssociations"),
      make_pair("EcalBarrelClusters",  "EcalBarrelClusterAssociations"),
      make_pair("LFHCALClusters",      "LFHCALClusterAssociations"),
      make_pair("EcalEndcapPClusters", "EcalEndcapPClusterAssociations")
    }
  };  // end Config

  public:

    CheckClustParEPProcessor() { SetTypeName(NAME_OF_THIS); }

    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

  private:

    // output histograms
    vector<TH1D*> m_vecHistEP;

};

#endif

// end ------------------------------------------------------------------------
