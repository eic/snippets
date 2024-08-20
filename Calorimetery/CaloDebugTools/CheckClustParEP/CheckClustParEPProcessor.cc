// ----------------------------------------------------------------------------
// 'CheckClustParEPProcessor.h'
// Derek Anderson
// 05.07.2024
//
// A JANA plugin to check the E/p of calorimeter
// clusters and their associated MCParticles
// ----------------------------------------------------------------------------

#define CHECKCLUSTPAREPPROCESSOR_CC

// plugin definition
#include "CheckClustParEPProcessor.h"

// The following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app -> Add(new CheckClustParEPProcessor);
  }
}



// public methods -------------------------------------------------------------

void CheckClustParEPProcessor::InitWithGlobalRootLock() {

  // create output directory
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc -> GetHistFile();
  rootfile -> mkdir(m_config.sPlugin.data()) -> cd();

  // hist binning
  const tuple<size_t, float, float> epBins = make_tuple(1000, 0., 10.);

  // hist names
  const vector<string> vecHistNames = {
    "hNHCalEP",
    "hNECalEP",
    "hBHCalEP",
    "hBECalEP",
    "hPHCalEP",
    "hPECalEP"
  };

  // turn on errors
  TH1::SetDefaultSumw2(true);

  // create hists
  for (const string& name : vecHistNames) {
    m_vecHistEP.push_back(
      new TH1D(name.data(), "", get<0>(epBins), get<1>(epBins), get<2>(epBins))
    );
  }
  return;

}  // end 'InitWithGlobalRootLock()'



void CheckClustParEPProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // grab mc particles
  const auto& particles = *(event -> GetCollection<edm4hep::MCParticle>(m_config.sMCPars.data()));

  // loop over calorimeters
  for (size_t iCalo = 0; iCalo < m_config.sClustAndAssocs.size(); iCalo++) {

    // grab collections
    const auto& clusters = *(event -> GetCollection<edm4eic::Cluster>(m_config.sClustAndAssocs.at(iCalo).first.data()));
    const auto& assocs   = *(event -> GetCollection<edm4eic::MCRecoClusterParticleAssociation>(m_config.sClustAndAssocs.at(iCalo).second.data()));

    // loop over clusters and associations
    for (auto cluster : clusters) {
      for (auto assoc : assocs) {

        // consider only current cluster
        const bool isSameClust = (assoc.getRecID() == cluster.getObjectID().index);
        if (!isSameClust) continue;

        // grab associated particle
        auto mcPar = assoc.getSim();

        // calculate e/p and fill hist
        const double ep = cluster.getEnergy() / edm4hep::utils::magnitude(mcPar.getMomentum());
        m_vecHistEP.at(iCalo) -> Fill(ep);

      } // end assoc loop
    }  // end cluster loop
  }  // end calo loop
  return;

}  // end 'ProcessSequential(std::shared_ptr<const JEvent>&)'



void CheckClustParEPProcessor::FinishWithGlobalRootLock() {

  /* nothing to do */

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
