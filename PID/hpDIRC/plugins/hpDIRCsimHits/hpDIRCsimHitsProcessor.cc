
//
// Template for this file generated with eicmkplugin.py
// Author: Nilanga Wickramaarachchi
//

#include "hpDIRCsimHitsProcessor.h"
#include "services/rootfile/RootFile_service.h"

// Include appropriate class headers. e.g.
#include <edm4hep/SimCalorimeterHitCollection.h>

// The following just makes this a JANA plugin
extern "C" {
    void InitPlugin(JApplication *app) {
        InitJANAPlugin(app);
        app->Add(new hpDIRCsimHitsProcessor);
    }
}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void hpDIRCsimHitsProcessor::InitWithGlobalRootLock(){
    auto rootfile_svc = GetApplication()->GetService<RootFile_service>();
    auto rootfile = rootfile_svc->GetHistFile();
    rootfile->mkdir("hpDIRCsimHits")->cd();

    // Create histograms here. e.g.
    // hEraw  = new TH1D("Eraw",  "BEMC hit energy (raw)",  100, 0, 0.075);
    // hEdigi = new TH2D("Edigi", "BEMC hit energy (digi) vs. raw", 200, 0, 2000.0, 100, 0, 0.075);
    hHitposition_direct = new TH2D("Hitposition_direct", "hit position Y vs. hit position X; hit position x (mm); hit position y (mm)", 130, 600, 990, 130, -190, 200);
    hnhits = new TH1D("nhits","; no. of hits per event", 400, 0, 2000);
    htime = new TH1D("leadtime","; photon time (ns)", 2000, 0, 100);
    hwavelength = new TH1D("wavelength","; wavelength of photon (nm)", 500, 0, 1000);
    hist_parent_pdg = new TH1I("parent_pdg","; parent PDG", 10000, -5000, 5000);
}

//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void hpDIRCsimHitsProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {
    // Data objects we will need from JANA e.g.
  //const auto &rawhits = *static_cast<const edm4hep::SimCalorimeterHitCollection*>(event->GetCollectionBase("EcalBarrelScFiHits"));

    // Fill histograms here. e.g.
    // for (auto hit : rawhits) hEraw->Fill(hit.getEnergy());

  Int_t count = 0;

  for(auto hit : simhits())
    {
      int parent_pdg = hit->getMCParticle().getParents(0).getPDG();
      bool is_from_primary_track = hit->getMCParticle().getParents(0).getGeneratorStatus();
      
      if(!is_from_primary_track) continue;
      hist_parent_pdg->Fill(parent_pdg);
      
      hHitposition_direct->Fill(hit->getPosition()[0], hit->getPosition()[1]);
      count = count+1;

      // wavelength calculation
      double hit_energy = hit->getEDep() * 1e9; // convert to eV
      double wavelength = (1.2398 / hit_energy) * 1000;
      hwavelength->Fill(wavelength);

      // photon time
      double leadtime = hit->getTime(); // in ns
      htime->Fill(leadtime);
      
    }
  hnhits->Fill(count);
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void hpDIRCsimHitsProcessor::FinishWithGlobalRootLock() {

    // Do any final calculations here.
}

