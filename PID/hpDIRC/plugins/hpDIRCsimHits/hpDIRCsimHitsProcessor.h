
//
// Template for this file generated with eicmkplugin.py
// Author: Nilanga Wickramaarachchi
//

#include <JANA/JEventProcessorSequentialRoot.h>
#include <TH2D.h>
#include <TFile.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/MCParticle.h>

class hpDIRCsimHitsProcessor: public JEventProcessorSequentialRoot {
private:
  PrefetchT<edm4hep::SimTrackerHit> simhits = {this, "DIRCBarHits"};
  
    // Declare histogram and tree pointers here. e.g.
    // TH1D* hEraw  = nullptr;
    // TH2D* hEdigi = nullptr ;
  TH2D* hHitposition_direct = nullptr;
  TH1D* hnhits = nullptr;
  TH1D* htime = nullptr;
  TH1D* hwavelength = nullptr;
  TH1I* hist_parent_pdg = nullptr;
  
public:
    hpDIRCsimHitsProcessor() { SetTypeName(NAME_OF_THIS); }

    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;
};
