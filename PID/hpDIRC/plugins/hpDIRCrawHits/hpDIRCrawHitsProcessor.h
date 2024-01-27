
//
// Template for this file generated with eicmkplugin.py
// Author: Nilanga Wickramaarachchi
//

#include <JANA/JEventProcessorSequentialRoot.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#define MAXHIT 50000

#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/MCRecoTrackerHitAssociation.h>
#include <edm4hep/MCParticle.h>

class hpDIRCrawHitsProcessor: public JEventProcessorSequentialRoot {
private:
  PrefetchT<edm4eic::RawTrackerHit> rawhits = {this, "DIRCRawHits"};
  PrefetchT<edm4eic::MCRecoTrackerHitAssociation> in_hit_assocs = {this, "DIRCRawHitsAssociations"};
  
    // Declare histogram and tree pointers here. e.g.
    // TH1D* hEraw  = nullptr;
    // TH2D* hEdigi = nullptr ;
  TH2D* hHitposition = nullptr;
  TH2D* hHitposition_after = nullptr;
  TH1D* hnhits = nullptr;
  TH1D* hnhits_qe = nullptr;
  TH1D* hnhits_qe_primary = nullptr;
  TH1I* hist_num_sim_hits = nullptr;
  TH1I* hist_parent_pdg = nullptr;

  TTree* dirctree;

  Int_t nhits;
  Int_t mcp_id[MAXHIT];
  Int_t pixel_id[MAXHIT];
  Double_t hit_pos_x[MAXHIT];
  Double_t hit_pos_y[MAXHIT];
  Double_t hit_pos_z[MAXHIT];
  Double_t hit_time[MAXHIT];

  double timeResolution = 1/16.0;

  double prism_width = 350.0;
  double lens_height = 50.0;
  double dirc_r_avg = 708.5; // --- r_min = 700 mm, r_max = 717 mm, r_avg = (r_min + r_max)/2
  double prism_x_edge = dirc_r_avg - lens_height/2;
  double prism_y_edge = prism_width/2;
  double NCol = 4, NRow = 6;
  double prism_x_dim = lens_height + (300 * tan(32*TMath::DegToRad()));
  double prism_y_dim = prism_width;
  double MCP_total_dim = 57.0;
  double MCP_active_dim = 53.0;
  double pixel_dim = MCP_active_dim / 16.0;

  double gap_x = (prism_x_dim - NCol*MCP_total_dim) / (NCol + 1);
  double gap_y = (prism_y_dim - NRow*MCP_total_dim) / (NRow + 1);

  double mcp_active_x_edge = prism_x_edge + gap_x + 2;
  double mcp_active_y_edge = prism_y_edge - gap_y - 2;
  
public:
    hpDIRCrawHitsProcessor() { SetTypeName(NAME_OF_THIS); }

    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;
};
