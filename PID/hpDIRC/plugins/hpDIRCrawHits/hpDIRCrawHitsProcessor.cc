
//
// Template for this file generated with eicmkplugin.py
// Author: Nilanga Wickramaarachchi
//

#include "hpDIRCrawHitsProcessor.h"
#include "services/rootfile/RootFile_service.h"

// Include appropriate class headers. e.g.
#include <edm4hep/SimCalorimeterHitCollection.h>

// The following just makes this a JANA plugin
extern "C" {
    void InitPlugin(JApplication *app) {
        InitJANAPlugin(app);
        app->Add(new hpDIRCrawHitsProcessor);
    }
}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void hpDIRCrawHitsProcessor::InitWithGlobalRootLock(){
    auto rootfile_svc = GetApplication()->GetService<RootFile_service>();
    auto rootfile = rootfile_svc->GetHistFile();
    rootfile->mkdir("hpDIRCrawHits")->cd();

    // Create histograms here. e.g.
    // hEraw  = new TH1D("Eraw",  "BEMC hit energy (raw)",  100, 0, 0.075);
    // hEdigi = new TH2D("Edigi", "BEMC hit energy (digi) vs. raw", 200, 0, 2000.0, 100, 0, 0.075);
    hHitposition = new TH2D("Hitposition", "hit position Y vs. hit position X; hit position x (mm); hit position y (mm)", 110, 683.5, 1047.875, 110, -175, 189.375);
    hHitposition_after = new TH2D("Hitposition_after", "hit position Y vs. hit position X; hit position x (mm); hit position y (mm)", 110, 683.5, 1047.875, 110, -175, 189.375);

    
    hnhits = new TH1D("nhits","; no. of hits per event", 500, 0, 1000);
    hnhits_qe = new TH1D("nhits_qe","; no. of hits per event", 500, 0, 1000);
    hnhits_qe_primary = new TH1D("nhits_qe_primary","; no. of hits per event", 500, 0, 1000);

    hist_num_sim_hits = new TH1I("num_sim_hits","; no. of sim hits for raw hit", 10, 0, 10);
    hist_parent_pdg = new TH1I("parent_pdg","; parent PDG", 10000, -5000, 5000);
    
    dirctree = new TTree("dirctree", "Tree for hpDIRC");

    /// Event level
    dirctree->Branch("nhits", &nhits, "nhits/I");

    /// Hit level
    dirctree->Branch("hit_pos_x", hit_pos_x, "hit_pos_x[nhits]/D");
    dirctree->Branch("hit_pos_y", hit_pos_y, "hit_pos_y[nhits]/D");
    dirctree->Branch("hit_pos_z", hit_pos_z, "hit_pos_z[nhits]/D");
    dirctree->Branch("mcp_id", mcp_id, "mcp_id[nhits]/I");
    dirctree->Branch("pixel_id", pixel_id, "pixel_id[nhits]/I");
    dirctree->Branch("hit_time", hit_time, "hit_time[nhits]/D");
}

//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void hpDIRCrawHitsProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

    // Fill histograms here. e.g.
    // for (auto hit : rawhits) hEraw->Fill(hit.getEnergy());
  Int_t count = 0;
  Int_t count_qe = 0;
  Int_t count_qe_primary = 0;
  
  for(auto hit : rawhits())
    {
      count_qe = count_qe+1;
      
      hit_time[count] = hit->getTimeStamp() * timeResolution;
      edm4hep::Vector3d hit_position;
      int parent_pdg;
      bool is_from_primary_track;
      
      for(const auto hit_assoc : in_hit_assocs()) {
            if(hit_assoc->getRawHit().isAvailable()) {
              if(hit_assoc->getRawHit().id() == hit->id()) {
                if(hit_assoc->simHits_size() > 0) {
                  hit_position = hit_assoc->getSimHits(0).getPosition();
		  hist_num_sim_hits->Fill(hit_assoc->simHits_size());

		  parent_pdg = hit_assoc->getSimHits(0).getMCParticle().getParents(0).getPDG();
		  is_from_primary_track = hit_assoc->getSimHits(0).getMCParticle().getParents(0).getGeneratorStatus();
		}
	      }
	    }
      }

      if(!is_from_primary_track) continue;
      count_qe_primary = count_qe_primary + 1;

      hHitposition->Fill(hit_position[0], hit_position[1]);
      hist_parent_pdg->Fill(parent_pdg);

      hit_pos_x[count] = hit_position[0];
      hit_pos_y[count] = hit_position[1];
      hit_pos_z[count] = hit_position[2];
            
      if(hit_pos_x[count] > (prism_x_edge + prism_x_dim)) continue;

      int mx = (hit_pos_x[count] - (prism_x_edge + gap_x)) / (MCP_total_dim + gap_x);
      int my = ((prism_y_edge - gap_y) - hit_pos_y[count]) / (MCP_total_dim + gap_y);
      mcp_id[count] = (6 * mx) + my;

      float hit_mcp_x_edge = mcp_active_x_edge + (MCP_active_dim + gap_x + 4)*mx;
      float hit_mcp_y_edge = mcp_active_y_edge - (MCP_active_dim + gap_y + 4)*my;

      if((hit_pos_x[count] - hit_mcp_x_edge) < 0) continue;
      if((hit_pos_y[count] - hit_mcp_y_edge) > 0) continue; 
      
      int pixel_x = (hit_pos_x[count] - hit_mcp_x_edge) / pixel_dim;
      int pixel_y = (hit_mcp_y_edge - hit_pos_y[count]) / pixel_dim;

      if(pixel_x > 15 || pixel_y > 15) continue;      
      
      if(((hit_pos_x[count] - hit_mcp_x_edge) > 0) && (fmod(hit_pos_x[count] - hit_mcp_x_edge, pixel_dim) == 0)) pixel_x -= 1;
      if(((hit_mcp_y_edge - hit_pos_y[count]) > 0) && (fmod(hit_mcp_y_edge - hit_pos_y[count], pixel_dim) == 0)) pixel_y -= 1;

      pixel_id[count] = (16 * pixel_x) + pixel_y;

      hHitposition_after->Fill(hit_position[0], hit_position[1]);
      count = count+1;
    }

  hnhits->Fill(count);
  hnhits_qe->Fill(count_qe);
  hnhits_qe_primary->Fill(count_qe_primary);

  nhits = count;
  dirctree->Fill();
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void hpDIRCrawHitsProcessor::FinishWithGlobalRootLock() {

    // Do any final calculations here.
}

