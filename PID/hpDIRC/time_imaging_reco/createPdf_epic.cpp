// Author: Nilanga Wickramaarachchi
// Based on https://github.com/rdom/eicdirc - by Roman Dzhygadlo

#include "prttools.cpp"

#include <iostream>

const int nch(24*256);
TH1F *hlef[nch], *hles[nch];

void createPdf_epic(TString in="eicrecon.root", int pid=321) // use pid=11 for e+
{
  for(int i=0; i<nch; i++){
    hlef[i] = new TH1F(Form("lef_%d",i),";LE time [ns]; entries/N [#]", 2000,0,100);
    hles[i] = new TH1F(Form("les_%d",i),";LE time [ns]; entries/N [#]", 2000,0,100);
  }

  TH1F* h_pip_nph = new TH1F("pip_nph",";multiplicity [#]; entries [#]", 220,0,220);

  prt_ch = new TChain("hpDIRCrawHits/dirctree");
  prt_ch->Add(in);
  int nEvents = prt_ch->GetEntries();
  std::cout << "Entries = " << nEvents << std::endl;

  const int arr_size = 10000;

  int hit_size = 0;
  int Particle_id = 0;

  prt_ch->SetBranchAddress("nhits", &hit_size);

  int mcp_num[arr_size], pixel_id[arr_size];
  Double_t lead_time[arr_size];
  
  prt_ch->SetBranchAddress("mcp_id", &mcp_num);
  prt_ch->SetBranchAddress("pixel_id", &pixel_id);
  prt_ch->SetBranchAddress("hit_time", &lead_time);
  
  int printstep = 2000;
  double time = 0;
 
  int pdg(0), totalf(0),totals(0), ch(0);
  double mean_nph_pip, mean_nph_pip_err, sigma_nph_pip, sigma_nph_pip_err;
  double nhits_cut;

  for (int ievent=0; ievent < nEvents; ievent++)
    {
      prt_ch->GetEntry(ievent);
    
      if(ievent < nEvents/2) Particle_id = 211;
      else Particle_id = 321;

      if(Particle_id != 211) continue;
      int numHits = hit_size;
      int nph_pip = 0;
      
      for(int i=0; i < numHits; i++) 
	{
	  nph_pip++;
	}
	
      h_pip_nph->Fill(nph_pip);
    }

  prt_canvasAdd("nph_fit",800,400);

  TFitResultPtr ptr = h_pip_nph->Fit("gaus","SQ","",0,220);
  mean_nph_pip = ptr->Parameter(1);
  mean_nph_pip_err = ptr->ParError(1);
  sigma_nph_pip = ptr->Parameter(2);
  sigma_nph_pip_err = ptr->ParError(2);
    
  nhits_cut = mean_nph_pip - 3*sigma_nph_pip;
  double nhits_cut_val = nhits_cut; // can be used to reject events with low no. of hits (tail on the left side in photon yield)

  //std::cout << "number of hits cut at " << nhits_cut << std::endl;

  prt_canvasGet("nph_fit")->Update();
  
  TPaveText *pt = new TPaveText();
  pt->AddText(Form("mean = %1.2f #pm %1.2f", mean_nph_pip, mean_nph_pip_err));
  pt->AddText(Form("#sigma = %1.2f #pm %1.2f", sigma_nph_pip, sigma_nph_pip_err));
  pt->AddText(Form("mean - 3#sigma = %1.2f", nhits_cut)); 

  pt->SetX1NDC(0.12);
  pt->SetX2NDC(0.24);
  pt->SetY1NDC(0.65);
  pt->SetY2NDC(0.85);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->SetLineColor(0);
  pt->Draw();

  prt_canvasSave(".",0);
  
 for (int ievent=0; ievent < nEvents; ievent++)
   { 
     prt_ch->GetEntry(ievent);

     if(ievent < nEvents/2) Particle_id = 211;
     else Particle_id = 321;
     pdg = Particle_id;

     //if((pdg == 211 && ievent < 2500) || (pdg == 11 && ievent < 27500)) continue; 
     if((pdg == 211 && ievent < 2500) || (pdg == 321 && ievent < 27500)) continue; 

     int nHits = hit_size;
     if(ievent%printstep==0 && ievent!=0) std::cout << "Event # " << ievent << " # hits "<< nHits << std::endl;

     if(nHits < 5) continue;
     //if(nHits < nhits_cut) continue;
     for(int i=0; i < nHits; i++)
      {
	if(pixel_id[i] < 0) continue;
	ch = mcp_num[i]*256 + pixel_id[i];
	time = lead_time[i] + gRandom->Gaus(0,0.1);
      
	if(pdg==pid){
	  totalf++;
	  hlef[ch]->Fill(time);
	}
	if(pdg==211){
	  totals++;
	  hles[ch]->Fill(time);
	}
      }
   }

 std::cout<<"#1 "<< totalf <<"  #2 "<<totals <<std::endl;
  
 if(totalf>=0 && totals>0) 
   {
    in.ReplaceAll(".root",".pdf.root");
    TFile efile(in,"RECREATE");

    for(int i=0; i<nch; i++){
      hlef[i]->Scale(1/(double)totalf);
      hles[i]->Scale(1/(double)totals);
      
      hlef[i]->SetName(Form("hf_%d",i));
      hles[i]->SetName(Form("hs_%d",i));
      hlef[i]->Write();
      hles[i]->Write();
      
      if(0){
	int nrow=4, ncol=6, p = i/256;
	int np =p%ncol*nrow + p/ncol;
	hles[i]->SetName(Form("mcp%dpix%d",np,i%256));
	prt_canvasAdd(Form("pdf_mcp%dpix%d",np,i%256),800,400);
      	prt_normalize(hlef[i],hles[i]);
      	hlef[i]->SetLineColor(2);
	hles[i]->SetLineColor(4);
	hles[i]->GetXaxis()->SetRangeUser(5, 20);
      	hles[i]->Draw("hist");	
      	hlef[i]->Draw("hist same");
	prt_canvasSave("data/pdfs_7",0,0,1);
      }
    }
    
    TTree *tree_nph = new TTree("nph_pip","nph_pip");
    tree_nph->Branch("nhits_cut",&nhits_cut_val,"nhits_cut/D");
    tree_nph->Fill();
    tree_nph->Write();
    
    efile.Write();
    efile.Close();
   }

}
