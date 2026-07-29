// Macro to read corryverkan outpufile event by event basis
// Shyam Kumar; INFN Bari, shyam.kumar@ba.infn.it

void Compare_Svtx(TString name="x"){
	
	gStyle->SetPalette(kRainBow);
	gStyle->SetTitleSize(0.045,"XY");	
	gStyle->SetTitleSize(0.04,"XY");	
	gStyle->SetLabelSize(0.04,"XY");	
	gStyle->SetTitleOffset(1.0,"XY");	
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
	gStyle->SetOptTitle(0);
	gStyle->SetGridColor(kBlack);     
	gStyle->SetGridWidth(2);        
	gStyle->SetGridStyle(2);


	TFile *f = TFile::Open("test.root");

	const int nplots = 1;

 	TCanvas *c[nplots]; 
	for (int i=0; i<nplots; ++i){
	c[i]= new TCanvas(Form("c%d",i),Form("c%d",i),0,52,1400,1000);
	c[i]->SetGrid();
	c[i]->SetMargin(0.10, 0.03 ,0.12,0.07);
	} 
	
	 c[0]->cd();
	 TH1D *hvtx_x_fit = (TH1D*)f->Get(Form("hRes_SV%s_Helixfit",name.Data()));
	 TH1D *hvtx_x_centre = (TH1D*)f->Get(Form("hRes_SV%s_HelixCalculation",name.Data()));
      
      hvtx_x_fit->SetLineColor(kRed);
      hvtx_x_fit->SetLineStyle(10);
      hvtx_x_fit->SetMaximum(hvtx_x_fit->GetMaximum()*1.20);
      hvtx_x_fit->SetLineWidth(2);
      hvtx_x_fit->Draw("hist");
      c[0]->Modified(); c[0]->Update();
      TPaveStats *stats = (TPaveStats*)	hvtx_x_fit->FindObject("stats"); // First need to update to make sure stats is there
      stats->SetTextColor(kRed);
      stats->SetY2NDC(stats->GetY2NDC());
      stats->SetY1NDC(stats->GetY2NDC()-0.15); 
      double y1 = stats->GetY1NDC();  
    
      hvtx_x_centre->SetLineColor(kBlue);
      hvtx_x_centre->SetLineStyle(9);
      hvtx_x_centre->SetLineWidth(2);
      hvtx_x_centre->Draw("sames");
      c[0]->Modified(); c[0]->Update();
      stats = (TPaveStats*)	hvtx_x_centre->FindObject("stats"); // First need to update to make sure stats is there
      stats->SetTextColor(kBlue);
      stats->SetY2NDC(y1);	
	    stats->SetY1NDC(y1-0.15);
	    y1 = stats->GetY1NDC();   
      c[0]->Modified();
	    c[0]->Update();
	 
    
      TLegend *lhit = new TLegend(0.15,0.70,0.35,0.92);
      lhit->SetHeader("Secondary Vertex (D^{0})","C");
      lhit->SetTextSize(0.035);
      lhit->SetBorderSize(0);
      lhit->AddEntry(hvtx_x_fit,"Chi2 Minimization");
      lhit->AddEntry(hvtx_x_centre,"Mid Point (Analytical)");
      lhit->Draw();
      c[0]->SaveAs(Form("cmp_secondary_vertex_%s_D0.png",name.Data()));

}   



















