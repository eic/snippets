// Macro to superimpose the pull distributions
// Shyam Kumar; INFN Bari, shyam.kumar@ba.infn.it

void compare_Pulls(TString name ="objectname", TString histtitle=""){
	
	gStyle->SetPalette(kRainBow);
	gStyle->SetTitleSize(0.045,"XY");	
	gStyle->SetTitleSize(0.04,"XY");	
	gStyle->SetLabelSize(0.04,"XY");	
	gStyle->SetTitleOffset(1.0,"XY");	
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
	gStyle->SetOptTitle(1);
	gStyle->SetGridColor(kBlack);     
	gStyle->SetGridWidth(2);        
	gStyle->SetGridStyle(2);


	TFile *f1 = TFile::Open("./Pull_distribution_vertex_test_D0.list.root");
	TFile *f2 = TFile::Open("./Pull_distribution_vertex_test_DIS.list.root");

	const int nplots = 1;
	double arr[2];
 
 	TCanvas *c[nplots]; 
	for (int i=0; i<nplots; ++i){
	c[i]= new TCanvas(Form("c%d",i),Form("c%d",i),0,52,1400,1000);
	c[i]->SetGrid();
	c[i]->SetMargin(0.10, 0.03 ,0.12,0.07);
	} 
	
	 c[0]->cd();
	 TH1D *h0 = (TH1D*)f1->Get(Form("%s",name.Data()));
	 TH1D *h1 = (TH1D*)f2->Get(Form("%s",name.Data()));
	 
	 arr[0] = h0->GetMaximum(); 	  arr[1] = h1->GetMaximum(); 	  
	 int size = sizeof(arr) / sizeof(arr[0]);
	 int max_val = *std::max_element(arr, arr + size);
      
      h0->SetLineColor(kBlue);
      h0->SetLineStyle(10);
      h0->SetMaximum(max_val*1.10);
      h0->SetTitle(Form(";%s; Entries (a.u.)",histtitle.Data()));
      h0->SetLineWidth(2);
      h0->Draw("hist");
      c[0]->Modified(); c[0]->Update();
      TPaveStats *stats = (TPaveStats*)	h0->FindObject("stats"); // First need to update to make sure stats is there
      stats->SetTextColor(kBlue);
      stats->SetY2NDC(stats->GetY2NDC());
      stats->SetY1NDC(stats->GetY2NDC()-0.15); 
      double y1 = stats->GetY1NDC();  
    
      h1->SetLineColor(kRed);
      h1->SetLineStyle(9);
      h1->SetLineWidth(2);
      h1->Draw("sames");
      c[0]->Modified(); c[0]->Update();
      stats = (TPaveStats*)	h1->FindObject("stats"); // First need to update to make sure stats is there
      stats->SetTextColor(kRed);
      stats->SetY2NDC(y1);	
	 stats->SetY1NDC(y1-0.15);
	 y1 = stats->GetY1NDC();   
      c[0]->Modified();
	 c[0]->Update();
	 
	 TLegend *lhit = new TLegend(0.18,0.70,0.38,0.92);
      lhit->SetHeader(Form("(10X100) Q2>1 GeV^{2} [25.04.1]"),"C");
      lhit->SetTextSize(0.035);
      lhit->SetBorderSize(0);
      lhit->AddEntry(h0,"D0 Sample");
      lhit->AddEntry(h1,"DIS Sample");
      lhit->Draw();
      c[0]->SaveAs(Form("%s_%s.png",name.Data(),histtitle.Data()));

}   



















