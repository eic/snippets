// Code to create the optimisation plot
// Shyam Kumar; INFN Bari; shyam055119@gmail.com

void Plot_D0()
{
gStyle->SetPalette(kRainBow);
gStyle->SetTitleSize(0.045,"XY");	
gStyle->SetTitleSize(0.04,"t");	
gStyle->SetLabelSize(0.045,"XY");	
gStyle->SetTitleOffset(1.10,"XY");	
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
gStyle->SetOptTitle(1);
gStyle->SetGridColor(kCyan);     
//gStyle->SetGridWidth(2);        
//gStyle->SetGridStyle(2); 

const Int_t nPtbins = 5;
const Int_t nYbins = 4;
double pt[nPtbins]={0.,1.,2.,5.,10.};  
double y[nYbins]={-3.0,-1.0,1.0,3.0};
    
gSystem->Exec("rm -rf Output && mkdir Output");

TCanvas *c1 = new TCanvas("c1","c1",0,52,1400,1000);
//c1->SetGrid();
c1->SetMargin(0.10, 0.12 ,0.1,0.06);

TFile *f = TFile::Open("test.root");
// Project 2D histograms
TH3D *h3D_Sig_All = (TH3D*)f->Get("h3InvMass_signal_all");
c1->cd();
h3D_Sig_All->GetZaxis()->SetRangeUser(1.6,2.0);
TH2D* h3D_Sig_All_XY = (TH2D*)h3D_Sig_All->Project3D("yx");
h3D_Sig_All_XY->Draw("colz");
c1->SaveAs("signal_all.png");
c1->Clear();

TH3D *h3D_Sig_DCA = (TH3D*)f->Get("h3InvMass_signal_DCA");
c1->cd();
h3D_Sig_DCA->GetZaxis()->SetRangeUser(1.6,2.0);
TH2D* h3D_Sig_DCA_XY = (TH2D*)h3D_Sig_DCA->Project3D("yx");
h3D_Sig_DCA_XY->Draw("colz");
c1->SaveAs("signal_DCA.png");
c1->Clear();

TH3D *h3D_Bkg_All = (TH3D*)f->Get("h3InvMass_bkg_all");
c1->cd();
h3D_Bkg_All->GetZaxis()->SetRangeUser(1.6,2.0);
TH2D* h3D_Bkg_All_XY = (TH2D*)h3D_Bkg_All->Project3D("yx");
h3D_Bkg_All_XY->Draw("colz");
c1->SaveAs("bkg_all.png");
c1->Clear();

TH3D *h3D_Bkg_DCA = (TH3D*)f->Get("h3InvMass_bkg_DCA");
c1->cd();
h3D_Bkg_DCA->GetZaxis()->SetRangeUser(1.6,2.0);
TH2D* h3D_Bkg_DCA_XY = (TH2D*)h3D_Bkg_DCA->Project3D("yx");
h3D_Bkg_DCA_XY->Draw("colz");
c1->SaveAs("bkg_dca.png");
c1->Clear();

// Project 1D histograms
TLegend* legend;
for (int i =0; i < nYbins-1; ++i){
	
	int bin_eta_min = h3D_Sig_All->GetYaxis()->FindBin(y[i]+0.0001); 
	int bin_eta_max = h3D_Sig_All->GetYaxis()->FindBin(y[i+1]-0.0001); 
   
   for (int j =0; j < nPtbins-1; ++j){
   
   c1->Clear();	
   int bin_pt_min = h3D_Sig_All->GetXaxis()->FindBin(pt[j]+0.0001); 
   int bin_pt_max = h3D_Sig_All->GetXaxis()->FindBin(pt[j+1]-0.0001);  
   
   TH1D *hSig_all = (TH1D*) h3D_Sig_All->ProjectionZ(Form("hSig_all_eta_%1.1f_%1.1f_pt_%1.1f_%1.1f",y[i],y[i+1],pt[j],pt[j+1]), bin_pt_min, bin_pt_max, bin_eta_min, bin_eta_max,"o");
   hSig_all->SetTitle(Form("hSig_all_eta_%1.1f_%1.1f_pt_%1.1f_%1.1fGeV/c",y[i],y[i+1],pt[j],pt[j+1]));
   hSig_all->SetLineColor(kBlack);
   
   TH1D *hSig_DCA = (TH1D*) h3D_Sig_DCA->ProjectionZ(Form("hSig_dca_eta_%1.1f_%1.1f_pt_%1.1f_%1.1f",y[i],y[i+1],pt[j],pt[j+1]), bin_pt_min, bin_pt_max, bin_eta_min, bin_eta_max,"o");
   hSig_DCA->SetTitle(Form("hSig_dca_eta_%1.1f_%1.1f_pt_%1.1f_%1.1fGeV/c",y[i],y[i+1],pt[j],pt[j+1]));
   hSig_DCA->SetLineColor(kRed);
   c1->cd();
   hSig_all->Draw("EP");
   hSig_DCA->Draw("EP-same");
   legend = new TLegend(0.2, 0.75, 0.4, 0.9); 
   legend->SetBorderSize(0);   
   legend->SetFillColor(0);    
   legend->SetTextSize(0.03);  
   legend->AddEntry(hSig_all, "All", "l"); 
   legend->AddEntry(hSig_DCA, "Topo Cuts", "l"); 
    // Draw the legend
    legend->Draw();
    
   c1->SaveAs(Form("Output/hSig_all_eta_%1.1f_%1.1f_pt_%1.1f_%1.1f.png",y[i],y[i+1],pt[j],pt[j+1]));   
  	
}	
	
}





}
