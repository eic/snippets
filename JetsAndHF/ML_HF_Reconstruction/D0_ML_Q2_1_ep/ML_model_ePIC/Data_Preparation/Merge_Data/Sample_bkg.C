// Sampling of Bkg after Preselection cuts
// Shyam Kumar, INFN Bari, Italy; shyam055119@gmail.com
Double_t pol2Fit(Double_t *x, Double_t *par);
void Sample_bkg(double ymin = -1.0, double ymax = 1.0, double ptmin = 1.0, double ptmax = 2.0, double nTotalEvents = 4.72162e+9){

 // Styles
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleOffset(1.2,"X");gStyle->SetTitleOffset(1.6,"Y");
   gStyle->SetTitleSize(.04,"X");gStyle->SetTitleSize(.04,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(1);
   gStyle->SetMarkerSize(1.3);
   gStyle->SetStatW(0.12);
   gStyle->SetStatH(0.12);
   gStyle->SetStatX(0.95);
   gStyle->SetStatY(0.95);
   
TFile *fout = new TFile(Form("merged_bkg_y_D0_%1.1f_%1.1f_pt_D0_%1.1f_%1.1f.root",ymin,ymax,ptmin,ptmax),"recreate");   

TCanvas *c1 = new TCanvas("c1","c1",0,52,1400,1000);
c1->SetGrid();
c1->SetMargin(0.10, 0.03 ,0.12,0.07);
TString std_cuts ="(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1";
// Reading Event histgram for DIS sample		
TFile *file_DISsample = TFile::Open("../../../DIS_Sample/test.root");
TH1D *hEvents_DISsample = (TH1D*)file_DISsample->Get("hEventStat");
printf("Number of Events DIS: = %f \n",hEvents_DISsample->GetBinContent(1));
// Luminosity correspond to 6.5M events 
double scaling_factor = nTotalEvents/hEvents_DISsample->GetBinContent(1); // applying scaling factor
double min_mass = 1.6, max_mass = 2.5;
Int_t nTotalPar = 3;
// Background from DIS sample after preselection
TFile *f_bkg_DISsample = TFile::Open("../Filtered_DISSample/BkgD0.root");	
TTree *t_bkg_DISsample = (TTree*)f_bkg_DISsample->Get("treeMLBkg");
if(ymax > 0.5){
t_bkg_DISsample->Draw("mass_D0>>h_Bkg_DISSample(200,1.6,2.5)",Form("(y_D0>%f && y_D0<%f) && (pt_D0>%f && pt_D0<%f)",ymin,ymax,ptmin,ptmax),"goff");
}
else{
t_bkg_DISsample->Draw("mass_D0>>h_Bkg_DISSample(200,1.6,2.5)",Form("(y_D0>%f && y_D0<%f) && (pt_D0>%f && pt_D0<%f) && %s",ymin,ymax,ptmin,ptmax,std_cuts.Data()),"goff");
}
TF1 *fitFunc = new TF1("fitFunc", pol2Fit, min_mass, max_mass, nTotalPar);
fitFunc->SetParNames("p0", "p1", "p2");

c1->cd();
c1->SetLogy();
TH1D *h0_DIS = (TH1D*)gDirectory->Get("h_Bkg_DISSample"); // access dynamically created histogram
h0_DIS->SetMarkerColor(kBlack);
h0_DIS->SetLineColor(kBlack);
h0_DIS->SetMarkerStyle(20);
h0_DIS->SetMarkerSize(2.0);
h0_DIS->SetTitle(Form("(y_D0>%1.1f && y_D0<%1.1f) && (pt_D0>%1.1f && pt_D0<%1.1f);Invmass m_{D^{0}} (GeV/c^{2}); Entries (a.u.)",ymin,ymax,ptmin,ptmax));
h0_DIS->SetTitleOffset(1.10,"XYZ");
h0_DIS->SetMaximum(1.0e+6);
h0_DIS->Draw("EP");

fitFunc->SetParameters(1., 1., 1.);
fitFunc->FixParameter(2, 0.);
// Perform fit
h0_DIS->Fit(fitFunc, "NR");
h0_DIS->Fit(fitFunc, "R");

TH1D *h0_DIS_Scaled = (TH1D*)h0_DIS->Clone("h_Bkg_DISSample_100");
h0_DIS_Scaled->SetMarkerColor(kGreen+2);
h0_DIS_Scaled->SetLineColor(kGreen+2);
h0_DIS_Scaled->SetMarkerStyle(47);
h0_DIS_Scaled->SetMarkerSize(2.0);
Int_t total_entries = scaling_factor*h0_DIS->GetEntries();
h0_DIS_Scaled->Reset();
for (int i =0; i<total_entries; ++i){

double random = fitFunc->GetRandom(); // resampling
h0_DIS_Scaled->Fill(random);
}
h0_DIS_Scaled->SetMaximum(h0_DIS_Scaled->GetMaximum()*2.0);
h0_DIS_Scaled->Draw("EPsames");
c1->Modified(); c1->Update(); // First need to update to make sure stats is there
TPaveStats *stats = (TPaveStats*)h0_DIS_Scaled->FindObject("stats");
stats->SetX1NDC(stats->GetX1NDC()-0.23); 
stats->SetX2NDC(stats->GetX2NDC()-0.23);
stats->SetTextColor(kGreen+2);	
c1->Modified();  c1->Update();

TLegend *lsample = new TLegend(0.13,0.80,0.33,0.92);
lsample->SetTextSize(0.03);
lsample->SetBorderSize(0);
lsample->AddEntry(h0_DIS,"Bkg_DIS_Sample");
lsample->AddEntry(h0_DIS_Scaled,"Bkg_DIS_Sample Resampled");
lsample->Draw();
c1->SaveAs(Form("Output/merged_bkg_y_D0_%1.1f_%1.1f_pt_D0_%1.1f_%1.1f.png",ymin,ymax,ptmin,ptmax));

fout->cd();
h0_DIS_Scaled->SetName("hData_D0");
h0_DIS_Scaled->Write();
fout->Close();


}

// Pure 2nd-order polynomial (pol2) fit function
Double_t pol2Fit(Double_t *x, Double_t *par) {
    Double_t a = par[0]; // Constant term
    Double_t b = par[1]; // Linear term
    Double_t c = par[2]; // Quadratic term

    return a + b * x[0] + c * x[0] * x[0];
}



