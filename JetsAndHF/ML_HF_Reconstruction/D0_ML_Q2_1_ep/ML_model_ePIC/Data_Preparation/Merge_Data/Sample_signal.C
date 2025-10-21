// Sampling of Signals after Preselection cuts
// Shyam Kumar, INFN Bari, Italy; shyam055119@gmail.com
Double_t studentTPeak(Double_t *x, Double_t *par);
Double_t pol2Fit(Double_t *x, Double_t *par);
void Sample_signal(double ymin = -1.0, double ymax = 1.0, double ptmin = 1.0, double ptmax = 2.0, double nEvents_D0Sample = 1.747e+9,double nTotalEvents = 4.72162e+9){

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
   
TFile *fout = new TFile(Form("merged_signal_y_D0_%1.1f_%1.1f_pt_D0_%1.1f_%1.1f.root",ymin,ymax,ptmin,ptmax),"recreate");   

TCanvas *c1 = new TCanvas("c1","c1",0,52,1400,1000);
c1->SetGrid();
c1->SetMargin(0.10, 0.03 ,0.12,0.07);
double min_mass = 1.6, max_mass = 2.5;
Int_t nTotalPar = 4;
// Signal from D0 sample after preselection
TFile *f_sig_D0sample = TFile::Open("../Filtered_D0Sample/SignalD0.root");	
TTree *t_sig_D0sample = (TTree*)f_sig_D0sample->Get("treeMLSig");
t_sig_D0sample->Draw("mass_D0>>h_D0Sample(200,1.6,2.5)",Form("(y_D0>%f && y_D0<%f) && (pt_D0>%f && pt_D0<%f)",ymin,ymax,ptmin,ptmax),"goff");
TH1D *h0_D0Sample = (TH1D*)gDirectory->Get("h_D0Sample"); // access dynamically created histogram
h0_D0Sample->SetMarkerColor(kGreen+2);

// Simulated events correspond 1.747e+9; but we need to 10fb^-1 = 4.72162e+9 
Int_t nexpected_D0 = (nTotalEvents/nEvents_D0Sample)*h0_D0Sample->GetEntries(); 


TF1 *fitFunc = new TF1("fitFunc", studentTPeak, min_mass, max_mass, nTotalPar);
fitFunc->SetParNames("Norm", "Mean", "Width", "nDOF");

c1->cd();
c1->SetLogy();
h0_D0Sample->SetMarkerColor(kBlack);
h0_D0Sample->SetLineColor(kBlack);
h0_D0Sample->SetMarkerStyle(20);
h0_D0Sample->SetMarkerSize(2.0);
h0_D0Sample->SetTitle(Form("(y_D0>%1.1f && y_D0<%1.1f) && (pt_D0>%1.1f && pt_D0<%1.1f);Invmass m_{D^{0}} (GeV/c^{2}); Entries (a.u.)",ymin,ymax,ptmin,ptmax));
h0_D0Sample->SetTitleOffset(1.10,"XYZ");
h0_D0Sample->SetMaximum(1.0e+6);
h0_D0Sample->Draw("EP");

fitFunc->SetParameters(100, 1.865, 0.01, 5);
fitFunc->SetParLimits(3, 1, 10); // nDOF reasonable range
fitFunc->SetParLimits(1, 1.84, 1.89); // Mean reasonable range
fitFunc->SetParLimits(2, 0.001, 0.1); // Sigma reasonable range
h0_D0Sample->SetMaximum(h0_D0Sample->GetMaximum()*2.0);
// Perform fit
h0_D0Sample->Fit(fitFunc, "NR");
h0_D0Sample->Fit(fitFunc, "R");

TH1D *h0_expect_D0Sample = (TH1D*)h0_D0Sample->Clone("hData_D0");
h0_expect_D0Sample->Reset();

for (int i =0; i<nexpected_D0; ++i){

double random = fitFunc->GetRandom(); // resampling
h0_expect_D0Sample->Fill(random);
}
h0_expect_D0Sample->SetMarkerColor(kGreen+2);
h0_expect_D0Sample->SetLineColor(kGreen+2);
h0_expect_D0Sample->Draw("EPsames");
c1->Modified(); c1->Update(); // First need to update to make sure stats is there
TPaveStats *stats = (TPaveStats*)h0_expect_D0Sample->FindObject("stats");
stats->SetX1NDC(stats->GetX1NDC()-0.23); 
stats->SetX2NDC(stats->GetX2NDC()-0.23);
stats->SetTextColor(kGreen+2);	
c1->Modified();  c1->Update();

TLegend *lsample = new TLegend(0.13,0.80,0.33,0.92);
lsample->SetTextSize(0.03);
lsample->SetBorderSize(0);
lsample->AddEntry(h0_D0Sample,"Signal_D0_Sample");
lsample->AddEntry(h0_expect_D0Sample,"Expected_D0_DIS");
lsample->Draw();
c1->SaveAs(Form("Output/signal_y_D0_%1.1f_%1.1f_pt_D0_%1.1f_%1.1f.png",ymin,ymax,ptmin,ptmax));

fout->cd();
h0_expect_D0Sample->SetName("hData_D0");
h0_expect_D0Sample->Write();
fout->Close();


}

// Student's t-distribution signal
Double_t studentTPeak(Double_t *x, Double_t *par) {
    Double_t norm = par[0];    // Height
    Double_t mean = par[1];    // Mean
    Double_t width = par[2];   // Width
    Double_t n = par[3];       // Degrees of freedom

    Double_t t = (x[0] - mean) / width;
    Double_t val = norm * TMath::Gamma((n+1)/2.0) /
                   (TMath::Sqrt(n*TMath::Pi()) * width * TMath::Gamma(n/2.0)) *
                   TMath::Power(1 + t*t/n, -(n+1)/2.0);
    return val;
}


// Pure 2nd-order polynomial (pol2) fit function
Double_t pol2Fit(Double_t *x, Double_t *par) {
    Double_t a = par[0]; // Constant term
    Double_t b = par[1]; // Linear term
    Double_t c = par[2]; // Quadratic term

    return a + b * x[0] + c * x[0] * x[0];
}



