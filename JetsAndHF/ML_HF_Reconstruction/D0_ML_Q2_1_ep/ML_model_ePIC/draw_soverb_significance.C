// Code to draw S/B and Significance as a function of BDT Threshold
// Shyam Kumar; INFN Bari; shyam055119@gmail.com; shyam.kumar@ba.infn.it 


#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>	

void draw_soverb_significance(float ymin = -1.0, float ymax = 1.0, float ptmin=1.0, float ptmax = 2.0){
	
	
  // Styles
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleOffset(1.2,"X");gStyle->SetTitleOffset(1.2,"Y");
   gStyle->SetTitleSize(.04,"X");gStyle->SetTitleSize(.04,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetMarkerSize(1.3);
   gStyle->SetStatW(0.12);
   gStyle->SetStatH(0.12);
   gStyle->SetStatX(0.95);
   gStyle->SetStatY(0.95);
	
   // Read the mass fit results without ML 
   TFile *file = TFile::Open(Form("./Data_Preparation/Merge_Data/DzeroMassPlots_Student/DzeroInvMassCan_pTbin_%1.1f_%1.1f_%1.1f_%1.1f.root",ymin,ymax,ptmin,ptmax));
   TH1F *hFitResult = (TH1F*)file->Get("hFitResult");
   double signal = hFitResult->GetBinContent(1);
   double background = hFitResult->GetBinContent(2);

   TCanvas *c = new TCanvas("c","c",1200,1000);


  TMultiGraph *mg = new TMultiGraph("mg",";BDT Threshold; Results");

   // Read the BDT efficiencies for signal and background 
   TGraph *greff_signal = new TGraph(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/bdt_efficiency_signal_vs_threshold.txt",ymin,ymax,ptmin,ptmax),"%lg %lg");
   // Read the BDT efficiencies for signal and background 
   TGraph *greff_bkg = new TGraph(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/bdt_efficiency_bkg_vs_threshold.txt",ymin,ymax,ptmin,ptmax),"%lg %lg");

	greff_signal->SetMarkerColor(kRed);
    greff_signal->SetLineColor(kRed);
    greff_signal->SetLineWidth(2);

    greff_bkg->SetMarkerColor(kBlack);
    greff_bkg->SetLineColor(kBlack);
    greff_bkg->SetLineWidth(2);

    TGraph *gr_SoverB = (TGraph*)greff_signal->Clone();
    TGraph *gr_Significance= (TGraph*)greff_signal->Clone();

    gr_SoverB->SetMarkerColor(kGreen+2);
    gr_SoverB->SetLineColor(kGreen+2);
    gr_SoverB->SetLineWidth(2);

    gr_Significance->SetMarkerColor(kBlue);
    gr_Significance->SetLineColor(kBlue);
    gr_Significance->SetLineWidth(2);

    // Calculation of S/B with BDT Cut
    for (int i=0; i<gr_SoverB->GetN(); ++i){
     double threshold, eff_signal, eff_kg;
     greff_signal->GetPoint(i,threshold,eff_signal);
     greff_bkg->GetPoint(i,threshold,eff_kg);
     
     double signal_witheff = signal*eff_signal;
     double bkg_witheff = background*eff_kg;

     if (bkg_witheff>0) gr_SoverB->SetPoint(i,threshold,signal_witheff/bkg_witheff);
     if (bkg_witheff>0) gr_Significance->SetPoint(i,threshold,signal_witheff/sqrt(signal_witheff+bkg_witheff));
   
    }
    TLegend *leg = new TLegend(.35,.20,.55,.40);
    leg->SetFillColor(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->AddEntry(greff_signal,"Signal Efficiency");
    leg->AddEntry(greff_bkg,"Background Efficiency");
    leg->AddEntry(gr_SoverB,"Signal/Background");
    leg->AddEntry(gr_Significance,"Significance");

    
    mg->Add(greff_signal);
    mg->Add(greff_bkg);
    mg->Add(gr_SoverB);
    mg->Add(gr_Significance);

    c->cd();
    c->SetGrid();
    c->SetTicks();
   // c->SetLogy();
    mg->Draw("ACP");
    leg->DrawClone("Same");
    c->SaveAs(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/final_results_sb_significance.png",ymin,ymax,ptmin,ptmax));
    c->SaveAs(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/final_results_sb_significance.root",ymin,ymax,ptmin,ptmax));

}
