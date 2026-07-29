// Code to draw S/B and Significance as a function of BDT Threshold
// Shyam Kumar; INFN Bari; shyam055119@gmail.com; shyam.kumar@ba.infn.it 


#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>	
void DrawFineGrid(TGraph* gr);
void Superimpose_BDT_efficiencies(float ymin = -1.0, float ymax = 1.0, float ptmin=1.0, float ptmax = 2.0){
	
	
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
	

 TCanvas *c = new TCanvas("c","c",1200,1000);
 c->SetMargin(0.12,0.02,0.1,0.02);
 TMultiGraph *mg = new TMultiGraph("mg",";BDT Threshold; Efficiencies");

   // Read the BDT efficiencies for signal and background 
   TGraph *greff_signal = new TGraph(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/bdt_efficiency_signal_vs_threshold.txt",ymin,ymax,ptmin,ptmax),"%lg %lg");
   // Read the BDT efficiencies for signal and background 
   TGraph *greff_bkg = new TGraph(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/bdt_efficiency_bkg_vs_threshold.txt",ymin,ymax,ptmin,ptmax),"%lg %lg");

    greff_signal->SetMarkerColor(kRed);
    greff_signal->SetLineColor(kRed);
    greff_signal->SetLineWidth(2);


    greff_bkg->SetMarkerColor(kBlue);
    greff_bkg->SetLineColor(kBlue);
    greff_bkg->SetLineWidth(2);

    TLegend *leg = new TLegend(.70,.80,.90,.95);
    leg->SetHeader(Form("%1.0f < y < %1.0f, %1.0f < p_{T} < %1.0f  (GeV/c)",ymin,ymax,ptmin,ptmax),"C");
    leg->SetFillColor(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->AddEntry(greff_signal,"Signal");
    leg->AddEntry(greff_bkg,"Background");
    
    c->cd();
    c->SetGrid();
    c->SetTicks();
    greff_signal->SetTitle(";BDT Threshold; Efficiencies");
    greff_signal->Draw("ACP");
    greff_bkg->Draw("same");
    DrawFineGrid(greff_signal);
 
    leg->DrawClone("Same");
    c->SaveAs(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/BDT_efficiencies_plot_graph.png",ymin,ymax,ptmin,ptmax));
    c->SaveAs(Form("./ML_Output_Optuna_%1.1f_%1.1f_%1.1f_%1.1f/BDT_efficiencies_plot_graph.root",ymin,ymax,ptmin,ptmax));

}

#include "TGraph.h"
#include "TLine.h"
#include "TStyle.h"
#include <cmath>

void DrawFineGrid(TGraph* gr)
{
    if (!gr || !gPad) return;
    if (!gr->GetHistogram()) gr->Draw("APL"); // ensure frame exists
    gPad->Update();

    // Fixed spacing for grid lines
    const double xStep = 0.02;
    const double yStep = 0.02;

    double xMin = gPad->GetUxmin(), xMax = gPad->GetUxmax();
    double yMin = gPad->GetUymin(), yMax = gPad->GetUymax();

    auto L = [&](double x1, double y1, double x2, double y2) {
        auto l = new TLine(x1, y1, x2, y2);
        l->SetLineStyle(3);       // dashed
        l->SetLineColor(kGray+1); // light gray
        l->SetLineWidth(1);
        l->Draw("same");
    };

    double startX = std::ceil(xMin / xStep);
    for (double k = startX; ; ++k) {
        double xx = k * xStep;
        if (xx > xMax + 1e-12) break;
        L(xx, yMin, xx, yMax);
    }

   double startY = std::ceil(yMin / yStep);
    for (double k = startY; ; ++k) {
        double yy = k * yStep;
        if (yy > yMax + 1e-12) break;
        L(xMin, yy, xMax, yy);
    }

    gPad->RedrawAxis();
}

