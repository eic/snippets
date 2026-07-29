void Plot_Performance(){

gStyle->SetPalette(kRainBow);
gStyle->SetTitleSize(0.045,"XY");	
gStyle->SetTitleSize(0.05,"t");	
gStyle->SetLabelSize(0.045,"XY");	
gStyle->SetTitleOffset(1.0,"XY");	
gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetOptTitle(1);
gStyle->SetGridColor(kBlack);  

TCanvas *c1 = new TCanvas("c1","c1",0,52,1400,1000);
c1->SetGrid();
c1->SetTicks();
c1->SetLogy();
c1->SetMargin(0.10, 0.05 ,0.1,0.06);

TGraph *plot_signal = new TGraph("./results_ml.txt","%lg %lg %*lg %*lg %*lg %*lg"); // Signal
TGraph *plot_bkg = new TGraph("./results_ml.txt","%lg %*lg %lg %*lg %*lg %*lg"); // Bkg
TGraph *plot_sbyb = new TGraph("./results_ml.txt","%lg %*lg %*lg %lg %*lg %*lg"); // S/B
TGraph *plot_signif = new TGraph("./results_ml.txt","%lg %*lg %*lg %*lg %lg %*lg"); // Significance
TGraph *plot_sigma = new TGraph("./results_ml.txt","%lg %*lg %*lg %*lg %*lg %lg"); // Sigma

TMultiGraph *mg = new TMultiGraph("mg",";BDT Theshold; Estimate");
TLegend *leg = new TLegend(.70,.20,.90,.40);
leg->SetFillColor(0);
leg->SetFillColor(0);
leg->SetTextSize(0.035);
leg->SetBorderSize(0);

plot_signal->SetMarkerColor(kBlue);
plot_signal->SetLineColor(kBlue);
plot_signal->SetMarkerSize(1.2);
plot_signal->SetMarkerStyle(20);
mg->Add(plot_signal);

plot_bkg->SetMarkerColor(kRed);
plot_bkg->SetLineColor(kRed);
plot_bkg->SetMarkerSize(1.2);
plot_bkg->SetMarkerStyle(21);
mg->Add(plot_bkg);

plot_sbyb->SetMarkerColor(kBlack);
plot_sbyb->SetLineColor(kBlack);
plot_sbyb->SetMarkerSize(1.2);
plot_sbyb->SetMarkerStyle(34);
mg->Add(plot_sbyb);

plot_signif->SetMarkerColor(kMagenta);
plot_signif->SetLineColor(kMagenta);
plot_signif->SetMarkerSize(1.2);
plot_signif->SetMarkerStyle(33);
mg->Add(plot_signif);

plot_sigma->SetMarkerColor(kGreen+2);
plot_sigma->SetLineColor(kGreen+2);
plot_sigma->SetMarkerSize(1.2);
plot_sigma->SetMarkerStyle(47);
mg->Add(plot_sigma);
mg->Draw("AP");

// Draw the Legend

leg->AddEntry(plot_signal,"Signal");
leg->AddEntry(plot_bkg,"Background");
leg->AddEntry(plot_sbyb,"Signal/Background");
leg->AddEntry(plot_signif,"Significance");
leg->AddEntry(plot_sigma,"Sigma (mass peak)");
leg->DrawClone("Same");

c1->SaveAs("Performance_BDT.png");
}
