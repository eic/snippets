// ============================================================================
//! \file    MakeJetValidationPlots.C
//! \authors Brian Page (bpage@bnl.gov),
//!          adapted by Derek Anderson (derek.murphy.anderson@protonmail.com)
// ----------------------------------------------------------------------------
//! \brief Adaption of the Jet Benchmark to run standalone for validation.
//!   This macro generates plots from the output file produced by
//!   `MakeJetValidationHists.C`.
//!
//! \usage In eic-shell:
//!     root -b -q MakeJetValidationPlots.C'(<output path>, \
//!                                          <output suffix>, \
//!                                          <input file>)'
// ============================================================================

#include <edm4eic/EDM4eicVersion.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <string>

#include "fmt/color.h"
#include "fmt/core.h"

///! Default output file path.
const std::string DefaultOutPath = ".";

///! Default output file suffix.
const std::string DefaultOutSuffix = "files26071.py8ncdis10x100q100t1000";

///! Default input file.
const std::string DefaultInFile = "hists.files26071.py8ncdis10x100q100t1000.root";

// ----------------------------------------------------------------------------
// Macro body
// ----------------------------------------------------------------------------
/*! Process input ROOT file to generate a set of reconstructed,
 *! generated jet plots and save them as PNG files.
 *!
 *! \param[out] results_path   Location to save PNGs to
 *! \param[out] results_suffix Suffix to append to PNGs
 *! \param[in]  input_file     Input ROOT file to use
 */
int MakeJetValidationPlots(
  const std::string& results_path = DefaultOutPath,
  const std::string& results_suffix = DefaultOutSuffix,
  const std::string& input_file = DefaultInFile
) {

  const bool PRINT = true;

  const int seabornRed = TColor::GetColor(213, 94, 0);

  // Seaborn Green: #009E73 -> (0, 158, 115)
  const int seabornGreen = TColor::GetColor(0, 158, 115);

  // Seaborn Blue: #56B4E9 -> (86, 180, 233)
  const int seabornBlue = TColor::GetColor(100, 149, 237);

  // Get input histograms
  TFile* in_file = new TFile(input_file.c_str(), "read");
  if (!in_file) {
    std::cerr << "PANIC: Couldn't open input file:\n  " << input_file << std::endl;
    return 1;
  }

  // Check if area histograms exist
  //   --> Using jet EDM if they do!
  bool useNewEDM = false;
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8,9,0)
  {
    TH1D* recoAreaHistTest = (TH1D*) in_file->Get("recoChargedJetAreaECut");
    if (recoAreaHistTest) {
      useNewEDM = true;
    }
  }
#endif
  if (PRINT) {
    std::cout << "INFO: Using new EDM? " << useNewEDM << std::endl;
  }

  // Reco
  TH1D *numRecoChargedJetsECutHist = (TH1D*) in_file->Get("numRecoChargedJetsECut");
  TH1D *recoChargedJetEHist = (TH1D*) in_file->Get("recoChargedJetE");
  TH1D *recoChargedJetEtaECutHist = (TH1D*) in_file->Get("recoChargedJetEtaECut");
  TH1D *recoChargedJetAreaECutHist = nullptr;
  TH2D *recoChargedJetEvsAreaHist = nullptr;
  if (useNewEDM) {
    recoChargedJetAreaECutHist = (TH1D*) in_file->Get("recoChargedJetAreaECut");
    recoChargedJetEvsAreaHist = (TH2D*) in_file->Get("recoChargedJetEvsArea");
  }
  TH2D *recoChargedJetEvsEtaHist = (TH2D*) in_file->Get("recoChargedJetEvsEta");
  TH2D *recoChargedJetPhiVsEtaECutHist = (TH2D*) in_file->Get("recoChargedJetPhiVsEtaECut");

  TH1D *numRecoChargedJetsECutNoElecHist = (TH1D*) in_file->Get("numRecoChargedJetsECutNoElec");
  TH1D *recoChargedJetENoElecHist = (TH1D*) in_file->Get("recoChargedJetENoElec");
  TH1D *recoChargedJetEtaECutNoElecHist = (TH1D*) in_file->Get("recoChargedJetEtaECutNoElec");
  TH1D *recoChargedJetAreaECutNoElecHist = nullptr;
  TH2D *recoChargedJetEvsAreaNoElecHist = nullptr;
  if (useNewEDM) {
    recoChargedJetAreaECutNoElecHist = (TH1D*) in_file->Get("recoHargedJetAreaECutNoElec");
    recoChargedJetEvsAreaNoElecHist = (TH2D*) in_file->Get("recoChargedJetEvsAreaNoElec");
  }
  TH2D *recoChargedJetEvsEtaNoElecHist = (TH2D*) in_file->Get("recoChargedJetEvsEtaNoElec");
  TH2D *recoChargedJetPhiVsEtaECutNoElecHist = (TH2D*) in_file->Get("recoChargedJetPhiVsEtaECutNoElec");

  TH1D *numRecoChargedJetPartsHist = (TH1D*) in_file->Get("numRecoChargedJetParts");
  TH1D *recoChargedJetPartPHist = (TH1D*) in_file->Get("recoChargedJetPartP");
  TH1D *recoChargedJetPartEtaHist = (TH1D*) in_file->Get("recoChargedJetPartEta");
  TH2D *recoChargedJetPartPvsEtaHist = (TH2D*) in_file->Get("recoChargedJetPartPvsEta");
  TH2D *recoChargedJetPartPhiVsEtaHist = (TH2D*) in_file->Get("recoChargedJetPartPhiVsEta");

  TH1D *numRecoChargedJetPartsNoElecHist = (TH1D*) in_file->Get("numRecoChargedJetPartsNoElec");
  TH1D *recoChargedJetPartPNoElecHist = (TH1D*) in_file->Get("recoChargedJetPartPNoElec");
  TH1D *recoChargedJetPartEtaNoElecHist = (TH1D*) in_file->Get("recoChargedJetPartEtaNoElec");
  TH2D *recoChargedJetPartPvsEtaNoElecHist = (TH2D*) in_file->Get("recoChargedJetPartPvsEtaNoElec");
  TH2D *recoChargedJetPartPhiVsEtaNoElecHist = (TH2D*) in_file->Get("recoChargedJetPartPhiVsEtaNoElec");

  TH1D *recoChargedJetPartPairwiseDeltaRHist = (TH1D*) in_file->Get("recoChargedJetPartPairwiseDeltaRHist");

  // Gen
  TH1D *numGenChargedJetsECutHist = (TH1D*) in_file->Get("numGenChargedJetsECut");
  TH1D *genChargedJetEHist = (TH1D*) in_file->Get("genChargedJetE");
  TH1D *genChargedJetEtaECutHist = (TH1D*) in_file->Get("genChargedJetEtaECut");
  TH1D *genChargedJetAreaECutHist = nullptr;
  TH2D *genChargedJetEvsAreaHist = nullptr;
  if (useNewEDM) {
    genChargedJetAreaECutHist = (TH1D*) in_file->Get("genChargedJetAreaECut");
    genChargedJetEvsAreaHist = (TH2D*) in_file->Get("genChargedJetEvsAreaHist");
  }
  TH2D *genChargedJetEvsEtaHist = (TH2D*) in_file->Get("genChargedJetEvsEta");
  TH2D *genChargedJetPhiVsEtaECutHist = (TH2D*) in_file->Get("genChargedJetPhiVsEtaECut");

  TH1D *numGenChargedJetsECutNoElecHist = (TH1D*) in_file->Get("numGenChargedJetsECutNoElec");
  TH1D *genChargedJetENoElecHist = (TH1D*) in_file->Get("genChargedJetENoElec");
  TH1D *genChargedJetEtaECutNoElecHist = (TH1D*) in_file->Get("genChargedJetEtaECutNoElec");
  TH1D *genChargedJetAreaECutNoElecHist = nullptr;
  TH2D *genChargedJetEvsAreaNoElecHist = nullptr;
  if (useNewEDM) {
    genChargedJetAreaECutNoElecHist = (TH1D*) in_file->Get("genChargedJetAreaECutNoElec");
    genChargedJetEvsAreaNoElecHist = (TH2D*) in_file->Get("genChargedJetEvsAreaNoElec");
  }
  TH2D *genChargedJetEvsEtaNoElecHist = (TH2D*) in_file->Get("genChargedJetEvsEtaNoElec");
  TH2D *genChargedJetPhiVsEtaECutNoElecHist = (TH2D*) in_file->Get("genChargedJetPhiVsEtaECutNoElec");

  TH1D *numGenChargedJetPartsHist = (TH1D*) in_file->Get("numGenChargedJetParts");
  TH1D *genChargedJetPartPHist = (TH1D*) in_file->Get("genChargedJetPartP");
  TH1D *genChargedJetPartEtaHist = (TH1D*) in_file->Get("genChargedJetPartEta");
  TH2D *genChargedJetPartPvsEtaHist = (TH2D*) in_file->Get("genChargedJetPartPvsEta");
  TH2D *genChargedJetPartPhiVsEtaHist = (TH2D*) in_file->Get("genChargedJetPartPhiVsEta");

  TH1D *numGenChargedJetPartsNoElecHist = (TH1D*) in_file->Get("numGenChargedJetPartsNoElec");
  TH1D *genChargedJetPartPNoElecHist = (TH1D*) in_file->Get("genChargedJetPartPNoElec");
  TH1D *genChargedJetPartEtaNoElecHist = (TH1D*) in_file->Get("genChargedJetPartEtaNoElec");
  TH2D *genChargedJetPartPvsEtaNoElecHist = (TH2D*) in_file->Get("genChargedJetPartPvsEtaNoElec");
  TH2D *genChargedJetPartPhiVsEtaNoElecHist = (TH2D*) in_file->Get("genChargedJetPartPhiVsEtaNoElec");

  TH1D *genChargedJetPartPairwiseDeltaRHist = (TH1D*) in_file->Get("genChargedJetPartPairwiseDeltaRHist");

  // Matched
  TH1D *matchJetDeltaRHist = (TH1D*) in_file->Get("matchJetDeltaR");
  TH1D *matchJetDeltaRBackHist = (TH1D*) in_file->Get("matchJetDeltaRBack");
  TH2D *recoVsGenChargedJetEtaHist = (TH2D*) in_file->Get("recoVsGenChargedJetEta");
  TH2D *recoVsGenChargedJetPhiHist = (TH2D*) in_file->Get("recoVsGenChargedJetPhi");
  TH2D *recoVsGenChargedJetAreaHist = nullptr;
  if (useNewEDM) {
    recoVsGenChargedJetAreaHist = (TH2D*) in_file->Get("recoVsGenChargedJetArea");
  }
  TH2D *recoVsGenChargedJetEHist = (TH2D*) in_file->Get("recoVsGenChargedJetE");
  TH2D *recoVsGenChargedJetENoDRHist = (TH2D*) in_file->Get("recoVsGenChargedJetENoDRHist");
  TH2D *recoVsGenChargedJetENoDupHist = (TH2D*) in_file->Get("recoVsGenChargedJetENoDup");

  TH2D *jetResVsEtaHist = (TH2D*) in_file->Get("jetResVsEta");
  TH2D *jetResVsEHist = (TH2D*) in_file->Get("jetResVsE");
  TH2D *jetResVsENegEtaHist = (TH2D*) in_file->Get("jetResVsENegEta");
  TH2D *jetResVsEMidEtaHist = (TH2D*) in_file->Get("jetResVsEMidEta");
  TH2D *jetResVsEPosEtaHist = (TH2D*) in_file->Get("jetResVsEPosEta");

  TH2D *jetResVsENegEtaNoDupHist = (TH2D*) in_file->Get("jetResVsENegEtaNoDup");
  TH2D *jetResVsEMidEtaNoDupHist = (TH2D*) in_file->Get("jetResVsEMidEtaNoDup");
  TH2D *jetResVsEPosEtaNoDupHist = (TH2D*) in_file->Get("jetResVsEPosEtaNoDup");

  gStyle->SetOptStat(0);
  ////////////////////////  Reconstructed Jets Plots  ////////////////////////
  // Reco Number
  TCanvas *c1 = new TCanvas("c1","Number Reco Jets",800,600);
  c1->Clear();
  c1->Divide(1,1);

  c1->cd(1);
  numRecoChargedJetsECutHist->Draw("HIST");
  numRecoChargedJetsECutNoElecHist->SetLineColor(seabornRed);
  numRecoChargedJetsECutNoElecHist->Draw("HISTSAME");
  numRecoChargedJetsECutHist->SetLineWidth(2); // Set line width to 2 (adjust as needed)
numRecoChargedJetsECutNoElecHist->SetLineWidth(2);
  numRecoChargedJetsECutHist->SetTitle("Reconstructed Jets per Event (|eta| < 2.5 && E > 5);Number");

TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
legend1->AddEntry(numRecoChargedJetsECutHist, "With Electrons", "l");
legend1->AddEntry(numRecoChargedJetsECutNoElecHist, "No Electrons", "l");
legend1->Draw();



  gPad->SetLogy();
  if(PRINT) c1->Print((results_path+"/numberRecoJets."+results_suffix+".png").c_str()); // Number of reconstructed jets per event with energy > 5 GeV and Abs(eta) < 2.5
   delete c1;

  // Reco Energy
  TCanvas *c2 = new TCanvas("c2","Reco Jet Energy",800,600);
  c2->Clear();
  c2->Divide(1,1);

  c2->cd(1);
  recoChargedJetEHist->Draw("HIST");
  recoChargedJetENoElecHist->SetLineColor(seabornRed);
  recoChargedJetENoElecHist->Draw("HISTSAME");

  recoChargedJetEHist->SetLineWidth(2);
  recoChargedJetENoElecHist->SetLineWidth(2);
  recoChargedJetEHist->SetTitle("Reconstructed Jet Energy (|eta| < 2.5);Energy [GeV]");

TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
legend2->AddEntry(recoChargedJetEHist, "With Electrons", "l");
legend2->AddEntry(recoChargedJetENoElecHist, "No Electrons", "l");
legend2->Draw();

  gPad->SetLogy();
  if(PRINT) c2->Print((results_path+"/recoJetEnergy."+results_suffix+".png").c_str()); // Energy spectrum of reconstructed jets with Abs(eta) < 2.5

    delete c2;
  // Reco Eta
  TCanvas *c3 = new TCanvas("c3","Reco Jet Eta",800,600);
  c3->Clear();
  c3->Divide(1,1);

  c3->cd(1);
  recoChargedJetEtaECutHist->Draw("HIST");
  recoChargedJetEtaECutNoElecHist->SetLineColor(seabornRed);
  recoChargedJetEtaECutNoElecHist->Draw("HISTSAME");

  recoChargedJetEtaECutHist->SetLineWidth(2);
  recoChargedJetEtaECutNoElecHist->SetLineWidth(2);
  recoChargedJetEtaECutHist->SetTitle("Reconstructed Jet Eta (E > 5);Eta");

//add legend 
TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
legend3->AddEntry(recoChargedJetEtaECutHist, "With Electrons", "l");
legend3->AddEntry(recoChargedJetEtaECutNoElecHist, "No Electrons", "l");
legend3->Draw();

  gPad->SetLogy();
  if(PRINT) c3->Print((results_path+"/recoJetEta."+results_suffix+".png").c_str()); // Eta spectrum of reconstructed jets with energy > 5 GeV
    delete c3;

  // Reco Area
  if(useNewEDM) {
    TCanvas *c3_1 = new TCanvas("c3_1","Reco Jet Area",800,600);
    c3_1->Clear();
    c3_1->Divide(1,1);

    c3_1->cd(1);
    recoChargedJetAreaECutHist->Draw("HIST");
    recoChargedJetAreaECutNoElecHist->SetLineColor(seabornRed);
    recoChargedJetAreaECutNoElecHist->Draw("HISTSAME");

    recoChargedJetAreaECutHist->SetLineWidth(2);
    recoChargedJetAreaECutNoElecHist->SetLineWidth(2);
    recoChargedJetAreaECutHist->SetTitle("Reconstructed Jet Area (E > 5);Area");

    //add legend
    TLegend *legend3_1 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
    legend3_1->AddEntry(recoChargedJetAreaECutHist, "With Electrons", "l");
    legend3_1->AddEntry(recoChargedJetAreaECutNoElecHist, "No Electrons", "l");
    legend3_1->Draw();

    gPad->SetLogy();
    if(PRINT) c3_1->Print((results_path+"/recoJetArea."+results_suffix+".png").c_str()); // Area spectrum of reconstructed jets with energy > 5 GeV
      delete c3_1;
   }

  // Reco E Vs Eta
  TCanvas *c4 = new TCanvas("c4","Reco Jet E Vs Eta",800,600);
  c4->Clear();
  c4->Divide(1,1);

  c4->cd(1);
  recoChargedJetEvsEtaHist->Draw("COLZ");
  recoChargedJetEvsEtaHist->SetTitle("Reconstructed Jet Energy Vs Eta;Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c4->Print((results_path+"/recoJetEnergyvsEta."+results_suffix+".png").c_str()); // Energy vs eta of reconstructed jets

  // Reco E Vs Area
  if (useNewEDM) {
    TCanvas *c4_1 = new TCanvas("c4_1","Reco Jet E Vs Area",800,600);
    c4_1->Clear();
    c4_1->Divide(1,1);

    c4_1->cd(1);
    recoChargedJetEvsAreaHist->Draw("COLZ");
    recoChargedJetEvsAreaHist->SetTitle("Reconstructed Jet Energy Vs Area;Area;Energy [GeV]");
    gPad->SetLogz();
    if(PRINT) c4_1->Print((results_path+"/recoJetEnergyvsArea."+results_suffix+".png").c_str()); // Energy vs area of reconstructed jets
  }

  // Reco Phi Vs Eta
  TCanvas *c5 = new TCanvas("c5","Reco Jet Phi Vs Eta",800,600);
  c5->Clear();
  c5->Divide(1,1);

  c5->cd(1);
  recoChargedJetPhiVsEtaECutHist->Draw("COLZ");
  recoChargedJetPhiVsEtaECutHist->SetTitle("Reconstructed Jet Phi Vs Eta (E > 5);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c5->Print((results_path+"/recoJetPhivsEta."+results_suffix+".png").c_str()); // Phi vs eta of reconstructed jets

  // Num Particles Per Reco Jet
  TCanvas *c6 = new TCanvas("c6","Number Constituents Per Reco Jet",800,600);
  c6->Clear();
  c6->Divide(1,1);

  c6->cd(1);
  numRecoChargedJetPartsHist->Draw("HIST");
  numRecoChargedJetPartsNoElecHist->SetLineColor(seabornRed);
  numRecoChargedJetPartsNoElecHist->Draw("HISTSAME");

  numRecoChargedJetPartsHist->SetLineWidth(2);
  numRecoChargedJetPartsNoElecHist->SetLineWidth(2);
  numRecoChargedJetPartsHist->SetTitle("Number of Constituents Per Reco Jet;Number of Constituents");

TLegend *legend6 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
legend6->AddEntry(numRecoChargedJetPartsHist, "With Electrons", "l");
legend6->AddEntry(numRecoChargedJetPartsNoElecHist, "No Electrons", "l");
legend6->Draw();

  gPad->SetLogy();
  if(PRINT) c6->Print((results_path+"/numConstituentsPerRecoJet."+results_suffix+".png").c_str()); // Number of constituents in reconstructed jets

  // Reco Part Energy
  TCanvas *c7 = new TCanvas("c7","Reco Jet Constituent Momentum",800,600);
  c7->Clear();
  c7->Divide(1,1);

  c7->cd(1);
  recoChargedJetPartPHist->Draw("HIST");
  recoChargedJetPartPNoElecHist->SetLineColor(seabornRed);
  recoChargedJetPartPNoElecHist->Draw("HISTSAME");

  recoChargedJetPartPHist->SetLineWidth(2);
  recoChargedJetPartPNoElecHist->SetLineWidth(2);
  recoChargedJetPartPHist->SetTitle("Reconstructed Jet Constituent Momentum;Momentum [GeV/c]");

  TLegend *legend7 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend7->AddEntry(recoChargedJetPartPHist, "With Electrons", "l");
  legend7->AddEntry(recoChargedJetPartPNoElecHist, "No Electrons", "l");
  legend7->Draw();

  gPad->SetLogy();
  if(PRINT) c7->Print((results_path+"/recoJetConstituentMomentum."+results_suffix+".png").c_str()); // Momentum of reconstructed jet constituents

  // Reco Part Eta
  TCanvas *c8 = new TCanvas("c8","Reco Jet Constituent Eta",800,600);
  c8->Clear();
  c8->Divide(1,1);

  c8->cd(1);
  recoChargedJetPartEtaHist->Draw("HIST");
  recoChargedJetPartEtaNoElecHist->SetLineColor(seabornRed);
  recoChargedJetPartEtaNoElecHist->Draw("HISTSAME");

  recoChargedJetPartEtaHist->SetLineWidth(2);
  recoChargedJetPartEtaNoElecHist->SetLineWidth(2);

  recoChargedJetPartEtaHist->SetTitle("Reconstructed Jet Constituent Eta;Eta");

  TLegend *legend8 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend8->AddEntry(recoChargedJetPartEtaHist, "With Electrons", "l");
  legend8->AddEntry(recoChargedJetPartEtaNoElecHist, "No Electrons", "l");
  legend8->Draw();

  gPad->SetLogy();
  if(PRINT) c8->Print((results_path+"/recoJetConstituentEta."+results_suffix+".png").c_str()); // Eta of reconstructed jet constituents

  // Reco Part P Vs Eta
  TCanvas *c9 = new TCanvas("c9","Reco Jet Constituent Momentum Vs Eta",800,600);
  c9->Clear();
  c9->Divide(1,1);

  c9->cd(1);
  recoChargedJetPartPvsEtaHist->Draw("COLZ");
  recoChargedJetPartPvsEtaHist->SetTitle("Reconstructed Jet Constituent Momentum Vs Eta;Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c9->Print((results_path+"/recoJetConstituentMomentumVsEta."+results_suffix+".png").c_str()); // Momentum vs eta of reconstructed jet constituents

  // Reco Part Phi Vs Eta
  TCanvas *c10 = new TCanvas("c10","Reco Jet Constituent Phi Vs Eta",800,600);
  c10->Clear();
  c10->Divide(1,1);

  c10->cd(1);
  recoChargedJetPartPhiVsEtaHist->Draw("COLZ");
  recoChargedJetPartPhiVsEtaHist->SetTitle("Reconstructed Jet Constituent Phi Vs Eta;Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c10->Print((results_path+"/recoJetConstituentPhiVsEta."+results_suffix+".png").c_str()); // Phi vs eta of reconstructed jet constituents

  // Reco Constituent Pairwise delta R
  TCanvas *c11 = new TCanvas("c11","Reco Jet Constituent Pairwise Delta R",800,600);
  c11->Clear();
  c11->Divide(1,1);

  c11->cd(1);
  recoChargedJetPartPairwiseDeltaRHist->Draw("COLZ");
  recoChargedJetPartPairwiseDeltaRHist->SetTitle("Pairwise Constituent Delta R;Delta R");
  recoChargedJetPartPairwiseDeltaRHist->GetXaxis()->SetRangeUser(0,0.5);
  gPad->SetLogy();
  if(PRINT) c11->Print((results_path+"/recoJetConstituentPairwiseDR."+results_suffix+".png").c_str()); // Distance between each pair of constituents in reconstructed jets

  // Reco E Vs Eta No Electron Jets
  TCanvas *c12 = new TCanvas("c12","Reco Jet E Vs Eta (No Electrons)",800,600);
  c12->Clear();
  c12->Divide(1,1);

  c12->cd(1);
  recoChargedJetEvsEtaNoElecHist->Draw("COLZ");
  recoChargedJetEvsEtaNoElecHist->SetTitle("Reconstructed Jet Energy Vs Eta (No Electrons);Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c12->Print((results_path+"/recoJetEnergyVsEtaNoElectron."+results_suffix+".png").c_str()); // Reconstructed jet energy - no jets containing electrons included

  // Reco Phi Vs Eta No Electron Jets
  TCanvas *c13 = new TCanvas("c13","Reco Jet Phi Vs Eta (No Electrons)",800,600);
  c13->Clear();
  c13->Divide(1,1);

  c13->cd(1);
  recoChargedJetPhiVsEtaECutNoElecHist->Draw("COLZ");
  recoChargedJetPhiVsEtaECutNoElecHist->SetTitle("Reconstructed Jet Phi Vs Eta (E > 5) (No Electrons);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c13->Print((results_path+"/recoJetPhiVsEtaNoElectron."+results_suffix+".png").c_str()); // Reconstructed Jet phi vs eta - no jets containing electrons included

  // Reco Part P Vs Eta No Electron Jets
  TCanvas *c14 = new TCanvas("c14","Reco Jet Constituent Momentum Vs Eta (No Electrons)",800,600);
  c14->Clear();
  c14->Divide(1,1);

  c14->cd(1);
  recoChargedJetPartPvsEtaNoElecHist->Draw("COLZ");
  recoChargedJetPartPvsEtaNoElecHist->SetTitle("Reconstructed Jet Constituent Momentum Vs Eta (No Electrons);Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c14->Print((results_path+"/recoJetConstituentMomentumVsEtaNoElectron."+results_suffix+".png").c_str()); // Reconstructed jet constituent momentum vs eta - no jets containing electrons included

  // Reco Part Phi Vs Eta No Electron Jets
  TCanvas *c15 = new TCanvas("c15","Reco Jet Constituent Phi Vs Eta (No Electrons)",800,600);
  c15->Clear();
  c15->Divide(1,1);

  c15->cd(1);
  recoChargedJetPartPhiVsEtaNoElecHist->Draw("COLZ");
  recoChargedJetPartPhiVsEtaNoElecHist->SetTitle("Reconstructed Jet Constituent Phi Vs Eta (No Electrons);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c15->Print((results_path+"/recoJetConstituentPhiVsEtaNoElectron."+results_suffix+".png").c_str()); // Reconstructed jet constituent phi vs eta - no jets containing electrons included

  
  ////////////////////////  Generated Jets Plots  ////////////////////////
  // Gen Number
  TCanvas *c16 = new TCanvas("c16","Number Gen Jets",800,600);
  c16->Clear();
  c16->Divide(1,1);

  c16->cd(1);
  numGenChargedJetsECutHist->Draw("HIST");
  numGenChargedJetsECutNoElecHist->SetLineColor(seabornRed);
  numGenChargedJetsECutNoElecHist->Draw("HISTSAME");

  numGenChargedJetsECutHist->SetLineWidth(2);
  numGenChargedJetsECutNoElecHist->SetLineWidth(2);

  numGenChargedJetsECutHist->SetTitle("Generator Jets per Event (|eta| < 2.5 && E > 5);Number");

  TLegend *legend16 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend16->AddEntry(numGenChargedJetsECutHist, "With Electrons", "l");
  legend16->AddEntry(numGenChargedJetsECutNoElecHist, "No Electrons", "l");
  legend16->Draw();

  gPad->SetLogy();
  if(PRINT) c16->Print((results_path+"/numberGenJets."+results_suffix+".png").c_str()); // Number of generator jets per event with energy > 5 GeV and Abs(eta) < 2.5

  // Gen Energy
  TCanvas *c17 = new TCanvas("c17","Gen Jet Energy",800,600);
  c17->Clear();
  c17->Divide(1,1);

  c17->cd(1);
  genChargedJetEHist->Draw("HIST");
  genChargedJetENoElecHist->SetLineColor(seabornRed);
  genChargedJetENoElecHist->Draw("HISTSAME");

  genChargedJetEHist->SetLineWidth(2);
  genChargedJetENoElecHist->SetLineWidth(2);

  genChargedJetEHist->SetTitle("Generator Jet Energy (|eta| < 2.5);Energy [GeV]");
  TLegend *legend17 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend17->AddEntry(genChargedJetEHist, "With Electrons", "l");
  legend17->AddEntry(genChargedJetENoElecHist, "No Electrons", "l");
  legend17->Draw();

  gPad->SetLogy();
  if(PRINT) c17->Print((results_path+"/genJetEnergy."+results_suffix+".png").c_str()); // Energy spectrum of generated jets with Abs(eta) < 2.5

  // Gen Eta
  TCanvas *c18 = new TCanvas("c18","Gen Jet Eta",800,600);
  c18->Clear();
  c18->Divide(1,1);

  c18->cd(1);
  genChargedJetEtaECutHist->Draw("HIST");
  genChargedJetEtaECutNoElecHist->SetLineColor(seabornRed);
  genChargedJetEtaECutNoElecHist->Draw("HISTSAME");

  genChargedJetEtaECutHist->SetLineWidth(2);
  genChargedJetEtaECutNoElecHist->SetLineWidth(2);

  genChargedJetEtaECutHist->SetTitle("Generator Jet Eta (E > 5);Eta");

  TLegend *legend18 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend18->AddEntry(genChargedJetEtaECutHist, "With Electrons", "l");
  legend18->AddEntry(genChargedJetEtaECutNoElecHist, "No Electrons", "l");
  legend18->Draw();

  gPad->SetLogy();
  if(PRINT) c18->Print((results_path+"/genJetEta."+results_suffix+".png").c_str()); // Eta spectrum of generator jets with energy > 5 GeV

  // Gen Area
  if(useNewEDM) {
    TCanvas *c18_1 = new TCanvas("c18_1","Gen Jet Area",800,600);
    c18_1->Clear();
    c18_1->Divide(1,1);

    c18_1->cd(1);
    genChargedJetAreaECutHist->Draw("HIST");
    genChargedJetAreaECutNoElecHist->SetLineColor(seabornRed);
    genChargedJetAreaECutNoElecHist->Draw("HISTSAME");

    genChargedJetAreaECutHist->SetLineWidth(2);
    genChargedJetAreaECutNoElecHist->SetLineWidth(2);

    genChargedJetAreaECutHist->SetTitle("Generator Jet Area (E > 5);Area");

    TLegend *legend18_1 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
    legend18_1->AddEntry(genChargedJetAreaECutHist, "With Electrons", "l");
    legend18_1->AddEntry(genChargedJetAreaECutNoElecHist, "No Electrons", "l");
    legend18_1->Draw();

    gPad->SetLogy();
    if(PRINT) c18_1->Print((results_path+"/genJetArea."+results_suffix+".png").c_str()); // Area spectrum of generator jets with energy > 5 GeV
  }

  // Gen E Vs Eta
  TCanvas *c19 = new TCanvas("c19","Gen Jet E Vs Eta",800,600);
  c19->Clear();
  c19->Divide(1,1);

  c19->cd(1);
  genChargedJetEvsEtaHist->Draw("COLZ");
  genChargedJetEvsEtaHist->SetTitle("Generator Jet Energy Vs Eta;Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c19->Print((results_path+"/genJetEnergyvsEta."+results_suffix+".png").c_str()); // Energy vs eta of generator jets

  // Gen E Vs Area
  if(useNewEDM) {
    TCanvas *c19_1 = new TCanvas("c19_1","Gen Jet E Vs Area",800,600);
    c19_1->Clear();
    c19_1->Divide(1,1);

    c19_1->cd(1);
    genChargedJetEvsAreaHist->Draw("COLZ");
    genChargedJetEvsAreaHist->SetTitle("Generator Jet Energy Vs Area;Area;Energy [GeV]");
    gPad->SetLogz();
    if(PRINT) c19_1->Print((results_path+"/genJetEnergyvsArea."+results_suffix+".png").c_str()); // Energy vs area of generator jets
  }

  // Gen Phi Vs Eta
  TCanvas *c20 = new TCanvas("c20","Gen Jet Phi Vs Eta",800,600);
  c20->Clear();
  c20->Divide(1,1);

  c20->cd(1);
  genChargedJetPhiVsEtaECutHist->Draw("COLZ");
  genChargedJetPhiVsEtaECutHist->SetTitle("Generator Jet Phi Vs Eta (E > 5);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c20->Print((results_path+"/genJetPhiVsEta."+results_suffix+".png").c_str()); // Phi vs eta of generator jets

  // Num Particles Per Gen Jet
  TCanvas *c21 = new TCanvas("c21","Number Constituents Per Gen Jet",800,600);
  c21->Clear();
  c21->Divide(1,1);

  c21->cd(1);
  numGenChargedJetPartsHist->Draw("HIST");
  numGenChargedJetPartsNoElecHist->SetLineColor(seabornRed);
  numGenChargedJetPartsNoElecHist->Draw("HISTSAME");

  numGenChargedJetPartsHist->SetLineWidth(2);
  numGenChargedJetPartsNoElecHist->SetLineWidth(2);

  numGenChargedJetPartsHist->SetTitle("Number of Constituents Per Gen Jet;Number of Constituents");

  TLegend *legend21 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend21->AddEntry(numGenChargedJetPartsHist, "With Electrons", "l");
  legend21->AddEntry(numGenChargedJetPartsNoElecHist, "No Electrons", "l");
  legend21->Draw();
  gPad->SetLogy();
  if(PRINT) c21->Print((results_path+"/numConstituentsPerGenJet."+results_suffix+".png").c_str()); // Number of constituents in generator jets

  // Gen Part Momentum
  TCanvas *c22 = new TCanvas("c22","Gen Jet Constituent Momentum",800,600);
  c22->Clear();
  c22->Divide(1,1);

  c22->cd(1);
  genChargedJetPartPHist->Draw("HIST");
  genChargedJetPartPNoElecHist->SetLineColor(seabornRed);
  genChargedJetPartPNoElecHist->Draw("HISTSAME");

  genChargedJetPartPHist->SetLineWidth(2);
  genChargedJetPartPNoElecHist->SetLineWidth(2);

  genChargedJetPartPHist->SetTitle("Generator Jet Constituent Momentum;Momentum [GeV/c]");

  TLegend *legend22 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend22->AddEntry(genChargedJetPartPHist, "With Electrons", "l");
  legend22->AddEntry(genChargedJetPartPNoElecHist, "No Electrons", "l");
  legend22->Draw();

  gPad->SetLogy();
  if(PRINT) c22->Print((results_path+"/genJetConstituentMomentum."+results_suffix+".png").c_str()); // Momentum of generator jet constituents

  // Gen Part Eta
  TCanvas *c23 = new TCanvas("c23","Gen Jet Constituent Eta",800,600);
  c23->Clear();
  c23->Divide(1,1);

  c23->cd(1);
  genChargedJetPartEtaHist->Draw("HIST");
  genChargedJetPartEtaNoElecHist->SetLineColor(seabornRed);
  genChargedJetPartEtaNoElecHist->Draw("HISTSAME");

  genChargedJetPartEtaHist->SetLineWidth(2);
  genChargedJetPartEtaNoElecHist->SetLineWidth(2);

  genChargedJetPartEtaHist->SetTitle("Generator Jet Constituent Eta;Eta");

  TLegend *legend23 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend23->AddEntry(genChargedJetPartEtaHist, "With Electrons", "l");
  legend23->AddEntry(genChargedJetPartEtaNoElecHist, "No Electrons", "l");
  legend23->Draw();

  gPad->SetLogy();
  if(PRINT) c23->Print((results_path+"/genJetConstituentEta."+results_suffix+".png").c_str()); // Eta of generator jet constituents

  // Gen Part P Vs Eta
  TCanvas *c24 = new TCanvas("c24","Gen Jet Constituent Momentum Vs Eta",800,600);
  c24->Clear();
  c24->Divide(1,1);

  c24->cd(1);
  genChargedJetPartPvsEtaHist->Draw("COLZ");
  genChargedJetPartPvsEtaHist->SetTitle("Generator Jet Constituent Momentum Vs Eta;Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c24->Print((results_path+"/genJetConstituentMomentumVsEta."+results_suffix+".png").c_str()); // Momentum vs eta of generator jet constituents

  // Gen Part Phi Vs Eta
  TCanvas *c25 = new TCanvas("c25","Gen Jet Constituent Phi Vs Eta",800,600);
  c25->Clear();
  c25->Divide(1,1);

  c25->cd(1);
  genChargedJetPartPhiVsEtaHist->Draw("COLZ");
  genChargedJetPartPhiVsEtaHist->SetTitle("Generator Jet Constituent Phi Vs Eta;Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c25->Print((results_path+"/genJetConstituentPhiVsEta."+results_suffix+".png").c_str()); // Phi vs eta of generator jet constituents

  // Gen Constituent Pairwise delta R
  TCanvas *c26 = new TCanvas("c26","Gen Jet Constituent Pairwise Delta R",800,600);
  c26->Clear();
  c26->Divide(1,1);

  c26->cd(1);
  genChargedJetPartPairwiseDeltaRHist->Draw("COLZ");
  genChargedJetPartPairwiseDeltaRHist->SetTitle("Generator Jet Pairwise Constituent Delta R;Delta R");
  genChargedJetPartPairwiseDeltaRHist->GetXaxis()->SetRangeUser(0,0.5);
  gPad->SetLogy();
  if(PRINT) c26->Print((results_path+"/genJetConstituentPairwiseDR."+results_suffix+".png").c_str()); // Distance between each pair of constituents in generator jets

  // Gen E Vs Eta No Electron Jets
  TCanvas *c27 = new TCanvas("c27","Gen Jet E Vs Eta (No Electrons)",800,600);
  c27->Clear();
  c27->Divide(1,1);

  c27->cd(1);
  genChargedJetEvsEtaNoElecHist->Draw("COLZ");
  genChargedJetEvsEtaNoElecHist->SetTitle("Generator Jet Energy Vs Eta (No Electrons);Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c27->Print((results_path+"/genJetEnergyVsEtaNoElectron."+results_suffix+".png").c_str()); // Generator jet energy vs eta - no jets containing electrons included

  // Gen Phi Vs Eta No Electron Jets
  TCanvas *c28 = new TCanvas("c28","Gen Jet Phi Vs Eta (No Electrons)",800,600);
  c28->Clear();
  c28->Divide(1,1);

  c28->cd(1);
  genChargedJetPhiVsEtaECutNoElecHist->Draw("COLZ");
  genChargedJetPhiVsEtaECutNoElecHist->SetTitle("Generator Jet Phi Vs Eta (E > 5) (No Electrons);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c28->Print((results_path+"/genJetPhiVsEtaNoElectron."+results_suffix+".png").c_str()); // Generator Jet phi vs eta - no jets containing electrons included

  // Gen Part P Vs Eta No Electron Jets
  TCanvas *c29 = new TCanvas("c29","Gen Jet Constituent Momentum Vs Eta (No Electrons)",800,600);
  c29->Clear();
  c29->Divide(1,1);

  c29->cd(1);
  genChargedJetPartPvsEtaNoElecHist->Draw("COLZ");
  genChargedJetPartPvsEtaNoElecHist->SetTitle("Generator Jet Constituent Momentum Vs Eta (No Electrons);Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c29->Print((results_path+"/genJetConstituentMomentumVsEtaNoElectron."+results_suffix+".png").c_str()); // Generator jet constituent momentum vs eta - no jets containing electrons included

  // Gen Part Phi Vs Eta No Electron Jets
  TCanvas *c30 = new TCanvas("c30","Gen Jet Constituent Phi Vs Eta (No Electrons)",800,600);
  c30->Clear();
  c30->Divide(1,1);

  c30->cd(1);
  genChargedJetPartPhiVsEtaNoElecHist->Draw("COLZ");
  genChargedJetPartPhiVsEtaNoElecHist->SetTitle("Generator Jet Constituent Phi Vs Eta (No Electrons);Eta;Phi");
  gPad->SetLogz();
  //c30->Print((results_path+"/recoJetEvsEta.png").c_str());
  if(PRINT) c30->Print((results_path+"/genJetConstituentPhiVsEtaNoElectron."+results_suffix+".png").c_str()); // Generator jet constituent phi vs eta - no jets containing electrons included

  
  ////////////////////////  Matched Jets Plots  ////////////////////////
  // Matched Delta R
  TCanvas *c31 = new TCanvas("c31","Gen - Reco Delta R",800,600);
  c31->Clear();
  c31->Divide(1,1);

  c31->cd(1);
  matchJetDeltaRHist->Draw("HIST");
  matchJetDeltaRBackHist->SetLineColor(seabornRed);
  //matchJetDeltaRBackHist->Draw("HISTSAME");
  matchJetDeltaRHist->SetTitle("Matched Gen - Reco Jet Delta R;Delta R");
  gPad->SetLogy();
  if(PRINT) c31->Print((results_path+"/genRecoJetDeltaR."+results_suffix+".png").c_str()); // Distance between closest generated and reconstructed jet pair

  // Matched Reco Vs Gen Eta
  TCanvas *c32 = new TCanvas("c32","Reco Vs Gen Eta",800,600);
  c32->Clear();
  c32->Divide(1,1);

  c32->cd(1);
  recoVsGenChargedJetEtaHist->Draw("COLZ");
  recoVsGenChargedJetEtaHist->SetTitle("Reconstructed Vs Generator Jet Eta;Gen Eta;Reco Eta");
  gPad->SetLogz();
  if(PRINT) c32->Print((results_path+"/matchedRecoVsGenJetEta."+results_suffix+".png").c_str()); // Matched Reconstructed Vs Generator Jet eta

  // Matched Reco Vs Gen Phi
  TCanvas *c33 = new TCanvas("c33","Reco Vs Gen Phi",800,600);
  c33->Clear();
  c33->Divide(1,1);

  c33->cd(1);
  recoVsGenChargedJetPhiHist->Draw("COLZ");
  recoVsGenChargedJetPhiHist->SetTitle("Reconstructed Vs Generator Jet Phi;Gen Phi;Reco Phi");
  gPad->SetLogz();
  if(PRINT) c33->Print((results_path+"/matchedRecoVsGenJetPhi."+results_suffix+".png").c_str()); // Matched reconstructed vs generator jet phi

  // Matched Reco Vs Gen Area
  if(useNewEDM) {
    TCanvas *c33_1 = new TCanvas("c33_1","Reco Vs Gen Area",800,600);
    c33_1->Clear();
    c33_1->Divide(1,1);

    c33_1->cd(1);
    recoVsGenChargedJetAreaHist->Draw("COLZ");
    recoVsGenChargedJetAreaHist->SetTitle("Reconstructed Vs Generator Jet Area;Gen Area;Reco Area");
    gPad->SetLogz();
    if(PRINT) c33_1->Print((results_path+"/matchedRecoVsGenJetArea."+results_suffix+".png").c_str()); // Matched reconstructed vs generator jet area
  }

  // Matched Reco Vs Gen Energy
  TCanvas *c34 = new TCanvas("c34","Reco Vs Gen Energy",800,600);
  c34->Clear();
  c34->Divide(1,1);

  TF1 *f1_34 = new TF1("f1_34","1.0*x + 0.0",1,100);
  TF1 *f2_34 = new TF1("f2_34","2.0*x + 0.0",1,100);
  TF1 *f3_34 = new TF1("f3_34","3.0*x + 0.0",1,100);

  c34->cd(1);
  recoVsGenChargedJetEHist->Draw("COLZ");
  recoVsGenChargedJetEHist->SetTitle("Reconstructed Vs Generator Jet Energy;Gen E;Reco E");
  f1_34->Draw("SAME");
  f2_34->Draw("SAME");
  f3_34->Draw("SAME");
  gPad->SetLogz();
  if(PRINT) c34->Print((results_path+"/matchedRecoVsGenJetEnergy."+results_suffix+".png").c_str()); // Matched reconstructed vs generator jet energy

  // Jet Res Vs Gen Eta
  TCanvas *c35 = new TCanvas("c35","Jet Res Vs Gen Eta",800,600);
  c35->Clear();
  c35->Divide(1,1);

  c35->cd(1);
  jetResVsEtaHist->Draw("COLZ");
  jetResVsEtaHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Eta;Gen Eta;Res");
  gPad->SetLogz();
  if(PRINT) c35->Print((results_path+"/matchedJetResolutionVsEta."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet eta

  // Jet Res Vs Gen E
  TCanvas *c36 = new TCanvas("c36","Jet Res Vs Gen E",800,600);
  c36->Clear();
  c36->Divide(1,1);

  c36->cd(1);
  jetResVsEHist->Draw("COLZ");
  jetResVsEHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c36->Print((results_path+"/matchedJetResolutionVsEnergy."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy

  // Jet Res Vs Gen E Neg Eta
  TCanvas *c37 = new TCanvas("c37","Jet Res Vs Gen E (-2.5 < eta < -1.0)",800,600);
  c37->Clear();
  c37->Divide(1,1);

  c37->cd(1);
  jetResVsENegEtaNoDupHist->Draw("COLZ");
  jetResVsENegEtaNoDupHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy (-2.5 < eta < -1.0) No Duplicate;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c37->Print((results_path+"/matchedJetResolutionVsEnergyNegEta."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy -2.5 < eta < -1.0

  // Jet Res Vs Gen E Mid Eta
  TCanvas *c38 = new TCanvas("c38","Jet Res Vs Gen E (-1.0 < eta < 1.0)",800,600);
  c38->Clear();
  c38->Divide(1,1);

  c38->cd(1);
  jetResVsEMidEtaNoDupHist->Draw("COLZ");
  jetResVsEMidEtaNoDupHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy (-1.0 < eta < 1.0) No Duplicate;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c38->Print((results_path+"/matchedJetResolutionVsEnergyMidEta."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy -1.0 < eta < 1.0
    delete c38;
  // Jet Res Vs Gen E Pos Eta
  TCanvas *c39 = new TCanvas("c39","Jet Res Vs Gen E (1.0 < eta < 2.5)",800,600);
  c39->Clear();
  c39->Divide(1,1);

  c39->cd(1);
  jetResVsEPosEtaNoDupHist->Draw("COLZ");
  jetResVsEPosEtaNoDupHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy (1.0 < eta < 2.5) No Duplicate;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c39->Print((results_path+"/matchedJetResolutionVsEnergyPosEta."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy 1.0 < eta < 2.5
  delete c39;

  
  // Generate Resolution Plots
  const int BINS = 20;
  double binCent[BINS];
  double jesVsENeg[BINS];
  double jesVsEMid[BINS];
  double jesVsEPos[BINS];
  double jerVsENeg[BINS];
  double jerVsEMid[BINS];
  double jerVsEPos[BINS];

  std::fill(std::begin(binCent), std::end(binCent), -999.);
  std::fill(std::begin(jesVsENeg), std::end(jesVsENeg), -999.);
  std::fill(std::begin(jesVsEMid), std::end(jesVsEMid), -999.);
  std::fill(std::begin(jesVsEPos), std::end(jesVsEPos), -999.);
  std::fill(std::begin(jerVsENeg), std::end(jerVsENeg), -999.);
  std::fill(std::begin(jerVsEMid), std::end(jerVsEMid), -999.);
  std::fill(std::begin(jerVsEPos), std::end(jerVsEPos), -999.);

  TH1D *pxA = jetResVsENegEtaNoDupHist->ProjectionX("pxA",1,10000);
  for(int i=0; i<BINS; i++)
    {
      binCent[i] = pxA->GetBinCenter(i+1);
    }

  TCanvas *c40 = new TCanvas("c40","Negative Rapidity Fit Results",800,600);
  c40->Clear();
  c40->Divide(5,4);

  TH1D *hA[20];
  for(int i=1; i<21; i++)
    {
      hA[i-1] = (TH1D *)jetResVsENegEtaNoDupHist->ProjectionY(Form("projYA_%d",i),i,i);

      TF1 *myGausA = new TF1("myGausA","gaus",-0.5,0.5);
      myGausA->SetParameters(hA[i-1]->GetMaximum(),0.0,0.01);

      c40->cd(i);
      //hA[i-1]->Draw("HIST");
      hA[i-1]->Fit("myGausA","B","",-0.5,0.5);
      hA[i-1]->GetXaxis()->SetRangeUser(-1,1);
      gPad->SetLogy();

      if(hA[i-1]->GetEntries() > 2)
	{
	  auto fA = hA[i-1]->GetFunction("myGausA");

	  if(fA->GetParError(2)/fA->GetParameter(2) < 0.5)
	    {
	      jesVsENeg[i-1] = fA->GetParameter(1);
	      jerVsENeg[i-1] = fA->GetParameter(2);
	    }
	  //cout << fA->GetParameter(0) << " " << fA->GetParameter(1) << " " << fA->GetParameter(2) << endl;
	  //cout << fA->GetParError(0) << " " << fA->GetParError(1) << " " << fA->GetParError(2) << endl;
	}
    }
  if(PRINT) c40->Print((results_path+"/matchedJetResolutionVsEnergyNegEtaFitSummary."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy -2.5 < eta < -1.0 fits
    delete c40;

  TCanvas *c41 = new TCanvas("c41","Mid Rapidity Fit Results",800,600);
  c41->Clear();
  c41->Divide(5,4);

  TH1D *hB[20];
  for(int i=1; i<21; i++)
    {
      hB[i-1] = (TH1D *)jetResVsEMidEtaNoDupHist->ProjectionY(Form("projYB_%d",i),i,i);

      TF1 *myGausB = new TF1("myGausB","gaus",-0.5,0.5);
      myGausB->SetParameters(hB[i-1]->GetMaximum(),0.0,0.01);

      c41->cd(i);
      //hA[i-1]->Draw("HIST");
      hB[i-1]->Fit("myGausB","B","",-0.5,0.5);
      hB[i-1]->GetXaxis()->SetRangeUser(-1,1);
      gPad->SetLogy();

      if(hB[i-1]->GetEntries() > 2)
	{
	  auto fB = hB[i-1]->GetFunction("myGausB");

	  if(fB->GetParError(2)/fB->GetParameter(2) < 0.5)
	    {
	      jesVsEMid[i-1] = fB->GetParameter(1);
	      jerVsEMid[i-1] = fB->GetParameter(2);
	    }
	  //cout << fB->GetParameter(0) << " " << fB->GetParameter(1) << " " << fB->GetParameter(2) << endl;
	  //cout << fB->GetParError(0) << " " << fB->GetParError(1) << " " << fB->GetParError(2) << endl;
	}
    }
  if(PRINT) c41->Print((results_path+"/matchedJetResolutionVsEnergyMidEtaFitSummary."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy -1.0 < eta < 1.0 fits
    delete c41;

  TCanvas *c42 = new TCanvas("c42","Positive Rapidity Fit Results",800,600);
  c42->Clear();
  c42->Divide(5,4);

  TH1D *hC[20];
  for(int i=1; i<21; i++)
    {
      hC[i-1] = (TH1D *)jetResVsEPosEtaNoDupHist->ProjectionY(Form("projYC_%d",i),i,i);

      TF1 *myGausC = new TF1("myGausC","gaus",-0.5,0.5);
      myGausC->SetParameters(hC[i-1]->GetMaximum(),0.0,0.01);

      c42->cd(i);
      //hA[i-1]->Draw("HIST");
      hC[i-1]->Fit("myGausC","B","",-0.5,0.5);
      hC[i-1]->GetXaxis()->SetRangeUser(-1,1);
      gPad->SetLogy();

      if(hC[i-1]->GetEntries() > 2)
	{
	  auto fC = hC[i-1]->GetFunction("myGausC");

	  if(fC->GetParError(2)/fC->GetParameter(2) < 0.5)
	    {
	      jesVsEPos[i-1] = fC->GetParameter(1);
	      jerVsEPos[i-1] = fC->GetParameter(2);
	    }
	  //cout << fC->GetParameter(0) << " " << fC->GetParameter(1) << " " << fC->GetParameter(2) << endl;
	  //cout << fC->GetParError(0) << " " << fC->GetParError(1) << " " << fC->GetParError(2) << endl;
	}
    }
  if(PRINT) c42->Print((results_path+"/matchedJetResolutionVsEnergyPosEtaFitSummary."+results_suffix+".png").c_str()); // Matched jet resolution vs generator jet energy 1.0 < eta < 2.5 fits
    delete c42;
  TCanvas *c43 = new TCanvas("c43","Positive JES/JER",800,600);
  c43->Clear();
  c43->Divide(1,1);

  TGraph *gJESvsENeg = new TGraph(BINS,binCent,jesVsENeg);
  TGraph *gJERvsENeg = new TGraph(BINS,binCent,jerVsENeg);

  TGraph *gJESvsEMid = new TGraph(BINS,binCent,jesVsEMid);
  TGraph *gJERvsEMid = new TGraph(BINS,binCent,jerVsEMid);

  TGraph *gJESvsEPos = new TGraph(BINS,binCent,jesVsEPos);
  TGraph *gJERvsEPos = new TGraph(BINS,binCent,jerVsEPos);

  TH2D *test43 = new TH2D("test43","Jet Energy Scale / Resolution Vs Eta;True Eta;JES/JER",1,0.,100.,1,-0.2,0.2);
  test43->Draw();

  c43->cd(1);
  gJERvsENeg->Draw("*");
  gJERvsENeg->SetMarkerStyle(21);
  gJERvsENeg->SetMarkerSize(1);
  gJERvsENeg->SetMarkerColor(seabornBlue);

  gJESvsENeg->Draw("*");
  gJESvsENeg->SetMarkerStyle(26);
  gJESvsENeg->SetMarkerSize(1);
  gJESvsENeg->SetMarkerColor(seabornBlue);

  gJERvsEMid->Draw("*");
  gJERvsEMid->SetMarkerStyle(21);
  gJERvsEMid->SetMarkerSize(1);
  gJERvsEMid->SetMarkerColor(seabornRed);

  gJESvsEMid->Draw("*");
  gJESvsEMid->SetMarkerStyle(26);
  gJESvsEMid->SetMarkerSize(1);
  gJESvsEMid->SetMarkerColor(seabornRed);

  gJERvsEPos->Draw("*");
  gJERvsEPos->SetMarkerStyle(21);
  gJERvsEPos->SetMarkerSize(1);
  gJERvsEPos->SetMarkerColor(seabornGreen);

  gJESvsEPos->Draw("*");
  gJESvsEPos->SetMarkerStyle(26);
  gJESvsEPos->SetMarkerSize(1);
  gJESvsEPos->SetMarkerColor(seabornGreen);

  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9); 
  legend->AddEntry(gJERvsENeg, "JER, (-2.5 < #eta < -1)", "p");
  legend->AddEntry(gJESvsENeg, "JES, (-2.5 < #eta < -1)","p");
  legend->AddEntry(gJERvsEMid, "JER, (-1 < #eta < 1)", "p");
  legend->AddEntry(gJESvsEMid, "JES, (-1 < #eta < 1)", "p");
  legend->AddEntry(gJERvsEPos, "JER, (1 < #eta < 2.5) ", "p");
  legend->AddEntry(gJESvsEPos, "JES, (1 < #eta < 2.5)", "p");
  legend->Draw();

  if(PRINT) c43->Print((results_path+"/matchedJetScaleResolutionSummary."+results_suffix+".png").c_str()); // Matched jet JER/JES summary
    delete c43;

  return 0;
}
