// File: fitD0withStudentT.C

#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <vector>

// Background mode flag: 0 = Polynomial (2nd order), 1 = Exponential
int bkgtype_mode = 0;
const Int_t nBins = 3; // n = 12

Double_t y_arr[nBins+1]  = {-3.0, -1.0,  1.0,  3.0};
Double_t pt_arr[nBins+2] = { 0.0, 1.0,  2.0,  5.0, 10.0};
double min_mass = 1.6, max_mass = 2.5;

// Polynomial or exponential background
Double_t background(Double_t *x, Double_t *par) {
  if (bkgtype_mode == 0) {
    // Polynomial: par[0] + par[1]*x + par[2]*x^2
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  } else {
    // Exponential: par[0] * exp(-x / par[1])
    return par[0] * TMath::Exp(-x[0]/par[1]);
  }
}

// Student's t-distribution signal
Double_t studentTPeak(Double_t *x, Double_t *par) {
  Double_t norm  = par[0];   // Height
  Double_t mean  = par[1];   // Mean
  Double_t width = par[2];   // Width
  Double_t n     = par[3];   // Degrees of freedom

  Double_t t  = (x[0] - mean) / width;
  Double_t val = norm * TMath::Gamma((n+1)/2.0) /
                 (TMath::Sqrt(n*TMath::Pi()) * width * TMath::Gamma(n/2.0)) *
                 TMath::Power(1 + t*t/n, -(n+1)/2.0);
  return val;
}

// Combined fit function: signal + background
Double_t fitFunction(Double_t *x, Double_t *par) {
  Double_t signal = studentTPeak(x, par);
  return signal + background(x, &par[4]);
}

void fitD0withStudentT() {

  // Styles
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(1.2,"X"); gStyle->SetTitleOffset(1.6,"Y");
  gStyle->SetTitleSize(.04,"X");   gStyle->SetTitleSize(.04,"Y");
  gStyle->SetLabelSize(.04,"X");   gStyle->SetLabelSize(.04,"Y");
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(1.3);
  gStyle->SetStatW(0.12);
  gStyle->SetStatH(0.12);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.95);

  gSystem->Exec("rm -rf DzeroMassPlots_Student/");
  gSystem->Exec("mkdir -p DzeroMassPlots_Student/");

  TCanvas *canvas[nBins][nBins];
  TH1D *hMassDzero[nBins][nBins]; // keep structure
  TH1D *hMassFit[nBins][nBins];   // keep structure
  TFile *file[nBins][nBins];
  TFile *fout[nBins][nBins];

  TH1F *hFitResult = new TH1F("hFitResult", "Event statistics", 2, 0, 2);
  hFitResult->GetXaxis()->SetBinLabel(1, "Signal");
  hFitResult->GetXaxis()->SetBinLabel(2, "Background");

  // Parameters count
  int nParSig   = 4;
  int nParBkg   = (bkgtype_mode == 0) ? 3 : 2;
  int nTotalPar = nParSig + nParBkg;

  for (Int_t iy=0; iy<nBins; ++iy) {
    for (Int_t ipt=0; ipt<nBins; ++ipt) {

      file[iy][ipt] = TFile::Open(Form("final_merged_%1.1f_%1.1f_%1.1f_%1.1f.root",
                                       y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]));
      if (!file[iy][ipt] || file[iy][ipt]->IsZombie()) {
        ::Error("fitD0withStudentT", "Cannot open input file for bin (y:%d,pt:%d). Skipping.", iy, ipt);
        continue;
      }

      fout[iy][ipt] = new TFile(Form("DzeroMassPlots_Student/DzeroInvMassCan_pTbin_%1.1f_%1.1f_%1.1f_%1.1f.root",
                                     y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]), "recreate");

      TF1 *fitFunc = new TF1("fitFunc", fitFunction, min_mass, max_mass, nTotalPar);
      fitFunc->SetParNames("Norm", "Mean", "Width", "nDOF");

      TH1F *hMass = (TH1F*) file[iy][ipt]->Get("hData_D0");
      if (!hMass) {
        ::Error("fitD0withStudentT", "Histogram hData_D0 not found. Skipping (y:%d, pt:%d).", iy, ipt);
        fout[iy][ipt]->Close();
        delete fout[iy][ipt];
        continue;
      }

      hMass->SetTitle(Form("%1.1f < #it{y} < %1.1f && %1.1f < #it{p}_{T} (GeV/#it{c}) < %1.1f",
                           y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]));
      hMass->GetXaxis()->SetTitle("#it{m}_{D^{0}} (GeV/#it{c}^{2})");
      hMass->GetXaxis()->CenterTitle();
      hMass->GetYaxis()->SetTitle("Entries");
      hMass->GetYaxis()->CenterTitle();
      hMass->SetStats(false);
      hMass->SetMarkerStyle(20);
      hMass->SetName(Form("hMass_%1.1f_%1.1f_%1.1f_%1.1f",
                          y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]));
      hMass->SetMinimum(0.);
      int maxbin = hMass->GetXaxis()->FindBin(1.864);
      hMass->SetMaximum(hMass->GetBinContent(maxbin)*1.50);
      hMass->GetXaxis()->SetRangeUser(min_mass, max_mass);
      hMass->SetMarkerColor(kBlack);
      hMass->SetMarkerSize(1.5);
      hMass->SetLineColor(kBlack);

      // Canvas
      canvas[iy][ipt] = new TCanvas(Form("DzeroInvMassCan_pTbin_%1.1f_%1.1f_%1.1f_%1.1f",
                                         y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]),
                                    Form("DzeroInvMassCan_pTbin_%1.1f_%1.1f_%1.1f_%1.1f",
                                         y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]),
                                    1200, 1000);
      canvas[iy][ipt]->SetMargin(0.13, 0.05 ,0.1,0.08);

      // Initial values & limits
      fitFunc->SetParameters(100, 1.865, 0.01, 5);
      fitFunc->SetParLimits(3, 1, 10);       // nDOF range
      fitFunc->SetParLimits(1, 1.84, 1.89);  // mean range
      fitFunc->SetParLimits(2, 0.001, 0.1);  // width range

      // Background initial guesses
      if (bkgtype_mode == 0) {
        fitFunc->SetParNames("Norm", "Mean", "Width", "nDOF", "Bkg0", "Bkg1", "Bkg2");
        fitFunc->SetParameters(100, 1.865, 0.01, 5, 100, 0, 0);
        // only valid for poly2 case (param index 6 exists)
        fitFunc->FixParameter(6, 0.0);
      } else {
        fitFunc->SetParNames("Norm", "Mean", "Width", "nDOF", "BkgAmp", "BkgTau");
        fitFunc->SetParameters(100, 1.865, 0.01, 5, 100, 0.1);
        // keep tau positive for stability
        fitFunc->SetParLimits(5, 1e-4, 1e3);
      }

      // Fit (store result to access covariance)
      hMass->Fit(fitFunc, "NR");
      TFitResultPtr fitres = hMass->Fit(fitFunc, "RS"); // store + quiet

      hMass->Draw("EP");

      // Background component
      TF1 *bkgFunc = new TF1("bkgFunc", background, min_mass, max_mass, nParBkg);
      bkgFunc->SetParameters(&fitFunc->GetParameters()[nParSig]);
      bkgFunc->SetLineColor(kGreen+2);
      bkgFunc->SetLineWidth(3);
      bkgFunc->Draw("same");

      // Signal component
      TF1 *sigFunc = new TF1("sigFunc", studentTPeak, min_mass, max_mass, nParSig);
      sigFunc->SetParameters(&fitFunc->GetParameters()[0]);
      sigFunc->SetLineColor(kRed);
      sigFunc->SetLineWidth(3);
      sigFunc->SetLineStyle(2);
      sigFunc->Draw("same");

      // Total fit
      fitFunc->SetLineColor(kBlue);
      fitFunc->SetLineWidth(3);
      fitFunc->Draw("same");

      // Legend
      TLegend *leg = new TLegend(0.60, 0.65, 0.90, 0.90);
      leg->AddEntry(hMass, "D^{0} Candidates", "lep");
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.025);
      leg->AddEntry((TObject*)nullptr, Form("Entries = %1.1f", hMass->GetEntries()), "");
      leg->AddEntry(fitFunc, "Total Fit", "l");
      leg->AddEntry(sigFunc, "Student-t Signal", "l");
      leg->AddEntry(bkgFunc, bkgtype_mode == 0 ? "Poly2 Background" : "Exponential Background", "l");

      // Fitted parameters
      Double_t mean  = fitFunc->GetParameter(1);
      Double_t width = fitFunc->GetParameter(2);
      Double_t nDOF  = fitFunc->GetParameter(3);

      // ±2σ window
      Double_t xLow  = mean - 2 * width;
      Double_t xHigh = mean + 2 * width;

      // Central values of integrals
      const Double_t binw = hMass->GetBinWidth(1);
      Double_t signalInt = sigFunc->Integral(xLow, xHigh) / binw;
      Double_t bkgInt    = bkgFunc->Integral(xLow, xHigh) / binw;

      // S/B and significance
      Double_t sobRatio     = (bkgInt > 0) ? (signalInt / bkgInt) : 0.0;
      Double_t significance = (signalInt + bkgInt > 0)
                                ? (signalInt / TMath::Sqrt(signalInt + bkgInt))
                                : 0.0;

      // χ2/ndf
      Double_t chi2    = fitFunc->GetChisquare();
      Int_t    ndf     = fitFunc->GetNDF();
      Double_t chi2ndf = (ndf > 0) ? (chi2 / ndf) : 0.0;

      // --- Propagate covariance to S and B (±2σ) ---
      Double_t sigErr = 0.0, bkgErr = 0.0, sobErr = 0.0;

      if (fitres.Get()) {
        const TMatrixDSym& covFull = fitres->GetCovarianceMatrix();

        // Flatten signal covariance (first 4 params: Norm..nDOF)
        std::vector<Double_t> covSigFlat(4*4);
        for (int r=0; r<4; ++r)
          for (int c=0; c<4; ++c)
            covSigFlat[r*4 + c] = covFull(r, c);

        // Flatten background covariance (next nParBkg params)
        std::vector<Double_t> covBkgFlat(nParBkg*nParBkg);
        for (int r=0; r<nParBkg; ++r)
          for (int c=0; c<nParBkg; ++c)
            covBkgFlat[r*nParBkg + c] = covFull(r+4, c+4);

        // Parameter buffers for these TF1s
        Double_t parsSig[4];        sigFunc->GetParameters(parsSig);
        std::vector<Double_t> parsBkg(nParBkg);
        bkgFunc->GetParameters(parsBkg.data());

        sigErr = sigFunc->IntegralError(xLow, xHigh, parsSig,        covSigFlat.data()) / binw;
        bkgErr = bkgFunc->IntegralError(xLow, xHigh, parsBkg.data(), covBkgFlat.data()) / binw;

        if (signalInt>0 && bkgInt>0) {
          sobErr = sobRatio * TMath::Sqrt( TMath::Power(sigErr/signalInt,2)
                                         + TMath::Power(bkgErr/bkgInt,  2) );
        }
      }

      // Legend entries with uncertainties
      leg->AddEntry((TObject*)nullptr, Form("S (2#sigma) = %.2f #pm %.2f", signalInt, sigErr), "");
      leg->AddEntry((TObject*)nullptr, Form("B (2#sigma) = %.2f #pm %.2f", bkgInt, bkgErr), "");
      leg->AddEntry((TObject*)nullptr, Form("S/B (2#sigma) = %.2f %s", sobRatio,
                                            (sobErr>0 ? Form("#pm %.2f", sobErr) : "")), "");
      leg->AddEntry((TObject*)nullptr, Form("Significance = %.2f", significance), "");
      leg->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f", chi2ndf), "");
      leg->Draw();

      // Save plot and summary
      canvas[iy][ipt]->SaveAs(Form("DzeroMassPlots_Student/DzeroInvMassCan_pTbin_%1.1f_%1.1f_%1.1f_%1.1f.png",
                                   y_arr[iy], y_arr[iy+1], pt_arr[ipt], pt_arr[ipt+1]));
      hFitResult->SetBinContent(1, signalInt);
      hFitResult->SetBinContent(2, bkgInt);

      fout[iy][ipt]->cd();
      hFitResult->Write();
      fout[iy][ipt]->Close();
      file[iy][ipt]->Close();
    }
  }
}

