double fit_slewing(double *x, double *par) {
  return par[0] + par[1] * pow(x[0] - par[2], par[3]);
}

double fit_landau(double *x, double *par) {
  return par[2] * TMath::Landau(x[0], par[0], par[1], kTRUE);
}

void timeWalkCor(TString filename) {

  TChain *ch = new TChain("events");
  ch->Add(filename);
  TTreeReader reader(ch);

  // CALOROCHits(p-going)
  TTreeReaderArray<int> pos_phase(reader, "EcalBarrelScFiPCALOROCHits.samplePhase");
  TTreeReaderArray<int> pos_stamp(reader, "EcalBarrelScFiPCALOROCHits.timeStamp");
  TTreeReaderArray<unsigned int> pos_begin(reader, "EcalBarrelScFiPCALOROCHits.bSamples_begin");
  TTreeReaderArray<unsigned int> pos_end(reader, "EcalBarrelScFiPCALOROCHits.bSamples_end");
  TTreeReaderArray<uint16_t> pos_lowADC( reader, "_EcalBarrelScFiPCALOROCHits_bSamples.lowGainADC");
  TTreeReaderArray<uint16_t> pos_TOA(reader, "_EcalBarrelScFiPCALOROCHits_bSamples.timeOfArrival");

  // CALOROCHits(e-going)
  TTreeReaderArray<int> neg_phase(reader, "EcalBarrelScFiNCALOROCHits.samplePhase");
  TTreeReaderArray<int> neg_stamp(reader, "EcalBarrelScFiNCALOROCHits.timeStamp");
  TTreeReaderArray<unsigned int> neg_begin(reader, "EcalBarrelScFiNCALOROCHits.bSamples_begin");
  TTreeReaderArray<unsigned int> neg_end(reader, "EcalBarrelScFiNCALOROCHits.bSamples_end");
  TTreeReaderArray<uint16_t> neg_lowADC(reader, "_EcalBarrelScFiNCALOROCHits_bSamples.lowGainADC");
  TTreeReaderArray<uint16_t> neg_TOA(reader, "_EcalBarrelScFiNCALOROCHits_bSamples.timeOfArrival");

  // Config variables
  const std::size_t n_samples = 7;
  const double capADC = 1024;
  const double dyRangeADC = 2500;
  const double capTOA = 1024;
  const double toa_thres = 7;
  const double timeWindow = 25;
  const double pulse_sigma = 10;

  // To perform time-walk correction
  TH2D *h2_ADC_dt = new TH2D("h2_ADC_dt", "", 80, 0, 400, 200, -25, 25);
  TF1 *f1_slewing = new TF1("f1_slewing", fit_slewing, 5, 500, 4);
  f1_slewing->SetNpx(10000);
  f1_slewing->SetParameters(-13.7915, 33.5238, 3.15088, -0.313885);

  // To reconstruct pulse shapes
  TGraph *gr_pulse = new TGraph();
  TF1 *f1_landau = new TF1("f1_landau", fit_landau, 0, 300, 3);

  // Build the pulse graph, fit it, and fill h2_ADC_dt for one side.
  auto processHits = [&](TTreeReaderArray<int> &phase_arr, TTreeReaderArray<int> &stamp_arr,
                         TTreeReaderArray<unsigned int> &begin_arr, TTreeReaderArray<unsigned int> &end_arr,
                         TTreeReaderArray<uint16_t> &lowADC, TTreeReaderArray<uint16_t> &TOA) {
    for (int i = 0; i < phase_arr.GetSize(); i++) {
      int begin = begin_arr[i];
      int end = end_arr[i];
      double phase = phase_arr[i];
      double stamp = stamp_arr[i];

      int idx_toa = -1;
      int idx_toa_rel = -1;
      double t_rec = 0;
      double adc[n_samples] = {0,};
      double adcSum = 0;

      for (int j = begin; j < end; j++) {
        adcSum += lowADC[j];
        adc[j - begin] = lowADC[j];

        if (TOA[j] > 0) {
          idx_toa = j;
          idx_toa_rel = j - begin;
        }
      }

      if (idx_toa > -1 && adcSum > toa_thres) {
        // Add (t_rec, toa_thres) point from TOA
        t_rec = (phase - TOA[idx_toa]) * (timeWindow / capTOA) +
                (stamp + idx_toa_rel) * timeWindow;
        gr_pulse->Set(0);
        gr_pulse->SetPoint(0, t_rec, toa_thres);

        // Add (t, Npe) points from ADC samples
        for (int j = idx_toa_rel; j < n_samples; j++) {
          if (adc[j] > 0)
            gr_pulse->SetPoint(gr_pulse->GetN(),
                               t_rec + TOA[idx_toa] * (timeWindow / capTOA) +
                                   timeWindow * (j - idx_toa_rel),
                               adc[j] * (dyRangeADC / capADC));
        }

        // Fit pulse shape
        f1_landau->SetParameters(t_rec + 20, pulse_sigma, adcSum * 55.37);
        f1_landau->SetParLimits(1, pulse_sigma - 2, pulse_sigma + 2);
        int status = gr_pulse->Fit("f1_landau", "RBQ0");
        if (status < 0)
          continue;

        double chi2 = f1_landau->GetChisquare();
        double npe_fit = f1_landau->GetMaximum();
        double t_fit = f1_landau->GetX(0.3 * npe_fit, 0, f1_landau->GetParameter(0));
        if (std::isnan(t_fit))
          t_fit = 0;

        if (adcSum > 0 && t_rec > 0 && t_fit > 0 && chi2 < 2) {
          h2_ADC_dt->Fill(adcSum, t_rec - t_fit);
        }
      }
    }
  };

  int nevent = 0;
  while (reader.Next()) {
    nevent++;
    if (nevent % 500 == 0)
      printf(">>> %d\n", nevent);

    processHits(pos_phase, pos_stamp, pos_begin, pos_end, pos_lowADC, pos_TOA);
    processHits(neg_phase, neg_stamp, neg_begin, neg_end, neg_lowADC, neg_TOA);
  }

  TCanvas *canv = new TCanvas("canv", "", 600, 600);
  h2_ADC_dt->Draw();
  h2_ADC_dt->Fit("f1_slewing", "R");
}
