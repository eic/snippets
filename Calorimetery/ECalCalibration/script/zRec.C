void zRec(TString filename) {

  TChain *ch = new TChain("events");
  ch->Add(filename);
  TTreeReader reader(ch);

  // True z position
  TTreeReaderArray<uint64_t> id(reader, "EcalBarrelScFiPAttenuatedHits.cellID");
  TTreeReaderArray<float> z_true(reader, "EcalBarrelScFiPAttenuatedHits.position.z");

  // CALOROCHits(p-going)
  TTreeReaderArray<uint64_t> pos_id(reader, "EcalBarrelScFiPCALOROCHits.cellID");
  TTreeReaderArray<int> pos_phase(reader, "EcalBarrelScFiPCALOROCHits.samplePhase");
  TTreeReaderArray<int> pos_stamp(reader, "EcalBarrelScFiPCALOROCHits.timeStamp");
  TTreeReaderArray<unsigned int> pos_begin(reader, "EcalBarrelScFiPCALOROCHits.bSamples_begin");
  TTreeReaderArray<unsigned int> pos_end(reader, "EcalBarrelScFiPCALOROCHits.bSamples_end");
  TTreeReaderArray<uint16_t> pos_lowADC(reader, "_EcalBarrelScFiPCALOROCHits_bSamples.lowGainADC");
  TTreeReaderArray<uint16_t> pos_TOA(reader, "_EcalBarrelScFiPCALOROCHits_bSamples.timeOfArrival");

  // CALOROCHits(e-going)
  TTreeReaderArray<uint64_t> neg_id(reader, "EcalBarrelScFiNCALOROCHits.cellID");
  TTreeReaderArray<int> neg_phase(reader, "EcalBarrelScFiNCALOROCHits.samplePhase");
  TTreeReaderArray<int> neg_stamp(reader, "EcalBarrelScFiNCALOROCHits.timeStamp");
  TTreeReaderArray<unsigned int> neg_begin(reader, "EcalBarrelScFiNCALOROCHits.bSamples_begin");
  TTreeReaderArray<unsigned int> neg_end(reader, "EcalBarrelScFiNCALOROCHits.bSamples_end");
  TTreeReaderArray<uint16_t> neg_lowADC(reader, "_EcalBarrelScFiNCALOROCHits_bSamples.lowGainADC");
  TTreeReaderArray<uint16_t> neg_TOA(reader, "_EcalBarrelScFiNCALOROCHits_bSamples.timeOfArrival");

  // Config variables
  const double capTOA = 1024;
  const double toa_thres = 7;
  const double timeWindow = 25;

  // z position of the e-going end
  const double z0_neg = -2637.5;
  // Time-walk correction parameters from timeWalkCor.C
  const double pars[4] = {-13.7915, 33.5238, 3.15088, -0.313885};

  // To extract z_true
  TH2D* h2_dt_z = new TH2D("h2_dt_z", "", 200, -20, 30, 100, 500, 4500);
  TF1* f1_dt_z = new TF1("f1_dt_z", "[0]*x+[1]", -15, 25);
  f1_dt_z->SetParameters(83.3221, 2219.58);
  std::unordered_map<uint64_t, double> id_t_pos, id_t_neg, id_z_true;

  // Reconstruct the time-walk-corrected time of each hit and store it by cellID.
  auto processHits = [&](std::unordered_map<uint64_t, double> &id_t,
          		 TTreeReaderArray<uint64_t> &id_arr, 
			 TTreeReaderArray<int> &phase_arr, TTreeReaderArray<int> &stamp_arr,
          		 TTreeReaderArray<unsigned int> &begin_arr, TTreeReaderArray<unsigned int> &end_arr,
          		 TTreeReaderArray<uint16_t> &lowADC, TTreeReaderArray<uint16_t> &TOA) {
    for (int i = 0; i < id_arr.GetSize(); i++) {
      uint64_t cellID = id_arr[i];
      int begin = begin_arr[i];
      int end = end_arr[i];
      double phase = phase_arr[i];
      double stamp = stamp_arr[i];

      int idx_toa = -1;
      int idx_toa_rel = -1;
      double adcSum = 0;

      for (int j = begin; j < end; j++) {
        adcSum += lowADC[j];

        if (TOA[j] > 0) {
          idx_toa = j;
          idx_toa_rel = j - begin;
        }
      }

      if (idx_toa > -1 && adcSum > toa_thres) {
        double t_rec = (phase - TOA[idx_toa]) * (timeWindow / capTOA) +
                       (stamp + idx_toa_rel) * timeWindow;
        id_t[cellID] = t_rec - (pars[1] * pow(adcSum - pars[2], pars[3]) + pars[0]);
      }
    }
  };

  int nevent = 0;
  while (reader.Next()) {
    nevent++;
    if (nevent % 500 == 0)
      printf(">>> %d\n", nevent);

    id_t_pos.clear();
    id_t_neg.clear();
    id_z_true.clear();

    for (int i = 0; i < id.GetSize(); i++)
      id_z_true[id[i]] = z_true[i];

    processHits(id_t_pos, pos_id, pos_phase, pos_stamp, pos_begin, pos_end,
                pos_lowADC, pos_TOA);
    processHits(id_t_neg, neg_id, neg_phase, neg_stamp, neg_begin, neg_end,
                neg_lowADC, neg_TOA);

    for (auto &[cellID, t_pos] : id_t_pos) {
      auto it = id_t_neg.find(cellID);
      if (it != id_t_neg.end())
        h2_dt_z->Fill(it->second - t_pos, id_z_true[cellID] - z0_neg);
    }
  }

  TCanvas *canv = new TCanvas("canv", "", 600, 600);
  h2_dt_z->Draw();
  h2_dt_z->Fit("f1_dt_z", "R");
}
