
// Example script to read jet branch, find associated jet constituents, and confirm that constituents return the jet kinematics
// Author: B. Page (bpage@bnl.gov)
//
// Usage: root -l -q jetReader_TTreeReader.C'("/path/to/eicrecon/output/file")'

void jetReader_TTreeReader(TString infile="")
{
  // Input
  TChain *mychain = new TChain("events");
  mychain->Add(infile);

  // Output
  TFile *ofile = TFile::Open("test.hist.root","RECREATE");

  // TTreeReader
  TTreeReader tree_reader(mychain);

  // Reco Jets
  TTreeReaderArray<int> recoType = {tree_reader, "ReconstructedChargedJets.type"};
  TTreeReaderArray<float> recoNRG = {tree_reader, "ReconstructedChargedJets.energy"};
  TTreeReaderArray<int> recoPDG = {tree_reader, "ReconstructedChargedJets.PDG"};
  TTreeReaderArray<float> recoMomX = {tree_reader, "ReconstructedChargedJets.momentum.x"};
  TTreeReaderArray<float> recoMomY = {tree_reader, "ReconstructedChargedJets.momentum.y"};
  TTreeReaderArray<float> recoMomZ = {tree_reader, "ReconstructedChargedJets.momentum.z"};
  TTreeReaderArray<float> recoM = {tree_reader, "ReconstructedChargedJets.mass"};
  TTreeReaderArray<unsigned int> partsBegin = {tree_reader, "ReconstructedChargedJets.particles_begin"};
  TTreeReaderArray<unsigned int> partsEnd = {tree_reader, "ReconstructedChargedJets.particles_end"};

  TTreeReaderArray<int> recoClustIndex = {tree_reader, "_ReconstructedChargedJets_clusters.index"};
  TTreeReaderArray<int> recoTrackIndex = {tree_reader, "_ReconstructedChargedJets_tracks.index"};
  TTreeReaderArray<int> recoPartIndex = {tree_reader, "_ReconstructedChargedJets_particles.index"};

  // Reconstructed Particles
  TTreeReaderArray<float> recoPartMomX = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> recoPartMomY = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> recoPartMomZ = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> recoPartM = {tree_reader, "ReconstructedChargedParticles.mass"};

  // Generated Jets
  TTreeReaderArray<int> genType = {tree_reader, "GeneratedChargedJets.type"};
  TTreeReaderArray<float> genNRG = {tree_reader, "GeneratedChargedJets.energy"};
  TTreeReaderArray<int> genPDG = {tree_reader, "GeneratedChargedJets.PDG"};
  TTreeReaderArray<float> genMomX = {tree_reader, "GeneratedChargedJets.momentum.x"};
  TTreeReaderArray<float> genMomY = {tree_reader, "GeneratedChargedJets.momentum.y"};
  TTreeReaderArray<float> genMomZ = {tree_reader, "GeneratedChargedJets.momentum.z"};
  TTreeReaderArray<float> genM = {tree_reader, "GeneratedChargedJets.mass"};
  TTreeReaderArray<unsigned int> genPartsBegin = {tree_reader, "GeneratedChargedJets.particles_begin"};
  TTreeReaderArray<unsigned int> genPartsEnd = {tree_reader, "GeneratedChargedJets.particles_end"};
  
  TTreeReaderArray<int> genPartIndex = {tree_reader, "_GeneratedChargedJets_particles.index"};
  
  // MC
  TTreeReaderArray<float> mcMomX = {tree_reader, "GeneratedParticles.momentum.x"};
  TTreeReaderArray<float> mcMomY = {tree_reader, "GeneratedParticles.momentum.y"};
  TTreeReaderArray<float> mcMomZ = {tree_reader, "GeneratedParticles.momentum.z"};
  TTreeReaderArray<float> mcM = {tree_reader, "GeneratedParticles.mass"};
  TTreeReaderArray<int> pdg = {tree_reader, "GeneratedParticles.PDG"};

  // Define Histograms
  // Reco
  TH1D *numRecoJetsEventHist = new TH1D("numRecoJetsEvent","",20,0.,20.);
  TH2D *recoJetEvsEtaHist = new TH2D("recoJetEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaHist = new TH2D("recoJetPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numRecoJetPartsHist = new TH1D("numRecoJetParts","",20,0.,20.);
  TH2D *recoJetPartEvsEtaHist = new TH2D("recoJetPartEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPartPhiVsEtaHist = new TH2D("recoJetPartPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *recoJetPartDeltaRHist = new TH1D("recoJetPartDeltaR","",5000,0.,5.);

  TH2D *recoJetEvsPartESumHist = new TH2D("recoJetEvsPartESum","",3000,0.,300.,3000,0.,300.);
  TH1D *recoJetEDiffHist = new TH1D("recoJetEDiff","",500,-10.,10.);

  TH2D *recoJetEvsEtaBadHist = new TH2D("recoJetEvsEtaBad","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaBadHist = new TH2D("recoJetPhiVsEtaBad","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  // Gen
  TH1D *numGenJetsEventHist = new TH1D("numGenJetsEvent","",20,0.,20.);
  TH2D *genJetEvsEtaHist = new TH2D("genJetEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaHist = new TH2D("genJetPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numGenJetPartsHist = new TH1D("numGenJetParts","",20,0.,20.);
  TH2D *genJetPartEvsEtaHist = new TH2D("genJetPartEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPartPhiVsEtaHist = new TH2D("genJetPartPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *genJetPartDeltaRHist = new TH1D("genJetPartDeltaR","",5000,0.,5.);

  TH2D *genJetEvsPartESumHist = new TH2D("genJetEvsPartESum","",3000,0.,300.,3000,0.,300.);
  TH1D *genJetEDiffHist = new TH1D("genJetEDiff","",500,-10.,10.);

  TH2D *genJetEvsEtaBadHist = new TH2D("genJetEvsEtaBad","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaBadHist = new TH2D("genJetPhiVsEtaBad","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());


  // Loop Through Events
  int NEVENTS = 0;
  while(tree_reader.Next()) {

    if(NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

    // Analyze Reonstructed Jets
    numRecoJetsEventHist->Fill(recoType.GetSize());
    for(unsigned int i=0; i<recoType.GetSize(); i++)
      {
	TVector3 jetMom(recoMomX[i],recoMomY[i],recoMomZ[i]);

	recoJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	recoJetPhiVsEtaHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	double esum = 0.0;
	for(unsigned int j=partsBegin[i]; j<partsEnd[i]; j++)
	  {
	    // partsbegin and partsEnd specify the entries from _ReconstructedChargedJets_particles.index that make up the jet
	    // _ReconstructedChargedJets_particles.index stores the ReconstructedChargedParticles index of the jet constituent
	    double mX = recoPartMomX[recoPartIndex[j]];
	    double mY = recoPartMomY[recoPartIndex[j]];
	    double mZ = recoPartMomZ[recoPartIndex[j]];
	    double mM = recoPartM[recoPartIndex[j]];

	    double tmpE = TMath::Sqrt(mX*mX + mY*mY + mZ*mZ + mM*mM);

	    esum += tmpE;

	    TVector3 partMom(mX,mY,mZ);

	    recoJetPartEvsEtaHist->Fill(partMom.PseudoRapidity(),tmpE);
	    recoJetPartPhiVsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Phi());

	    // Distance between jet axis and particle
	    double dEta = jetMom.PseudoRapidity() - partMom.PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - partMom.Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

	    recoJetPartDeltaRHist->Fill(dR);
	  }

	numRecoJetPartsHist->Fill(partsEnd[i] - partsBegin[i]);
	recoJetEvsPartESumHist->Fill(recoNRG[i],esum);
	recoJetEDiffHist->Fill(recoNRG[i]-esum);
	
	if(TMath::Abs(esum - recoNRG[i]) > 0.00001)
	  {
	    recoJetEvsEtaBadHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	    recoJetPhiVsEtaBadHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	  }
      }


    // Analyze Generated Jets
    numRecoJetsEventHist->Fill(genType.GetSize());
    for(unsigned int i=0; i<genType.GetSize(); i++)
      {
	TVector3 jetMom(genMomX[i],genMomY[i],genMomZ[i]);

	genJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
	genJetPhiVsEtaHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	double esumG = 0.0;
	for(unsigned int j=genPartsBegin[i]; j<genPartsEnd[i]; j++)
	  {
	    double mX = mcMomX[genPartIndex[j]];
	    double mY = mcMomY[genPartIndex[j]];
	    double mZ = mcMomZ[genPartIndex[j]];
	    double mM = mcM[genPartIndex[j]];

	    double tmpE = TMath::Sqrt(mX*mX + mY*mY + mZ*mZ + mM*mM);

	    esumG += tmpE;

	    TVector3 partMom(mX,mY,mZ);

	    genJetPartEvsEtaHist->Fill(partMom.PseudoRapidity(),tmpE);
	    genJetPartPhiVsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Phi());

	    // Distance between jet axis and particle
	    double dEta = jetMom.PseudoRapidity() - partMom.PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - partMom.Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

	    genJetPartDeltaRHist->Fill(dR);
	  }

	numGenJetPartsHist->Fill(genPartsEnd[i] - genPartsBegin[i]);
	genJetEvsPartESumHist->Fill(genNRG[i],esumG);
	genJetEDiffHist->Fill(genNRG[i]-esumG);
	
	if(TMath::Abs(esumG - genNRG[i]) > 0.00001)
	  {
	    genJetEvsEtaBadHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
	    genJetPhiVsEtaBadHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	  }
      }

    NEVENTS++;
  }

  ofile->Write();
  ofile->Close();

}
