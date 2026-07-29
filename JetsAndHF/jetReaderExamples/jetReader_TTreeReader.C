
// Example script to read jet branch, find associated jet constituents, and confirm that constituents return the jet kinematics
// Author: B. Page (bpage@bnl.gov)
//
// Usage: root -l -q jetReader_TTreeReader.C'("/path/to/eicrecon/output/file")'

void jetReader_TTreeReader(TString infile="root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.05.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.18*.eicrecon.edm4eic.root") // 1950
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

  TTreeReaderArray<int> recoPartIndex = {tree_reader, "_ReconstructedChargedJets_particles.index"};

  // Reconstructed Particles
  TTreeReaderArray<float> recoPartMomX = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> recoPartMomY = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> recoPartMomZ = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> recoPartM = {tree_reader, "ReconstructedChargedParticles.mass"};
  TTreeReaderArray<int> recoPartPDG = {tree_reader, "ReconstructedChargedParticles.PDG"};
  TTreeReaderArray<float> recoPartNRG = {tree_reader, "ReconstructedChargedParticles.energy"};

  // Uncomment the following two lines if using an eic-shell version older than 25.12.0-stable
  // Refer https://chat.epic-eic.org/main/pl/n9mbqtf4fiyptjz9q13f7m8xne

  //TTreeReaderArray<unsigned int> recoPartAssocRec = {tree_reader, "ReconstructedChargedParticleAssociations.recID"}; // Reco <-> MCParticle
  //TTreeReaderArray<unsigned int> recoPartAssocSim = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};

  // updated code after eic-shell --version 25.12.0-stable
  TTreeReaderArray<int> recoPartAssocRec = {tree_reader, "_ReconstructedChargedParticleAssociations_rec.index"}; // Reco <-> MCParticle
  TTreeReaderArray<int> recoPartAssocSim = {tree_reader, "_ReconstructedChargedParticleAssociations_sim.index"};

  TTreeReaderArray<float> recoPartAssocWeight = {tree_reader, "ReconstructedChargedParticleAssociations.weight"};

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
  TTreeReaderArray<float> mcE = {tree_reader, "GeneratedParticles.energy"};
  TTreeReaderArray<int> pdg = {tree_reader, "GeneratedParticles.PDG"};

  TTreeReaderArray<int> mcGenStat = {tree_reader, "MCParticles.generatorStatus"};
  TTreeReaderArray<double> mcMomXPart = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomYPart = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomZPart = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcMPart = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> pdgMCPart = {tree_reader, "MCParticles.PDG"};


  // Define Histograms
  // Reco
  TH1D *numRecoJetsEventHist = new TH1D("numRecoJetsEvent","",20,0.,20.);
  TH1D *numRecoJetsNoElecEventHist = new TH1D("numRecoJetsNoElecEvent","",20,0.,20.);

  TH1D *recoJetEHist = new TH1D("recoJetE","",300,0.,300.);
  TH2D *recoJetEvsEtaHist = new TH2D("recoJetEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaHist = new TH2D("recoJetPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *recoJetENoElecHist = new TH1D("recoJetENoElec","",300,0.,300.);
  TH2D *recoJetEvsEtaNoElecHist = new TH2D("recoJetEvsEtaNoElec","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaNoElecHist = new TH2D("recoJetPhiVsEtaNoElec","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *recoJetEElecHist = new TH1D("recoJetEElec","",300,0.,300.);
  TH2D *recoJetEvsEtaElecHist = new TH2D("recoJetEvsEtaElec","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaElecHist = new TH2D("recoJetPhiVsEtaElec","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *constituentRecoPHist = new TH1D("constituentRecoP","",500,0.,50.);
  TH2D *constituentRecoPVsEtaHist = new TH2D("constituentRecoPVsEta","",100,-5.,5.,500,0.,50.);
  TH2D *constituentRecoPhiVsEtaHist = new TH2D("constituentRecoPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *constituentRecoTrueDeltaRHist = new TH1D("constituentRecoTrueDeltaR","",5000,0.,5.);

  TH2D *constituentRecoVsTruePHist = new TH2D("constituentRecoVsTrueP","",500,0.,50.,500,0.,50.);
  TH2D *constituentRecoVsTrueEtaHist = new TH2D("constituentRecoVsTrueEta","",100,-5.,5.,100,-5.,5.);
  TH2D *constituentRecoVsTruePhiHist = new TH2D("constituentRecoVsTruePhi","",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());

  TH1D *constituentPResHist = new TH1D("constituentPRes","",2000,-10.,10.);

  TH2D *constituentRecoPVsEtaNoPIDHist = new TH2D("constituentRecoPVsEtaNoPID","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaNoPIDHist = new TH2D("constituentAssocPartPVsEtaNoPID","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPDGElecHist = new TH1D("constituentPDGElec","",3000,0.,3000.);

  TH1D *constituentPResElecHist = new TH1D("constituentPResElec","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaElecHist = new TH2D("constituentRecoPVsEtaElec","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaElecHist = new TH2D("constituentAssocPartPVsEtaElec","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResNoPIDElecHist = new TH1D("constituentPResNoPIDElec","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaNoPIDElecHist = new TH2D("constituentRecoPVsEtaNoPIDElec","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaNoPIDElecHist = new TH2D("constituentAssocPartPVsEtaNoPIDElec","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResBadPIDElecHist = new TH1D("constituentPResBadPIDElec","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaBadPIDElecHist = new TH2D("constituentRecoPVsEtaBadPIDElec","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaBadPIDElecHist = new TH2D("constituentAssocPartPVsEtaBadPIDElec","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPDGPionHist = new TH1D("constituentPDGPion","",3000,0.,3000.);

  TH1D *constituentPResPionHist = new TH1D("constituentPResPion","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaPionHist = new TH2D("constituentRecoPVsEtaPion","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaPionHist = new TH2D("constituentAssocPartPVsEtaPion","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResNoPIDPionHist = new TH1D("constituentPResNoPIDPion","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaNoPIDPionHist = new TH2D("constituentRecoPVsEtaNoPIDPion","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaNoPIDPionHist = new TH2D("constituentAssocPartPVsEtaNoPIDPion","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResBadPIDPionHist = new TH1D("constituentPResBadPIDPion","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaBadPIDPionHist = new TH2D("constituentRecoPVsEtaBadPIDPion","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaBadPIDPionHist = new TH2D("constituentAssocPartPVsEtaBadPIDPion","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResKaonHist = new TH1D("constituentPResKaon","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaKaonHist = new TH2D("constituentRecoPVsEtaKaon","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaKaonHist = new TH2D("constituentAssocPartPVsEtaKaon","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResNoPIDKaonHist = new TH1D("constituentPResNoPIDKaon","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaNoPIDKaonHist = new TH2D("constituentRecoPVsEtaNoPIDKaon","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaNoPIDKaonHist = new TH2D("constituentAssocPartPVsEtaNoPIDKaon","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResBadPIDKaonHist = new TH1D("constituentPResBadPIDKaon","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaBadPIDKaonHist = new TH2D("constituentRecoPVsEtaBadPIDKaon","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaBadPIDKaonHist = new TH2D("constituentAssocPartPVsEtaBadPIDKaon","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResProtonHist = new TH1D("constituentPResProton","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaProtonHist = new TH2D("constituentRecoPVsEtaProton","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaProtonHist = new TH2D("constituentAssocPartPVsEtaProton","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResNoPIDProtonHist = new TH1D("constituentPResNoPIDProton","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaNoPIDProtonHist = new TH2D("constituentRecoPVsEtaNoPIDProton","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaNoPIDProtonHist = new TH2D("constituentAssocPartPVsEtaNoPIDProton","",100,-5.,5.,500,0.,50.);

  TH1D *constituentPResBadPIDProtonHist = new TH1D("constituentPResBadPIDProton","",2000,-10.,10.);
  TH2D *constituentRecoPVsEtaBadPIDProtonHist = new TH2D("constituentRecoPVsEtaBadPIDProton","",100,-5.,5.,500,0.,50.);
  TH2D *constituentAssocPartPVsEtaBadPIDProtonHist = new TH2D("constituentAssocPartPVsEtaBadPIDProton","",100,-5.,5.,500,0.,50.);

  TH1D *numRecoJetPartsHist = new TH1D("numRecoJetParts","",20,0.,20.);
  TH2D *recoJetEvsPartESumHist = new TH2D("recoJetEvsPartESum","",3000,0.,300.,3000,0.,300.);
  TH1D *recoJetEDiffHist = new TH1D("recoJetEDiff","",20000,-0.001,0.001);

  TH2D *recoJetEvsEtaBadHist = new TH2D("recoJetEvsEtaBad","",100,-5.,5.,300,0.,300.);
  TH2D *recoJetPhiVsEtaBadHist = new TH2D("recoJetPhiVsEtaBad","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  // Gen
  TH1D *numGenJetsEventHist = new TH1D("numGenJetsEvent","",20,0.,20.);
  TH1D *numGenJetsNoElecEventHist = new TH1D("numGenJetsNoElecEvent","",20,0.,20.);

  TH1D *genJetEHist = new TH1D("genJetE","",300,0.,300.);
  TH2D *genJetEvsEtaHist = new TH2D("genJetEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaHist = new TH2D("genJetPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *genJetENoElecHist = new TH1D("genJetENoElec","",300,0.,300.);
  TH2D *genJetEvsEtaNoElecHist = new TH2D("genJetEvsEtaNoElec","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaNoElecHist = new TH2D("genJetPhiVsEtaNoElec","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *genJetEElecHist = new TH1D("genJetEElec","",300,0.,300.);
  TH2D *genJetEvsEtaElecHist = new TH2D("genJetEvsEtaElec","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaElecHist = new TH2D("genJetPhiVsEtaElec","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numGenJetPartsHist = new TH1D("numGenJetParts","",20,0.,20.);
  TH2D *genJetPartEvsEtaHist = new TH2D("genJetPartEvsEta","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPartPhiVsEtaHist = new TH2D("genJetPartPhiVsEta","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());

  TH2D *genJetEvsPartESumHist = new TH2D("genJetEvsPartESum","",3000,0.,300.,3000,0.,300.);
  TH1D *genJetEDiffHist = new TH1D("genJetEDiff","",20000,-0.001,0.001);

  TH2D *genJetEvsEtaBadHist = new TH2D("genJetEvsEtaBad","",100,-5.,5.,300,0.,300.);
  TH2D *genJetPhiVsEtaBadHist = new TH2D("genJetPhiVsEtaBad","",100,-5.,5.,100,-TMath::Pi(),TMath::Pi());


  // Loop Through Events
  int NEVENTS = 0;
  while(tree_reader.Next()) {

    if(NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

    // Analyze Reonstructed Jets
    int numRecoJetsNoElec = 0;
    numRecoJetsEventHist->Fill(recoType.GetSize());
    for(unsigned int i=0; i<recoType.GetSize(); i++)
      {
	TVector3 jetMom(recoMomX[i],recoMomY[i],recoMomZ[i]);

	if(TMath::Abs(jetMom.PseudoRapidity()) > 2.5)
	  continue;

	recoJetEHist->Fill(recoNRG[i]);
	recoJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	recoJetPhiVsEtaHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	
	// Check if Jet Contains an Electron - Use Particle Matching to Find True PID
	bool noElectron = true;
	for(unsigned int m=partsBegin[i]; m<partsEnd[i]; m++) // Loop over jet constituents
	  {
	    int elecIndex = -1;
	    double elecIndexWeight = -1.0;
	    int chargePartIndex = recoPartIndex[m]; // ReconstructedChargedParticle Index for m'th Jet Component
	    for(unsigned int n=0; n<recoPartAssocRec.GetSize(); n++) // Loop Over All ReconstructedChargedParticleAssociations
	      {
		if(recoPartAssocRec[n] == chargePartIndex) // Select Entry Matching the ReconstructedChargedParticle Index
		  {
		    if(recoPartAssocWeight[n] > elecIndexWeight) // Find Particle with Greatest Weight = Contributed Most Hits to Track
		      {
			elecIndex = recoPartAssocSim[n]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
			elecIndexWeight = recoPartAssocWeight[n];
		      }
		  }
	      }

	    if(pdgMCPart[elecIndex] == 11) // Test if Matched Particle is an Electron
	      noElectron = false;
	  }
	
	if(noElectron)
	  {
	    recoJetENoElecHist->Fill(recoNRG[i]);
	    recoJetEvsEtaNoElecHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	    recoJetPhiVsEtaNoElecHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	    numRecoJetsNoElec++;
	  }

	if(!noElectron)
	  {
	    recoJetEElecHist->Fill(recoNRG[i]);
	    recoJetEvsEtaElecHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	    recoJetPhiVsEtaElecHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	  }
	

	// Look at Constituents
	double esum = 0.0;
	double eRes = -999.;
	for(unsigned int j=partsBegin[i]; j<partsEnd[i]; j++)
	  {
	    // partsbegin and partsEnd specify the entries from _ReconstructedChargedJets_particles.index that make up the jet
	    // _ReconstructedChargedJets_particles.index stores the ReconstructedChargedParticles index of the jet constituent

	    TVector3 recoPartMom(recoPartMomX[recoPartIndex[j]],recoPartMomY[recoPartIndex[j]],recoPartMomZ[recoPartIndex[j]]);
	    double mM = recoPartM[recoPartIndex[j]];
	    double mE = recoPartNRG[recoPartIndex[j]];
	    int mPDG = recoPartPDG[recoPartIndex[j]];

	    esum += mE;

	    constituentRecoPHist->Fill(recoPartMom.Mag());
	    constituentRecoPVsEtaHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
	    constituentRecoPhiVsEtaHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Phi());

	    // Find Associated MC Particle (Association with Largest Weight)
	    int simIndex = -1;
	    double simIndexWeight = -1.0;
	    int chargePartIndex = recoPartIndex[j];
	    for(unsigned int n=0; n<recoPartAssocRec.GetSize(); n++)
	      {
		if(recoPartAssocRec[n] == chargePartIndex)
		  {
		    if(recoPartAssocWeight[n] > simIndexWeight)
		      {
			simIndex = recoPartAssocSim[n];
			simIndexWeight = recoPartAssocWeight[n];
		      }
		  }
	      }

	    // Define Matching Truth Particle
	    TVector3 genPartMom(mcMomXPart[simIndex],mcMomYPart[simIndex],mcMomZPart[simIndex]);
	    double mPartM = mcMPart[simIndex];
	    int mPartPDG = pdgMCPart[simIndex];
	    int mPartGenStat = mcGenStat[simIndex];

	    double momRes = (recoPartMom.Mag() - genPartMom.Mag())/genPartMom.Mag();
	    constituentPResHist->Fill(momRes);

	    double dEta = recoPartMom.PseudoRapidity() - genPartMom.PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(recoPartMom.Phi() - genPartMom.Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

	    constituentRecoTrueDeltaRHist->Fill(dR);
	    
	    constituentRecoVsTruePHist->Fill(genPartMom.Mag(),recoPartMom.Mag());
	    constituentRecoVsTrueEtaHist->Fill(genPartMom.PseudoRapidity(),recoPartMom.PseudoRapidity());
	    constituentRecoVsTruePhiHist->Fill(genPartMom.Phi(),recoPartMom.Phi());

	    // Look at PID
	    if(mPDG == 0)
	      {
		constituentRecoPVsEtaNoPIDHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		constituentAssocPartPVsEtaNoPIDHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
	      }

	    if(TMath::Abs(mPartPDG) == 11)
	      {
		constituentPDGElecHist->Fill(TMath::Abs(mPDG));

		constituentPResElecHist->Fill(momRes);
		constituentRecoPVsEtaElecHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		constituentAssocPartPVsEtaElecHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());

		if(mPDG == 0)
		  {
		    constituentPResNoPIDElecHist->Fill(momRes);
		    constituentRecoPVsEtaNoPIDElecHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaNoPIDElecHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }

		if(mPDG != mPartPDG)
		  {
		    constituentPResBadPIDElecHist->Fill(momRes);
		    constituentRecoPVsEtaBadPIDElecHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaBadPIDElecHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }
	      }
	    if(TMath::Abs(mPartPDG) == 211)
	      {
		constituentPDGPionHist->Fill(TMath::Abs(mPDG));

		constituentPResPionHist->Fill(momRes);
		constituentRecoPVsEtaPionHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		constituentAssocPartPVsEtaPionHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());

		if(mPDG == 0)
		  {
		    constituentPResNoPIDPionHist->Fill(momRes);
		    constituentRecoPVsEtaNoPIDPionHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaNoPIDPionHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }

		if(mPDG != mPartPDG)
		  {
		    constituentPResBadPIDPionHist->Fill(momRes);
		    constituentRecoPVsEtaBadPIDPionHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaBadPIDPionHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }
	      }
	    if(TMath::Abs(mPartPDG) == 321)
	      {
		constituentPResKaonHist->Fill(momRes);
		constituentRecoPVsEtaKaonHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		constituentAssocPartPVsEtaKaonHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());

		if(mPDG == 0)
		  {
		    constituentPResNoPIDKaonHist->Fill(momRes);
		    constituentRecoPVsEtaNoPIDKaonHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaNoPIDKaonHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }

		if(mPDG != mPartPDG)
		  {
		    constituentPResBadPIDKaonHist->Fill(momRes);
		    constituentRecoPVsEtaBadPIDKaonHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaBadPIDKaonHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }
	      }
	    if(TMath::Abs(mPartPDG) == 2212)
	      {
		constituentPResProtonHist->Fill(momRes);
		constituentRecoPVsEtaProtonHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		constituentAssocPartPVsEtaProtonHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());

		if(mPDG == 0)
		  {
		    constituentPResNoPIDProtonHist->Fill(momRes);
		    constituentRecoPVsEtaNoPIDProtonHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaNoPIDProtonHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }

		if(mPDG != mPartPDG)
		  {
		    constituentPResBadPIDProtonHist->Fill(momRes);
		    constituentRecoPVsEtaBadPIDProtonHist->Fill(recoPartMom.PseudoRapidity(),recoPartMom.Mag());
		    constituentAssocPartPVsEtaBadPIDProtonHist->Fill(genPartMom.PseudoRapidity(),genPartMom.Mag());
		  }
	      }
	  }

	numRecoJetsNoElecEventHist->Fill(numRecoJetsNoElec);
	numRecoJetPartsHist->Fill(partsEnd[i] - partsBegin[i]);
	recoJetEvsPartESumHist->Fill(recoNRG[i],esum);
	recoJetEDiffHist->Fill(recoNRG[i]-esum);
	
	if(TMath::Abs(esum - recoNRG[i]) > 0.000001)
	  {
	    recoJetEvsEtaBadHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
	    recoJetPhiVsEtaBadHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	  }
      }


    // Analyze Generated Jets
    int numGenJetsNoElec = 0;
    numGenJetsEventHist->Fill(genType.GetSize());
    for(unsigned int i=0; i<genType.GetSize(); i++)
      {
	TVector3 jetMom(genMomX[i],genMomY[i],genMomZ[i]);

	if(TMath::Abs(jetMom.PseudoRapidity()) > 2.5)
	  continue;

	genJetEHist->Fill(genNRG[i]);
	genJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
	genJetPhiVsEtaHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	double esumG = 0.0;
	bool noGenElectron = true;
	for(unsigned int j=genPartsBegin[i]; j<genPartsEnd[i]; j++)
	  {
	    double mX = mcMomX[genPartIndex[j]];
	    double mY = mcMomY[genPartIndex[j]];
	    double mZ = mcMomZ[genPartIndex[j]];
	    double mM = mcM[genPartIndex[j]];
	    int mPDG = pdg[genPartIndex[j]];

	    double tmpE = TMath::Sqrt(mX*mX + mY*mY + mZ*mZ + mM*mM);

	    esumG += tmpE;

	    if(mPDG == 11)
	      noGenElectron = false;

	    TVector3 partMom(mX,mY,mZ);

	    genJetPartEvsEtaHist->Fill(partMom.PseudoRapidity(),tmpE);
	    genJetPartPhiVsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Phi());
	  }

	if(noGenElectron)
	  {
	    genJetENoElecHist->Fill(genNRG[i]);
	    genJetEvsEtaNoElecHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
	    genJetPhiVsEtaNoElecHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	    numGenJetsNoElec++;
	  }

	if(!noGenElectron)
	  {
	    genJetEElecHist->Fill(genNRG[i]);
	    genJetEvsEtaElecHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
	    genJetPhiVsEtaElecHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());
	  }

	numGenJetsNoElecEventHist->Fill(numGenJetsNoElec);
	numGenJetPartsHist->Fill(genPartsEnd[i] - genPartsBegin[i]);
	genJetEvsPartESumHist->Fill(genNRG[i],esumG);
	genJetEDiffHist->Fill(genNRG[i]-esumG);
	
	if(TMath::Abs(esumG - genNRG[i]) > 0.000001)
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
