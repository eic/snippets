// ============================================================================
//! \file    JetValidation.C
//! \authors Brian Page (bpage@bnl.gov),
//!          adapted by Derek Anderson (derek.murphy.anderson@protonmail.com)
// ----------------------------------------------------------------------------
//! \brief Adaption of the Jet Benchmark to run
//!   standalone to generate plots for validation.
//!
//! \usage In eic-shell:
//!     root -b -q JetValidation.C'(<input file list>, \
//!                                 <n files to read>, \
//!                                 <output path>)'
// ============================================================================

#include <edm4eic/EDM4eicVersion.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TVector3.h>
#include <fstream>
#include <string>
#include <vector>

#include "fmt/color.h"
#include "fmt/core.h"

///! Default input file list
const std::string DefaultInFileList = "filelists/files26060.py8ncdis10x100q100t1000.list";

///! Default no. of files
const std::size_t DefaultNFiles = 1000;

///! Default output file path
const std::string DefaultOutPath = ".";

// ----------------------------------------------------------------------------
// Does a branch exist?
// ----------------------------------------------------------------------------
/*! Checks if a branch exists in a TTree.
 *!
 *! \param[in] tree The tree to check
 *! \param[in] branch The branch to check for
 */
bool branchExists(TTree* tree, const std::string& branch) {
  bool exists = false;
  if(tree->GetBranch(branch.c_str())) {
    exists = true;
  }
  return exists;
}

// ----------------------------------------------------------------------------
// Macro body
// ----------------------------------------------------------------------------
/*! Process input files to generate a set of reconstructed,
 *! generated jet distributions and save them as PNGs.
 *!
 *! \param[in]  filelist     Input filelist to use
 *! \param[in]  n_files      Number of files to read
 *! \param[out] results_path Location to save PNGs to
 */
int JetValidation(
  const std::string& filelist = DefaultInFileList,
  const std::size_t n_files = DefaultNFiles,
  const std::string& results_path = DefaultOutPath
) {

  const bool PRINT = true;

  std::ifstream mylist(filelist);
  if (!mylist.is_open()) {
    std::cerr << "PANIC: Couldn't open fileslist '" << filelist << "'!" << std::endl;
    return 1;
  }

  // Load input files
  std::string file;
  std::size_t i_file = 0;
  std::vector<std::string> rec_files;
  while (std::getline(mylist, file)) {
    rec_files.push_back(file);
    ++i_file;
    if (i_file == n_files) {
      break;
    }
  }

  // Input
  TChain *mychain = new TChain("events");
  for (const auto& rec_file : rec_files) {
    mychain->Add(rec_file.c_str());
  }

  const int seabornRed = TColor::GetColor(213, 94, 0);

  // Seaborn Green: #009E73 -> (0, 158, 115)
  const int seabornGreen = TColor::GetColor(0, 158, 115);

  // Seaborn Blue: #56B4E9 -> (86, 180, 233)
  const int seabornBlue = TColor::GetColor(100, 149, 237);

  // TTreeReader
  TTreeReader tree_reader(mychain);

  // Set Delta R Cut
  float DELTARCUT = 0.05;

  // Check if area branch exists
  //   --> Using jet EDM if it does!
  bool useNewEDM = false;
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8,9,0)
  {
    auto file = TFile::Open(rec_files.front().c_str(), "READ");
    auto tree = file->Get<TTree>("events");
    bool hasRecoArea = branchExists(tree, "ReconstructedChargedJets.area");
    bool hasGenArea = branchExists(tree, "GeneratedChargedJets.area");
    useNewEDM = hasRecoArea && hasGenArea;
  }
#endif

  // Set area branches to dummy values if not using jet EDM
  //   --> These branches won't be used if not using
  //       jet EDM
  std::string recoAreaBranch = "ReconstructedChargedJets.energy";
  std::string genAreaBranch = "GeneratedChargedJets.energy";
  if(useNewEDM) {
    recoAreaBranch = "ReconstructedChargedJets.area";
    genAreaBranch = "GeneratedChargedJets.area";
  }

  // Use updated constituent branch names if using jet EDM
  std::string recoCstsBeginBranch = "ReconstructedChargedJets.particles_begin";
  std::string recoCstsEndBranch = "ReconstructedChargedJets.particles_end";
  std::string recoCstIndexBranch = "_ReconstructedChargedJets_particles.index";
  if (useNewEDM) {
    recoCstsBeginBranch = "ReconstructedChargedJets.constituents_begin";
    recoCstsEndBranch = "ReconstructedChargedJets.constituents_end";
    recoCstIndexBranch = "_ReconstructedChargedJets_constituents.index";
  }

  std::string genCstsBeginBranch = "GeneratedChargedJets.particles_begin";
  std::string genCstsEndBranch = "GeneratedChargedJets.particles_end";
  std::string genCstIndexBranch = "_GeneratedChargedJets_particles.index";
  if (useNewEDM) {
    genCstsBeginBranch = "GeneratedChargedJets.constituents_begin";
    genCstsEndBranch = "GeneratedChargedJets.constituents_end";
    genCstIndexBranch = "_GeneratedChargedJets_constituents.index";
  }

  // Reco Jets
  TTreeReaderArray<int> recoType = {tree_reader, "ReconstructedChargedJets.type"};
  TTreeReaderArray<float> recoNRG = {tree_reader, "ReconstructedChargedJets.energy"};
  TTreeReaderArray<float> recoMomX = {tree_reader, "ReconstructedChargedJets.momentum.x"};
  TTreeReaderArray<float> recoMomY = {tree_reader, "ReconstructedChargedJets.momentum.y"};
  TTreeReaderArray<float> recoMomZ = {tree_reader, "ReconstructedChargedJets.momentum.z"};
  TTreeReaderArray<float> recoArea = {tree_reader, recoAreaBranch.c_str()};
  TTreeReaderArray<unsigned int> recoCstsBegin = {tree_reader, recoCstsBeginBranch.c_str()};
  TTreeReaderArray<unsigned int> recoCstsEnd = {tree_reader, recoCstsEndBranch.c_str()};

  TTreeReaderArray<int> recoCstIndex = {tree_reader, recoCstIndexBranch.c_str()};

  // Reconstructed Particles
  TTreeReaderArray<float> recoPartMomX = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> recoPartMomY = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> recoPartMomZ = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> recoPartM = {tree_reader, "ReconstructedChargedParticles.mass"};
  TTreeReaderArray<int> recoPartPDG = {tree_reader, "ReconstructedChargedParticles.PDG"};
  TTreeReaderArray<float> recoPartNRG = {tree_reader, "ReconstructedChargedParticles.energy"};

  TTreeReaderArray<int> recoPartAssocRec = {tree_reader, "_ReconstructedChargedParticleAssociations_rec.index"}; // Reco <-> MCParticle
  TTreeReaderArray<int> recoPartAssocSim = {tree_reader, "_ReconstructedChargedParticleAssociations_sim.index"};
  TTreeReaderArray<float> recoPartAssocWeight = {tree_reader, "ReconstructedChargedParticleAssociations.weight"};

  // Generated Jets
  TTreeReaderArray<int> genType = {tree_reader, "GeneratedChargedJets.type"};
  TTreeReaderArray<float> genNRG = {tree_reader, "GeneratedChargedJets.energy"};
  TTreeReaderArray<float> genMomX = {tree_reader, "GeneratedChargedJets.momentum.x"};
  TTreeReaderArray<float> genMomY = {tree_reader, "GeneratedChargedJets.momentum.y"};
  TTreeReaderArray<float> genMomZ = {tree_reader, "GeneratedChargedJets.momentum.z"};
  TTreeReaderArray<float> genArea = {tree_reader, genAreaBranch.c_str()};
  TTreeReaderArray<unsigned int> genCstsBegin = {tree_reader, genCstsBeginBranch.c_str()};
  TTreeReaderArray<unsigned int> genCstsEnd = {tree_reader, genCstsEndBranch.c_str()};

  TTreeReaderArray<int> genPartIndex = {tree_reader, genCstIndexBranch.c_str()};
  //TTreeReaderArray<int> genChargedIndex = {tree_reader, "GeneratedChargedParticles_objIdx.index"};
  
  // MC
  //TTreeReaderArray<int> mcGenStat = {tree_reader, "MCParticles.generatorStatus"};
  TTreeReaderArray<float> mcMomX = {tree_reader, "GeneratedParticles.momentum.x"};
  TTreeReaderArray<float> mcMomY = {tree_reader, "GeneratedParticles.momentum.y"};
  TTreeReaderArray<float> mcMomZ = {tree_reader, "GeneratedParticles.momentum.z"};
  TTreeReaderArray<float> mcM = {tree_reader, "GeneratedParticles.mass"};
  TTreeReaderArray<int> pdg = {tree_reader, "GeneratedParticles.PDG"};

  TTreeReaderArray<int> mcGenStat = {tree_reader, "MCParticles.generatorStatus"};
  TTreeReaderArray<double> mcMomXPart = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomYPart = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomZPart = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcMPart = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> pdgMCPart = {tree_reader, "MCParticles.PDG"};

  // Define Histograms
  TH1D *counter = new TH1D("counter","",10,0.,10.);

  
  // Reco
  TH1D *numRecoChargedJetsECutHist = new TH1D("numRecoChargedJetsECut","",20,0.,20.);
  TH1D *recoChargedJetEHist = new TH1D("recoChargedJetE","",300,0.,300.);
  TH1D *recoChargedJetEtaECutHist = new TH1D("recoChargedJetEtaECut","",60,-3.,3.);
  TH1D *recoChargedJetAreaECutHist = nullptr;
  TH2D *recoChargedJetEvsAreaHist = nullptr;
  if (useNewEDM) {
    recoChargedJetAreaECutHist = new TH1D("recoChargedJetAreaECut","",250,0.,5.);
    recoChargedJetEvsAreaHist = new TH2D("recoChargedJetEvsArea","",250,0.,5.,300,0.,300.);
  }
  TH2D *recoChargedJetEvsEtaHist = new TH2D("recoChargedJetEvsEta","",60,-3.,3.,300,0.,300.);
  TH2D *recoChargedJetPhiVsEtaECutHist = new TH2D("recoChargedJetPhiVsEtaECut","",60,-3.,3.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numRecoChargedJetsECutNoElecHist = new TH1D("numRecoChargedJetsECutNoElec","",20,0.,20.);
  TH1D *recoChargedJetENoElecHist = new TH1D("recoChargedJetENoElec","",300,0.,300.);
  TH1D *recoChargedJetEtaECutNoElecHist = new TH1D("recoChargedJetEtaECutNoElec","",60,-3.,3.);
  TH1D *recoChargedJetAreaECutNoElecHist = nullptr;
  TH2D *recoChargedJetEvsAreaNoElecHist = nullptr;
  if (useNewEDM) {
    recoChargedJetAreaECutNoElecHist = new TH1D("recoHargedJetAreaECutNoElec","",250,0.,5.);
    recoChargedJetEvsAreaNoElecHist = new TH2D("recoChargedJetEvsAreaNoElec","",250,0.,5.,300,0.,300.);
  }
  TH2D *recoChargedJetEvsEtaNoElecHist = new TH2D("recoChargedJetEvsEtaNoElec","",60,-3.,3.,300,0.,300.);
  TH2D *recoChargedJetPhiVsEtaECutNoElecHist = new TH2D("recoChargedJetPhiVsEtaECutNoElec","",60,-3.,3.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numRecoChargedJetPartsHist = new TH1D("numRecoChargedJetParts","",20,0.,20.);
  TH1D *recoChargedJetPartPHist = new TH1D("recoChargedJetPartP","",500,0.,100.);
  TH1D *recoChargedJetPartEtaHist = new TH1D("recoChargedJetPartEta","",80,-4.,4.);
  TH2D *recoChargedJetPartPvsEtaHist = new TH2D("recoChargedJetPartPvsEta","",80,-4.,4.,500,0.,100.);
  TH2D *recoChargedJetPartPhiVsEtaHist = new TH2D("recoChargedJetPartPhiVsEta","",80,-4.,4.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numRecoChargedJetPartsNoElecHist = new TH1D("numRecoChargedJetPartsNoElec","",20,0.,20.);
  TH1D *recoChargedJetPartPNoElecHist = new TH1D("recoChargedJetPartPNoElec","",500,0.,100.);
  TH1D *recoChargedJetPartEtaNoElecHist = new TH1D("recoChargedJetPartEtaNoElec","",80,-4.,4.);
  TH2D *recoChargedJetPartPvsEtaNoElecHist = new TH2D("recoChargedJetPartPvsEtaNoElec","",80,-4.,4.,500,0.,100.);
  TH2D *recoChargedJetPartPhiVsEtaNoElecHist = new TH2D("recoChargedJetPartPhiVsEtaNoElec","",80,-4.,4.,100,-TMath::Pi(),TMath::Pi());

  TH1D *recoChargedJetPartPairwiseDeltaRHist = new TH1D("recoChargedJetPartPairwiseDeltaRHist","",5000,0.,5.);

  // Gen
  TH1D *numGenChargedJetsECutHist = new TH1D("numGenChargedJetsECut","",20,0.,20.);
  TH1D *genChargedJetEHist = new TH1D("genChargedJetE","",300,0.,300.);
  TH1D *genChargedJetEtaECutHist = new TH1D("genChargedJetEtaECut","",60,-3.,3.);
  TH1D *genChargedJetAreaECutHist = nullptr;
  TH2D *genChargedJetEvsAreaHist = nullptr;
  if (useNewEDM) {
    genChargedJetAreaECutHist = new TH1D("genChargedJetAreaECut","",250,0.,5.);
    genChargedJetEvsAreaHist = new TH2D("genChargedJetEvsAreaHist","",250,0.,5.,300,0.,300.);
  }
  TH2D *genChargedJetEvsEtaHist = new TH2D("genChargedJetEvsEta","",60,-3.,3.,300,0.,300.);
  TH2D *genChargedJetPhiVsEtaECutHist = new TH2D("genChargedJetPhiVsEtaECut","",60,-3.,3.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numGenChargedJetsECutNoElecHist = new TH1D("numGenChargedJetsECutNoElec","",20,0.,20.);
  TH1D *genChargedJetENoElecHist = new TH1D("genChargedJetENoElec","",300,0.,300.);
  TH1D *genChargedJetEtaECutNoElecHist = new TH1D("genChargedJetEtaECutNoElec","",60,-3.,3.);
  TH1D *genChargedJetAreaECutNoElecHist = nullptr;
  TH2D *genChargedJetEvsAreaNoElecHist = nullptr;
  if (useNewEDM) {
    genChargedJetAreaECutNoElecHist = new TH1D("genChargedJetAreaECutNoElec","",250,0.,5.);
    genChargedJetEvsAreaNoElecHist = new TH2D("genChargedJetEvsAreaNoElec","",250,0.,5.,300,0.,300.);
  }
  TH2D *genChargedJetEvsEtaNoElecHist = new TH2D("genChargedJetEvsEtaNoElec","",60,-3.,3.,300,0.,300.);
  TH2D *genChargedJetPhiVsEtaECutNoElecHist = new TH2D("genChargedJetPhiVsEtaECutNoElec","",60,-3.,3.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numGenChargedJetPartsHist = new TH1D("numGenChargedJetParts","",20,0.,20.);
  TH1D *genChargedJetPartPHist = new TH1D("genChargedJetPartP","",500,0.,100.);
  TH1D *genChargedJetPartEtaHist = new TH1D("genChargedJetPartEta","",80,-4.,4.);
  TH2D *genChargedJetPartPvsEtaHist = new TH2D("genChargedJetPartPvsEta","",80,-4.,4.,500,0.,100.);
  TH2D *genChargedJetPartPhiVsEtaHist = new TH2D("genChargedJetPartPhiVsEta","",80,-4.,4.,100,-TMath::Pi(),TMath::Pi());

  TH1D *numGenChargedJetPartsNoElecHist = new TH1D("numGenChargedJetPartsNoElec","",20,0.,20.);
  TH1D *genChargedJetPartPNoElecHist = new TH1D("genChargedJetPartPNoElec","",500,0.,100.);
  TH1D *genChargedJetPartEtaNoElecHist = new TH1D("genChargedJetPartEtaNoElec","",80,-4.,4.);
  TH2D *genChargedJetPartPvsEtaNoElecHist = new TH2D("genChargedJetPartPvsEtaNoElec","",80,-4.,4.,500,0.,100.);
  TH2D *genChargedJetPartPhiVsEtaNoElecHist = new TH2D("genChargedJetPartPhiVsEtaNoElec","",80,-4.,4.,100,-TMath::Pi(),TMath::Pi());

  TH1D *genChargedJetPartPairwiseDeltaRHist = new TH1D("genChargedJetPartPairwiseDeltaRHist","",5000,0.,5.);

  // Matched
  TH1D *matchJetDeltaRHist = new TH1D("matchJetDeltaR","",5000,0.,5.);
  TH1D *matchJetDeltaRBackHist = new TH1D("matchJetDeltaRBack","",5000,0.,5.);
  TH2D *recoVsGenChargedJetEtaHist = new TH2D("recoVsGenChargedJetEta","",80,-4.,4.,80,-4.,4.);
  TH2D *recoVsGenChargedJetPhiHist = new TH2D("recoVsGenChargedJetPhi","",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  TH2D *recoVsGenChargedJetAreaHist = nullptr;
  if (useNewEDM) {
    recoVsGenChargedJetAreaHist = new TH2D("recoVsGenChargedJetArea","",250,0.,5.,250,0.,5.);
  }
  TH2D *recoVsGenChargedJetEHist = new TH2D("recoVsGenChargedJetE","",100,0.,100.,100,0.,100.);
  TH2D *recoVsGenChargedJetENoDRHist = new TH2D("recoVsGenChargedJetENoDRHist","",100,0.,100.,100,0.,100.);
  TH2D *recoVsGenChargedJetENoDupHist = new TH2D("recoVsGenChargedJetENoDup","",100,0.,100.,100,0.,100.);

  TH2D *jetResVsEtaHist = new TH2D("jetResVsEta","",80,-4.,4.,10000,-10.,10.);
  TH2D *jetResVsEHist = new TH2D("jetResVsE","",100,0.,100.,10000,-10.,10.);
  TH2D *jetResVsENegEtaHist = new TH2D("jetResVsENegEta","",20,0.,100.,10000,-10.,10.);
  TH2D *jetResVsEMidEtaHist = new TH2D("jetResVsEMidEta","",20,0.,100.,10000,-10.,10.);
  TH2D *jetResVsEPosEtaHist = new TH2D("jetResVsEPosEta","",20,0.,100.,10000,-10.,10.);

  TH2D *jetResVsENegEtaNoDupHist = new TH2D("jetResVsENegEtaNoDup","",20,0.,100.,10000,-10.,10.);
  TH2D *jetResVsEMidEtaNoDupHist = new TH2D("jetResVsEMidEtaNoDup","",20,0.,100.,10000,-10.,10.);
  TH2D *jetResVsEPosEtaNoDupHist = new TH2D("jetResVsEPosEtaNoDup","",20,0.,100.,10000,-10.,10.);


  // Loop Through Events
  int NEVENTS = 0;
  while(tree_reader.Next()) {

    if(NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

    counter->Fill(0);

    //////////////////////////////////////////////////////////////////////////
    //////////////////////  Analyze Reconstructed Jets  //////////////////////
    //////////////////////////////////////////////////////////////////////////
    int numRecoChargedJets = 0;
    int numRecoChargedJetsNoElec = 0;
    for(unsigned int i=0; i<recoType.GetSize(); i++)
      {
	TVector3 jetMom(recoMomX[i],recoMomY[i],recoMomZ[i]);

	counter->Fill(3);

	// Place eta cut to avoid edges of tracking acceptance
	if(TMath::Abs(jetMom.PseudoRapidity()) > 2.5) continue;

	// Place a minimum energy condition for several plots
	bool ECut = recoNRG[i] > 5.0;

	if(ECut) numRecoChargedJets++; 

	recoChargedJetEHist->Fill(recoNRG[i]);
	if(ECut) recoChargedJetEtaECutHist->Fill(jetMom.PseudoRapidity());
	recoChargedJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
        if(useNewEDM) {
          if(ECut) recoChargedJetAreaECutHist->Fill(recoArea[i]);
          recoChargedJetEvsAreaHist->Fill(recoArea[i],recoNRG[i]);
        }
	if(ECut) recoChargedJetPhiVsEtaECutHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	// Find Jets with Electrons
	bool noElectron = true;
	for(unsigned int m=recoCstsBegin[i]; m<recoCstsEnd[i]; m++) // Loop over jet constituents
	  {
	    int elecIndex = -1;
	    double elecIndexWeight = -1.0;
	    int chargePartIndex = recoCstIndex[m]; // ReconstructedChargedParticle Index for m'th Jet Component
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
	
	if(ECut)
	  {
	    for(unsigned int j=recoCstsBegin[i]; j<recoCstsEnd[i]; j++)
	      {
		// recoCstsBegin and recoCstsEnd specify the entries from _ReconstructedChargedJets_particles.index that make up the jet
		// _ReconstructedChargedJets_particles.index stores the ReconstructedChargedParticles index of the jet constituent
		double mX = recoPartMomX[recoCstIndex[j]];
		double mY = recoPartMomY[recoCstIndex[j]];
		double mZ = recoPartMomZ[recoCstIndex[j]];
		double mM = recoPartM[recoCstIndex[j]];
		//double tmpE = TMath::Sqrt(mX*mX + mY*mY + mZ*mZ + mM*mM);
		
		TVector3 partMom(mX,mY,mZ);
		
		recoChargedJetPartPHist->Fill(partMom.Mag());
		recoChargedJetPartEtaHist->Fill(partMom.PseudoRapidity());
		recoChargedJetPartPvsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Mag());
		recoChargedJetPartPhiVsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Phi());

		if(noElectron)
		  {
		    recoChargedJetPartPNoElecHist->Fill(partMom.Mag());
		    recoChargedJetPartEtaNoElecHist->Fill(partMom.PseudoRapidity());
		    recoChargedJetPartPvsEtaNoElecHist->Fill(partMom.PseudoRapidity(),partMom.Mag());
		    recoChargedJetPartPhiVsEtaNoElecHist->Fill(partMom.PseudoRapidity(),partMom.Phi());
		  }

		// Pairwise Distance Between Constituents
		if(j<(recoCstsEnd[i]-1))
		  {
		    for(unsigned int k=j+1; k<recoCstsEnd[i]; k++)
		      {
			double mXB = recoPartMomX[recoCstIndex[k]];
			double mYB = recoPartMomY[recoCstIndex[k]];
			double mZB = recoPartMomZ[recoCstIndex[k]];

			TVector3 partMomB(mXB,mYB,mZB);

			double dEta = partMom.PseudoRapidity() - partMomB.PseudoRapidity();
			double dPhi = TVector2::Phi_mpi_pi(partMom.Phi() - partMomB.Phi());
			double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

			recoChargedJetPartPairwiseDeltaRHist->Fill(dR);
		      }
		  }
	      }
	    numRecoChargedJetPartsHist->Fill(recoCstsEnd[i] - recoCstsBegin[i]);
	    if(noElectron) numRecoChargedJetPartsNoElecHist->Fill(recoCstsEnd[i] - recoCstsBegin[i]);
	  }

	// No Electrons
	if(noElectron)
	  {
	    recoChargedJetENoElecHist->Fill(recoNRG[i]);
	    if(ECut) recoChargedJetEtaECutNoElecHist->Fill(jetMom.PseudoRapidity());
	    recoChargedJetEvsEtaNoElecHist->Fill(jetMom.PseudoRapidity(),recoNRG[i]);
            if(useNewEDM) {
              if(ECut) recoChargedJetAreaECutNoElecHist->Fill(recoArea[i]);
              recoChargedJetEvsAreaNoElecHist->Fill(recoArea[i],recoNRG[i]);
            }
	    if(ECut) recoChargedJetPhiVsEtaECutNoElecHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	    if(ECut) numRecoChargedJetsNoElec++; 
	  }
      }
    numRecoChargedJetsECutHist->Fill(numRecoChargedJets);
    numRecoChargedJetsECutNoElecHist->Fill(numRecoChargedJetsNoElec);

    //////////////////////////////////////////////////////////////////////////
    ////////////////////////  Analyze Generator Jets  ////////////////////////
    //////////////////////////////////////////////////////////////////////////
    int numGenChargedJets = 0;
    int numGenChargedJetsNoElec = 0;
    for(unsigned int i=0; i<genType.GetSize(); i++)
      {
	TVector3 jetMom(genMomX[i],genMomY[i],genMomZ[i]);

	counter->Fill(4);

	// Place eta cut to avoid edges of tracking acceptance
	if(TMath::Abs(jetMom.PseudoRapidity()) > 2.5) continue;

	// Place a minimum energy condition for several plots
	bool ECut = genNRG[i] > 5.0;

	if(ECut) numGenChargedJets++; 

	genChargedJetEHist->Fill(genNRG[i]);
	if(ECut) genChargedJetEtaECutHist->Fill(jetMom.PseudoRapidity());
	genChargedJetEvsEtaHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
        if(useNewEDM) {
          if(ECut) genChargedJetAreaECutHist->Fill(genArea[i]);
          genChargedJetEvsAreaHist->Fill(genArea[i],genNRG[i]);
        }
	if(ECut) genChargedJetPhiVsEtaECutHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	// Find Jets with Electrons
	bool noElectron = true;
	for(unsigned int m=genCstsBegin[i]; m<genCstsEnd[i]; m++)
	  {
	    if(pdg[genPartIndex[m]] == 11)
	      noElectron = false;
	  }

	if(ECut)
	  {
	    for(unsigned int j=genCstsBegin[i]; j<genCstsEnd[i]; j++)
	      {
		// genCstsBegin and genCstsEnd specify the entries from _GeneratedChargedJets_particles.index that make up the jet
		// _GeneratedChargedJets_particles.index stores the GeneratedChargedParticles index of the jet constituent
		double mX = mcMomX[genPartIndex[j]];
		double mY = mcMomY[genPartIndex[j]];
		double mZ = mcMomZ[genPartIndex[j]];
		double mM = mcM[genPartIndex[j]];
		//double tmpE = TMath::Sqrt(mX*mX + mY*mY + mZ*mZ + mM*mM);
		
		TVector3 partMom(mX,mY,mZ);
		
		genChargedJetPartPHist->Fill(partMom.Mag());
		genChargedJetPartEtaHist->Fill(partMom.PseudoRapidity());
		genChargedJetPartPvsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Mag());
		genChargedJetPartPhiVsEtaHist->Fill(partMom.PseudoRapidity(),partMom.Phi());

		if(noElectron)
		  {
		    genChargedJetPartPNoElecHist->Fill(partMom.Mag());
		    genChargedJetPartEtaNoElecHist->Fill(partMom.PseudoRapidity());
		    genChargedJetPartPvsEtaNoElecHist->Fill(partMom.PseudoRapidity(),partMom.Mag());
		    genChargedJetPartPhiVsEtaNoElecHist->Fill(partMom.PseudoRapidity(),partMom.Phi());
		  }

		// Pairwise Distance Between Constituents
		if(j<(genCstsEnd[i]-1))
		  {
		    for(unsigned int k=j+1; k<genCstsEnd[i]; k++)
		      {
			double mXB = mcMomX[genPartIndex[k]];
			double mYB = mcMomY[genPartIndex[k]];
			double mZB = mcMomZ[genPartIndex[k]];

			TVector3 partMomB(mXB,mYB,mZB);

			double dEta = partMom.PseudoRapidity() - partMomB.PseudoRapidity();
			double dPhi = TVector2::Phi_mpi_pi(partMom.Phi() - partMomB.Phi());
			double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

			genChargedJetPartPairwiseDeltaRHist->Fill(dR);
		      }
		  }
	      }
	    numGenChargedJetPartsHist->Fill(genCstsEnd[i] - genCstsBegin[i]);
	    if(noElectron) numGenChargedJetPartsNoElecHist->Fill(genCstsEnd[i] - genCstsBegin[i]);
	  }

	// No Electrons
	if(noElectron)
	  {
	    genChargedJetENoElecHist->Fill(genNRG[i]);
	    if(ECut) genChargedJetEtaECutNoElecHist->Fill(jetMom.PseudoRapidity());
	    genChargedJetEvsEtaNoElecHist->Fill(jetMom.PseudoRapidity(),genNRG[i]);
            if (useNewEDM) {
              if(ECut) genChargedJetAreaECutNoElecHist->Fill(genArea[i]);
              genChargedJetEvsAreaHist->Fill(genArea[i],genNRG[i]);
            }
	    if(ECut) genChargedJetPhiVsEtaECutNoElecHist->Fill(jetMom.PseudoRapidity(),jetMom.Phi());

	    if(ECut) numGenChargedJetsNoElec++; 
	  }
      }
    numGenChargedJetsECutHist->Fill(numGenChargedJets);
    numGenChargedJetsECutNoElecHist->Fill(numGenChargedJetsNoElec);

    
    //////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Matched Jets  /////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    for(unsigned int i=0; i<genType.GetSize(); i++)
      {
	TVector3 jetMom(genMomX[i],genMomY[i],genMomZ[i]);

	// Place eta cut to avoid edges of tracking acceptance
	//if(TMath::Abs(jetMom.PseudoRapidity()) > 2.5) continue;

	// Place a minimum energy condition
	//if(genNRG[i] < 5.0) continue;
	
	// Don't Look at Electron Jets
	bool hasElectron = false;
	// Find Jets with Electrons
	for(unsigned int m=genCstsBegin[i]; m<genCstsEnd[i]; m++)
	  {
	    if(pdg[genPartIndex[m]] == 11)
	      hasElectron = true;
	  }
	//if(hasElectron) continue;

	// Find Matching Reconstructed Jet
	double minDeltaR = 999.;
	int minIndex = -1;
	for(unsigned int j=0; j<recoType.GetSize(); j++)
	  {
	    TVector3 recoMom(recoMomX[j],recoMomY[j],recoMomZ[j]);

	    double dEta = jetMom.PseudoRapidity() - recoMom.PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(jetMom.Phi() - recoMom.Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

	    if(dR < minDeltaR)
	      {
		minDeltaR = dR;
		minIndex = j;
	      }
	  }

	// Do Backwards Match
	double minDeltaRBack = 999.;
	double minIndexBack = -1;
	if(minIndex > -1)
	  {
	    TVector3 recoMatchMom(recoMomX[minIndex],recoMomY[minIndex],recoMomZ[minIndex]);
	    for(unsigned int j=0; j<genType.GetSize(); j++)
	      {
		TVector3 genMom(genMomX[j],genMomY[j],genMomZ[j]);
		
		double dEta = recoMatchMom.PseudoRapidity() - genMom.PseudoRapidity();
		double dPhi = TVector2::Phi_mpi_pi(recoMatchMom.Phi() - genMom.Phi());
		double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
		
		if(dR < minDeltaRBack)
		  {
		    minDeltaRBack = dR;
		    minIndexBack = j;
		  }
	      }
	  }

	// Look at Best Match
	if(genNRG[i] > 5.0 && TMath::Abs(jetMom.PseudoRapidity()) < 2.5 && minIndex > -1 && !hasElectron) matchJetDeltaRHist->Fill(minDeltaR);
	if(genNRG[i] > 5.0 && TMath::Abs(jetMom.PseudoRapidity()) < 2.5 && minIndex > -1) matchJetDeltaRBackHist->Fill(minDeltaR);
	if(minIndex > -1 && genNRG[i] > 5.0 && TMath::Abs(jetMom.PseudoRapidity()) < 2.5 && !hasElectron)
	  {
	    TVector3 recoMatchMom(recoMomX[minIndex],recoMomY[minIndex],recoMomZ[minIndex]);

	    recoVsGenChargedJetENoDRHist->Fill(genNRG[i],recoNRG[minIndex]);

	    if(minDeltaR < DELTARCUT)
	      {
		recoVsGenChargedJetEtaHist->Fill(jetMom.PseudoRapidity(),recoMatchMom.PseudoRapidity());
		recoVsGenChargedJetPhiHist->Fill(jetMom.Phi(),recoMatchMom.Phi());
                if(useNewEDM) {
                  recoVsGenChargedJetAreaHist->Fill(genArea[i],recoArea[minIndex]);
                }
		recoVsGenChargedJetEHist->Fill(genNRG[i],recoNRG[minIndex]);
		
		double jetERes = (recoNRG[minIndex] - genNRG[i])/genNRG[i];
		
		jetResVsEtaHist->Fill(jetMom.PseudoRapidity(),jetERes);
		jetResVsEHist->Fill(genNRG[i],jetERes);
		if(jetMom.PseudoRapidity() > -2.5 && jetMom.PseudoRapidity() < -1.0)
		  jetResVsENegEtaHist->Fill(genNRG[i],jetERes);
		if(jetMom.PseudoRapidity() > -1.0 && jetMom.PseudoRapidity() < 1.0)
		  jetResVsEMidEtaHist->Fill(genNRG[i],jetERes);
		if(jetMom.PseudoRapidity() > 1.0 && jetMom.PseudoRapidity() < 2.5)
		  jetResVsEPosEtaHist->Fill(genNRG[i],jetERes);
		
		// Check for Duplicate Tracks
		bool noDuplicate = true;
		for(unsigned int j=recoCstsBegin[minIndex]; j<recoCstsEnd[minIndex]; j++)
		  {
		    double mX = recoPartMomX[recoCstIndex[j]];
		    double mY = recoPartMomY[recoCstIndex[j]];
		    double mZ = recoPartMomZ[recoCstIndex[j]];
		    double mM = recoPartM[recoCstIndex[j]];
		    double tmpE = TMath::Sqrt(mX*mX + mY*mY + mZ*mZ + mM*mM);
		    
		    TVector3 partMom(mX,mY,mZ);
		    
		    // Pairwise Distance Between Constituents
		    if(j<(recoCstsEnd[minIndex]-1))
		      {
			for(unsigned int k=j+1; k<recoCstsEnd[minIndex]; k++)
			  {
			    double mXB = recoPartMomX[recoCstIndex[k]];
			    double mYB = recoPartMomY[recoCstIndex[k]];
			    double mZB = recoPartMomZ[recoCstIndex[k]];
			    
			    TVector3 partMomB(mXB,mYB,mZB);
			    
			    double dEta = partMom.PseudoRapidity() - partMomB.PseudoRapidity();
			    double dPhi = TVector2::Phi_mpi_pi(partMom.Phi() - partMomB.Phi());
			    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
			    
			    if(dR < 0.02) noDuplicate = false;
			  }
		      }
		  }

		if(noDuplicate)
		  {
		    recoVsGenChargedJetENoDupHist->Fill(genNRG[i],recoNRG[minIndex]);

		    if(jetMom.PseudoRapidity() > -2.5 && jetMom.PseudoRapidity() < -1.0)
		      jetResVsENegEtaNoDupHist->Fill(genNRG[i],jetERes);
		    if(jetMom.PseudoRapidity() > -1.0 && jetMom.PseudoRapidity() < 1.0)
		      jetResVsEMidEtaNoDupHist->Fill(genNRG[i],jetERes);
		    if(jetMom.PseudoRapidity() > 1.0 && jetMom.PseudoRapidity() < 2.5)
		      jetResVsEPosEtaNoDupHist->Fill(genNRG[i],jetERes);
		  }
	      }
	  }
      }

    NEVENTS++;
  }

  
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
  if(PRINT) c1->Print((results_path+"/numberRecoJets.png").c_str()); // Number of reconstructed jets per event with energy > 5 GeV and Abs(eta) < 2.5
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
  if(PRINT) c2->Print((results_path+"/recoJetEnergy.png").c_str()); // Energy spectrum of reconstructed jets with Abs(eta) < 2.5

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
  if(PRINT) c3->Print((results_path+"/recoJetEta.png").c_str()); // Eta spectrum of reconstructed jets with energy > 5 GeV
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
    if(PRINT) c3_1->Print((results_path+"/recoJetArea.png").c_str()); // Area spectrum of reconstructed jets with energy > 5 GeV
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
  if(PRINT) c4->Print((results_path+"/recoJetEnergyvsEta.png").c_str()); // Energy vs eta of reconstructed jets

  // Reco E Vs Area
  if (useNewEDM) {
    TCanvas *c4_1 = new TCanvas("c4_1","Reco Jet E Vs Area",800,600);
    c4_1->Clear();
    c4_1->Divide(1,1);

    c4_1->cd(1);
    recoChargedJetEvsAreaHist->Draw("COLZ");
    recoChargedJetEvsAreaHist->SetTitle("Reconstructed Jet Energy Vs Area;Area;Energy [GeV]");
    gPad->SetLogz();
    if(PRINT) c4_1->Print((results_path+"/recoJetEnergyvsArea.png").c_str()); // Energy vs area of reconstructed jets
  }

  // Reco Phi Vs Eta
  TCanvas *c5 = new TCanvas("c5","Reco Jet Phi Vs Eta",800,600);
  c5->Clear();
  c5->Divide(1,1);

  c5->cd(1);
  recoChargedJetPhiVsEtaECutHist->Draw("COLZ");
  recoChargedJetPhiVsEtaECutHist->SetTitle("Reconstructed Jet Phi Vs Eta (E > 5);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c5->Print((results_path+"/recoJetPhivsEta.png").c_str()); // Phi vs eta of reconstructed jets

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
  if(PRINT) c6->Print((results_path+"/numConstituentsPerRecoJet.png").c_str()); // Number of constituents in reconstructed jets

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
  if(PRINT) c7->Print((results_path+"/recoJetConstituentMomentum.png").c_str()); // Momentum of reconstructed jet constituents

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
  if(PRINT) c8->Print((results_path+"/recoJetConstituentEta.png").c_str()); // Eta of reconstructed jet constituents

  // Reco Part P Vs Eta
  TCanvas *c9 = new TCanvas("c9","Reco Jet Constituent Momentum Vs Eta",800,600);
  c9->Clear();
  c9->Divide(1,1);

  c9->cd(1);
  recoChargedJetPartPvsEtaHist->Draw("COLZ");
  recoChargedJetPartPvsEtaHist->SetTitle("Reconstructed Jet Constituent Momentum Vs Eta;Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c9->Print((results_path+"/recoJetConstituentMomentumVsEta.png").c_str()); // Momentum vs eta of reconstructed jet constituents

  // Reco Part Phi Vs Eta
  TCanvas *c10 = new TCanvas("c10","Reco Jet Constituent Phi Vs Eta",800,600);
  c10->Clear();
  c10->Divide(1,1);

  c10->cd(1);
  recoChargedJetPartPhiVsEtaHist->Draw("COLZ");
  recoChargedJetPartPhiVsEtaHist->SetTitle("Reconstructed Jet Constituent Phi Vs Eta;Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c10->Print((results_path+"/recoJetConstituentPhiVsEta.png").c_str()); // Phi vs eta of reconstructed jet constituents

  // Reco Constituent Pairwise delta R
  TCanvas *c11 = new TCanvas("c11","Reco Jet Constituent Pairwise Delta R",800,600);
  c11->Clear();
  c11->Divide(1,1);

  c11->cd(1);
  recoChargedJetPartPairwiseDeltaRHist->Draw("COLZ");
  recoChargedJetPartPairwiseDeltaRHist->SetTitle("Pairwise Constituent Delta R;Delta R");
  recoChargedJetPartPairwiseDeltaRHist->GetXaxis()->SetRangeUser(0,0.5);
  gPad->SetLogy();
  if(PRINT) c11->Print((results_path+"/recoJetConstituentPairwiseDR.png").c_str()); // Distance between each pair of constituents in reconstructed jets

  // Reco E Vs Eta No Electron Jets
  TCanvas *c12 = new TCanvas("c12","Reco Jet E Vs Eta (No Electrons)",800,600);
  c12->Clear();
  c12->Divide(1,1);

  c12->cd(1);
  recoChargedJetEvsEtaNoElecHist->Draw("COLZ");
  recoChargedJetEvsEtaNoElecHist->SetTitle("Reconstructed Jet Energy Vs Eta (No Electrons);Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c12->Print((results_path+"/recoJetEnergyVsEtaNoElectron.png").c_str()); // Reconstructed jet energy - no jets containing electrons included

  // Reco Phi Vs Eta No Electron Jets
  TCanvas *c13 = new TCanvas("c13","Reco Jet Phi Vs Eta (No Electrons)",800,600);
  c13->Clear();
  c13->Divide(1,1);

  c13->cd(1);
  recoChargedJetPhiVsEtaECutNoElecHist->Draw("COLZ");
  recoChargedJetPhiVsEtaECutNoElecHist->SetTitle("Reconstructed Jet Phi Vs Eta (E > 5) (No Electrons);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c13->Print((results_path+"/recoJetPhiVsEtaNoElectron.png").c_str()); // Reconstructed Jet phi vs eta - no jets containing electrons included

  // Reco Part P Vs Eta No Electron Jets
  TCanvas *c14 = new TCanvas("c14","Reco Jet Constituent Momentum Vs Eta (No Electrons)",800,600);
  c14->Clear();
  c14->Divide(1,1);

  c14->cd(1);
  recoChargedJetPartPvsEtaNoElecHist->Draw("COLZ");
  recoChargedJetPartPvsEtaNoElecHist->SetTitle("Reconstructed Jet Constituent Momentum Vs Eta (No Electrons);Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c14->Print((results_path+"/recoJetConstituentMomentumVsEtaNoElectron.png").c_str()); // Reconstructed jet constituent momentum vs eta - no jets containing electrons included

  // Reco Part Phi Vs Eta No Electron Jets
  TCanvas *c15 = new TCanvas("c15","Reco Jet Constituent Phi Vs Eta (No Electrons)",800,600);
  c15->Clear();
  c15->Divide(1,1);

  c15->cd(1);
  recoChargedJetPartPhiVsEtaNoElecHist->Draw("COLZ");
  recoChargedJetPartPhiVsEtaNoElecHist->SetTitle("Reconstructed Jet Constituent Phi Vs Eta (No Electrons);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c15->Print((results_path+"/recoJetConstituentPhiVsEtaNoElectron.png").c_str()); // Reconstructed jet constituent phi vs eta - no jets containing electrons included

  
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
  if(PRINT) c16->Print((results_path+"/numberGenJets.png").c_str()); // Number of generator jets per event with energy > 5 GeV and Abs(eta) < 2.5

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
  if(PRINT) c17->Print((results_path+"/genJetEnergy.png").c_str()); // Energy spectrum of generated jets with Abs(eta) < 2.5

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
  if(PRINT) c18->Print((results_path+"/genJetEta.png").c_str()); // Eta spectrum of generator jets with energy > 5 GeV

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
    if(PRINT) c18_1->Print((results_path+"/genJetArea.png").c_str()); // Area spectrum of generator jets with energy > 5 GeV
  }

  // Gen E Vs Eta
  TCanvas *c19 = new TCanvas("c19","Gen Jet E Vs Eta",800,600);
  c19->Clear();
  c19->Divide(1,1);

  c19->cd(1);
  genChargedJetEvsEtaHist->Draw("COLZ");
  genChargedJetEvsEtaHist->SetTitle("Generator Jet Energy Vs Eta;Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c19->Print((results_path+"/genJetEnergyvsEta.png").c_str()); // Energy vs eta of generator jets

  // Gen E Vs Area
  if(useNewEDM) {
    TCanvas *c19_1 = new TCanvas("c19_1","Gen Jet E Vs Area",800,600);
    c19_1->Clear();
    c19_1->Divide(1,1);

    c19_1->cd(1);
    genChargedJetEvsAreaHist->Draw("COLZ");
    genChargedJetEvsAreaHist->SetTitle("Generator Jet Energy Vs Area;Area;Energy [GeV]");
    gPad->SetLogz();
    if(PRINT) c19_1->Print((results_path+"/genJetEnergyvsArea.png").c_str()); // Energy vs area of generator jets
  }

  // Gen Phi Vs Eta
  TCanvas *c20 = new TCanvas("c20","Gen Jet Phi Vs Eta",800,600);
  c20->Clear();
  c20->Divide(1,1);

  c20->cd(1);
  genChargedJetPhiVsEtaECutHist->Draw("COLZ");
  genChargedJetPhiVsEtaECutHist->SetTitle("Generator Jet Phi Vs Eta (E > 5);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c20->Print((results_path+"/genJetPhiVsEta.png").c_str()); // Phi vs eta of generator jets

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
  if(PRINT) c21->Print((results_path+"/numConstituentsPerGenJet.png").c_str()); // Number of constituents in generator jets

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
  if(PRINT) c22->Print((results_path+"/genJetConstituentMomentum.png").c_str()); // Momentum of generator jet constituents

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
  if(PRINT) c23->Print((results_path+"/genJetConstituentEta.png").c_str()); // Eta of generator jet constituents

  // Gen Part P Vs Eta
  TCanvas *c24 = new TCanvas("c24","Gen Jet Constituent Momentum Vs Eta",800,600);
  c24->Clear();
  c24->Divide(1,1);

  c24->cd(1);
  genChargedJetPartPvsEtaHist->Draw("COLZ");
  genChargedJetPartPvsEtaHist->SetTitle("Generator Jet Constituent Momentum Vs Eta;Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c24->Print((results_path+"/genJetConstituentMomentumVsEta.png").c_str()); // Momentum vs eta of generator jet constituents

  // Gen Part Phi Vs Eta
  TCanvas *c25 = new TCanvas("c25","Gen Jet Constituent Phi Vs Eta",800,600);
  c25->Clear();
  c25->Divide(1,1);

  c25->cd(1);
  genChargedJetPartPhiVsEtaHist->Draw("COLZ");
  genChargedJetPartPhiVsEtaHist->SetTitle("Generator Jet Constituent Phi Vs Eta;Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c25->Print((results_path+"/genJetConstituentPhiVsEta.png").c_str()); // Phi vs eta of generator jet constituents

  // Gen Constituent Pairwise delta R
  TCanvas *c26 = new TCanvas("c26","Gen Jet Constituent Pairwise Delta R",800,600);
  c26->Clear();
  c26->Divide(1,1);

  c26->cd(1);
  genChargedJetPartPairwiseDeltaRHist->Draw("COLZ");
  genChargedJetPartPairwiseDeltaRHist->SetTitle("Generator Jet Pairwise Constituent Delta R;Delta R");
  genChargedJetPartPairwiseDeltaRHist->GetXaxis()->SetRangeUser(0,0.5);
  gPad->SetLogy();
  if(PRINT) c26->Print((results_path+"/genJetConstituentPairwiseDR.png").c_str()); // Distance between each pair of constituents in generator jets

  // Gen E Vs Eta No Electron Jets
  TCanvas *c27 = new TCanvas("c27","Gen Jet E Vs Eta (No Electrons)",800,600);
  c27->Clear();
  c27->Divide(1,1);

  c27->cd(1);
  genChargedJetEvsEtaNoElecHist->Draw("COLZ");
  genChargedJetEvsEtaNoElecHist->SetTitle("Generator Jet Energy Vs Eta (No Electrons);Eta;Energy [GeV]");
  gPad->SetLogz();
  if(PRINT) c27->Print((results_path+"/genJetEnergyVsEtaNoElectron.png").c_str()); // Generator jet energy vs eta - no jets containing electrons included

  // Gen Phi Vs Eta No Electron Jets
  TCanvas *c28 = new TCanvas("c28","Gen Jet Phi Vs Eta (No Electrons)",800,600);
  c28->Clear();
  c28->Divide(1,1);

  c28->cd(1);
  genChargedJetPhiVsEtaECutNoElecHist->Draw("COLZ");
  genChargedJetPhiVsEtaECutNoElecHist->SetTitle("Generator Jet Phi Vs Eta (E > 5) (No Electrons);Eta;Phi");
  gPad->SetLogz();
  if(PRINT) c28->Print((results_path+"/genJetPhiVsEtaNoElectron.png").c_str()); // Generator Jet phi vs eta - no jets containing electrons included

  // Gen Part P Vs Eta No Electron Jets
  TCanvas *c29 = new TCanvas("c29","Gen Jet Constituent Momentum Vs Eta (No Electrons)",800,600);
  c29->Clear();
  c29->Divide(1,1);

  c29->cd(1);
  genChargedJetPartPvsEtaNoElecHist->Draw("COLZ");
  genChargedJetPartPvsEtaNoElecHist->SetTitle("Generator Jet Constituent Momentum Vs Eta (No Electrons);Eta;Momentum [GeV/c]");
  gPad->SetLogz();
  if(PRINT) c29->Print((results_path+"/genJetConstituentMomentumVsEtaNoElectron.png").c_str()); // Generator jet constituent momentum vs eta - no jets containing electrons included

  // Gen Part Phi Vs Eta No Electron Jets
  TCanvas *c30 = new TCanvas("c30","Gen Jet Constituent Phi Vs Eta (No Electrons)",800,600);
  c30->Clear();
  c30->Divide(1,1);

  c30->cd(1);
  genChargedJetPartPhiVsEtaNoElecHist->Draw("COLZ");
  genChargedJetPartPhiVsEtaNoElecHist->SetTitle("Generator Jet Constituent Phi Vs Eta (No Electrons);Eta;Phi");
  gPad->SetLogz();
  //c30->Print((results_path+"/recoJetEvsEta.png").c_str());
  if(PRINT) c30->Print((results_path+"/genJetConstituentPhiVsEtaNoElectron.png").c_str()); // Generator jet constituent phi vs eta - no jets containing electrons included

  
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
  if(PRINT) c31->Print((results_path+"/genRecoJetDeltaR.png").c_str()); // Distance between closest generated and reconstructed jet pair

  // Matched Reco Vs Gen Eta
  TCanvas *c32 = new TCanvas("c32","Reco Vs Gen Eta",800,600);
  c32->Clear();
  c32->Divide(1,1);

  c32->cd(1);
  recoVsGenChargedJetEtaHist->Draw("COLZ");
  recoVsGenChargedJetEtaHist->SetTitle("Reconstructed Vs Generator Jet Eta;Gen Eta;Reco Eta");
  gPad->SetLogz();
  if(PRINT) c32->Print((results_path+"/matchedRecoVsGenJetEta.png").c_str()); // Matched Reconstructed Vs Generator Jet eta

  // Matched Reco Vs Gen Phi
  TCanvas *c33 = new TCanvas("c33","Reco Vs Gen Phi",800,600);
  c33->Clear();
  c33->Divide(1,1);

  c33->cd(1);
  recoVsGenChargedJetPhiHist->Draw("COLZ");
  recoVsGenChargedJetPhiHist->SetTitle("Reconstructed Vs Generator Jet Phi;Gen Phi;Reco Phi");
  gPad->SetLogz();
  if(PRINT) c33->Print((results_path+"/matchedRecoVsGenJetPhi.png").c_str()); // Matched reconstructed vs generator jet phi

  // Matched Reco Vs Gen Area
  if(useNewEDM) {
    TCanvas *c33_1 = new TCanvas("c33_1","Reco Vs Gen Area",800,600);
    c33_1->Clear();
    c33_1->Divide(1,1);

    c33_1->cd(1);
    recoVsGenChargedJetAreaHist->Draw("COLZ");
    recoVsGenChargedJetAreaHist->SetTitle("Reconstructed Vs Generator Jet Area;Gen Area;Reco Area");
    gPad->SetLogz();
    if(PRINT) c33_1->Print((results_path+"/matchedRecoVsGenJetArea.png").c_str()); // Matched reconstructed vs generator jet area
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
  if(PRINT) c34->Print((results_path+"/matchedRecoVsGenJetEnergy.png").c_str()); // Matched reconstructed vs generator jet energy

  // Jet Res Vs Gen Eta
  TCanvas *c35 = new TCanvas("c35","Jet Res Vs Gen Eta",800,600);
  c35->Clear();
  c35->Divide(1,1);

  c35->cd(1);
  jetResVsEtaHist->Draw("COLZ");
  jetResVsEtaHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Eta;Gen Eta;Res");
  gPad->SetLogz();
  if(PRINT) c35->Print((results_path+"/matchedJetResolutionVsEta.png").c_str()); // Matched jet resolution vs generator jet eta

  // Jet Res Vs Gen E
  TCanvas *c36 = new TCanvas("c36","Jet Res Vs Gen E",800,600);
  c36->Clear();
  c36->Divide(1,1);

  c36->cd(1);
  jetResVsEHist->Draw("COLZ");
  jetResVsEHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c36->Print((results_path+"/matchedJetResolutionVsEnergy.png").c_str()); // Matched jet resolution vs generator jet energy

  // Jet Res Vs Gen E Neg Eta
  TCanvas *c37 = new TCanvas("c37","Jet Res Vs Gen E (-2.5 < eta < -1.0)",800,600);
  c37->Clear();
  c37->Divide(1,1);

  c37->cd(1);
  jetResVsENegEtaNoDupHist->Draw("COLZ");
  jetResVsENegEtaNoDupHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy (-2.5 < eta < -1.0) No Duplicate;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c37->Print((results_path+"/matchedJetResolutionVsEnergyNegEta.png").c_str()); // Matched jet resolution vs generator jet energy -2.5 < eta < -1.0

  // Jet Res Vs Gen E Mid Eta
  TCanvas *c38 = new TCanvas("c38","Jet Res Vs Gen E (-1.0 < eta < 1.0)",800,600);
  c38->Clear();
  c38->Divide(1,1);

  c38->cd(1);
  jetResVsEMidEtaNoDupHist->Draw("COLZ");
  jetResVsEMidEtaNoDupHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy (-1.0 < eta < 1.0) No Duplicate;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c38->Print((results_path+"/matchedJetResolutionVsEnergyMidEta.png").c_str()); // Matched jet resolution vs generator jet energy -1.0 < eta < 1.0
    delete c38;
  // Jet Res Vs Gen E Pos Eta
  TCanvas *c39 = new TCanvas("c39","Jet Res Vs Gen E (1.0 < eta < 2.5)",800,600);
  c39->Clear();
  c39->Divide(1,1);

  c39->cd(1);
  jetResVsEPosEtaNoDupHist->Draw("COLZ");
  jetResVsEPosEtaNoDupHist->SetTitle("(Reco - Gen)/Gen Jet Energy Vs Gen Energy (1.0 < eta < 2.5) No Duplicate;Gen E;Res");
  gPad->SetLogz();
  if(PRINT) c39->Print((results_path+"/matchedJetResolutionVsEnergyPosEta.png").c_str()); // Matched jet resolution vs generator jet energy 1.0 < eta < 2.5
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
  if(PRINT) c40->Print((results_path+"/matchedJetResolutionVsEnergyNegEtaFitSummary.png").c_str()); // Matched jet resolution vs generator jet energy -2.5 < eta < -1.0 fits
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
  if(PRINT) c41->Print((results_path+"/matchedJetResolutionVsEnergyMidEtaFitSummary.png").c_str()); // Matched jet resolution vs generator jet energy -1.0 < eta < 1.0 fits
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
  if(PRINT) c42->Print((results_path+"/matchedJetResolutionVsEnergyPosEtaFitSummary.png").c_str()); // Matched jet resolution vs generator jet energy 1.0 < eta < 2.5 fits
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

  if(PRINT) c43->Print((results_path+"/matchedJetScaleResolutionSummary.png").c_str()); // Matched jet JER/JES summary
    delete c43;

delete mychain;

  return 0;
}
