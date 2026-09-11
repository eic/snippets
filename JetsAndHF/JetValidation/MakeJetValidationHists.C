// ============================================================================
//! \file    MakeJetValidationHists.C
//! \authors Brian Page (bpage@bnl.gov),
//!          adapted by Derek Anderson (derek.murphy.anderson@protonmail.com)
// ----------------------------------------------------------------------------
//! \brief Adaption of the Jet Benchmark to run standalone for validation.
//!   This macro generates histograms from eicrecon output.
//!
//! \usage In eic-shell:
//!     root -b -q MakeJetValidationHists.C'(<output path>, \
//!                                          <output suffix>, \
//!                                          <input file list>, \
//!                                          <n files to read>, \
//!                                          <n events to process>)'
// ============================================================================

#include <edm4eic/EDM4eicVersion.h>
#include <TChain.h>
#include <TFile.h>
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

///! Default output file path.
const std::string DefaultOutPath = ".";

///! Default output file suffix.
const std::string DefaultOutSuffix = "files26071.py8ncdis10x100q100t1000";

///! Default input file list.
const std::string DefaultInFileList = "filelists/files26071.py8ncdis10x100q100t1000.list";

///! Default no. of files to process.
///! -1 means process all files.
const std::int32_t DefaultNFiles = -1;

///! Default total no. of events to process.
///! -1 means process all events.
const std::int32_t DefaultNEvents = -1;

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
 *! generated jet histograms and save them to a ROOT file.
 *!
 *! \param[out] results_path   Location to save output file to
 *! \param[out] results_suffix Suffix to append to output file
 *! \param[in]  filelist       Input filelist to use
 *! \param[in]  n_files        Number of files to read
 *! \param[in]  n_events       Number of total events to process
 */
int MakeJetValidationHists(
  const std::string& results_path = DefaultOutPath,
  const std::string& results_suffix = DefaultOutSuffix,
  const std::string& filelist = DefaultInFileList,
  const std::int32_t n_files = DefaultNFiles,
  const std::int32_t n_events = DefaultNEvents
) {

  const bool PRINT = true;

  std::ifstream mylist(filelist);
  if (!mylist.is_open()) {
    std::cerr << "PANIC: Couldn't open fileslist '" << filelist << "'!" << std::endl;
    return 1;
  }

  // Load input files
  std::string file;
  std::int32_t i_file = 0;
  std::vector<std::string> rec_files;
  while (std::getline(mylist, file)) {
    rec_files.push_back(file);
    ++i_file;
    if ((n_files > -1) && (i_file == n_files)) {
      break;
    }
  }

  // Input
  TChain *mychain = new TChain("events");
  for (const auto& rec_file : rec_files) {
    mychain->Add(rec_file.c_str());
  }

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
  if (PRINT) {
    std::cout << "INFO: Using new EDM? " << useNewEDM << std::endl;
  }

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

  // open output file
  TFile* out_file = new TFile((results_path+"/hists."+results_suffix+".root").c_str(), "recreate");
  std::cout << "INFO: writing histograms to " << out_file->GetName() << std::endl;

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

    if ((n_events > -1) && (NEVENTS == n_events)) {
      break;
    }
    if (NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

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

  // Save and close histograms
  out_file->cd();
  out_file->Write();
  out_file->Close();

  delete mychain;
  return 0;
}
