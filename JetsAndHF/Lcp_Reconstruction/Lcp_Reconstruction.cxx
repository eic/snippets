// Λc⁺ → p K⁻ π⁺ reconstruction with (6.28 ± 0.32)%
// root [1] TDatabasePDG *pdgDB = TDatabasePDG::Instance();
//root [2] TParticlePDG *particle = pdgDB->GetParticle(4122);
// TParticlePDG *particle = pdgDB->GetParticle("Lambda_c+");
//root [3] particle->Print()  select all and (ctrl+shift+P) reindent all the code 
// Mass:  2.28646 ± 0.00014 GeV/c²
// Shyam Kumar; INFN Bari, Italy
// shyam.kumar@ba.infn.it; shyam055119@gmail.com 

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"
#endif
#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

using namespace std;

const double gPionMass = 0.13957;
const double gKaonMass = 0.493677;
const double gProtonMass = 0.938272;

const double twoPi = 2.*3.1415927;
const double eMass = 0.000511;
const int pdg_Lcp = 4122, pdg_p = 2212, pdg_k = 321, pdg_pi = 211;
bool debug = false;


const double bField = -1.7; // Tesla

TVector3 GetDCAToPrimaryVertex(const int index, TVector3 vtx);

TLorentzVector GetDCADaughters(const int index1, const int index2, const int index3, TVector3 vtx,
  float *dcaDaughters, float &cosTheta, float &cosTheta_xy, float &decayLength, float &V0DcaToVtx, float &sigma_vtx);

TTreeReaderArray<float> *rcMomPx2;
TTreeReaderArray<float> *rcMomPy2;
TTreeReaderArray<float> *rcMomPz2;
TTreeReaderArray<float> *rcCharge2;

TTreeReaderArray<float> *rcTrkLoca2;
TTreeReaderArray<float> *rcTrkLocb2;
TTreeReaderArray<float> *rcTrkTheta2;
TTreeReaderArray<float> *rcTrkPhi2;

int main(int argc, char **argv)
{
  if(argc!=3 && argc!=1) return 0;

  TString listname;
  TString outname;

  if(argc==1)
  {
    listname  = "test.list";
    outname = "test.root";
  }

  if(argc==3)
  {
    listname = argv[1];
    outname = argv[2];
  }  

  TChain *chain = new TChain("events");

  int nfiles = 0;
  char filename[512];
  ifstream *inputstream = new ifstream;
  inputstream->open(listname.Data());
  if(!inputstream)
  {
    printf("[e] Cannot open file list: %s\n", listname.Data());
  }
  while(inputstream->good())
  {
    inputstream->getline(filename, 512);
    if(inputstream->good())
    {
     TFile *ftmp = TFile::Open(filename, "read");
     if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) 
     {
       printf("[e] Could you open file: %s\n", filename);
     } 
     else
     {
       cout<<"[i] Add "<<nfiles<<"th file: "<<filename<<endl;
       chain->Add(filename);
       nfiles++;
     }
   }
 }
 inputstream->close();
 printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());

 TH1F *hEventStat = new TH1F("hEventStat", "Event statistics", 7, 0, 7);
 hEventStat->GetXaxis()->SetBinLabel(1, "MC events");
 hEventStat->GetXaxis()->SetBinLabel(2, "#Lambda_{c}^{+}");
 hEventStat->GetXaxis()->SetBinLabel(3, "#Lambda_{c}^{+} -> p + K^{-} + #pi^{+}");
 hEventStat->GetXaxis()->SetBinLabel(4, "Reco Signal #Lambda_{c}^{+}/#Lambda_{c}^{-}");
 hEventStat->GetXaxis()->SetBinLabel(5, "Reco Signal #Lambda_{c}^{+}");
 hEventStat->GetXaxis()->SetBinLabel(6, "Reco Signal #Lambda_{c}^{-}");
 hEventStat->GetXaxis()->SetBinLabel(7, "Reco Bkg #Lambda_{c}^{+}");


 TH1F *hMcMult = new TH1F("hMcMult", "MC multiplicity (|#eta| < 3.5);N_{MC}", 50, 0, 50);

 TH1F *hMcVtxX = new TH1F("hMcVtxX", "x position of MC vertex;x (mm)", 500, -5.0, 5.0);
 TH1F *hMcVtxY = new TH1F("hMcVtxY", "y position of MC vertex;y (mm)", 500, -5.0, 5.0);
 TH1F *hMcVtxZ = new TH1F("hMcVtxZ", "z position of MC vertex;z (mm)", 800, -200, 200);

 TH1F *hRecVtxX = new TH1F("hRecVtxX", "x position of Rec vertex;x (mm)", 500, -5.0, 5.0);
 TH1F *hRecVtxY = new TH1F("hRecVtxY", "y position of Rec vertex;y (mm)", 500, -5.0, 5.0);
 TH1F *hRecVtxZ = new TH1F("hRecVtxZ", "z position of Rec vertex;z (mm)", 800, -200, 200);

 TH2F *hLcpDecayVxVy = new TH2F("hLcpDecayVxVy", "#Lambda_{c}^{+} decay vertex to primary vertex;#Deltav_{x} (mm);#Deltav_{y} (mm)", 400, -1-0.0025, 1-0.0025, 400, -1-0.0025, 1-0.0025);
 TH2F *hLcpDecayVrVz = new TH2F("hLcpDecayVrVz", "#Lambda_{c}^{+} decay vertex to primary vertex;#Deltav_{z} (mm);#Deltav_{r} (mm)", 100, -2, 2, 100, -0.2, 1.8);

 TH2F *hMCLcpPtRap = new TH2F("hMCLcpPtRap", "MC #Lambda_{c}^{+};y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);

 TH2F *hMcPPtEta = new TH2F("hMcPPtEta", "MC P from #Lambda_{c}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
 TH2F *hMcPPtEtaReco = new TH2F("hMcPPtEtaReco", "RC P from #Lambda_{c}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

 TH2F *hMcKPtEta = new TH2F("hMcKPtEta", "MC K from #Lambda_{c}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
 TH2F *hMcKPtEtaReco = new TH2F("hMcKPtEtaReco", "RC K from #Lambda_{c}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

 TH2F *hMcPiPtEta = new TH2F("hMcPiPtEta", "MC #pi from #Lambda_{c}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
 TH2F *hMcPiPtEtaReco = new TH2F("hMcPiPtEtaReco", "RC #pi from #Lambda_{c}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

 TH1F *hNRecoVtx = new TH1F("hNRecoVtx", "Number of reconstructed vertices;N", 10, 0, 10);

 const char* part_name[3] = {"P", "K", "Pi"};
 const char* part_title[3] = {"P", "K", "#pi"};
 TH3F *hRcSecPartLocaToRCVtx[3];
 TH3F *hRcSecPartLocbToRCVtx[3];
 TH3F *hRcPrimPartLocaToRCVtx[3];
 TH3F *hRcPrimPartLocbToRCVtx[3];
 for(int i=0; i<3; i++)
 {
  hRcSecPartLocaToRCVtx[i] = new TH3F(Form("hRcSec%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for #Lambda_{c}^{+} decayed %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
  hRcSecPartLocbToRCVtx[i] = new TH3F(Form("hRcSec%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for #Lambda_{c}^{+} decayed %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
  hRcPrimPartLocaToRCVtx[i] = new TH3F(Form("hRcPrim%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
  hRcPrimPartLocbToRCVtx[i] = new TH3F(Form("hRcPrim%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
}

const char* axis_name[3] = {"x", "y", "z"};
const int nDimDca = 4;
const int nBinsDca[nDimDca] = {50, 20, 500, 50};
const double minBinDca[nDimDca] = {0, -5, -1+0.002, 0};
const double maxBinDca[nDimDca] = {5, 5, 1+0.002, 50};
THnSparseF *hPrimTrkDcaToRCVtx[3][3];
for(int i=0; i<3; i++)
{
  for(int j=0; j<3; j++)
  {
   hPrimTrkDcaToRCVtx[i][j] = new THnSparseF(Form("hPrim%sDca%sToRCVtx",part_name[i],axis_name[j]), Form("DCA_{%s} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{%s} (mm);N_{MC}",axis_name[j],part_title[i],axis_name[j]), nDimDca, nBinsDca, minBinDca, maxBinDca);
 }
}

TH3F *h3PairDca12[2], *h3PairDca23[2], *h3PairDca13[2];
TH3F *h3PairCosTheta[2];
TH3F *h3PairDca[2];
TH3F *h3PairDecayLength[2];
const char* pair_name[2] = {"signal", "bkg"};
const char* pair_title[2] = {"Signal", "Background"};
for(int i=0; i<2; i++)
{
  h3PairDca12[i] = new TH3F(Form("h3PairDca12_%s", pair_name[i]), Form("%s pair DCA_{12};p_{T} (GeV/c);#eta;DCA_{12} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
  h3PairDca23[i] = new TH3F(Form("h3PairDca23_%s", pair_name[i]), Form("%s pair DCA_{23};p_{T} (GeV/c);#eta;DCA_{23} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
  h3PairDca13[i] = new TH3F(Form("h3PairDca13_%s", pair_name[i]), Form("%s pair DCA_{13};p_{T} (GeV/c);#eta;DCA_{13} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
  h3PairCosTheta[i] = new TH3F(Form("h3PairCosTheta_%s", pair_name[i]), Form("%s pair cos(#theta);p_{T} (GeV/c);#eta;cos(#theta)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, -1, 1);

  h3PairDca[i] = new TH3F(Form("h3PairDca_%s", pair_name[i]), Form("%s pair DCA;p_{T} (GeV/c);#eta;DCA_{pair} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);

  h3PairDecayLength[i] = new TH3F(Form("h3PairDecayLength_%s", pair_name[i]), Form("%s pair decay length;p_{T} (GeV/c);#eta;L (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
}

  // Invariant mass
const char* cut_name[2] = {"all", "DCA"};
TH3F *h3InvMass[2][2];
for(int i=0; i<2; i++)
{
  for(int j=0; j<2; j++)
  {
   h3InvMass[i][j] = new TH3F(Form("h3InvMass_%s_%s", pair_name[i], cut_name[j]), "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 2.0, 3.0);
 }
}

TTreeReader treereader(chain);
  // MC  
TTreeReaderArray<int> mcPartGenStatus = {treereader, "MCParticles.generatorStatus"};
TTreeReaderArray<int> mcPartPdg = {treereader, "MCParticles.PDG"};
TTreeReaderArray<float> mcPartCharge = {treereader, "MCParticles.charge"};
TTreeReaderArray<unsigned int> mcPartParent_begin = {treereader, "MCParticles.parents_begin"};
TTreeReaderArray<unsigned int> mcPartParent_end = {treereader, "MCParticles.parents_end"};
TTreeReaderArray<int> mcPartParent_index = {treereader, "_MCParticles_parents.index"};
TTreeReaderArray<unsigned int> mcPartDaughter_begin = {treereader, "MCParticles.daughters_begin"};
TTreeReaderArray<unsigned int> mcPartDaughter_end = {treereader, "MCParticles.daughters_end"};
TTreeReaderArray<int> mcPartDaughter_index = {treereader, "_MCParticles_daughters.index"};
TTreeReaderArray<double> mcPartMass = {treereader, "MCParticles.mass"};
TTreeReaderArray<double> mcPartVx = {treereader, "MCParticles.vertex.x"};
TTreeReaderArray<double> mcPartVy = {treereader, "MCParticles.vertex.y"};
TTreeReaderArray<double> mcPartVz = {treereader, "MCParticles.vertex.z"};
TTreeReaderArray<float> mcMomPx = {treereader, "MCParticles.momentum.x"};
TTreeReaderArray<float> mcMomPy = {treereader, "MCParticles.momentum.y"};
TTreeReaderArray<float> mcMomPz = {treereader, "MCParticles.momentum.z"};
TTreeReaderArray<double> mcEndPointX = {treereader, "MCParticles.endpoint.x"};
TTreeReaderArray<double> mcEndPointY = {treereader, "MCParticles.endpoint.y"};
TTreeReaderArray<double> mcEndPointZ = {treereader, "MCParticles.endpoint.z"};

TTreeReaderArray<unsigned int> assocChSimID = {treereader, "ReconstructedChargedParticleAssociations.simID"};
TTreeReaderArray<unsigned int> assocChRecID = {treereader, "ReconstructedChargedParticleAssociations.recID"};
TTreeReaderArray<float> assocWeight = {treereader, "ReconstructedChargedParticleAssociations.weight"};

TTreeReaderArray<float> rcMomPx = {treereader, "ReconstructedChargedParticles.momentum.x"};
TTreeReaderArray<float> rcMomPy = {treereader, "ReconstructedChargedParticles.momentum.y"};
TTreeReaderArray<float> rcMomPz = {treereader, "ReconstructedChargedParticles.momentum.z"};
TTreeReaderArray<float> rcPosx = {treereader, "ReconstructedChargedParticles.referencePoint.x"};
TTreeReaderArray<float> rcPosy = {treereader, "ReconstructedChargedParticles.referencePoint.y"};
TTreeReaderArray<float> rcPosz = {treereader, "ReconstructedChargedParticles.referencePoint.z"};
TTreeReaderArray<float> rcCharge = {treereader, "ReconstructedChargedParticles.charge"};
TTreeReaderArray<int>   rcPdg = {treereader, "ReconstructedChargedParticles.PDG"};

TTreeReaderArray<float> rcTrkLoca = {treereader, "CentralCKFTrackParameters.loc.a"};
TTreeReaderArray<float> rcTrkLocb = {treereader, "CentralCKFTrackParameters.loc.b"};
TTreeReaderArray<float> rcTrkqOverP = {treereader, "CentralCKFTrackParameters.qOverP"};
TTreeReaderArray<float> rcTrkTheta = {treereader, "CentralCKFTrackParameters.theta"};
TTreeReaderArray<float> rcTrkPhi = {treereader, "CentralCKFTrackParameters.phi"};

rcMomPx2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.x"};
rcMomPy2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.y"};
rcMomPz2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.z"};
rcCharge2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.charge"};

rcTrkLoca2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.loc.a"};
rcTrkLocb2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.loc.b"};
rcTrkTheta2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.theta"};
rcTrkPhi2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.phi"};

TTreeReaderArray<float> CTVx = {treereader, "CentralTrackVertices.position.x"};
TTreeReaderArray<float> CTVy = {treereader, "CentralTrackVertices.position.y"};
TTreeReaderArray<float> CTVz = {treereader, "CentralTrackVertices.position.z"};
TTreeReaderArray<int> CTVndf = {treereader, "CentralTrackVertices.ndf"};
TTreeReaderArray<float> CTVchi2 = {treereader, "CentralTrackVertices.chi2"};
TTreeReaderArray<float> CTVerr_xx = {treereader, "CentralTrackVertices.positionError.xx"};
TTreeReaderArray<float> CTVerr_yy = {treereader, "CentralTrackVertices.positionError.yy"};
TTreeReaderArray<float> CTVerr_zz = {treereader, "CentralTrackVertices.positionError.zz"};

TTreeReaderArray<int> prim_vtx_index = {treereader, "PrimaryVertices_objIdx.index"};

TTreeReaderArray<unsigned int> vtxAssocPart_begin = {treereader, "CentralTrackVertices.associatedParticles_begin"};
TTreeReaderArray<unsigned int> vtxAssocPart_end = {treereader, "CentralTrackVertices.associatedParticles_end"};
TTreeReaderArray<int> vtxAssocPart_index = {treereader, "_CentralTrackVertices_associatedParticles.index"};

  // File to create efficiency of D0 meson
TFile *file_gen = new TFile("SignalLcpGen.root", "RECREATE");
TTree *tree_gen = new TTree("treeMLSigGen", "treeMLSigGen"); 

   // Define variables to store in the Ntuple
float pt_Lcp_gen, y_Lcp_gen;
tree_gen->Branch("pt_Lcp_gen", &pt_Lcp_gen, "pt_Lcp_gen/F");  
tree_gen->Branch("y_Lcp_gen", &y_Lcp_gen, "y_Lcp_gen/F"); 


    // Create a ROOT file to store the Ntuple
TFile *file_signal = new TFile("SignalLcp.root", "RECREATE");
TTree *tree_sig = new TTree("treeMLSig", "treeMLSig"); 

  // Define variables to store in the Ntuple
float d0_p_sig, d0_k_sig, d0_pi_sig, d0xy_p_sig, d0xy_k_sig, d0xy_pi_sig, sum_d0xy_sig, dca_12_sig, dca_Lcp_sig, decay_length_sig;
float costheta_sig, costhetaxy_sig, pt_Lcp_sig, y_Lcp_sig, mass_Lcp_sig, sigma_vtx_sig, mult_sig;

    // Link the variables to the TTree branches
tree_sig->Branch("d0_p", &d0_p_sig, "d0_p/F");
tree_sig->Branch("d0_k", &d0_k_sig, "d0_k/F");
tree_sig->Branch("d0_pi", &d0_pi_sig, "d0_pi/F"); 
tree_sig->Branch("d0xy_p", &d0xy_p_sig, "d0xy_p/F");
tree_sig->Branch("d0xy_k", &d0xy_k_sig, "d0xy_k/F");
tree_sig->Branch("d0xy_pi", &d0xy_pi_sig, "d0xy_pi/F");  
tree_sig->Branch("sum_d0xy", &sum_d0xy_sig, "sum_d0xy/F");        
tree_sig->Branch("dca_12", &dca_12_sig, "dca_12/F");
tree_sig->Branch("dca_Lcp", &dca_Lcp_sig, "dca_Lcp/F");
tree_sig->Branch("pt_Lcp", &pt_Lcp_sig, "pt_Lcp/F");  
tree_sig->Branch("y_Lcp", &y_Lcp_sig, "y_Lcp/F");  
tree_sig->Branch("mass_Lcp", &mass_Lcp_sig, "mass_Lcp/F");              
tree_sig->Branch("decay_length", &decay_length_sig, "decay_length/F");   
tree_sig->Branch("costheta", &costheta_sig, "costheta/F"); 
tree_sig->Branch("costheta_xy", &costhetaxy_sig, "costheta_xy/F"); 
tree_sig->Branch("sigma_vtx", &sigma_vtx_sig, "sigma_vtx/F"); 
tree_sig->Branch("mult", &mult_sig, "mult/F");             

TFile *file_bkg = new TFile("BkgLcp.root", "RECREATE");
TTree *tree_bkg = new TTree("treeMLBkg", "treeMLBkg"); 

  // Define variables to store in the Ntuple
float d0_p_bkg, d0_k_bkg, d0_pi_bkg, d0xy_p_bkg, d0xy_k_bkg, d0xy_pi_bkg, sum_d0xy_bkg, dca_12_bkg, dca_Lcp_bkg, decay_length_bkg;
float costheta_bkg, costhetaxy_bkg, pt_Lcp_bkg, y_Lcp_bkg, mass_Lcp_bkg, sigma_vtx_bkg, mult_bkg;
  // Link the variables to the TTree branches
tree_bkg->Branch("d0_p", &d0_p_bkg, "d0_p/F");
tree_bkg->Branch("d0_k", &d0_k_bkg, "d0_k/F");
tree_bkg->Branch("d0_pi", &d0_pi_bkg, "d0_pi/F");
tree_bkg->Branch("d0xy_p", &d0xy_p_bkg, "d0xy_p/F");
tree_bkg->Branch("d0xy_k", &d0xy_k_bkg, "d0xy_k/F");
tree_bkg->Branch("d0xy_pi", &d0xy_pi_bkg, "d0xy_pi/F");  
tree_bkg->Branch("sum_d0xy", &sum_d0xy_bkg, "sum_d0xy/F");            
tree_bkg->Branch("dca_12", &dca_12_bkg, "dca_12/F");
tree_bkg->Branch("dca_Lcp", &dca_Lcp_bkg, "dca_Lcp/F");
tree_bkg->Branch("pt_Lcp", &pt_Lcp_bkg, "pt_Lcp/F");  
tree_bkg->Branch("y_Lcp", &y_Lcp_bkg, "y_Lcp/F");  
tree_bkg->Branch("mass_Lcp", &mass_Lcp_bkg, "mass_Lcp/F");              
tree_bkg->Branch("decay_length", &decay_length_bkg, "decay_length/F");   
tree_bkg->Branch("costheta", &costheta_bkg, "costheta/F");  
tree_bkg->Branch("costheta_xy", &costhetaxy_bkg, "costheta_xy/F"); 
tree_bkg->Branch("sigma_vtx", &sigma_vtx_bkg, "sigma_vtx/F"); 
tree_bkg->Branch("mult", &mult_bkg, "mult/F");                    


int nevents = 0;
int count = 0;

while(treereader.Next())
{
  if(nevents%1000==0) printf("\nEvent No.-----> %d\n",nevents);
  // find MC primary vertex
  int nMCPart = mcPartMass.GetSize();

  TVector3 vertex_mc(-999., -999., -999.);
  for(int imc=0; imc<nMCPart; imc++)
  {
   // Status 4 describes the incoming electron/proton and the end point of electron or proton is MC Vertex  
   if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11)
   {
     vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
     break;
   }
 }
 hEventStat->Fill(0.5);
 hMcVtxX->Fill(vertex_mc.x());
 hMcVtxY->Fill(vertex_mc.y());
 hMcVtxZ->Fill(vertex_mc.z());


      // get RC primary vertex
 TVector3 vertex_rc(-999., -999., -999.);
 if(prim_vtx_index.GetSize()>0)
 {
   int rc_vtx_index = prim_vtx_index[0];
   vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
 }

 hRecVtxX->Fill(vertex_rc.x());
 hRecVtxY->Fill(vertex_rc.y());
 hRecVtxZ->Fill(vertex_rc.z());

      // map MC and RC particles
 int nAssoc = assocChRecID.GetSize();
      map<int, int> assoc_map_mc_to_rc;  // MC --> RC Associations
      map<int, int> assoc_map_rc_to_mc;  // RC --> MC Associations

      for(unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
      {
	  // loop over the association to find the matched MC particle
	  // with largest weight
       double max_weight = 0;
       int matched_mc_index = -1;
       for(int j=0; j<nAssoc; j++)
       {
         if(assocChRecID[j] != rc_index) continue;
         if(assocWeight[j] > max_weight)
         {
          max_weight = assocWeight[j];
          matched_mc_index = assocChSimID[j];
        }
      }

	  // build the map
      assoc_map_mc_to_rc[matched_mc_index] = rc_index;
      assoc_map_rc_to_mc[rc_index] = matched_mc_index;
    }
        // MC --> RC associations (kv; key,value for map)
   // std::cout << "MC to RC map:\n";
   // for (const auto& kv : assoc_map_mc_to_rc) {
   //     std::cout << "MC index " << kv.first << " ---> RC index " << kv.second << "\n";
   // }

    // RC --> MC associations
   // std::cout << "\nRC to MC map:\n";
   // for (const auto& kv : assoc_map_rc_to_mc) {
    //    std::cout << "RC index " << kv.first << " --->  MC index " << kv.second << "\n";
   // }


      // Loop over primary particles
    int nMcPart = 0;
    for(int imc=0; imc<nMCPart; imc++)
    {
     if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0)
     {
       double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
       if(dist < 1e-4)
       {
		  // count primary charged particles within |eta| < 3.5
        TVector3 mc_mom(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);
        double mcEta = mc_mom.PseudoRapidity();
        if(fabs(mcEta) < 3.5) nMcPart++;
      }
    }
  }
  // Primary charged multiplicity at generated level with acceptance cut
  hMcMult->Fill(nMcPart);


  for(int imc=0; imc<nMCPart; imc++)
  {
   if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0)
   {
     double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
     if(dist < 1e-4)
     {		  
		  // check if the MC particle is reconstructed
      int rc_index = -1;
      if(assoc_map_mc_to_rc.find(imc) != assoc_map_mc_to_rc.end()) rc_index = assoc_map_mc_to_rc[imc];

      if(rc_index>=0)
      {
        TVector3 dcaToVtx = GetDCAToPrimaryVertex(rc_index, vertex_rc);

        int ip = -1;
        if(fabs(mcPartPdg[imc]) == pdg_p) ip = 0;
        if(fabs(mcPartPdg[imc]) == pdg_k) ip = 1;
        if(fabs(mcPartPdg[imc]) == pdg_pi) ip = 2;
        if(ip>=0)
        {
         TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
         if(ip<3)
         {
           hRcPrimPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Perp());
           hRcPrimPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
         }

         double fill1[] = {mom.Pt(), mom.Eta(), dcaToVtx.x(), nMcPart*1.};
         double fill2[] = {mom.Pt(), mom.Eta(), dcaToVtx.y(), nMcPart*1.};
         double fill3[] = {mom.Pt(), mom.Eta(), dcaToVtx.z(), nMcPart*1.};
         hPrimTrkDcaToRCVtx[ip][0]->Fill(fill1);
         hPrimTrkDcaToRCVtx[ip][1]->Fill(fill2);
         hPrimTrkDcaToRCVtx[ip][2]->Fill(fill3);
       }
     }
   }
 }
}

 // Look for Λc⁺ → p K⁻ π⁺
bool hasLc = false;
vector<int> mc_index_Lcp_p;
vector<int> mc_index_Lcp_k;
vector<int> mc_index_Lcp_pi;
mc_index_Lcp_p.clear();      
mc_index_Lcp_k.clear();
mc_index_Lcp_pi.clear();

for(int imc=0; imc<nMCPart; imc++)
{
 if(fabs(mcPartPdg[imc]) != pdg_Lcp) continue;
 hEventStat->Fill(1.5);

  // printf("MC Particle (Momentum) = (%f, %f, %f) \n",mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);		  	
 int nDuaghters = mcPartDaughter_end[imc]-mcPartDaughter_begin[imc];
 if(nDuaghters!=3) continue;

	  // find Lc+ that decay into pKPi
 bool is_pkpi_decay = false;	  
 int daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
 int daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
 int daug_index_3 = mcPartDaughter_index[mcPartDaughter_begin[imc]+2];	  
 int daug_pdg_1 = mcPartPdg[daug_index_1];
 int daug_pdg_2 = mcPartPdg[daug_index_2];
 int daug_pdg_3 = mcPartPdg[daug_index_3];

 if( (fabs(daug_pdg_1)==pdg_p && fabs(daug_pdg_2)==pdg_k && fabs(daug_pdg_3)==pdg_pi) || (fabs(daug_pdg_1)==pdg_p && fabs(daug_pdg_2)==pdg_pi && fabs(daug_pdg_3)==pdg_k)  || 
   (fabs(daug_pdg_1)==pdg_k && fabs(daug_pdg_2)==pdg_p && fabs(daug_pdg_3)==pdg_pi) || (fabs(daug_pdg_1)==pdg_k && fabs(daug_pdg_2)==pdg_pi && fabs(daug_pdg_3)==pdg_p) ||
   (fabs(daug_pdg_1)==pdg_pi && fabs(daug_pdg_2)==pdg_k && fabs(daug_pdg_3)==pdg_p) || (fabs(daug_pdg_1)==pdg_pi && fabs(daug_pdg_2)==pdg_p && fabs(daug_pdg_3)==pdg_k) )
 {
   is_pkpi_decay = true;
 }
 if(!is_pkpi_decay) continue;
   // Efficiency for Lc+
 TLorentzVector mc_part_gen;
 mc_part_gen.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
 // Fill the momentum at Gen level Aug3, 2025
 float mcRap_gen = mc_part_gen.Rapidity();
 float mcPt_gen = mc_part_gen.Pt();  
 pt_Lcp_gen = mcPt_gen;
 y_Lcp_gen = mcRap_gen;
 tree_gen->Fill();

 hEventStat->Fill(2.5);

 if((fabs(daug_pdg_1)==pdg_p && fabs(daug_pdg_2)==pdg_k && fabs(daug_pdg_3)==pdg_pi))
 {
   mc_index_Lcp_p.push_back(daug_index_1);
   mc_index_Lcp_k.push_back(daug_index_2);
   mc_index_Lcp_pi.push_back(daug_index_3);
 }
 else if((fabs(daug_pdg_1)==pdg_p && fabs(daug_pdg_2)==pdg_pi && fabs(daug_pdg_3)==pdg_k))
 {
   mc_index_Lcp_p.push_back(daug_index_1);
   mc_index_Lcp_k.push_back(daug_index_3);
   mc_index_Lcp_pi.push_back(daug_index_2);
 }
 else if((fabs(daug_pdg_1)==pdg_k && fabs(daug_pdg_2)==pdg_p && fabs(daug_pdg_3)==pdg_pi))
 {
   mc_index_Lcp_p.push_back(daug_index_2);
   mc_index_Lcp_k.push_back(daug_index_1);
   mc_index_Lcp_pi.push_back(daug_index_3);
 }
 else if((fabs(daug_pdg_1)==pdg_k && fabs(daug_pdg_2)==pdg_pi && fabs(daug_pdg_3)==pdg_p))
 {
   mc_index_Lcp_p.push_back(daug_index_3);
   mc_index_Lcp_k.push_back(daug_index_1);
   mc_index_Lcp_pi.push_back(daug_index_2);
 }
 else if((fabs(daug_pdg_1)==pdg_pi && fabs(daug_pdg_2)==pdg_k && fabs(daug_pdg_3)==pdg_p))
 {
   mc_index_Lcp_p.push_back(daug_index_3);
   mc_index_Lcp_k.push_back(daug_index_2);
   mc_index_Lcp_pi.push_back(daug_index_1);
 }	    	    	    	    
   else 
   {
     mc_index_Lcp_p.push_back(daug_index_2);
     mc_index_Lcp_k.push_back(daug_index_3);
     mc_index_Lcp_pi.push_back(daug_index_1);
   }
   hasLc = true;

	  // Lcp kinematics
   TLorentzVector mc_mom_vec;
   mc_mom_vec.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
//printf("Mass  = %f \n",mcPartMass[imc]);

   double mcRap = mc_mom_vec.Rapidity();
   double mcPt = mc_mom_vec.Pt();
   hMCLcpPtRap->Fill(mcRap, mcPt); 

	  // decay daughter kinematics
   for(int ip = 0; ip<3; ip++)
   {
     int mc_part_index;
     if(ip==0) mc_part_index = mc_index_Lcp_p[mc_index_Lcp_p.size()-1];
     if(ip==1) mc_part_index = mc_index_Lcp_k[mc_index_Lcp_k.size()-1];
     if(ip==2) mc_part_index = mc_index_Lcp_pi[mc_index_Lcp_pi.size()-1];

     TLorentzVector mc_part_vec;
     mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
     if(ip==0) hMcPPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
     if(ip==1) hMcKPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
     if(ip==2) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());	      

     int rc_part_index = -1;
     if(assoc_map_mc_to_rc.find(mc_part_index) != assoc_map_mc_to_rc.end()) rc_part_index = assoc_map_mc_to_rc[mc_part_index];
     if (debug) printf("Rec: ProngNo., MC Index, Reco Index = (%d, %d, %d) \n",ip,mc_part_index,rc_part_index);
     if(rc_part_index>=0)
     {
      TVector3 dcaToVtx = GetDCAToPrimaryVertex(rc_part_index, vertex_rc);
      TVector3 mom(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
      hRcSecPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
      hRcSecPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());

		  //printf("Sec %d: (%2.4f, %2.4f, %2.4f), mcStartPoint = (%2.4f, %2.4f, %2.4f)\n", rc_part_index, pos.x(), pos.y(), pos.z(), mcPartVx[mc_part_index], mcPartVy[mc_part_index], mcPartVz[mc_part_index]);
    }
  }
}

      // Get reconstructed protons, kaons and pions
      hNRecoVtx->Fill(CTVx.GetSize());
      const int pid_mode = 0; // 0 - truth; 1 - realistic
      vector<unsigned int> p_index;
      vector<unsigned int> k_index;      
      vector<unsigned int> pi_index;
      p_index.clear();
      k_index.clear();      
      pi_index.clear();

      for(unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
      {	  
       if(pid_mode==0)
       {
         int iSimPartID = -1;
         if(assoc_map_rc_to_mc.find(rc_index) != assoc_map_rc_to_mc.end()) iSimPartID = assoc_map_rc_to_mc[rc_index];
         if(iSimPartID>=0)
         {
          if(fabs(mcPartPdg[iSimPartID]) == pdg_p) p_index.push_back(rc_index);	
          if(fabs(mcPartPdg[iSimPartID]) == pdg_k) k_index.push_back(rc_index);
          if(fabs(mcPartPdg[iSimPartID]) == pdg_pi) pi_index.push_back(rc_index);
        }
      }
      else if(pid_mode==1)
      {
       if(fabs(rcPdg[rc_index]) == pdg_p) p_index.push_back(rc_index);
       if(fabs(rcPdg[rc_index]) == pdg_k) k_index.push_back(rc_index);
       if(fabs(rcPdg[rc_index]) == pdg_pi) pi_index.push_back(rc_index);	      
     }
   }

  // printf("Proton, Kaon, and Pion vector sizes: = (%d, %d, %d) \n",p_index.size(), k_index.size(),pi_index.size()); 
      // Combinatorics proton, kaon, and pion 
   for(unsigned int i=0; i<p_index.size(); i++)
   {
     TVector3 dcaToVtx_p = GetDCAToPrimaryVertex(p_index[i], vertex_rc);
     int q_proton = rcCharge[p_index[i]];

     for(unsigned int j=0; j<k_index.size(); j++)
     {
      TVector3 dcaToVtx_k = GetDCAToPrimaryVertex(k_index[j], vertex_rc);
      int q_kaon = rcCharge[k_index[j]];

      for(unsigned int k=0; k<pi_index.size(); k++)
      {
       TVector3 dcaToVtx_pi = GetDCAToPrimaryVertex(pi_index[k], vertex_rc);
       int q_pion = rcCharge[pi_index[k]];

       if ((q_proton == +1 && q_kaon == -1 && q_pion == +1) || (q_proton == -1 && q_kaon == +1 && q_pion == -1))       
       {
        bool is_Lcp_pkpi = false;

        int  mc_index_p = -1, mc_index_k = -1, mc_index_pi = -1;
        if(assoc_map_rc_to_mc.find(p_index[i]) != assoc_map_rc_to_mc.end()) mc_index_p = assoc_map_rc_to_mc[p_index[i]];
        if(assoc_map_rc_to_mc.find(k_index[j])  != assoc_map_rc_to_mc.end()) mc_index_k  = assoc_map_rc_to_mc[k_index[j]];
        if(assoc_map_rc_to_mc.find(pi_index[k])  != assoc_map_rc_to_mc.end()) mc_index_pi  = assoc_map_rc_to_mc[pi_index[k]];

        for(unsigned int idaugh=0; idaugh<mc_index_Lcp_pi.size(); idaugh++)
        {
          if(mc_index_p==mc_index_Lcp_p[idaugh] && mc_index_k==mc_index_Lcp_k[idaugh] && mc_index_pi==mc_index_Lcp_pi[idaugh])
          {
           is_Lcp_pkpi = true;
           break;
         }
       }

       float dcaDaughters[3], cosTheta, decayLength, V0DcaToVtx, cosTheta_xy, sigma_vtx;
       TLorentzVector parent = GetDCADaughters(p_index[i], k_index[j], pi_index[k], vertex_rc, dcaDaughters, cosTheta, cosTheta_xy, decayLength, V0DcaToVtx, sigma_vtx);


       if(is_Lcp_pkpi)
       {
        hEventStat->Fill(3.5);
        if (q_proton == 1 && q_kaon == -1 && q_pion == 1)
        hEventStat->Fill(4.5);   // Λc⁺
      else if (q_proton == -1 && q_kaon == 1 && q_pion == -1)
        hEventStat->Fill(5.5);   // Λc⁻

      h3PairDca12[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[0]);
      h3PairDca23[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[1]); 
      h3PairDca13[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[2]);		      
      h3PairCosTheta[0]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
      h3PairDca[0]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
      h3PairDecayLength[0]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
		      //printf("Signal: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
      h3InvMass[0][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());

                // Toplogical Variables for Signal
      d0_p_sig = dcaToVtx_p.Mag();
      d0_k_sig = dcaToVtx_k.Mag();
      d0_pi_sig = dcaToVtx_pi.Mag();		      
      d0xy_p_sig = dcaToVtx_p.Perp();
      d0xy_k_sig = dcaToVtx_k.Perp();
      d0xy_pi_sig = dcaToVtx_pi.Perp();		      
      sum_d0xy_sig = sqrt(d0xy_p_sig*d0xy_p_sig+d0xy_k_sig*d0xy_k_sig+d0xy_pi_sig*d0xy_pi_sig);		      
      dca_12_sig = *min_element(dcaDaughters, dcaDaughters + 3);
      dca_Lcp_sig = V0DcaToVtx;
      decay_length_sig = decayLength;
      costheta_sig = cosTheta;
      costhetaxy_sig = cosTheta_xy;
      pt_Lcp_sig = parent.Pt();
      y_Lcp_sig = parent.Rapidity();
      mass_Lcp_sig = parent.M();
      sigma_vtx_sig = sigma_vtx;
      mult_sig = nMcPart;
      tree_sig->Fill();  

    }
    else
    {
      hEventStat->Fill(6.5);
      h3PairDca12[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[0]);
      h3PairDca23[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[1]); 
      h3PairDca13[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[2]);			      
      h3PairCosTheta[1]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
      h3PairDca[1]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
      h3PairDecayLength[1]->Fill(parent.Pt(), parent.Rapidity(), decayLength);

		      //printf("Bkg: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
      h3InvMass[1][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());

		      // Toplogical Variables for Bkg
      d0_p_bkg = dcaToVtx_p.Mag();
      d0_k_bkg = dcaToVtx_k.Mag();
      d0_pi_bkg = dcaToVtx_pi.Mag();		      
      d0xy_p_bkg = dcaToVtx_p.Perp();
      d0xy_k_bkg = dcaToVtx_k.Perp();
      d0xy_pi_bkg = dcaToVtx_pi.Perp();		      
      sum_d0xy_bkg = sqrt(d0xy_p_bkg*d0xy_p_bkg+d0xy_k_bkg*d0xy_k_bkg+d0xy_pi_bkg*d0xy_pi_bkg);	      
      dca_12_bkg = *min_element(dcaDaughters, dcaDaughters + 3);
      dca_Lcp_bkg = V0DcaToVtx;
      decay_length_bkg = decayLength;
      costheta_bkg = cosTheta;
      costhetaxy_bkg = cosTheta_xy;
      pt_Lcp_bkg = parent.Pt();
      y_Lcp_bkg = parent.Rapidity();
      mass_Lcp_bkg = parent.M();
      sigma_vtx_bkg = sigma_vtx;
      mult_bkg = nMcPart;
      tree_bkg->Fill();
    } 

  }
}
}
}

nevents++;
}


file_gen->cd();  
tree_gen->Write();
file_gen->Close();

file_signal->cd();  
tree_sig->Write();
file_signal->Close();  

file_bkg->cd();  
tree_bkg->Write();
file_bkg->Close(); 

TFile *outfile = new TFile(outname.Data(), "recreate");
hEventStat->SetMarkerSize(2);
hEventStat->Write("");
hMcMult->Write();
hMcVtxX->Write();
hMcVtxY->Write();
hMcVtxZ->Write();
hRecVtxX->Write();
hRecVtxY->Write();
hRecVtxZ->Write();

hLcpDecayVxVy->Write();
hLcpDecayVrVz->Write();

hMCLcpPtRap->Write();
hMcPPtEta->Write();
hMcPPtEtaReco->Write();
hMcPiPtEta->Write();
hMcPiPtEtaReco->Write();
hMcKPtEta->Write();
hMcKPtEtaReco->Write();

hNRecoVtx->Write();

for(int ip=0; ip<3; ip++)
{
  hRcSecPartLocaToRCVtx[ip]->Write();
  hRcSecPartLocbToRCVtx[ip]->Write();
  hRcPrimPartLocaToRCVtx[ip]->Write();
  hRcPrimPartLocbToRCVtx[ip]->Write();
}

for(int i=0; i<3; i++)
{
  for(int j=0; j<3; j++)
  {
   hPrimTrkDcaToRCVtx[i][j]->Write();
 }
}

for(int i=0; i<2; i++)
{
  h3PairDca12[i]->Write();
  h3PairDca23[i]->Write();
  h3PairDca13[i]->Write();            
  h3PairCosTheta[i]->Write();
  h3PairDca[i]->Write();
  h3PairDecayLength[i]->Write();
}

for(int i=0; i<2; i++)
{
  for(int j=0; j<2; j++)
  {
   h3InvMass[i][j]->Write();
 }
}


outfile->Close();

}

//======================================
TVector3 GetDCAToPrimaryVertex(const int index, TVector3 vtx)
{
  //printf("check %d: (%2.4f, %2.4f, %2.4f) =? (%2.4f, %2.4f, %2.4f)\n", index, mom.x(), mom.y(), mom.z(), rcMomPx2->At(index), rcMomPy2->At(index), rcMomPz2->At(index));

  TVector3 pos(rcTrkLoca2->At(index) * sin(rcTrkPhi2->At(index)) * -1 * millimeter, rcTrkLoca2->At(index) * cos(rcTrkPhi2->At(index)) * millimeter, rcTrkLocb2->At(index) * millimeter);
  TVector3 mom(rcMomPx2->At(index), rcMomPy2->At(index), rcMomPz2->At(index));

  StPhysicalHelix pHelix(mom, pos, bField * tesla, rcCharge2->At(index));

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);
  
  pHelix.moveOrigin(pHelix.pathLength(vtx_tmp));
  TVector3 dcaToVtx = pHelix.origin() - vtx_tmp;

  dcaToVtx.SetXYZ(dcaToVtx.x()/millimeter, dcaToVtx.y()/millimeter, dcaToVtx.z()/millimeter);
  
  return dcaToVtx;
}

//======================================
// Local to global conversion 
// x = - l0 Sin (phi), y = l0 Cos (phi), z = l1 
TLorentzVector GetDCADaughters(const int index1, const int index2, const int index3, TVector3 vtx,
  float *dcaDaughters, float &cosTheta, float &cosTheta_xy, float &decayLength, float &V0DcaToVtx, float &sigma_vtx)
{
  // -- get helix
  TVector3 pos1(rcTrkLoca2->At(index1) * sin(rcTrkPhi2->At(index1)) * -1 * millimeter, rcTrkLoca2->At(index1) * cos(rcTrkPhi2->At(index1)) * millimeter, rcTrkLocb2->At(index1) * millimeter);
  TVector3 pos2(rcTrkLoca2->At(index2) * sin(rcTrkPhi2->At(index2)) * -1 * millimeter, rcTrkLoca2->At(index2) * cos(rcTrkPhi2->At(index2)) * millimeter, rcTrkLocb2->At(index2) * millimeter);
  TVector3 pos3(rcTrkLoca2->At(index3) * sin(rcTrkPhi2->At(index3)) * -1 * millimeter, rcTrkLoca2->At(index3) * cos(rcTrkPhi2->At(index3)) * millimeter, rcTrkLocb2->At(index3) * millimeter);
  
  TVector3 mom1(rcMomPx2->At(index1), rcMomPy2->At(index1), rcMomPz2->At(index1));
  TVector3 mom2(rcMomPx2->At(index2), rcMomPy2->At(index2), rcMomPz2->At(index2));
  TVector3 mom3(rcMomPx2->At(index3), rcMomPy2->At(index3), rcMomPz2->At(index3));

  float charge1 = rcCharge2->At(index1);
  float charge2 = rcCharge2->At(index2);
  float charge3 = rcCharge2->At(index3);  
  
  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, charge1);
  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla, charge2);
  StPhysicalHelix p3Helix(mom3, pos3, bField * tesla, charge3);

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);

  // -- Full calculation for pathlengths
  // Point1
  pair<double, double> const ss_12 = p1Helix.pathLengths(p2Helix); 
  pair<double, double> const ss_13 = p1Helix.pathLengths(p3Helix);
  // Point2 
  pair<double, double> const ss_23 = p2Helix.pathLengths(p3Helix); 
  pair<double, double> const ss_21 = p2Helix.pathLengths(p1Helix); 
  // Point3
  pair<double, double> const ss_31 = p3Helix.pathLengths(p1Helix); 
  pair<double, double> const ss_32 = p3Helix.pathLengths(p2Helix);  
  
  // Point1 w.r.t. helices 2 and 3, later average will be taken 
  TVector3 const P1P2 = p1Helix.at(ss_12.first); TVector3 const P1P3 = p1Helix.at(ss_13.first); 
  // Point2 w.r.t. helices 1 and 3, later average will be taken 
  TVector3 const P2P1 = p2Helix.at(ss_21.first); TVector3 const P2P3 = p2Helix.at(ss_23.first); 
  // Point3 w.r.t. helices 1 and 2, later average will be taken 
  TVector3 const P3P1 = p3Helix.at(ss_31.first); TVector3 const P3P2 = p3Helix.at(ss_32.first); 

 // printf("P1P2 = (%2.4f, %2.4f, %2.4f) \t P1P3 = (%2.4f, %2.4f, %2.4f)\n", P1P2.x(), P1P2.y(), P1P2.z(), P1P3.x(), P1P3.y(), P1P3.z());
 // printf("P2P1 = (%2.4f, %2.4f, %2.4f) \t P2P3 = (%2.4f, %2.4f, %2.4f)\n", P2P1.x(), P2P1.y(), P2P1.z(), P2P3.x(), P2P3.y(), P2P3.z());
 // printf("P3P1 = (%2.4f, %2.4f, %2.4f) \t P3P2 = (%2.4f, %2.4f, %2.4f)\n", P3P1.x(), P3P1.y(), P3P1.z(), P3P2.x(), P3P2.y(), P3P2.z());
  // Evaluate the DCA Points
  // Between track1 and track2
  TVector3 const p1 = 0.5*(P1P2+P1P3);
  TVector3 const p2 = 0.5*(P2P1+P2P3);
  TVector3 const p3 = 0.5*(P3P1+P3P2);


  // cout << ss.first << "  " << ss.second << endl;
  // printf("p1AtDcaToP2 origin = (%2.4f, %2.4f, %2.4f)\n", p1AtDcaToP2.x(), p1AtDcaToP2.y(), p1AtDcaToP2.z());
  // printf("p2AtDcaToP1 origin = (%2.4f, %2.4f, %2.4f)\n", p2AtDcaToP1.x(), p2AtDcaToP1.y(), p2AtDcaToP1.z());
  
  // -- calculate DCA between daughters 
  // Three dcaDaughters; dcaDaughters_12, dcaDaughters_23, dcaDaughters_13
  dcaDaughters[0] = (p1 - p2).Mag()/millimeter;
  dcaDaughters[1] = (p2 - p3).Mag()/millimeter;
  dcaDaughters[2] = (p3 - p1).Mag()/millimeter; 

  // -- calculate Momentum of each daughters at the DCA point
  TVector3 const p1MomAtDcap2 = p1Helix.momentumAt(ss_12.first,  bField * tesla);
  TVector3 const p1MomAtDcap3 = p1Helix.momentumAt(ss_13.first,  bField * tesla);
  TVector3 const p2MomAtDca3 = p2Helix.momentumAt(ss_23.first, bField * tesla);
  TVector3 const p2MomAtDca1 = p2Helix.momentumAt(ss_21.first, bField * tesla);
  TVector3 const p3MomAtDcap1 = p3Helix.momentumAt(ss_31.first, bField * tesla); 
  TVector3 const p3MomAtDcap2 = p3Helix.momentumAt(ss_32.first, bField * tesla); 

  // printf("p1MomAtDcap2 = (%2.4f, %2.4f, %2.4f) \t p1MomAtDcap3 = (%2.4f, %2.4f, %2.4f)\n", p1MomAtDcap2.x(), p1MomAtDcap2.y(), p1MomAtDcap2.z(), p1MomAtDcap3.x(), p1MomAtDcap3.y(), p1MomAtDcap3.z());
  // printf("p2MomAtDca3 = (%2.4f, %2.4f, %2.4f) \t p2MomAtDca1 = (%2.4f, %2.4f, %2.4f)\n", p2MomAtDca3.x(), p2MomAtDca3.y(), p2MomAtDca3.z(), p2MomAtDca1.x(), p2MomAtDca1.y(), p2MomAtDca1.z());
  // printf("p3MomAtDcap1 = (%2.4f, %2.4f, %2.4f) \t p3MomAtDcap2 = (%2.4f, %2.4f, %2.4f)\n", p3MomAtDcap1.x(), p3MomAtDcap1.y(), p3MomAtDcap1.z(), p3MomAtDcap2.x(), p3MomAtDcap2.y(), p3MomAtDcap2.z());

  
  TVector3 const p1MomAtDca = 0.5*(p1MomAtDcap2+p1MomAtDcap3);  
  TVector3 const p2MomAtDca = 0.5*(p2MomAtDca3+p2MomAtDca1); 
  TVector3 const p3MomAtDca = 0.5*(p3MomAtDcap1+p3MomAtDcap2); 
  
  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+gProtonMass*gProtonMass));
  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+gKaonMass*gKaonMass));
  TLorentzVector p3FourMom(p3MomAtDca, sqrt(p3MomAtDca.Mag2()+gPionMass*gPionMass));  
  
  TLorentzVector parent = p1FourMom + p2FourMom + p3FourMom;

  // -- calculate cosThetaStar
  // TLorentzVector pairFourMomReverse(-parent.Px(), -parent.Py(), -parent.Pz(), parent.E());
  // TLorentzVector p1FourMomStar = p1FourMom;
  // p1FourMomStar.Boost(pairFourMomReverse.Vect());
  // float cosThetaStar = std::cos(p1FourMomStar.Vect().Angle(parent.Vect()));

  //Approach1: Assign the decay vertex the centroid coordinate of the triangle
   TVector3 decayVertex = 1./3*(p1 + p2+ p3);
 
  //Approach2: Assign the decay vertex as the circumcenter of the triangle (p1, p2, p3)
  // Compute side vectors of triangle; a = p1, b = p2, c = p3
    // Side vectors
  TVector3 ab = p2 - p1; 
  TVector3 ac = p3 - p1;

    // Cross product of ab and ac
  TVector3 abXac = ab.Cross(ac);

  double ab2 = ab.Mag2();
  double ac2 = ac.Mag2();
  double abXac2 = abXac.Mag2();

    // Vector from point a to the circhumspere center
  TVector3 a_to_CircumsphereCenter = (abXac.Cross(ab) * ac2 + ac.Cross(abXac) * ab2) * (1.0 / (2.0 * abXac2));

    // Magnitude of the above vector is radius
  double radius = a_to_CircumsphereCenter.Mag();

    // Coordinate of the Circumsphere center is the decayvertex in 3D (A point at equal distance of three PCA)
 
  //TVector3 decayVertex = p1 + a_to_CircumsphereCenter; // uncomment it if you want to define as circhumspere center

  // printf("dist_p1_to_decayvertex = (%1.4f) \t dist_p2_to_decayvertex = (%1.4f) \t dist_p3_to_decayvertex = (%1.4f) \n", (p1-decayVertex).Mag(), (p2-decayVertex).Mag(), (p3-decayVertex).Mag());

    // Confirm by evaluating the distance with vertices

  sigma_vtx = sqrt((p1-decayVertex).Mag2()+(p2-decayVertex).Mag2()+(p3-decayVertex).Mag2())/millimeter;

  // -- calculate pointing angle and decay length with respect to primary vertex
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  TVector3 vtxToV0 = decayVertex - vtx_tmp;
  TVector3 vtxToV0_xy(vtxToV0.x(), vtxToV0.y(), 0.);
  TVector3 parent_xy(parent.Vect().x(),parent.Vect().y(),0.);
  float pointingAngle = vtxToV0.Angle(parent.Vect());
  float pointingAngle_xy = vtxToV0_xy.Angle(parent_xy);
  cosTheta = std::cos(pointingAngle);
  cosTheta_xy = std::cos(pointingAngle_xy);  
  decayLength = vtxToV0.Mag()/millimeter;

  // -- calculate V0 DCA to primary vertex
  V0DcaToVtx = decayLength * std::sin(pointingAngle);

  //TVector3 dcaToVtx = GetDCAToPrimaryVertex(parent.Vect(), decayVertex, 0, vtx);
  //V0DcaToVtx = dcaToVtx.Mag();

  return parent;
}

