// D0-tagged Jets

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
#include "TMinuit.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "Math/Minimizer.h"
#endif


#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"

using namespace std;

StPhysicalHelix* gHelix1 = nullptr;
StPhysicalHelix* gHelix2 = nullptr;

const double gPionMass = 0.13957;
const double gKaonMass = 0.493677;

const double twoPi = 2.*3.1415927;
const double eMass = 0.000511;

const double bField = -1.7; // Tesla

//JetClustering Parameters:
double minCstPt            = 0.2 ;				 // minimum pT of objects
double maxCstPt            = 100.;  			 // maximum pT of objects
double minJetPt            = 1.0 ;  			 // minimum jet pT
double ghostMaxRap         = 3.5;   			 // maximum rapidity of ghosts
double ghostArea           = 0.01; 				 // area per ghost
int numGhostRepeat         = 1;                              // reuse count of ghosts
int removeelectrons        = 1;
int nhitcut                = 0;

TVector3 getDcaToVtx(const int index, TVector3 vtx);
void fcnVertexFit(int& npar, double* grad, double& fval, double* par, int iflag);
void getDecayVertex_Chi2fit(const int index1, const int index2, double &s1, double &s2, TVector3 &vertex, double &chi2_ndf, double * parFitErr);

TLorentzVector getPairParent(const int index1, const int index2, TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &cosTheta_xy, float &decayLength, float &V0DcaToVtx, float &sigma_vtx, TVector3 &decayVertex,TVector3 &decayVertex_ana, double &chi2_ndf,  double * parFitErr);

TTreeReaderArray<float> *rcMomPx2;
TTreeReaderArray<float> *rcMomPy2;
TTreeReaderArray<float> *rcMomPz2;
TTreeReaderArray<float> *rcCharge2;

TTreeReaderArray<float> *rcTrkLoca2;
TTreeReaderArray<float> *rcTrkLocb2;
TTreeReaderArray<float> *rcTrkTheta2;
TTreeReaderArray<float> *rcTrkPhi2;
TTreeReaderArray<std::array<float, 21>> *rcTrkCov;



// Define function for Chi2 minimization between two helices	
struct Chi2Minimization {    
     StPhysicalHelix fhelix1, fhelix2;
     std::array<float, 21> fcov1, fcov2; // full covariance matrix
     
     Chi2Minimization(StPhysicalHelix helix1, StPhysicalHelix helix2, std::array<float, 21> cov1, std::array<float, 21> cov2) : fhelix1(helix1),fhelix2(helix2), fcov1(cov1), fcov2(cov2) {}
    // Implementation of the function to be minimized
    double operator() (const double *par) {
    double x = par[0];
    double y = par[1];
    double z = par[2];
    double s1 = par[3];
    double s2 = par[4];   
    double f = 0;
    TVector3 vertex(x, y, z);
    TVector3 p1 = fhelix1.at(s1);
    TVector3 p2 = fhelix2.at(s2);
    TVector3 mom1 = fhelix1.momentumAt(s1,  bField * tesla);
    TVector3 mom2 = fhelix2.momentumAt(s2,  bField * tesla);
    
   // x= −l0 sinϕ , y=l0 cosϕ , z=l
   // Recalculate l0 at PCA for error propagation
     float l0_track1 = p1.Pt(); float l1_track1 = p1.Z(); double phi_track1 = mom1.Phi();
     float l0_track2 = p2.Pt(); float l1_track2 = p2.Z(); double phi_track2 = mom2.Phi();
     // Track1: σx^2​=sin^2ϕ⋅σℓ0​^2​+ℓ0^2​cos^2ϕ⋅σϕ^2​+2⋅ℓ0​sinϕcosϕ⋅Cov(ℓ0​,ϕ)
     float sigx1_2 = sin(phi_track1)*sin(phi_track1)*fcov1[0] + l0_track1*l0_track1*cos(phi_track1)*cos(phi_track1)*fcov1[5]+ 2.0*l0_track1*sin(phi_track1)*cos(phi_track1)*fcov1[3]; 
     float sigx2_2 = sin(phi_track2)*sin(phi_track2)*fcov2[0] + l0_track2*l0_track2*cos(phi_track2)*cos(phi_track2)*fcov2[5]+ 2.0*l0_track2*sin(phi_track2)*cos(phi_track2)*fcov2[3];
  
     // σy^2​=cos^2ϕ⋅σℓ0^​2​+ℓ0^2​sin^2ϕ⋅σϕ^2​−2⋅ℓ0​sinϕcosϕ⋅Cov(ℓ0​,ϕ)
     float sigy1_2 = cos(phi_track1)*cos(phi_track1)*fcov1[0] + l0_track1*l0_track1*sin(phi_track1)*sin(phi_track1)*fcov1[5]-2.0*l0_track1*sin(phi_track1)*cos(phi_track1)*fcov1[3]; 
     float sigy2_2 = cos(phi_track2)*cos(phi_track2)*fcov2[0] + l0_track2*l0_track2*sin(phi_track2)*sin(phi_track2)*fcov2[5]-2.0*l0_track2*sin(phi_track2)*cos(phi_track2)*fcov2[3]; 
     
     // σz^2​
     float sigz1_2 = fcov1[2];
     float sigz2_2 = fcov2[2];  
     // convert to mm
    double d1_x = 10.*(vertex - p1).X(); double d2_x = 10.*(vertex - p2).X(); 
    double d1_y = 10.*(vertex - p1).Y(); double d2_y = 10.*(vertex - p2).Y();
    double d1_z = 10.*(vertex - p1).Z(); double d2_z = 10.*(vertex - p2).Z();
    
     f = d1_x*d1_x/sigx1_2 + d2_x*d2_x/sigx2_2 + d1_y*d1_y/sigy1_2 + d2_y*d2_y/sigy2_2 + d1_z*d1_z/sigz1_2+ d2_z*d2_z/sigz2_2; // chi2
      return f;
	}
	};

int main(int argc, char **argv)
{
	
  TString listname;
  TString outname;
  TString collname;
  float R_value;
  TString signal, bkg, gen;
  TString signalmc; 

  if(argc==1)
    {
      listname  = "test.list";
      outname = "test.root";
      collname = "ep";   
      R_value = 1.0;
      signal = "signal.root";
      bkg = "bkg.root";
      gen = "gen.root";
      signalmc = "signalmc.root";
    }

  if(argc>=8)
    {
      listname = argv[1];
      outname = argv[2];
      collname = argv[3];      
      R_value = std::atof(argv[4]);       
      signal = argv[5]; 
      bkg = argv[6];
      gen = argv[7];
      signalmc = "signalmc.root";
    }

  TChain *chain = new TChain("events");

  int nfiles = 0;
  char filename[512];
  ifstream *inputstream = new ifstream;
  inputstream->open(listname.Data());
  if (!inputstream->is_open())
  {
  printf("[e] Cannot open file list: %s\n", listname.Data());
  return 0; // or handle as needed
  }

  while (inputstream->good())
  {
  inputstream->getline(filename, 512);
  if (inputstream->good())
  {
    TFile *ftmp = TFile::Open(filename, "READ");
    if (!ftmp || !ftmp->IsOpen() || !ftmp->GetNkeys())
    {
      printf("[e] Skipping bad file: %s\n", filename);
      if (ftmp) { ftmp->Close(); delete ftmp; }
      continue; 
    }
    cout << "[i] Add " << nfiles << "th file: " << filename << endl;
    chain->Add(filename);
    nfiles++;

    ftmp->Close(); // cleanup
    delete ftmp;
  }
}

  inputstream->close();

  printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());

  while (inputstream->good())
  {
  inputstream->getline(filename, 512);
  if (inputstream->good())
  {
    TFile *ftmp = TFile::Open(filename, "READ");
    if (!ftmp || !ftmp->IsOpen() || !ftmp->GetNkeys())
    {
      printf("[e] Skipping bad file: %s\n", filename);
      if (ftmp) { ftmp->Close(); delete ftmp; }
      continue; 
    }
    cout << "[i] Add " << nfiles << "th file: " << filename << endl;
    chain->Add(filename);
    nfiles++;

    ftmp->Close(); // cleanup
    delete ftmp;
  }
}

  inputstream->close();

  printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());

  TH1F *hEventStat = new TH1F("hEventStat", "Event statistics", 7, 0, 7);
  hEventStat->GetXaxis()->SetBinLabel(1, "MC events");
  hEventStat->GetXaxis()->SetBinLabel(2, "D0");
  hEventStat->GetXaxis()->SetBinLabel(3, "D0 -> pi+K");
  hEventStat->GetXaxis()->SetBinLabel(4, "Reco D0");
  hEventStat->GetXaxis()->SetBinLabel(5, "Reco Signal D0");
  hEventStat->GetXaxis()->SetBinLabel(6, "Reco Signal D0bar");
  hEventStat->GetXaxis()->SetBinLabel(7, "Reco Bkg D0");

  TH1F *hMcMult = new TH1F("hMcMult", "MC multiplicity (|#eta| < 3.5);N_{MC}", 50, 0, 50);

  // Secondary vertex with chi2 fit method
  TH3F *hRes_SVx_Helixfit = new TH3F("hRes_SVx_Helixfit", "Fit method: Residual of SV (X); p_{T} (GeV/c); y ; SVx_{rec}-SVx_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVy_Helixfit = new TH3F("hRes_SVy_Helixfit", "Fit method: Residual of SV (Y); p_{T} (GeV/c); y ; SVy_{rec}-SVy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVz_Helixfit = new TH3F("hRes_SVz_Helixfit", "Fit method: Residual of SV (Z); p_{T} (GeV/c); y ; SVz_{rec}-SVz_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0);   
 
  // Secondary vertex with analytical method
  TH3F *hRes_SVx_Helixana = new TH3F("hRes_SVx_Helixana", "Analytical method: Residual of SV (X); p_{T} (GeV/c); y ; SVx_{rec}-SVx_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVy_Helixana = new TH3F("hRes_SVy_Helixana", "Analytical method: Residual of SV (Y); p_{T} (GeV/c); y ; SVy_{rec}-SVy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVz_Helixana = new TH3F("hRes_SVz_Helixana", "Analytical method: Residual of SV (Z); p_{T} (GeV/c); y ; SVz_{rec}-SVz_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0);   
 
  // Secondary vertex  (XY) fit method
  TH3F *hRes_SVxy_Helixfit = new TH3F("hRes_SVxy_Helixfit", "Fit method: Residual of SV (XY); p_{T} (GeV/c); y ; SVxy_{rec}-SVxy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 

  // Secondary vertex  (XY) analytical method
  TH3F *hRes_SVxy_Helixana = new TH3F("hRes_SVxy_Helixana", "Analytical method: Residual of SV (XY); p_{T} (GeV/c); y ; SVxy_{rec}-SVxy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
    
  TH3F *hRes_SVx_Helixfit_pull = new TH3F("hRes_SVx_Helixfit_pull", "Fit method: Pull of SV (X); p_{T} (GeV/c); y ; SVx_{rec}-SVx_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVy_Helixfit_pull = new TH3F("hRes_SVy_Helixfit_pull", "Fit method: Pull of SV (Y); p_{T} (GeV/c); y ; SVy_{rec}-SVy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVz_Helixfit_pull = new TH3F("hRes_SVz_Helixfit_pull", "Fit method: Pull of SV (Z); p_{T} (GeV/c); y ; SVz_{rec}-SVz_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  
  TH1F *hchi2_vtx = new TH1F("hchi2_vtx", "Helix Calculation: Chi2/ndf; #chi^{2}/ndf; Entries (a.u.)", 1000, 0.0, 50.0); 
  TH1F *hchi2_vtx_sig = new TH1F("hchi2_vtx_sig", "Helix Calculation: Chi2/ndf; #chi^{2}/ndf; Entries (a.u.)", 1000, 0.0, 50.0); 
  TH1F *hchi2_vtx_bkg = new TH1F("hchi2_vtx_bkg", "Helix Calculation: Chi2/ndf; #chi^{2}/ndf; Entries (a.u.)", 1000, 0.0, 50.0);   

  TH1F *hMcVtxX = new TH1F("hMcVtxX", "x position of MC vertex;x (mm)", 100, -5.05, 4.95);
  TH1F *hMcVtxY = new TH1F("hMcVtxY", "y position of MC vertex;y (mm)", 500, -5.01, 4.99);
  TH1F *hMcVtxZ = new TH1F("hMcVtxZ", "z position of MC vertex;z (mm)", 400, -200, 200);
  
  TH1F *hPullVtxX = new TH1F("hPullVtxX", "Pull x position of MC vertex;(Vx_{rec}-Vx_{mc})/#sigma_{vx}", 100, -5.05, 4.95);
  TH1F *hPullVtxY = new TH1F("hPullVtxY", "Pull y position of MC vertex;(Vy_{rec}-Vy_{mc})/#sigma_{vy}", 500, -5.01, 4.99);
  TH1F *hPullVtxZ = new TH1F("hPullVtxZ", "Pull z position of MC vertex;(Vz_{rec}-Vz_{mc})/#sigma_{vz}", 400, -200, 200);

  TH2F *hD0DecayVxVy = new TH2F("hD0DecayVxVy", "D^{0} decay vertex to primary vertex;#Deltav_{x} (mm);#Deltav_{y} (mm)", 400, -1-0.0025, 1-0.0025, 400, -1-0.0025, 1-0.0025);
  TH2F *hD0DecayVrVz = new TH2F("hD0DecayVrVz", "D^{0} decay vertex to primary vertex;#Deltav_{z} (mm);#Deltav_{r} (mm)", 100, -2, 2, 100, -0.2, 1.8);

  TH2F *hMCD0PtRap = new TH2F("hMCD0PtRap", "MC D^{0};y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcPiPtEta = new TH2F("hMcPiPtEta", "MC #pi from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcPiPtEtaReco = new TH2F("hMcPiPtEtaReco", "RC #pi from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcKPtEta = new TH2F("hMcKPtEta", "MC K from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcKPtEtaReco = new TH2F("hMcKPtEtaReco", "RC K from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH1F *hNRecoVtx = new TH1F("hNRecoVtx", "Number of reconstructed vertices;N", 10, 0, 10);

  const char* part_name[3] = {"Pi", "K", "P"};
  const char* part_title[3] = {"#pi", "K", "P"};
  TH3F *hRcSecPartLocaToRCVtx[2];
  TH3F *hRcSecPartLocbToRCVtx[2];
  TH3F *hRcPrimPartLocaToRCVtx[2];
  TH3F *hRcPrimPartLocbToRCVtx[2];
  for(int i=0; i<2; i++)
    {
      hRcSecPartLocaToRCVtx[i] = new TH3F(Form("hRcSec%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcSecPartLocbToRCVtx[i] = new TH3F(Form("hRcSec%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
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

  TH3F *h3PairDca12[2];
  TH3F *h3PairCosTheta[2];
  TH3F *h3PairDca[2];
  TH3F *h3PairDecayLength[2];
  const char* pair_name[2] = {"signal", "bkg"};
  const char* pair_title[2] = {"Signal", "Background"};
  for(int i=0; i<2; i++)
    {
      h3PairDca12[i] = new TH3F(Form("h3PairDca12_%s", pair_name[i]), Form("%s pair DCA_{12};p_{T} (GeV/c);#eta;DCA_{12} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);

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
	  h3InvMass[i][j] = new TH3F(Form("h3InvMass_%s_%s", pair_name[i], cut_name[j]), "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 1.6, 2.0);
	}
    }
  
  // Fragmentaion variable  
  TH3F *h3sig_z = new TH3F("h3sig_z","Signal_D0;Z;y;M_D0(GeV/c^{2})",1000, -10, 10, 20, -5, 5, 100, 1.6, 2.0);    
  TH3F *h3bkg_z = new TH3F("h3bkg_z","Bkg_D0;Z;y;M_D0(GeV/c^{2})",1000, -10, 10, 20, -5, 5, 100, 1.6, 2.0);
  TH1F *hreco_eta = new TH1F("hreco_eta", "Eta distribution of reco particles ;eta", 500, -5,+5);

  TH1F *hmc_eta   = new TH1F("hmc_eta",   "MC Particle #eta;#eta;Number of Particles", 500, -5, 5);
  TH1F *hmc_eta_e = new TH1F("hmc_eta_e", "Electron #eta;#eta;Number of Particles",       500, -5, 5);
  TH1F *hmc_eta_pi = new TH1F("hmc_eta_pi", "Pion #eta;#eta;Number of Particles",         500, -5, 5);
  TH1F *hmc_eta_k  = new TH1F("hmc_eta_k",  "Kaon #eta;#eta;Number of Particles",         500, -5, 5);
  TH1F *hmc_eta_p  = new TH1F("hmc_eta_p",  "Proton #eta;#eta;Number of Particles",       500, -5, 5);

  TH1F *hreco_eta_e  = new TH1F("hreco_eta_e",  "Electron #eta;#eta;Number of Particles", 500, -5, 5);
  TH1F *hreco_eta_pi = new TH1F("hreco_eta_pi", "Pion #eta;#eta;Number of Particles",     500, -5, 5);
  TH1F *hreco_eta_k  = new TH1F("hreco_eta_k",  "Kaon #eta;#eta;Number of Particles",     500, -5, 5);
  TH1F *hreco_eta_p  = new TH1F("hreco_eta_p",  "Proton #eta;#eta;Number of Particles",   500, -5, 5);
  
  
  TTreeReader treereader(chain);
  // MC Particles 
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
  TTreeReaderArray<double> mcMomPx = {treereader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomPy = {treereader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomPz = {treereader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcEndPointX = {treereader, "MCParticles.endpoint.x"};
  TTreeReaderArray<double> mcEndPointY = {treereader, "MCParticles.endpoint.y"};
  TTreeReaderArray<double> mcEndPointZ = {treereader, "MCParticles.endpoint.z"};
  // Reconstructed Particles
  TTreeReaderArray<int> assocChSimID = {treereader, "_ReconstructedChargedParticleAssociations_sim.index"};
  TTreeReaderArray<int> assocChRecID = {treereader, "_ReconstructedChargedParticleAssociations_rec.index"};
  TTreeReaderArray<float> assocWeight = {treereader, "ReconstructedChargedParticleAssociations.weight"};
 
  TTreeReaderArray<float> rcMomPx = {treereader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> rcMomPy = {treereader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> rcMomPz = {treereader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> rcPosx = {treereader, "ReconstructedChargedParticles.referencePoint.x"};
  TTreeReaderArray<float> rcPosy = {treereader, "ReconstructedChargedParticles.referencePoint.y"};
  TTreeReaderArray<float> rcPosz = {treereader, "ReconstructedChargedParticles.referencePoint.z"};
  TTreeReaderArray<float> rcCharge = {treereader, "ReconstructedChargedParticles.charge"};
  TTreeReaderArray<int>   rcPdg = {treereader, "ReconstructedChargedParticles.PDG"};
  TTreeReaderArray<float> TrkRecoE = {treereader, "ReconstructedChargedParticles.energy"};  
  TTreeReaderArray<float> TrkRecoM = {treereader, "ReconstructedChargedParticles.mass"};
  
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
  rcTrkCov = new TTreeReaderArray<std::array<float, 21>>{treereader, "CentralCKFTrackParameters.covariance.covariance[21]"};

  TTreeReaderArray<float> CTVx = {treereader, "CentralTrackVertices.position.x"};
  TTreeReaderArray<float> CTVy = {treereader, "CentralTrackVertices.position.y"};
  TTreeReaderArray<float> CTVz = {treereader, "CentralTrackVertices.position.z"};
  TTreeReaderArray<int> CTVndf = {treereader, "CentralTrackVertices.ndf"};
  TTreeReaderArray<float> CTVchi2 = {treereader, "CentralTrackVertices.chi2"};
  TTreeReaderArray<float> CTVerr_xx = {treereader, "CentralTrackVertices.positionError.xx"};
  TTreeReaderArray<float> CTVerr_yy = {treereader, "CentralTrackVertices.positionError.yy"};
  TTreeReaderArray<float> CTVerr_zz = {treereader, "CentralTrackVertices.positionError.zz"};
  TTreeReaderArray<float> incQ2 = {treereader, "InclusiveKinematicsTruth.Q2"};
  TTreeReaderArray<float> incxB = {treereader, "InclusiveKinematicsTruth.x"};

  TTreeReaderArray<int> prim_vtx_index = {treereader, "PrimaryVertices_objIdx.index"};

  TTreeReaderArray<unsigned int> vtxAssocPart_begin = {treereader, "CentralTrackVertices.associatedParticles_begin"};
  TTreeReaderArray<unsigned int> vtxAssocPart_end = {treereader, "CentralTrackVertices.associatedParticles_end"};
  TTreeReaderArray<int> vtxAssocPart_index = {treereader, "_CentralTrackVertices_associatedParticles.index"};


  TTreeReaderArray<float> EvtQ2 ={treereader, "InclusiveKinematicsElectron.Q2"};
  TTreeReaderArray<float> Evtx =  {treereader, "InclusiveKinematicsElectron.x"};
  TTreeReaderArray<float> EvtQ2Gen = {treereader, "InclusiveKinematicsTruth.Q2"};
  TTreeReaderArray<float> EvtxGen = {treereader, "InclusiveKinematicsTruth.x"};

  TTreeReaderArray<int>  ScatElecRecoId = {treereader, "_InclusiveKinematicsElectron_scat.index"};
  TTreeReaderArray<int> ScatElecGenId = {treereader, "MCScatteredElectrons_objIdx.index"};
   
  TTreeReaderArray<int> TrkPartAssocRec = {treereader, "_ReconstructedChargedParticleAssociations_rec.index"};
  TTreeReaderArray<int> TrkPartAssocSim = {treereader, "_ReconstructedChargedParticleAssociations_sim.index"};
  TTreeReaderArray<unsigned int> TrkRecoNhits = {treereader, "CentralCKFTrajectories.nMeasurements"};

  // Gen particles
  TTreeReaderArray<float>  TrkGenE = {treereader, "GeneratedParticles.energy"};
  TTreeReaderArray<float> TrkGenPx = {treereader, "GeneratedParticles.momentum.x"};
  TTreeReaderArray<float> TrkGenPy = {treereader, "GeneratedParticles.momentum.y"};
  TTreeReaderArray<float> TrkGenPz = {treereader, "GeneratedParticles.momentum.z"};
  TTreeReaderArray<float> TrkGenM = {treereader, "GeneratedParticles.mass"};
  TTreeReaderArray<int> TrkGenPDG = {treereader, "GeneratedParticles.PDG"};
  TTreeReaderArray<float> TrkGenCharge = {treereader, "GeneratedParticles.charge"};
   

  
    // Create a ROOT file to store the Ntuple
   TFile *file_signal = new TFile(signal.Data(), "RECREATE");
   TTree *tree_sig = new TTree("treeMLSig", "treeMLSig"); 

  // Define variables to store in the Ntuple
  float d0_pi_sig, d0_k_sig, d0xy_pi_sig, d0xy_k_sig, sum_d0xy_sig, dca_12_sig, dca_D0_sig, decay_length_sig, xB_sig, Q2_sig;
  float costheta_sig, costhetaxy_sig, pt_D0_sig, y_D0_sig, mass_D0_sig, sigma_vtx_sig, mult_sig, signif_d0xy_pi_sig, signif_d0xy_k_sig, chi2_dca_sig,z_sig,dr_sig,angle_sig,etajet_sig,pTjet_sig, D0eta_sig, pt_parent_sig;
  
    // Link the variables to the TTree branches
  tree_sig->Branch("d0_pi", &d0_pi_sig, "d0_pi/F");
  tree_sig->Branch("d0_k", &d0_k_sig, "d0_k/F"); 
  tree_sig->Branch("d0xy_pi", &d0xy_pi_sig, "d0xy_pi/F");
  tree_sig->Branch("d0xy_k", &d0xy_k_sig, "d0xy_k/F");
  tree_sig->Branch("sum_d0xy", &sum_d0xy_sig, "sum_d0xy/F");        
  tree_sig->Branch("dca_12", &dca_12_sig, "dca_12/F");
  tree_sig->Branch("dca_D0", &dca_D0_sig, "dca_D0/F");
  tree_sig->Branch("pt_D0", &pt_D0_sig, "pt_D0/F");  
  tree_sig->Branch("y_D0", &y_D0_sig, "y_D0/F");  
  tree_sig->Branch("mass_D0", &mass_D0_sig, "mass_D0/F");              
  tree_sig->Branch("decay_length", &decay_length_sig, "decay_length/F");   
  tree_sig->Branch("costheta", &costheta_sig, "costheta/F"); 
  tree_sig->Branch("costheta_xy", &costhetaxy_sig, "costheta_xy/F"); 
  tree_sig->Branch("sigma_vtx", &sigma_vtx_sig, "sigma_vtx/F"); 
  tree_sig->Branch("mult", &mult_sig, "mult/F"); 
  tree_sig->Branch("signif_d0xy_pi", &signif_d0xy_pi_sig, "signif_d0xy_pi/F");               
  tree_sig->Branch("signif_d0xy_k", &signif_d0xy_k_sig, "signif_d0xy_k/F"); 
  tree_sig->Branch("chi2_dca", &chi2_dca_sig, "chi2_dca/F");
  tree_sig->Branch("xB", &xB_sig, "xB/F");
  tree_sig->Branch("Q2", &Q2_sig, "Q2/F");
  tree_sig->Branch("z", &z_sig, "z/F");
  tree_sig->Branch("dr", &dr_sig, "dr/F");
  tree_sig->Branch("angle", &angle_sig, "angle/F");
  tree_sig->Branch("etajet",&etajet_sig,"etajet/F");
  tree_sig->Branch("pTjet",&pTjet_sig,"pTjet/F");
  tree_sig->Branch("etaD0",&D0eta_sig,"etaD0/F");
  tree_sig->Branch("pt_parent",&pt_parent_sig,"pt_parent/F");

  
  TFile *file_bkg = new TFile(bkg.Data(), "RECREATE");
  TTree *tree_bkg = new TTree("treeMLBkg", "treeMLBkg"); 
  
  // Define variables to store in the Ntuple
  float d0_pi_bkg, d0_k_bkg, d0xy_pi_bkg, d0xy_k_bkg, sum_d0xy_bkg, dca_12_bkg, dca_D0_bkg, decay_length_bkg, xB_bkg, Q2_bkg; 
  float costheta_bkg, costhetaxy_bkg, pt_D0_bkg, y_D0_bkg, mass_D0_bkg, sigma_vtx_bkg, mult_bkg, signif_d0xy_pi_bkg, signif_d0xy_k_bkg, chi2_dca_bkg,z_bkg,dr_bkg,angle_bkg, etajet_bkg, pTjet_bkg, D0eta_bkg,pt_parent_bkg;
  
  // Link the variables to the TTree branches
  tree_bkg->Branch("d0_pi", &d0_pi_bkg, "d0_pi/F");
  tree_bkg->Branch("d0_k", &d0_k_bkg, "d0_k/F");
  tree_bkg->Branch("d0xy_pi", &d0xy_pi_bkg, "d0xy_pi/F");
  tree_bkg->Branch("d0xy_k", &d0xy_k_bkg, "d0xy_k/F");
  tree_bkg->Branch("sum_d0xy", &sum_d0xy_bkg, "sum_d0xy/F");            
  tree_bkg->Branch("dca_12", &dca_12_bkg, "dca_12/F");
  tree_bkg->Branch("dca_D0", &dca_D0_bkg, "dca_D0/F");
  tree_bkg->Branch("pt_D0", &pt_D0_bkg, "pt_D0/F");  
  tree_bkg->Branch("y_D0", &y_D0_bkg, "y_D0/F");  
  tree_bkg->Branch("mass_D0", &mass_D0_bkg, "mass_D0/F");              
  tree_bkg->Branch("decay_length", &decay_length_bkg, "decay_length/F");   
  tree_bkg->Branch("costheta", &costheta_bkg, "costheta/F");  
  tree_bkg->Branch("costheta_xy", &costhetaxy_bkg, "costheta_xy/F"); 
  tree_bkg->Branch("sigma_vtx", &sigma_vtx_bkg, "sigma_vtx/F"); 
  tree_bkg->Branch("mult", &mult_bkg, "mult/F");
  tree_bkg->Branch("signif_d0xy_pi", &signif_d0xy_pi_bkg, "signif_d0xy_pi/F");  
  tree_bkg->Branch("signif_d0xy_k", &signif_d0xy_k_bkg, "signif_d0xy_k/F");  
  tree_bkg->Branch("chi2_dca", &chi2_dca_bkg, "chi2_dca/F");   
  tree_bkg->Branch("xB", &xB_bkg, "xB/F"); 
  tree_bkg->Branch("Q2", &Q2_bkg, "Q2/F"); 
  tree_bkg->Branch("z", &z_bkg, "z/F");
  tree_bkg->Branch("dr", &dr_bkg, "dr/F");
  tree_bkg->Branch("angle", &angle_bkg, "angle/F");
  tree_bkg->Branch("etajet",&etajet_bkg,"etajet/F");
  tree_bkg->Branch("pTjet",&pTjet_bkg,"pTjet/F");
  tree_bkg->Branch("etaD0",&D0eta_bkg,"etaD0/F");
  tree_bkg->Branch("pt_parent",&pt_parent_bkg,"pt_parent/F");


  // Generate the tree with true D0 Properties
  TFile *fout_mcgen = new TFile(signalmc.Data(),"RECREATE");
  TTree *tree_D0 = new TTree("D0Tree","D0 Meson Kinematics");

  float d0_px, d0_py, d0_pz;
  float d0_pt, d0_eta, d0_y, d0_mass;

  tree_D0->Branch("px",&d0_px,"px/F");
  tree_D0->Branch("py",&d0_py,"py/F");
  tree_D0->Branch("pz",&d0_pz,"pz/F");
  tree_D0->Branch("pt",&d0_pt,"pt/F");
  tree_D0->Branch("eta",&d0_eta,"eta/F");
  tree_D0->Branch("rapidity",&d0_y,"rapidity/F");
  tree_D0->Branch("mass",&d0_mass,"mass/F");

  //ROOT file for Gen Jets 
  TFile *file_gen = new TFile(gen.Data(), "RECREATE");
  TTree *tree_gen_sig = new TTree("GenTree_sig","GenTree_sig");//Gen Jet Tree sig
  TTree *tree_gen_bkg = new TTree("GenTree_bkg","GenTree_bkg");
  float pt_gjet_sig,eta_gjet_sig,y_gjet_sig,genZ_sig, gdr_eta_sig, gdr_y_sig,gD0jetAngle_sig;
  float pt_gjet_bkg,eta_gjet_bkg,y_gjet_bkg,genZ_bkg, gdr_eta_bkg, gdr_y_bkg,gD0jetAngle_bkg;
    
  tree_gen_sig->Branch("pt_gjet",&pt_gjet_sig,"pt_gjet/F");
  tree_gen_sig->Branch("eta_gjet",&eta_gjet_sig,"eta_gjet/F");
  tree_gen_sig->Branch("y_gjet",&y_gjet_sig,"y_gjet/F");
  tree_gen_sig->Branch("Z", &genZ_sig,"Z/F");
  tree_gen_sig->Branch("DeltaR_inEtaPhi", &gdr_eta_sig,"DeltaR_inEtaPhi/F");
  tree_gen_sig->Branch("DeltaR_inYPhi", &gdr_y_sig,"DeltaR_inYPhi/F");
  tree_gen_sig->Branch("Angle_bw_D0andJet", &gD0jetAngle_sig,"Angle_bw_D0andJet/F");
  
  tree_gen_bkg->Branch("pt_gjet",&pt_gjet_bkg,"pt_gjet/F");
  tree_gen_bkg->Branch("eta_gjet",&eta_gjet_bkg,"eta_gjet/F");
  tree_gen_bkg->Branch("y_gjet",&y_gjet_bkg,"y_gjet/F");
  tree_gen_bkg->Branch("Z", &genZ_bkg,"Z/F");
  tree_gen_bkg->Branch("DeltaR_inEtaPhi", &gdr_eta_bkg,"DeltaR_inEtaPhi/F");
  tree_gen_bkg->Branch("DeltaR_inYPhi", &gdr_y_bkg,"DeltaR_inYPhi/F");
  tree_gen_bkg->Branch("Angle_bw_D0andJet", &gD0jetAngle_bkg,"Angle_bw_D0andJet/F");
  
//Variables for Jet Clustering	
	  int NEVENTS = 0;
	  int EVETMULTRECO = 0; 
	  int EVETMULTGEN = 0;
	  int ScatteredERecId = 0;
	  int ScatteredEGenId = 0;
	  float EventQ2 = -999;
	  float Eventx = -999;
	  float EventQ2Gen = -999;
	  float EventxGen = -999;

//vectors to store jet clustering outputs
// Vertex
	std::vector<float> Vertex_x;
	std::vector<float> Vertex_y;
	std::vector<float> Vertex_z;
	std::vector<int>   Vertex_ndf;
	std::vector<float> Vertex_chi2;
	std::vector<float> VertexErr_xx;
	std::vector<float> VertexErr_yy;
	std::vector<float> VertexErr_zz;

	// Reco Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> RecoJet_pt;
	std::vector<float> RecoJet_eta;
	std::vector<float> RecoJet_rapidity;
	std::vector<float> RecoJet_phi;
	std::vector<float> RecoJet_E;
	std::vector<float> RecoJet_M;
	std::vector<bool>  RecoJet_hasElectron;
	std::vector<float> RecoJet_maxPtPart_pt;
        std::vector<float> RecoJet_rapi;
	// Reco jet constituents
	std::vector<std::vector<float>> RecoJet_constituent_pt;
	std::vector<std::vector<float>> RecoJet_constituent_eta;
	std::vector<std::vector<float>> RecoJet_constituent_phi;
	std::vector<std::vector<int>> RecoJet_constituent_nhits;
	std::vector<std::vector<int>> RecoJet_constituent_pdgid;
	std::vector<std::vector<int>> RecoJet_constituent_pdgidTruth;
	std::vector<std::vector<int>> RecoJet_constituent_idx;
        std::vector<std::vector<float>> RecoJet_constituent_energy;
        std::vector<std::vector<float>> RecoJet_constituent_rapi;
	
	// Gen Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> GenJet_pt;
	std::vector<float> GenJet_eta;
	std::vector<float> GenJet_phi;
	std::vector<float> GenJet_E;
	std::vector<float> GenJet_M;
	std::vector<bool> GenJet_hasElectron;
	std::vector<bool> GenJet_hasNeutral;
	std::vector<float> GenJet_maxPtPart_pt;
        std::vector<float> GenJet_rapi;
	
	// Gen jet constituents
	std::vector<std::vector<float>> GenJet_constituent_pt;
	std::vector<std::vector<float>> GenJet_constituent_eta;
	std::vector<std::vector<float>> GenJet_constituent_phi;
	std::vector<std::vector<int>> GenJet_constituent_pdgid;
	std::vector<std::vector<float>> GenJet_constituent_rapi;
	std::vector<std::vector<int>> GenJet_constituent_idx;
	std::vector<std::vector<float>> GenJet_constituent_energy;
	
  int nevents = 0;
  int mult_charged = 0;
  while(treereader.Next())
    {
      if(nevents%1000==0) printf("\n[i] New event %d\n",nevents);
      float Q2_mc = -1000.;
      float xB_mc = -1000.;
      if (incQ2.GetSize()>0){
      Q2_mc = incQ2[0];
      xB_mc = incxB[0];
      }
      // find MC primary vertex
      int nMCPart = mcPartMass.GetSize();

      TVector3 vertex_mc(-999., -999., -999.);
      for(int imc=0; imc<nMCPart; imc++)
	{
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
      TVector3 err_vertex_rc(-999., -999., -999.);
      if(prim_vtx_index.GetSize()>0)
	{
	  int rc_vtx_index = prim_vtx_index[0];
	  vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
	  err_vertex_rc.SetXYZ(sqrt(CTVerr_xx[rc_vtx_index]), sqrt(CTVerr_yy[rc_vtx_index]), sqrt(CTVerr_zz[rc_vtx_index]));
	}
	hPullVtxX->Fill((vertex_rc.x()-vertex_mc.x())/err_vertex_rc.x()); 
	hPullVtxY->Fill((vertex_rc.y()-vertex_mc.y())/err_vertex_rc.y()); 
	hPullVtxZ->Fill((vertex_rc.z()-vertex_mc.z())/err_vertex_rc.z()); 	

      // map MC and RC particles
      int nAssoc = assocChRecID.GetSize();
      map<int, int> assoc_map_to_rc;
      map<int, int> assoc_map_to_mc;
      
      for(unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{
	  // reco level eta using Real PID
	  TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
	  hreco_eta->Fill(mom.Eta());

	  if(abs(rcPdg[rc_index]) == 11) hreco_eta_e->Fill(mom.Eta());
	  if(abs(rcPdg[rc_index]) == 211) hreco_eta_pi->Fill(mom.Eta());
	  if(abs(rcPdg[rc_index]) == 321 ) hreco_eta_k->Fill(mom.Eta());
	  if(abs(rcPdg[rc_index]) == 2212) hreco_eta_p->Fill(mom.Eta());
	  
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
	  assoc_map_to_rc[matched_mc_index] = rc_index;
	  assoc_map_to_mc[rc_index] = matched_mc_index;
	}

      // Loop over primary particles
      int nMcPart = 0;
      for(int imc=0; imc<nMCPart; imc++)
	{
	  // mc level eta using Truth PID
	  TVector3 mom(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);
          hmc_eta->Fill(mom.Eta());

          if(abs(mcPartPdg[imc]) == 11) hmc_eta_e->Fill(mom.Eta());
          if(abs(mcPartPdg[imc]) == 211) hmc_eta_pi->Fill(mom.Eta());
          if(abs(mcPartPdg[imc]) == 321 ) hmc_eta_k->Fill(mom.Eta());
          if(abs(mcPartPdg[imc]) == 2212) hmc_eta_p->Fill(mom.Eta());

	  
	  if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0)
	    {
	      double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
	      if(dist < 1e-4)
		{
		  // count charged particles within |eta| < 3.5
		  TVector3 mc_mom(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);
		  double mcEta = mc_mom.PseudoRapidity();
		  if(fabs(mcEta) < 3.5) nMcPart++;
		}
	    }
	}
      mult_charged = nMcPart;
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
		  if(assoc_map_to_rc.find(imc) != assoc_map_to_rc.end()) rc_index = assoc_map_to_rc[imc];

		  if(rc_index>=0)
		    {
		      TVector3 dcaToVtx = getDcaToVtx(rc_index, vertex_rc);
		      
		      int ip = -1;
		      if(fabs(mcPartPdg[imc]) == 211) ip = 0;
		      if(fabs(mcPartPdg[imc]) == 321) ip = 1;
		      if(fabs(mcPartPdg[imc]) == 2212) ip = 2;
		      if(ip>=0)
			{
			  TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
			  if(ip<2)
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
      
      // Filling some values: 
      // For Vertex
      Vertex_x.clear();
      Vertex_y.clear();
      Vertex_z.clear();
      Vertex_ndf.clear();
      Vertex_chi2.clear();
      VertexErr_xx.clear();
      VertexErr_yy.clear();
      VertexErr_zz.clear();
      
      for (unsigned int ivtx = 0; ivtx < CTVx.GetSize(); ++ivtx) {
	Vertex_x.push_back(CTVx[ivtx]);
	Vertex_y.push_back(CTVy[ivtx]);
	Vertex_z.push_back(CTVz[ivtx]);
	Vertex_ndf.push_back(CTVndf[ivtx]);
	Vertex_chi2.push_back(CTVchi2[ivtx]);
	VertexErr_xx.push_back(CTVerr_xx[ivtx]);
	VertexErr_yy.push_back(CTVerr_yy[ivtx]);
	VertexErr_zz.push_back(CTVerr_zz[ivtx]);        
        }
        
      //  For Scattered electron
      
      if(ScatElecRecoId.GetSize() > 0){
	int iscatR = ScatElecRecoId[0];
	ScatteredERecId = iscatR;    
      }else{ScatteredERecId = -999;}
      
      
      if(ScatElecGenId.GetSize() > 0){
	int iscatS = ScatElecGenId[0];
	ScatteredEGenId = iscatS;    
      }else{ScatteredEGenId = -999;}
		
      // other reco and gen quantities
      if(EvtQ2.GetSize() > 0){
	float q2valueR = EvtQ2[0];
	EventQ2 = q2valueR;    
      }else{EventQ2 = -999;}
      
      if(Evtx.GetSize() > 0){
	float xvalueR = Evtx[0];
	Eventx = xvalueR;    
      }else{Eventx = -999;}
      
      if(EvtQ2Gen.GetSize() > 0){
	float q2valueG = EvtQ2Gen[0];
	EventQ2Gen = q2valueG;    
      }else{EventQ2Gen = -999;}
		
      if(EvtxGen.GetSize() > 0){
	float xvalueG = EvtxGen[0];
	EventxGen = xvalueG;    
      }else{EventxGen = -999;}
      
      
      // look for D0
      bool hasD0 = false;
      vector<int> mc_index_D0_pi;
      vector<int> mc_index_D0_k;
      mc_index_D0_pi.clear();
      mc_index_D0_k.clear();
      
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(fabs(mcPartPdg[imc]) != 421) continue;
	  hEventStat->Fill(1.5);
   //----------------------------------- Fill D0 gen tree with properties-----------------------------------------
    TLorentzVector d0_vec;

    d0_px = mcMomPx[imc];
    d0_py = mcMomPy[imc];
    d0_pz = mcMomPz[imc];
    d0_mass = mcPartMass[imc];

    d0_vec.SetXYZM(d0_px, d0_py, d0_pz, d0_mass);

    d0_pt  = d0_vec.Pt();
    d0_eta = d0_vec.Eta();
    d0_y   = d0_vec.Rapidity();
    tree_D0->Fill();
    //--------------------------------------------------------------------------------------------------------------------
	  
	  int nDuaghters = mcPartDaughter_end[imc]-mcPartDaughter_begin[imc];
	  if(nDuaghters!=2) continue;

	  // find D0 that decay into pi+K
	  bool is_pik_decay = false;	  
	  int daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
	  int daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
	  int daug_pdg_1 = mcPartPdg[daug_index_1];
	  int daug_pdg_2 = mcPartPdg[daug_index_2];
	  if( (fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==211) || (fabs(daug_pdg_1)==211 && fabs(daug_pdg_2)==321) )
	    {
	      is_pik_decay = true;
	      
	    }
	  if(!is_pik_decay) continue;
	  if(fabs(daug_pdg_1)==211)
	    {
	      mc_index_D0_pi.push_back(daug_index_1);
	      mc_index_D0_k.push_back(daug_index_2);
	    }
	  else
	    {
	      mc_index_D0_pi.push_back(daug_index_2);
	      mc_index_D0_k.push_back(daug_index_1);
	    }
	  hasD0 = true;
	  hEventStat->Fill(2.5);

	  // D0 kinematics
	  TLorentzVector mc_mom_vec;
	  mc_mom_vec.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
	  double mcRap = mc_mom_vec.Rapidity();
	  double mcPt = mc_mom_vec.Pt();
	  hMCD0PtRap->Fill(mcRap, mcPt);

	  // decay dauther kinematics
	  for(int ip = 0; ip<2; ip++)
	    {
	      int mc_part_index;
	      if(ip==0) mc_part_index = mc_index_D0_pi[mc_index_D0_pi.size()-1];
	      if(ip==1) mc_part_index = mc_index_D0_k[mc_index_D0_k.size()-1];
	      
	      TLorentzVector mc_part_vec;
	      mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
	      if(ip==0) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
	      if(ip==1) hMcKPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
		  
	      int rc_part_index = -1;
	      if(assoc_map_to_rc.find(mc_part_index) != assoc_map_to_rc.end()) rc_part_index = assoc_map_to_rc[mc_part_index];
	      if(rc_part_index>=0)
		{
		  TVector3 dcaToVtx = getDcaToVtx(rc_part_index, vertex_rc);
		  
		  TVector3 mom(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
		  hRcSecPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
		  hRcSecPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
		  
		  //printf("Sec %d: (%2.4f, %2.4f, %2.4f), mcStartPoint = (%2.4f, %2.4f, %2.4f)\n", rc_part_index, pos.x(), pos.y(), pos.z(), mcPartVx[mc_part_index], mcPartVy[mc_part_index], mcPartVz[mc_part_index]);
		}
	    }
	}

      // Get reconstructed pions and kaons
      hNRecoVtx->Fill(CTVx.GetSize());
      const int pid_mode = 1; // 0 - truth; 1 - realistic
      vector<unsigned int> pi_index;
      vector<unsigned int> k_index;
      pi_index.clear();
      k_index.clear();
      for(unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{	  
	  if(pid_mode==0)
	    {
	      int iSimPartID = -1;
	      if(assoc_map_to_mc.find(rc_index) != assoc_map_to_mc.end()) iSimPartID = assoc_map_to_mc[rc_index];
	      if(iSimPartID>=0)
		{
		  if(fabs(mcPartPdg[iSimPartID]) == 211) pi_index.push_back(rc_index);
		  if(fabs(mcPartPdg[iSimPartID]) == 321) k_index.push_back(rc_index);
		}
	    }
	  else if(pid_mode==1)
	    {
	      if(fabs(rcPdg[rc_index]) == 211) pi_index.push_back(rc_index);
	      if(fabs(rcPdg[rc_index]) == 321) k_index.push_back(rc_index);
	    }
	}

      // pair pion and kaon
      for(unsigned int i=0; i<pi_index.size(); i++)
	{

	  //Check for pion using Truth level info
	  /*int mcIdx_pii = -1;
	  if(assoc_map_to_mc.find(pi_index[i]) != assoc_map_to_mc.end()) mcIdx_pii = assoc_map_to_mc[pi_index[i]];
	  
	  if(mcIdx_pii >=0 ){ 
	    if(abs(mcPartPdg[mcIdx_pii]) != 211 ) continue ;
	  }
	  */
	  bool not_pi = false;
	  TVector3 dcaToVtx = getDcaToVtx(pi_index[i], vertex_rc);
	  std::array<float, 21>& cov_pion = rcTrkCov->At(pi_index[i]);
	  int q_pion = rcCharge[pi_index[i]];
	   
	  for(unsigned int j=0; j<k_index.size(); j++)
	    {
	      bool not_k = false;
	      //Check for kaon using Truth level info
	      /*int mcIdx_kk = -1;
	      if(assoc_map_to_mc.find(k_index[i]) != assoc_map_to_mc.end()) mcIdx_kk = assoc_map_to_mc[k_index[i]];

	      if(mcIdx_kk>=0){
		if(abs(mcPartPdg[mcIdx_kk]) != 321) continue ;
	      }
	      */
	      TVector3 dcaToVtx2 = getDcaToVtx(k_index[j], vertex_rc);
	      std::array<float, 21>& cov_kaon = rcTrkCov->At(k_index[j]); 
	      int q_kaon = rcCharge[k_index[j]];
	      
	      if(rcCharge[pi_index[i]]*rcCharge[k_index[j]]<0)
		{
		  //printf("[i] Check pair (%d, %d)\n", pi_index[i], k_index[j]);
		  // -- only look at unlike-sign pi+k pair
		  bool is_D0_pik = false;
		  int mc_index_pi = -1, mc_index_k = -1;
		  if(assoc_map_to_mc.find(pi_index[i]) != assoc_map_to_mc.end()) mc_index_pi = assoc_map_to_mc[pi_index[i]];
		  if(assoc_map_to_mc.find(k_index[j])  != assoc_map_to_mc.end()) mc_index_k  = assoc_map_to_mc[k_index[j]];

		  for(unsigned int k=0; k<mc_index_D0_pi.size(); k++)
		    {
		      if(mc_index_pi==mc_index_D0_pi[k] && mc_index_k==mc_index_D0_k[k])
			{
			  is_D0_pik = true;
			  break;
			}
      		    }

		  float dcaDaughters, cosTheta, decayLength, V0DcaToVtx, cosTheta_xy, sigma_vtx;
		  TVector3 decayVertex, decayVertex_ana; 
		  double chi2_ndf;
		  double err_Par[5];  // or whatever size is appropriate
		  //double* err_Par = errParArray;
		  TLorentzVector parent = getPairParent(pi_index[i], k_index[j], vertex_rc, dcaDaughters, cosTheta, cosTheta_xy, decayLength, V0DcaToVtx, sigma_vtx,decayVertex,decayVertex_ana,chi2_ndf, err_Par);
		  hchi2_vtx->Fill(chi2_ndf);

		  // Storing rc index of Pion and Kaon
         	  int pii_index = pi_index[i]; // rc index of pion
		  int kk_index =  k_index[j];  // rc index of kaon

	          //Signal D0
		  if(is_D0_pik) 
		    {
		   TVector3 MCVertex_Kaon(mcPartVx[mc_index_k], mcPartVy[mc_index_k], mcPartVz[mc_index_k]);
		   TVector3 MCVertex_Pion(mcPartVx[mc_index_pi], mcPartVy[mc_index_pi], mcPartVz[mc_index_pi]);
		   
		  // printf("Signal MC Vertex Kaon = (%f, %f, %f)\n",MCVertex_Kaon.X(), MCVertex_Kaon.Y(), MCVertex_Kaon.Z());	 	
		  // printf("Signal MC Vertex Pion = (%f, %f, %f)\n",MCVertex_Pion.X(), MCVertex_Pion.Y(), MCVertex_Pion.Z());
		      hRes_SVx_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.X()-MCVertex_Kaon.X()*0.1)*10);
		      hRes_SVy_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.Y()-MCVertex_Kaon.Y()*0.1)*10);
		      hRes_SVz_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.Z()-MCVertex_Kaon.Z()*0.1)*10);	
		      
		      hRes_SVx_Helixana->Fill(parent.Pt(), parent.Rapidity(), (decayVertex_ana.X()-MCVertex_Kaon.X()*0.1)*10);
		      hRes_SVy_Helixana->Fill(parent.Pt(), parent.Rapidity(), (decayVertex_ana.Y()-MCVertex_Kaon.Y()*0.1)*10);
		      hRes_SVz_Helixana->Fill(parent.Pt(), parent.Rapidity(), (decayVertex_ana.Z()-MCVertex_Kaon.Z()*0.1)*10);
		      
		      hRes_SVxy_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.Perp()-MCVertex_Kaon.Perp()*0.1)*10);
		      hRes_SVxy_Helixana->Fill(parent.Pt(), parent.Rapidity(), (decayVertex_ana.Perp()-MCVertex_Kaon.Perp()*0.1)*10);
		      
		      hRes_SVx_Helixfit_pull->Fill(parent.Pt(), parent.Rapidity(), ((decayVertex.X()-MCVertex_Kaon.X()*0.1))/(err_Par[0]));
		      hRes_SVy_Helixfit_pull->Fill(parent.Pt(), parent.Rapidity(), ((decayVertex.Y()-MCVertex_Kaon.Y()*0.1))/(err_Par[1]));
		      hRes_SVz_Helixfit_pull->Fill(parent.Pt(), parent.Rapidity(), ((decayVertex.Z()-MCVertex_Kaon.Z()*0.1))/(err_Par[2]));	
		      
		      hchi2_vtx_sig->Fill(chi2_ndf);		      
		            		   
		      hEventStat->Fill(3.5);
		      if (q_kaon == -1 && q_pion == 1)
		      hEventStat->Fill(4.5);   // D0
		      else if (q_kaon == 1 && q_pion == -1)
		      hEventStat->Fill(5.5);   // D0bar
		      h3PairDca12[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters);
		      h3PairCosTheta[0]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
		      h3PairDca[0]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
		      h3PairDecayLength[0]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
		      //printf("Signal: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
		      h3InvMass[0][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
               
           // Toplogical Variables for Signal
		      d0_pi_sig = dcaToVtx.Mag();
		      d0_k_sig = dcaToVtx2.Mag();
		      d0xy_pi_sig = dcaToVtx.Perp();
		      signif_d0xy_pi_sig = d0xy_pi_sig/sqrt(cov_pion[0]);
		      d0xy_k_sig = dcaToVtx2.Perp();
		      signif_d0xy_k_sig = d0xy_k_sig/sqrt(cov_kaon[0]);
		      sum_d0xy_sig = sqrt(d0xy_pi_sig*d0xy_pi_sig+d0xy_k_sig*d0xy_k_sig);		      
		      dca_12_sig = dcaDaughters;
		      dca_D0_sig = V0DcaToVtx;
		      decay_length_sig = decayLength;
		      costheta_sig = cosTheta;
		      costhetaxy_sig = cosTheta_xy;
		      pt_D0_sig = parent.Pt();
		      y_D0_sig = parent.Rapidity();
		      mass_D0_sig = parent.M();
		      sigma_vtx_sig = sigma_vtx;
		      mult_sig = mult_charged;
		      chi2_dca_sig = chi2_ndf;
		      xB_sig = xB_mc;
		      Q2_sig = Q2_mc;        
		        
//Jet Clustering for the signal D0 (signal pik pair)   
// Build pseudojets (using reco-particles)

		      std::vector<fastjet::PseudoJet> particles_reco; // Reco PseudoJet type Vector
		      int c1 =0;                                                           
		      
		      for (unsigned int i = 0; i < rcMomPx.GetSize(); ++i) { 
			
			bool is_PiorK = false;
			if(i == pii_index || i == kk_index){  //looking for pi and k daughters of signal D0  
			  is_PiorK = true;
			}
 
			if(is_PiorK) continue;
 
			TVector3 mom(rcMomPx[i], rcMomPy[i],rcMomPz[i]);
			if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
			if ( TrkRecoNhits[i] < nhitcut ) continue;
			
			if ( removeelectrons == 1 ){
			  // Find electron
			  int chargePartIndex = i; 
			  int elecIndex = -1;
			  float elecIndexWeight = -1.0;
			  for(unsigned int itrkass = 0; itrkass < assocWeight.GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
			    
			    if( assocWeight.GetSize() > 0 ){
			      if(TrkPartAssocRec[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
				if(assocWeight[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
				  elecIndex = TrkPartAssocSim[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
				  elecIndexWeight = assocWeight[itrkass];
				}
			      }
			    }
			  }

			  if(mcPartPdg[elecIndex] == 11) continue;
			}
			
			if ( removeelectrons == 2 && i == ScatteredERecId ) continue; 

			fastjet::PseudoJet p(rcMomPx[i], rcMomPy[i], rcMomPz[i], TrkRecoE[i]);
			p.set_user_index(i);
			particles_reco.push_back(p);
		      }

		      //----------Introducing D0 parent into recotype PseudoJets---------------------	
		      float dcaDaughters, cosTheta, decayLength, V0DcaToVtx, cosTheta_xy, sigma_vtx;
		      TVector3 decayVertex, decayVertex_ana; 
		      double chi2_ndf;
		      double err_Par[5];  // or whatever size is appropriate
		      //double* err_Par = errParArray;
		      TLorentzVector parent = getPairParent(pii_index,kk_index, vertex_rc, dcaDaughters, cosTheta, cosTheta_xy, decayLength, V0DcaToVtx, sigma_vtx,decayVertex,decayVertex_ana,chi2_ndf, err_Par);
		      pt_parent_sig = parent.Pt();
  
		      if ( parent.Pt() < minCstPt || parent.Pt() > maxCstPt ){ // if D0 parent don't fit in Jet Constituent Pt limits => introduce pi-k daughters  to pseudojets                                                //we should not do this, it would give jets with no D0!!!
    
			for (int l=0;l<2;l++){ 
			  int idx;
			  if (l==0) idx = pii_index;
			  if(l==1)  idx = kk_index;
			
			
			  TVector3 mom(rcMomPx[idx], rcMomPy[idx],rcMomPz[idx]);
			  if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
			  if ( TrkRecoNhits[i] < nhitcut ) continue;
			  
			  fastjet::PseudoJet p(rcMomPx[idx], rcMomPy[idx], rcMomPz[idx], TrkRecoE[idx]);
			  p.set_user_index(idx);
			  particles_reco.push_back(p);
			  
			} 
		      }
		      
		      if ( parent.Pt() > minCstPt && parent.Pt() < maxCstPt ){
			fastjet::PseudoJet p(parent.Px(), parent.Py(),parent.Pz(), parent.E());  
			p.set_user_index(-99999);                                 // setting D0 particle index = negative integer (-99999) 
			particles_reco.push_back(p);
		      }
		      
	  // Define algorithm
	fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
	fastjet::RecombinationScheme scheme = fastjet::E_scheme;
	  // Jet definition
	fastjet::JetDefinition jet_def(algo, R_value, scheme);
	fastjet::GhostedAreaSpec ghost_spec(ghostMaxRap, numGhostRepeat, ghostArea);
	fastjet::AreaType atype = fastjet::active_area;
	fastjet::AreaDefinition area_def(atype, ghost_spec);
                        
	  // Clear all vectors for the Reco Jet Clustering for new Signal D0 
	  RecoJet_pt.clear();
	  RecoJet_eta.clear();
	  RecoJet_rapidity.clear();
	  RecoJet_phi.clear();
	  RecoJet_E.clear();
	  RecoJet_M.clear();
	  RecoJet_rapi.clear();
	  RecoJet_hasElectron.clear();
	  RecoJet_maxPtPart_pt.clear();
	  RecoJet_constituent_pt.clear(); 
	  RecoJet_constituent_eta.clear();
	  RecoJet_constituent_phi.clear(); 
	  RecoJet_constituent_nhits.clear();
	  RecoJet_constituent_pdgid.clear();
	  RecoJet_constituent_pdgidTruth.clear();
	  RecoJet_constituent_idx.clear();
          RecoJet_constituent_energy.clear();
	  RecoJet_constituent_rapi.clear();
	  
	  EVETMULTRECO = (int)rcMomPx.GetSize();		

	  //Reco Jet clustering(Sig)
	  fastjet::ClusterSequenceArea cs_reco(particles_reco, jet_def, area_def);
	  std::vector<fastjet::PseudoJet> jets_reco = fastjet::sorted_by_pt(cs_reco.inclusive_jets(minJetPt));
	  for (auto &jet : jets_reco) {
	    RecoJet_pt.push_back(jet.pt());
	    RecoJet_eta.push_back(jet.eta());
	    RecoJet_rapidity.push_back(jet.rapidity());
	    RecoJet_phi.push_back(jet.phi_std());
	    RecoJet_E.push_back(jet.e());
	    RecoJet_M.push_back(jet.m());
	    RecoJet_rapi.push_back(jet.rap());
		

	    bool hasElectron = false;
	    float maxPtReco = -1.0;
	    std::vector<float> cpt, ceta, cphi, cenergy,crapi;
	    std::vector<int> chits, cpdgid, cpdgidtruth,cindex;
	    cpt.clear(); ceta.clear(); cphi.clear(); chits.clear(); cpdgid.clear(); cpdgidtruth.clear(), cindex.clear(), cenergy.clear(),crapi.clear();
	    
	    fastjet::PseudoJet D0_parent;

	    for (auto &c : jet.constituents()) {
	      int idx = c.user_index();
	      
        //Look for the Signal D0 constituent of the Jet 
	      bool is_D0_parent = false;
	      
	      if(idx == -99999){ //check for the signal D0
	      
		for(const auto& p : particles_reco ){ // loop over all the jet particles
		
		  if (p.user_index() == idx) {
		    D0_parent = p;
		    is_D0_parent = true;
		    break;
		  }
		}
                      
		cpt.push_back(D0_parent.pt());
		ceta.push_back(D0_parent.eta());
		cphi.push_back(D0_parent.phi());
		cenergy.push_back(D0_parent.E());
		chits.push_back(0);           // 0 hits for neutral D0 
		cpdgid.push_back(421);        // PID = 421 for D0
		cpdgidtruth.push_back(421);
		cindex.push_back(idx);
		crapi.push_back(D0_parent.rap());
		                                                          // Here, hasElectron = false
		if (D0_parent.pt() > maxPtReco) maxPtReco = D0_parent.pt();   
	      }
	      if (is_D0_parent) continue; // skip the rest of loop for signal D0 particle jet
		    
	      TVector3 v3(rcMomPx[idx], rcMomPy[idx], rcMomPz[idx]);
	      float rcEnergy = TrkRecoE[idx];
	      TLorentzVector lv(rcMomPx[idx],  rcMomPy[idx], rcMomPz[idx], rcEnergy);
	      
	      cpt.push_back(v3.Pt());
	      ceta.push_back(v3.Eta());
	      cphi.push_back(v3.Phi());
	      chits.push_back(TrkRecoNhits[idx]);
	      cindex.push_back(idx);
	      cenergy.push_back(rcEnergy);
	      crapi.push_back(lv.Rapidity());
	      if (v3.Pt() > maxPtReco) maxPtReco = v3.Pt();
	      
	      // Find electron
	      int chargePartIndex = idx; 
	      int elecIndex = -1;
	      float elecIndexWeight = -1.0;
	      for(unsigned int itrkass = 0; itrkass < assocWeight.GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations

		if( assocWeight.GetSize() > 0 ){
		  if(TrkPartAssocRec[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
		    if(assocWeight[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
		      elecIndex = TrkPartAssocSim[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
		      elecIndexWeight = assocWeight[itrkass];
		    }
		  }
		}
	      }

	      if(mcPartPdg[elecIndex] == 11){
		hasElectron = true;
	      }
		
	      cpdgid.push_back(rcPdg[idx]);
	      cpdgidtruth.push_back(mcPartPdg[elecIndex]);
	    }
            
	    RecoJet_constituent_pt.push_back(cpt);  
	    RecoJet_constituent_eta.push_back(ceta);
	    RecoJet_constituent_phi.push_back(cphi);
	    RecoJet_constituent_nhits.push_back(chits);
	    RecoJet_constituent_pdgid.push_back(cpdgid);
	    RecoJet_constituent_pdgidTruth.push_back(cpdgidtruth);
	    RecoJet_hasElectron.push_back(hasElectron);
	    RecoJet_maxPtPart_pt.push_back(maxPtReco);
	    RecoJet_constituent_idx.push_back(cindex);
            RecoJet_constituent_energy.push_back(cenergy);
            RecoJet_constituent_rapi.push_back(crapi);
	    
	  } 

// Reading Reco Jet(Signal):

	for(unsigned int i=0; i<RecoJet_pt.size(); i++){  //loop over reco jets  
	  
	  int  jetParticles = RecoJet_constituent_idx[i].size();  // Total no. of constituents of ith jet 
     
	  float jetPt  = RecoJet_pt[i];
	  float jetEta = RecoJet_eta[i];
	  float jetRapidity = RecoJet_rapi[i];
	  float jetPhi = RecoJet_phi[i];
	  float jetE = RecoJet_E[i];
     
	  TLorentzVector LvJet;
	  LvJet.SetPtEtaPhiE(jetPt,jetEta, jetPhi, jetE); // TLorentzVector of Jet 
	  TVector3 jetMom = LvJet.Vect(); // Jet Momentum Vector
       	  bool c2= false;
     
	  for(unsigned int j=0; j<RecoJet_constituent_idx[i].size(); j++){ //Loop over all constituent particles of the ith Jet

	    int particle_idx = RecoJet_constituent_idx[i][j]; // index of jth particle of ith jet
	    bool is_D0 = false;
            
	    //check if constituent particle is D0
	    if(particle_idx == -99999){        
	      is_D0 = true;
	      c2 = true;
	    }
	 
	    if(is_D0){  // if Jet constituent is the signal D0

	      float D0Pt  = RecoJet_constituent_pt[i][j];
	      float D0Eta = RecoJet_constituent_eta[i][j];
	      float D0Phi = RecoJet_constituent_phi[i][j];
	      float D0E = RecoJet_constituent_energy[i][j];
		  
	      TLorentzVector lvD0; 
	      lvD0.SetPtEtaPhiE(D0Pt, D0Eta, D0Phi, D0E);
	      
	      TVector3 D0_Mom;
	      D0_Mom.SetPtEtaPhi(D0Pt, D0Eta, D0Phi);
              float D0Mass = lvD0.M();
	      float D0Rapidity = lvD0.Rapidity();
	
	      float sig_z = float(jetMom.Dot(D0_Mom))/float(jetMom.Dot(jetMom)) ;  // Fragmentation Variable z for the sinal D0
              z_sig =sig_z;
              h3sig_z->Fill(sig_z,D0Rapidity,D0Mass);
              
                     
	      float dPhi = jetPhi - D0Phi;
              dPhi = TVector2::Phi_mpi_pi(dPhi);
              float dEta = jetEta - D0Eta;
              float deltaR = TMath::Sqrt(dEta*dEta + dPhi*dPhi); // sig_D0 distance from JetAxis in eta-Phi space
              
              float angleRad = jetMom.Angle(D0_Mom); // angle bw sig_D0 momentum vector and jetMomentum vector
              float angleDeg = angleRad * 180.0 / TMath::Pi(); 
	                    
              angle_sig	= angleDeg;
              dr_sig = deltaR;
              etajet_sig = jetEta;
	      pTjet_sig = jetPt;
	      D0eta_sig = D0Eta;
	      
	      break;//leaving the JetConstituent loop because we found the only signal D0 present
	   } 	 
	}
	  if(c2) break; // leaving the Jet loop because we found the only signal D0 present
      }
                  
      tree_sig->Fill();
       
		    }

		  //Background:   
		  else
		    {

		      hEventStat->Fill(6.5);   
		      hchi2_vtx_bkg->Fill(chi2_ndf);		            
		      h3PairDca12[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters);
		      h3PairCosTheta[1]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
		      h3PairDca[1]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
		      h3PairDecayLength[1]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
			  
		      //printf("Bkg: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
		      h3InvMass[1][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());

		      // Toplogical Variables for Bkg
		      d0_pi_bkg = dcaToVtx.Mag();
		      d0_k_bkg = dcaToVtx2.Mag();
		      d0xy_pi_bkg = dcaToVtx.Perp();
		      signif_d0xy_pi_bkg = d0xy_pi_bkg/sqrt(cov_pion[0]);
		      d0xy_k_bkg = dcaToVtx2.Perp();
		      signif_d0xy_k_bkg = d0xy_k_bkg/sqrt(cov_kaon[0]);
		      sum_d0xy_bkg = sqrt(d0xy_pi_bkg*d0xy_pi_bkg+d0xy_k_bkg*d0xy_k_bkg);		      
		      dca_12_bkg = dcaDaughters;
		      dca_D0_bkg = V0DcaToVtx;
		      decay_length_bkg = decayLength;
		      costheta_bkg = cosTheta;
		      costhetaxy_bkg = cosTheta_xy;
		      pt_D0_bkg = parent.Pt();
		      y_D0_bkg = parent.Rapidity();
		      mass_D0_bkg = parent.M();
		      sigma_vtx_bkg = sigma_vtx;
		      mult_bkg = mult_charged;
                      chi2_dca_bkg = chi2_ndf;
    		      xB_bkg = xB_mc;          
    		      Q2_bkg = Q2_mc; 
		      
//Jet Clustering for the Bkg D0 (Bkg pi-k pair)      
//Build pseudojets (recotype)


	std::vector<fastjet::PseudoJet> particles_reco;                                                      
        
        for (unsigned int i = 0; i < rcMomPx.GetSize(); ++i) { 
                                                                             
	  
	TVector3 mom(rcMomPx[i], rcMomPy[i],rcMomPz[i]);
	if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
	if ( TrkRecoNhits[i] < nhitcut ) continue;
	
	if ( removeelectrons == 1 ){ // ******
	  // Find electron
	  int chargePartIndex = i; 
	  int elecIndex = -1;
	  float elecIndexWeight = -1.0;
	  for(unsigned int itrkass = 0; itrkass < assocWeight.GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations

	    if( assocWeight.GetSize() > 0 ){
	    if(TrkPartAssocRec[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
	    if(assocWeight[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
		  elecIndex = TrkPartAssocSim[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
		  elecIndexWeight = assocWeight[itrkass];
		}
	      }
	    }
	  }
	  if(mcPartPdg[elecIndex] == 11){

	    if( i == pii_index || i == kk_index ){
	      if(i == pii_index) not_pi = true;
	      if(i == kk_index) not_k = true;
	      break;
	    }
	    continue;
	  }
	}
	
	if ( removeelectrons == 2 && i == ScatteredERecId ) continue; // *****

	fastjet::PseudoJet p(rcMomPx[i], rcMomPy[i], rcMomPz[i], TrkRecoE[i]);
	p.set_user_index(i);
	particles_reco.push_back(p);
        }

	if(not_k) continue;  // skip to next rc k
	if(not_pi) goto next_pi; // get out of k loop and move to next pi
	   
	//Introducing Signal D0 parent into recotype PseudoJets	
	float dcaDaughters, cosTheta, decayLength, V0DcaToVtx, cosTheta_xy, sigma_vtx;
	TVector3 decayVertex, decayVertex_ana; 
	double chi2_ndf;
	double err_Par[5];  // or whatever size is appropriate
	    //double* err_Par = errParArray;
        TLorentzVector parent = getPairParent(pii_index,kk_index, vertex_rc, dcaDaughters, cosTheta, cosTheta_xy, decayLength, V0DcaToVtx, sigma_vtx,decayVertex,decayVertex_ana,chi2_ndf, err_Par);
	pt_parent_bkg = parent.Pt();

//if bkg D0 Pt < minCstPt or > maxCstPt, then instead of D0 parent, introducing its pi-k daughters to recotype pseudoJet	
	if ( parent.Pt() < minCstPt || parent.Pt() > maxCstPt ){ //we should not do this as this will may give jet with no D0
        
	  for (int l=0;l<2;l++){ 
	    int idx;
            
            if (l==0){
              idx = pii_index;
            }
	  
	  if(l==1){
	    idx = kk_index;
	  }
	  
        TVector3 mom(rcMomPx[idx], rcMomPy[idx],rcMomPz[idx]);
        if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
        if ( TrkRecoNhits[i] < nhitcut ) continue;

        fastjet::PseudoJet p(rcMomPx[idx], rcMomPy[idx], rcMomPz[idx], TrkRecoE[idx]);
        p.set_user_index(idx);
        particles_reco.push_back(p);
	
	} 
	}

//If bkg D0 Pt is b/w minCstPt and maxCstPt
	if ( parent.Pt() > minCstPt && parent.Pt() < maxCstPt ){ 
	    fastjet::PseudoJet p(parent.Px(), parent.Py(),parent.Pz(), parent.E());  
	    p.set_user_index(-99999);// setting signal D0 index = negative integer (-99999) 
	    particles_reco.push_back(p);
        }
		      
	  // Define algorithm
	fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
	fastjet::RecombinationScheme scheme = fastjet::E_scheme;
	  // Jet definition
	fastjet::JetDefinition jet_def(algo, R_value, scheme);
	fastjet::GhostedAreaSpec ghost_spec(ghostMaxRap, numGhostRepeat, ghostArea);
	fastjet::AreaType atype = fastjet::active_area;
	fastjet::AreaDefinition area_def(atype, ghost_spec);
                        
	  // Clear all vectors for the new bkg D0
	  RecoJet_pt.clear();
	  RecoJet_eta.clear();
	  RecoJet_rapidity.clear();
	  RecoJet_phi.clear();
	  RecoJet_E.clear();
	  RecoJet_M.clear();
	  RecoJet_rapi.clear();
	  RecoJet_hasElectron.clear();
	  RecoJet_maxPtPart_pt.clear();
	  RecoJet_constituent_pt.clear(); 
	  RecoJet_constituent_eta.clear();
	  RecoJet_constituent_phi.clear(); 
	  RecoJet_constituent_nhits.clear();
	  RecoJet_constituent_pdgid.clear();
	  RecoJet_constituent_pdgidTruth.clear();
	  RecoJet_constituent_idx.clear();
          RecoJet_constituent_energy.clear();
	  RecoJet_constituent_rapi.clear();
	  
	  EVETMULTRECO = (int)rcMomPx.GetSize();  
			
	  //Reco Jet clustering 
	  fastjet::ClusterSequenceArea cs_reco(particles_reco, jet_def, area_def);
	  std::vector<fastjet::PseudoJet> jets_reco = fastjet::sorted_by_pt(cs_reco.inclusive_jets(minJetPt));
	  for (auto &jet : jets_reco) {
	    RecoJet_pt.push_back(jet.pt());
	    RecoJet_eta.push_back(jet.eta());
	    RecoJet_rapidity.push_back(jet.rapidity());
	    RecoJet_phi.push_back(jet.phi_std());
	    RecoJet_E.push_back(jet.e());
	    RecoJet_M.push_back(jet.m());
	    RecoJet_rapi.push_back(jet.rap());	

	    bool hasElectron = false;
	    float maxPtReco = -1.0;
	    std::vector<float> cpt, ceta, cphi, cenergy,crapi;
	    std::vector<int> chits, cpdgid, cpdgidtruth,cindex;
	    cpt.clear(); ceta.clear(); cphi.clear(); chits.clear(); cpdgid.clear(); cpdgidtruth.clear(), cindex.clear(), cenergy.clear(),crapi.clear();

	    fastjet::PseudoJet D0_parent;

	    for (auto &c : jet.constituents()) {
	      int idx = c.user_index();
	      
        //Look for the bkg D0 particle of the Jet 
	      bool is_D0_parent = false;
	      if(idx == -99999){ // check for bkg D0 particle
		for(const auto& p : particles_reco ){ // loop over all the particles of the jet
		  if (p.user_index() == idx) {
		    
		    D0_parent = p;
		    is_D0_parent = true;
		    break;
		  }
		}
                      
		cpt.push_back(D0_parent.pt());
		ceta.push_back(D0_parent.eta());
		cphi.push_back(D0_parent.phi());
		cenergy.push_back(D0_parent.E());
		chits.push_back(0);           // 0 hits for neutral D0 
		cpdgid.push_back(421);        // PID = 421 for D0
		cpdgidtruth.push_back(421);
		cindex.push_back(idx);
		crapi.push_back(D0_parent.rap());
		                                                          // Here, hasElectron = false
		if (D0_parent.pt() > maxPtReco) maxPtReco = D0_parent.pt();   
	      }
	      if (is_D0_parent) continue;// skip the rest of the loop for bkg D0
		    
	      TVector3 v3(rcMomPx[idx], rcMomPy[idx], rcMomPz[idx]);
	      float rcEnergy = TrkRecoE[idx];
	      TLorentzVector lv(rcMomPx[idx],  rcMomPy[idx], rcMomPz[idx], rcEnergy);

	      cpt.push_back(v3.Pt());
	      ceta.push_back(v3.Eta());
	      cphi.push_back(v3.Phi());
	      chits.push_back(TrkRecoNhits[idx]);
	      cindex.push_back(idx);
	      cenergy.push_back(rcEnergy);
	      crapi.push_back(lv.Rapidity());
	      if (v3.Pt() > maxPtReco) maxPtReco = v3.Pt();
	      
	      // Find electron
	      int chargePartIndex = idx; 
	      int elecIndex = -1;
	      float elecIndexWeight = -1.0;
	      for(unsigned int itrkass = 0; itrkass < assocWeight.GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations

		if( assocWeight.GetSize() > 0 ){
		  if(TrkPartAssocRec[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
		    if(assocWeight[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
		      elecIndex = TrkPartAssocSim[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
		      elecIndexWeight = assocWeight[itrkass];
		    }
		  }
		}
	      }

	      if(mcPartPdg[elecIndex] == 11){
                hasElectron = true;
              }

	      cpdgid.push_back(rcPdg[idx]);
	      cpdgidtruth.push_back(mcPartPdg[elecIndex]);
	    }

	    RecoJet_constituent_pt.push_back(cpt);  
	    RecoJet_constituent_eta.push_back(ceta);
	    RecoJet_constituent_phi.push_back(cphi);
	    RecoJet_constituent_nhits.push_back(chits);
	    RecoJet_constituent_pdgid.push_back(cpdgid);
	    RecoJet_constituent_pdgidTruth.push_back(cpdgidtruth);
	    RecoJet_hasElectron.push_back(hasElectron);
	    RecoJet_maxPtPart_pt.push_back(maxPtReco);
	    RecoJet_constituent_idx.push_back(cindex);
            RecoJet_constituent_energy.push_back(cenergy);
	    RecoJet_constituent_rapi.push_back(crapi);

	  } 
 

//Reading Reco Jet(BKG):
	for(unsigned int i=0; i<RecoJet_pt.size(); i++){  //loop over reco Jets  
	  
	  int  jetParticles = RecoJet_constituent_idx[i].size();  // Total no. of particles of the ith Jet
     
	  float jetPt  = RecoJet_pt[i];
	  float jetEta = RecoJet_eta[i];
	  float jetRapidity = RecoJet_rapi[i];
	  float jetPhi = RecoJet_phi[i];
	  float jetE = RecoJet_E[i];
     
	  TLorentzVector LvJet;
	  LvJet.SetPtEtaPhiE(jetPt,jetEta, jetPhi, jetE); // Jet TLorentzVector
	  TVector3 jetMom = LvJet.Vect(); // Jet Momentum Vector
       	  bool c2= false;
     
	  for(unsigned int j=0; j<RecoJet_constituent_idx[i].size(); j++){   //Loop over all particles of the ith Jet

	    int particle_idx = RecoJet_constituent_idx[i][j]; // index of jth particle of ith jet
	    bool is_D0 = false;
            
	    //check if particle is the bkg D0 
	    if(particle_idx == -99999){        
	      is_D0 = true;
	      c2 = true;
	    }
	 
	    if(is_D0){  // if Jet constituent is the bkg D0

	      float D0Pt  = RecoJet_constituent_pt[i][j];
	      float D0Eta = RecoJet_constituent_eta[i][j];
	      float D0Phi = RecoJet_constituent_phi[i][j];
	      float D0E = RecoJet_constituent_energy[i][j];
		  
	      TLorentzVector lvD0;
	      lvD0.SetPtEtaPhiE(D0Pt, D0Eta, D0Phi, D0E);
	      
	      TVector3 D0_Mom;
	      D0_Mom.SetPtEtaPhi(D0Pt, D0Eta, D0Phi);
              float D0Mass = lvD0.M();
	      float D0Rapidity = lvD0.Rapidity();
	
	      float bkg_z = float(jetMom.Dot(D0_Mom))/float(jetMom.Dot(jetMom)) ;  // Fragmentation Variable z for bkg D0

              z_bkg =bkg_z;
              h3bkg_z->Fill(bkg_z,D0Rapidity,D0Mass);
                            
              float dPhi = jetPhi - D0Phi;
              dPhi = TVector2::Phi_mpi_pi(dPhi);
              float dEta = jetEta - D0Eta;
              float deltaR = TMath::Sqrt(dEta*dEta + dPhi*dPhi); // bkg_D0 distance from JetAxis in eta-Phi space
              
              float angleRad = jetMom.Angle(D0_Mom); // angle bw bkg_D0 Momentum vector and jetMomentum vector
              float angleDeg = angleRad * 180.0 / TMath::Pi();
        
              angle_bkg = angleDeg;
              dr_bkg = deltaR;
              pTjet_bkg = jetPt;
              etajet_bkg= jetEta;
	      D0eta_bkg = D0Eta;

	      break; // leave the jet constituent loop once we find the only bkg D0 present
	    } 	 
	  }
	  if(c2) break; // leave the jet constituent loop once we find the only bkg D0 present 
	}
	
       tree_bkg->Fill();
		    }

		  if(dcaToVtx.Perp() >= 0.02 && dcaToVtx2.Perp() >= 0.02 &&
		     dcaDaughters < 0.07 && cosTheta > 0.95 && decayLength > 0.05 && V0DcaToVtx < 0.1)
		    {
		      if(is_D0_pik)
			{
			  h3InvMass[0][1]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
			}
		      else
			{
			  h3InvMass[1][1]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
			}
		    }
		}// if(unlike charge rc pi and k) 
	    }//for(rc k)
	next_pi:
	  continue;
	}//for(rc pi)

      //--------------------------------------------------GenJetClustering--------------------------------------------------------------------------


      vector<unsigned int> mc_pi_index;
      vector<unsigned int> mc_k_index;
      mc_pi_index.clear();
      mc_k_index.clear();
      
      // get all MC pi and K
      for(unsigned int mc_index=0; mc_index<mcPartMass.GetSize(); mc_index++)
	{	  
	  if(abs(mcPartPdg[mc_index]) == 211) mc_pi_index.push_back(mc_index);
	  if(abs(mcPartPdg[mc_index]) == 321) mc_k_index.push_back(mc_index);
	  
	}
		  
      // pair pion and kaon
      for(unsigned int i=0; i<mc_pi_index.size(); i++)//loop over all mc pi
	{
	  for(unsigned int j=0; j<mc_k_index.size(); j++)// loop over all mc k
	    {
	      if(mcPartCharge[mc_pi_index[i]]*mcPartCharge[mc_k_index[j]]<0)// only unlike charged mc pi and mc k
		{
		  bool is_D0_pik = false;
		  for(unsigned int k=0; k<mc_index_D0_pi.size(); k++)
		    {
		      if(mc_pi_index[i]==mc_index_D0_pi[k] && mc_k_index[j]==mc_index_D0_k[k])// checking for if mc pi-k pair coming from a D0 decay
			{
			  is_D0_pik = true;
			  break;
			}
      		    }

		  // Storing mc index of pi and k of the pair
		  int pii_index = mc_pi_index[i]; // mc index of pion
		  int kk_index =  mc_k_index[j];  // mc index of kaon

	          //Signal mc pi-k pair
		  if(is_D0_pik) 
		    {  
		      
		      // PseudoJet Building
		      std::vector<fastjet::PseudoJet> particles_gen;       
		      for (unsigned int i = 0; i < mcPartMass.GetSize(); ++i) {
			
			TVector3 mom(mcMomPx[i], mcMomPy[i], mcMomPz[i]);
			if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
			if ( mcPartCharge[i] == 0 ) continue;
			if ( removeelectrons == 1 && mcPartPdg[i] == 11 ) continue;
			if ( removeelectrons == 2 && i == ScatteredEGenId ) continue;
		    
			if(i == pii_index || i == kk_index) continue; //skipping the pi-k pair in pseudojet making                  	                
			
			float E =sqrt(mcMomPx[i]*mcMomPx[i] + mcMomPy[i]*mcMomPy[i] + mcMomPz[i]*mcMomPz[i] + mcPartMass[i]*mcPartMass[i]);//energyof mc
			fastjet::PseudoJet p(mcMomPx[i], mcMomPy[i], mcMomPz[i], E);   
			p.set_user_index(i);
			particles_gen.push_back(p);
		      }
		  
		      TLorentzVector lv_pii, lv_kk; //making TLorentzVector of the pair pion and kaon
		      lv_pii.SetXYZM(mcMomPx[pii_index], mcMomPy[pii_index], mcMomPz[pii_index], mcPartMass[pii_index]); 
		      lv_kk.SetXYZM(mcMomPx[kk_index], mcMomPy[kk_index], mcMomPz[kk_index], mcPartMass[kk_index]);

		      // calculating lorentz vector of the parent D0 of the pi-k pair
		      TLorentzVector lv_D0 = lv_pii + lv_kk;	
	
		      if ( lv_D0.Pt() < minCstPt || lv_D0.Pt() > maxCstPt ){// if Pt of parent D0 is not in require Pt limits---> Introducing pi and k of the pair to pseudoJet making
			
			for (int l=0;l<2;l++){
			  int idx;
			  
			  if (l==0){
			    idx = pii_index;
			  }
			  
			  if(l==1){
			    idx = kk_index;
			  }

			  TVector3 mom(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx]);
			  if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;

			  float E = sqrt(mcMomPx[idx]*mcMomPx[idx] + mcMomPy[idx]*mcMomPy[idx] + mcMomPz[idx]*mcMomPz[idx] + mcPartMass[idx]*mcPartMass[idx]);
       
			  fastjet::PseudoJet p(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx], E);	 
			  p.set_user_index(idx);
			  particles_gen.push_back(p);
			  
			}
		      }
        
                                                                                                    
		      if ( lv_D0.Pt() > minCstPt && lv_D0.Pt() < maxCstPt ){// introducing parent D0 to pseudojets
			
			fastjet::PseudoJet p(lv_D0.Px(), lv_D0.Py(),lv_D0.Pz(), lv_D0.E());
			p.set_user_index(-99999);// setting signal parent D0 index = negative integer (-99999)                                          
			particles_gen.push_back(p);
		      }
		
	                	    
		      // Define algorithm
		      fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
		      fastjet::RecombinationScheme scheme = fastjet::E_scheme;
		      // Jet definition
		      fastjet::JetDefinition jet_def(algo, R_value, scheme);
		      fastjet::GhostedAreaSpec ghost_spec(ghostMaxRap, numGhostRepeat, ghostArea);
		      fastjet::AreaType atype = fastjet::active_area;
		      fastjet::AreaDefinition area_def(atype, ghost_spec);      

		      //clearing the vectors
		      GenJet_pt.clear();
		      GenJet_eta.clear();
		      GenJet_phi.clear();
		      GenJet_E.clear();
		      GenJet_M.clear();
		      GenJet_rapi.clear();
		      GenJet_hasElectron.clear();
		      GenJet_hasNeutral.clear();
		      GenJet_maxPtPart_pt.clear();
		      GenJet_constituent_pt.clear(); 
		      GenJet_constituent_eta.clear();
		      GenJet_constituent_phi.clear(); 
		      GenJet_constituent_pdgid.clear();
		      GenJet_constituent_rapi.clear();
		      GenJet_constituent_idx.clear();
		      GenJet_constituent_energy.clear();
		      
		      EVETMULTGEN = (int)TrkGenPx.GetSize();
		      
		      //Gen clustering 
		      fastjet:: ClusterSequenceArea cs_gen(particles_gen, jet_def, area_def);
		      std::vector<fastjet::PseudoJet> jets_gen = fastjet::sorted_by_pt(cs_gen.inclusive_jets(minJetPt));
		      for (auto &jet : jets_gen) { //loop over genJets
			GenJet_pt.push_back(jet.pt());
			GenJet_eta.push_back(jet.eta());
			GenJet_phi.push_back(jet.phi_std());
			GenJet_E.push_back(jet.e());
			GenJet_M.push_back(jet.m());
			GenJet_rapi.push_back(jet.rap());

	
			bool hasGenElectron = false;
			bool hasGenNeutral = false;
			float maxPtGen = -1.0;
			std::vector<float> gpt, geta, gphi,grapi,genergy;
			std::vector<int> gpdgid,gidx;
			gpt.clear(); geta.clear(); gphi.clear(); gpdgid.clear();gidx.clear();grapi.clear();genergy.clear();
			
			fastjet::PseudoJet D0_parent;
			
			for (auto &c : jet.constituents()) {//loop over genJetConstituents
			  int idx = c.user_index();
			  
			  bool is_D0_parent = false;
			  if(idx == -99999){     //looking for D0 parent                                                                                        
			    for(const auto& p : particles_gen ){// loop over all the pseudojets                                                         

			      if (p.user_index() == idx) {
		  	
				D0_parent = p;
				is_D0_parent = true;
				break;
			      }
			    }
		
			    gpt.push_back(D0_parent.pt());
			    geta.push_back(D0_parent.eta());
			    gphi.push_back(D0_parent.phi());
			    genergy.push_back(D0_parent.E());                                             
			    gpdgid.push_back(421);                                                                                   
			    gidx.push_back(idx);
			    grapi.push_back(D0_parent.rap());
			    
			    if (D0_parent.pt() > maxPtGen) maxPtGen = D0_parent.pt();
			    hasGenElectron = false;
			    hasGenNeutral = true;
			  }
              
			  if (is_D0_parent) continue;// skip the rest of the loop for the parent D0
	  
	  
			  TVector3 gv(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx]);
			  TLorentzVector glv;
			  glv.SetXYZM(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx], mcPartMass[idx]);	

			  gpt.push_back(gv.Pt());
			  geta.push_back(gv.Eta());
			  gphi.push_back(gv.Phi());
			  gpdgid.push_back(mcPartPdg[idx]);
			  gidx.push_back(idx);
			  grapi.push_back(glv.Rapidity());
			  genergy.push_back(glv.E());
			  
			  if (gv.Pt() > maxPtGen) maxPtGen = gv.Pt();
			  
			  if (mcPartPdg[idx] == 11) hasGenElectron = true;
			  if (mcPartCharge[idx] == 0) hasGenNeutral = true;
			}
			
			GenJet_constituent_pt.push_back(gpt);
			GenJet_constituent_eta.push_back(geta);
			GenJet_constituent_phi.push_back(gphi);
			GenJet_hasElectron.push_back(hasGenElectron);
			GenJet_hasNeutral.push_back(hasGenNeutral);
			GenJet_maxPtPart_pt.push_back(maxPtGen);
			GenJet_constituent_pdgid.push_back(gpdgid);
			GenJet_constituent_rapi.push_back(grapi);
			GenJet_constituent_idx.push_back(gidx);
			GenJet_constituent_energy.push_back(genergy);
		      }
		      
		      //Gen Jet Reading

		      for(int i=0;i<GenJet_eta.size();i++){//loop over all the genJets

			int  jetmulti = GenJet_constituent_idx[i].size();  // Total no. of particles of the ith genJet
     
			float jetPt  = GenJet_pt[i];
			float jetEta = GenJet_eta[i];
			float jetRapi = GenJet_rapi[i];
			float jetPhi = GenJet_phi[i];
			float jetE = GenJet_E[i];
			
			TLorentzVector lvJet;
			lvJet.SetPtEtaPhiE(jetPt,jetEta, jetPhi, jetE); // Jet TLorentzVector	  
			TVector3 jetMom = lvJet.Vect(); // Jet Momentum Vector
			bool c2 = false;
			
			for(unsigned int j=0; j<GenJet_constituent_idx[i].size(); j++){   //Loop over all particles of the ith Jet

			  int genIdx = GenJet_constituent_idx[i][j]; // index of jth particle of ith jet
			  bool isD0 = false;
			 
			  
			  //check if jet constituent is the D0 parent
			  if( genIdx == -99999 ){        
			    isD0 = true;
			    c2 = true;
			  }
	 
			  if(isD0){  // if Jet constituent is the D0 parent

			    float D0Pt  = GenJet_constituent_pt[i][j];
			    float D0Eta = GenJet_constituent_eta[i][j];
			    float D0Phi = GenJet_constituent_phi[i][j];
			    float D0E = GenJet_constituent_energy[i][j];
			    float D0Rapi = GenJet_constituent_rapi[i][j];
	      
			    TLorentzVector lvD0;
			    lvD0.SetPtEtaPhiE(D0Pt, D0Eta, D0Phi, D0E);
			    
			    TVector3 D0_Mom;
			    D0_Mom.SetPtEtaPhi(D0Pt, D0Eta, D0Phi);
			    
			    float gz = float(jetMom.Dot(D0_Mom))/float(jetMom.Dot(jetMom)) ;  // Fragmentation Variable z 
			    cout<<"genZ = "<<gz<<endl;
                            
			    float dPhi = jetPhi - D0Phi;
			    dPhi = TVector2::Phi_mpi_pi(dPhi);
			    
			    float dEta = jetEta - D0Eta;
			    float drapi = jetRapi - D0Rapi;
	      
			    float deltaR_eta = TMath::Sqrt(dEta*dEta + dPhi*dPhi); // D0 distance from JetAxis in eta-Phi space
			    float deltaR_y = TMath::Sqrt(drapi*drapi + dPhi*dPhi); // D0 distance from JetAxis in y-Phi space 
	      
			    float angleRad = jetMom.Angle(D0_Mom); // angle bw parent D0 Momentum vector and jetMomentum vector
			    float angleDeg = angleRad * 180.0 / TMath::Pi();

			    //storing D0 parameters to branches of gen signal tree
			    genZ_sig = gz;
			    gdr_eta_sig = deltaR_eta;
			    gdr_y_sig = deltaR_y;
			    gD0jetAngle_sig = angleDeg;
			    pt_gjet_sig = jetPt;
			    eta_gjet_sig = jetEta;
			    y_gjet_sig = jetRapi;
	      
			    break;// break the jetConstituent loop immediately the D0 parent is found
			  }
			}
			if(c2) break;// break the jet loop immediately when the D0 parent is found
		      }
		      tree_gen_sig->Fill(); // tree filled only fill once, in the that jet case which have D0 in it
		    }


		  //gen BKG D0
		  if(!is_D0_pik) 
		    {  

		    		      
		      // PseudoJet Building
		      std::vector<fastjet::PseudoJet> particles_gen;       
		      for (unsigned int i = 0; i < mcPartMass.GetSize(); ++i) {
			
			TVector3 mom(mcMomPx[i], mcMomPy[i], mcMomPz[i]);
			if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;
			if ( mcPartCharge[i] == 0 ) continue;
			if ( removeelectrons == 1 && mcPartPdg[i] == 11 ) continue;
			if ( removeelectrons == 2 && i == ScatteredEGenId ) continue;
		    
			if(i == pii_index || i == kk_index) continue; //skipping the pi-k pair in pseudojet making                  	                
			
			float E =sqrt(mcMomPx[i]*mcMomPx[i] + mcMomPy[i]*mcMomPy[i] + mcMomPz[i]*mcMomPz[i] + mcPartMass[i]*mcPartMass[i]);//energyof mc
			fastjet::PseudoJet p(mcMomPx[i], mcMomPy[i], mcMomPz[i], E);   
			p.set_user_index(i);
			particles_gen.push_back(p);
		      }
		  
		      TLorentzVector lv_pii, lv_kk; //making TLorentzVector of the pair pion and kaon
		      lv_pii.SetXYZM(mcMomPx[pii_index], mcMomPy[pii_index], mcMomPz[pii_index], mcPartMass[pii_index]); 
		      lv_kk.SetXYZM(mcMomPx[kk_index], mcMomPy[kk_index], mcMomPz[kk_index], mcPartMass[kk_index]);

		      // calculating lorentz vector of the parent D0 of the pi-k pair
		      TLorentzVector lv_D0 = lv_pii + lv_kk;	
	
		      if ( lv_D0.Pt() < minCstPt || lv_D0.Pt() > maxCstPt ){// if Pt of parent D0 is not in require Pt limits---> Introducing pi and k of the pair to pseudoJet making
			
			for (int l=0;l<2;l++){
			  int idx;
			  
			  if (l==0){
			    idx = pii_index;
			  }
			  
			  if(l==1){
			    idx = kk_index;
			  }

			  TVector3 mom(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx]);
			  if ( mom.Pt() < minCstPt || mom.Pt() > maxCstPt ) continue;

			  float E = sqrt(mcMomPx[idx]*mcMomPx[idx] + mcMomPy[idx]*mcMomPy[idx] + mcMomPz[idx]*mcMomPz[idx] + mcPartMass[idx]*mcPartMass[idx]);
       
			  fastjet::PseudoJet p(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx], E);	 
			  p.set_user_index(idx);
			  particles_gen.push_back(p);
			  
			}
		      }
        
                                                                                                    
		      if ( lv_D0.Pt() > minCstPt && lv_D0.Pt() < maxCstPt ){// introducing parent D0 to pseudojets
			
			fastjet::PseudoJet p(lv_D0.Px(), lv_D0.Py(),lv_D0.Pz(), lv_D0.E());
			p.set_user_index(-99999);// setting bkg  parent D0 index = negative integer (-99999)                                          
			particles_gen.push_back(p);
		      }
		
	                	    
		      // Define algorithm
		      fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
		      fastjet::RecombinationScheme scheme = fastjet::E_scheme;
		      // Jet definition
		      fastjet::JetDefinition jet_def(algo, R_value, scheme);
		      fastjet::GhostedAreaSpec ghost_spec(ghostMaxRap, numGhostRepeat, ghostArea);
		      fastjet::AreaType atype = fastjet::active_area;
		      fastjet::AreaDefinition area_def(atype, ghost_spec);      

		      //clearing the vectors
		      GenJet_pt.clear();
		      GenJet_eta.clear();
		      GenJet_phi.clear();
		      GenJet_E.clear();
		      GenJet_M.clear();
		      GenJet_rapi.clear();
		      GenJet_hasElectron.clear();
		      GenJet_hasNeutral.clear();
		      GenJet_maxPtPart_pt.clear();
		      GenJet_constituent_pt.clear(); 
		      GenJet_constituent_eta.clear();
		      GenJet_constituent_phi.clear(); 
		      GenJet_constituent_pdgid.clear();
		      GenJet_constituent_rapi.clear();
		      GenJet_constituent_idx.clear();
		      GenJet_constituent_energy.clear();
		      
		      EVETMULTGEN = (int)TrkGenPx.GetSize();
		      
		      //Gen clustering 
		      fastjet:: ClusterSequenceArea cs_gen(particles_gen, jet_def, area_def);
		      std::vector<fastjet::PseudoJet> jets_gen = fastjet::sorted_by_pt(cs_gen.inclusive_jets(minJetPt));
		      for (auto &jet : jets_gen) { //loop over genJets
			GenJet_pt.push_back(jet.pt());
			GenJet_eta.push_back(jet.eta());
			GenJet_phi.push_back(jet.phi_std());
			GenJet_E.push_back(jet.e());
			GenJet_M.push_back(jet.m());
			GenJet_rapi.push_back(jet.rap());

	
			bool hasGenElectron = false;
			bool hasGenNeutral = false;
			float maxPtGen = -1.0;
			std::vector<float> gpt, geta, gphi,grapi,genergy;
			std::vector<int> gpdgid,gidx;
			gpt.clear(); geta.clear(); gphi.clear(); gpdgid.clear();gidx.clear();grapi.clear();genergy.clear();
			
			fastjet::PseudoJet D0_parent;
			
			for (auto &c : jet.constituents()) {//loop over genJetConstituents
			  int idx = c.user_index();
			  
			  bool is_D0_parent = false;
			  if(idx == -99999){     //looking for D0 parent                                                                                        
			    for(const auto& p : particles_gen ){// loop over all the pseudojets                                                         

			      if (p.user_index() == idx) {
		  	
				D0_parent = p;
				is_D0_parent = true;
				break;
			      }
			    }
		
			    gpt.push_back(D0_parent.pt());
			    geta.push_back(D0_parent.eta());
			    gphi.push_back(D0_parent.phi());
			    genergy.push_back(D0_parent.E());                                             
			    gpdgid.push_back(421);                                                                                   
			    gidx.push_back(idx);
			    grapi.push_back(D0_parent.rap());
			    
			    if (D0_parent.pt() > maxPtGen) maxPtGen = D0_parent.pt();
			    hasGenElectron = false;
			    hasGenNeutral = true;
			  }
              
			  if (is_D0_parent) continue;// skip the rest of the loop for the parent D0
	  
	  
			  TVector3 gv(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx]);
			  TLorentzVector glv;
			  glv.SetXYZM(mcMomPx[idx], mcMomPy[idx], mcMomPz[idx], mcPartMass[idx]);	

			  gpt.push_back(gv.Pt());
			  geta.push_back(gv.Eta());
			  gphi.push_back(gv.Phi());
			  gpdgid.push_back(mcPartPdg[idx]);
			  gidx.push_back(idx);
			  grapi.push_back(glv.Rapidity());
			  genergy.push_back(glv.E());
			  
			  if (gv.Pt() > maxPtGen) maxPtGen = gv.Pt();
			  
			  if (mcPartPdg[idx] == 11) hasGenElectron = true;
			  if (mcPartCharge[idx] == 0) hasGenNeutral = true;
			}
			
			GenJet_constituent_pt.push_back(gpt);
			GenJet_constituent_eta.push_back(geta);
			GenJet_constituent_phi.push_back(gphi);
			GenJet_hasElectron.push_back(hasGenElectron);
			GenJet_hasNeutral.push_back(hasGenNeutral);
			GenJet_maxPtPart_pt.push_back(maxPtGen);
			GenJet_constituent_pdgid.push_back(gpdgid);
			GenJet_constituent_rapi.push_back(grapi);
			GenJet_constituent_idx.push_back(gidx);
			GenJet_constituent_energy.push_back(genergy);
		      }
		      
		      //Gen Bkg Jet Reading

		      for(int i=0;i<GenJet_eta.size();i++){//loop over all the genJets

			int  jetmulti = GenJet_constituent_idx[i].size();  // Total no. of particles of the ith genJet
     
			float jetPt  = GenJet_pt[i];
			float jetEta = GenJet_eta[i];
			float jetRapi = GenJet_rapi[i];
			float jetPhi = GenJet_phi[i];
			float jetE = GenJet_E[i];
			
			TLorentzVector lvJet;
			lvJet.SetPtEtaPhiE(jetPt,jetEta, jetPhi, jetE); // Jet TLorentzVector	  
			TVector3 jetMom = lvJet.Vect(); // Jet Momentum Vector
			bool c2 = false;
			
			for(unsigned int j=0; j<GenJet_constituent_idx[i].size(); j++){   //Loop over all particles of the ith Jet

			  int genIdx = GenJet_constituent_idx[i][j]; // index of jth particle of ith jet
			  bool isD0 = false;
			 
			  
			  //check if jet constituent is the D0 parent
			  if( genIdx == -99999 ){        
			    isD0 = true;
			    c2 = true;
			  }
	 
			  if(isD0){  // if Jet constituent is the D0 parent

			    float D0Pt  = GenJet_constituent_pt[i][j];
			    float D0Eta = GenJet_constituent_eta[i][j];
			    float D0Phi = GenJet_constituent_phi[i][j];
			    float D0E = GenJet_constituent_energy[i][j];
			    float D0Rapi = GenJet_constituent_rapi[i][j];
	      
			    TLorentzVector lvD0;
			    lvD0.SetPtEtaPhiE(D0Pt, D0Eta, D0Phi, D0E);
			    
			    TVector3 D0_Mom;
			    D0_Mom.SetPtEtaPhi(D0Pt, D0Eta, D0Phi);
			    
			    float gz = float(jetMom.Dot(D0_Mom))/float(jetMom.Dot(jetMom)) ;  // Fragmentation Variable z 
			    cout<<"genZ = "<<gz<<endl;
                            
			    float dPhi = jetPhi - D0Phi;
			    dPhi = TVector2::Phi_mpi_pi(dPhi);
			    
			    float dEta = jetEta - D0Eta;
			    float drapi = jetRapi - D0Rapi;
	      
			    float deltaR_eta = TMath::Sqrt(dEta*dEta + dPhi*dPhi); // D0 distance from JetAxis in eta-Phi space
			    float deltaR_y = TMath::Sqrt(drapi*drapi + dPhi*dPhi); // D0 distance from JetAxis in y-Phi space 
	      
			    float angleRad = jetMom.Angle(D0_Mom); // angle bw parent D0 Momentum vector and jetMomentum vector
			    float angleDeg = angleRad * 180.0 / TMath::Pi();

			    //storing D0 parameters to branches of gen signal tree
			    genZ_bkg = gz;
			    gdr_eta_bkg = deltaR_eta;
			    gdr_y_bkg = deltaR_y;
			    gD0jetAngle_bkg = angleDeg;
			    pt_gjet_bkg = jetPt;
			    eta_gjet_bkg = jetEta;
			    y_gjet_bkg = jetRapi;
	      
			    break;// break the jetConstituent loop immediately the D0 parent is found
			  }
			}
			if(c2) break;// break the jet loop immediately when the D0 parent is found
		      }
		      tree_gen_bkg->Fill(); // tree filled only fill once, in the that jet case which have D0 in it
		    }

		}
	    }
	}

      
      //--------------------------------------------------GenJetClustering-----------------------------------------------------------------------------
      cout<<nevents<<" Completed !!"<<endl;
      nevents++;
    } 

  file_signal->cd();  
  tree_sig->Write();
  file_signal->Close();  
  
  file_bkg->cd();  
  tree_bkg->Write();
  file_bkg->Close();

  file_gen->cd();
  tree_gen_sig->Write();
  tree_gen_bkg->Write();
  file_gen->Close();
  
  fout_mcgen->cd();
  tree_D0->Write();
  fout_mcgen->Close();
  
  TFile *outfile = new TFile(outname.Data(), "recreate");

  hEventStat->Write();
  hMcMult->Write();
  hMcVtxX->Write();
  hMcVtxY->Write();
  hMcVtxZ->Write();
  
  hPullVtxX->Write();
  hPullVtxY->Write();  
  hPullVtxZ->Write();
  
  hRes_SVx_Helixfit->Write();
  hRes_SVy_Helixfit->Write();
  hRes_SVz_Helixfit->Write(); 
  hRes_SVxy_Helixfit->Write();  
 
  hRes_SVx_Helixana->Write();
  hRes_SVy_Helixana->Write();
  hRes_SVz_Helixana->Write();
  hRes_SVxy_Helixana->Write();    
      
  hchi2_vtx->Write();
  hchi2_vtx_sig->Write();
  hchi2_vtx_bkg->Write();  
  hRes_SVx_Helixfit_pull->Write();
  hRes_SVy_Helixfit_pull->Write(); 
  hRes_SVz_Helixfit_pull->Write(); 
    
  hD0DecayVxVy->Write();
  hD0DecayVrVz->Write();
  
  hMCD0PtRap->Write();

  hMcPiPtEta->Write();
  hMcPiPtEtaReco->Write();
  hMcKPtEta->Write();
  hMcKPtEtaReco->Write();
  
  hNRecoVtx->Write();
  h3sig_z->Write();
  h3bkg_z->Write();

  hmc_eta->Write();
  hmc_eta_e->Write();
  hmc_eta_pi->Write();
  hmc_eta_k->Write();
  hmc_eta_p->Write();

  hreco_eta->Write();
  hreco_eta_e->Write();
  hreco_eta_pi->Write();
  hreco_eta_k->Write();
  hreco_eta_p->Write();

  

  for(int ip=0; ip<2; ip++)
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
TVector3 getDcaToVtx(const int index, TVector3 vtx)
{
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
TLorentzVector getPairParent(const int index1, const int index2, TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &cosTheta_xy, float &decayLength, float &V0DcaToVtx, float &sigma_vtx, TVector3 &decayVertex, TVector3 &decayVertex_ana, double &chi2_ndf, double * parFitErr)
{
  // -- get helix
  TVector3 pos1(rcTrkLoca2->At(index1) * sin(rcTrkPhi2->At(index1)) * -1 * millimeter, rcTrkLoca2->At(index1) * cos(rcTrkPhi2->At(index1)) * millimeter, rcTrkLocb2->At(index1) * millimeter);
  TVector3 mom1(rcMomPx2->At(index1), rcMomPy2->At(index1), rcMomPz2->At(index1));

  TVector3 pos2(rcTrkLoca2->At(index2) * sin(rcTrkPhi2->At(index2)) * -1 * millimeter, rcTrkLoca2->At(index2) * cos(rcTrkPhi2->At(index2)) * millimeter, rcTrkLocb2->At(index2) * millimeter);
  TVector3 mom2(rcMomPx2->At(index2), rcMomPy2->At(index2), rcMomPz2->At(index2));

  float charge1 = rcCharge2->At(index1);
  float charge2 = rcCharge2->At(index2);
  
  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, charge1);
  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla, charge2);

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);
  
  double s1, s2;
  getDecayVertex_Chi2fit(index1, index2, s1, s2, decayVertex, chi2_ndf, parFitErr);
 
  TVector3 const p1AtDcaToP2 = p1Helix.at(s1);
  TVector3 const p2AtDcaToP1 = p2Helix.at(s2);
  // printf("p1AtDcaToP2 origin = (%2.4f, %2.4f, %2.4f)\n", p1AtDcaToP2.x(), p1AtDcaToP2.y(), p1AtDcaToP2.z());
  // printf("p2AtDcaToP1 origin = (%2.4f, %2.4f, %2.4f)\n", p2AtDcaToP1.x(), p2AtDcaToP1.y(), p2AtDcaToP1.z());
  
  // -- calculate DCA of particle1 to particle2 at their DCA
  dcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag()/millimeter;
	
  // -- calculate Lorentz vector of particle1-particle2 pair
  TVector3 const p1MomAtDca = p1Helix.momentumAt(s1,  bField * tesla);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(s2, bField * tesla);
  
  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+gPionMass*gPionMass));
  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+gKaonMass*gKaonMass));
  
  TLorentzVector parent = p1FourMom + p2FourMom;

  // -- calculate decay vertex (secondary or tertiary)
  decayVertex_ana = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
  sigma_vtx = sqrt((p1AtDcaToP2-decayVertex).Mag2()+(p2AtDcaToP1-decayVertex).Mag2())/millimeter;
	
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
    
  //TVector3 dcaToVtx = getDcaToVtx(parent.Vect(), decayVertex, 0, vtx);
  //V0DcaToVtx = dcaToVtx.Mag();
  return parent;
}



void getDecayVertex_Chi2fit(const int index1, const int index2, double &s1, double &s2, TVector3 &vertex, double &chi2_ndf, double *parFitErr)
{
    TVector3 pos1(rcTrkLoca2->At(index1) * sin(rcTrkPhi2->At(index1)) * -1 * millimeter,
                  rcTrkLoca2->At(index1) * cos(rcTrkPhi2->At(index1)) * millimeter,
                  rcTrkLocb2->At(index1) * millimeter);

    TVector3 mom1(rcMomPx2->At(index1), rcMomPy2->At(index1), rcMomPz2->At(index1));

    TVector3 pos2(rcTrkLoca2->At(index2) * sin(rcTrkPhi2->At(index2)) * -1 * millimeter,
                  rcTrkLoca2->At(index2) * cos(rcTrkPhi2->At(index2)) * millimeter,
                  rcTrkLocb2->At(index2) * millimeter);

    TVector3 mom2(rcMomPx2->At(index2), rcMomPy2->At(index2), rcMomPz2->At(index2));


    float charge1 = rcCharge2->At(index1);
    float charge2 = rcCharge2->At(index2);


    StPhysicalHelix helix1(mom1, pos1, bField * tesla, charge1);
    StPhysicalHelix helix2(mom2, pos2, bField * tesla, charge2);
   
    std::array<float, 21>& cov_track1 = rcTrkCov->At(index1);
    std::array<float, 21>& cov_track2 = rcTrkCov->At(index2);  
    
    pair<double, double> const ss = helix1.pathLengths(helix2);
    TVector3 const p1_init = helix1.at(ss.first);
    TVector3 const p2_init = helix2.at(ss.second); 
    TVector3 const mid_point = 0.5*(p1_init+p2_init); 

      // Perform Minimization
    const Int_t nPar = 5;
    Chi2Minimization d2Function(helix1,helix2,cov_track1,cov_track2);
    ROOT::Math::Functor fcn(d2Function,nPar); // 5 parameters
    ROOT::Fit::Fitter fitter;

    double pStart[nPar] = {mid_point.X(),mid_point.Y(),mid_point.Z(),ss.first,ss.second};
    fitter.SetFCN(fcn, pStart,nPar,1);
    
    fitter.Config().ParSettings(0).SetName("x0");
    fitter.Config().ParSettings(0).SetStepSize(0.01);
   // fitter.Config().ParSettings(0).SetLimits(-1., 1.);    
    // No limits for x, y, z

    fitter.Config().ParSettings(1).SetName("y0");
    fitter.Config().ParSettings(1).SetStepSize(0.01);
    //fitter.Config().ParSettings(1).SetLimits(-1., 1.);

    fitter.Config().ParSettings(2).SetName("z0");
    fitter.Config().ParSettings(2).SetStepSize(0.01);
   // fitter.Config().ParSettings(2).SetLimits(-10., 10.);    

    fitter.Config().ParSettings(3).SetName("s1");
    fitter.Config().ParSettings(3).SetValue(0.0);
    fitter.Config().ParSettings(3).SetStepSize(0.01);
    //fitter.Config().ParSettings(3).SetLimits(-1., 1.);

    fitter.Config().ParSettings(4).SetName("s2");
    fitter.Config().ParSettings(4).SetValue(0.0);
    fitter.Config().ParSettings(4).SetStepSize(0.01);
    //fitter.Config().ParSettings(4).SetLimits(-1., 1.);      
    
    fitter.Config().MinimizerOptions().SetMaxIterations(10000);	
    // do the fit 

    Bool_t ok = fitter.FitFCN();
    if (!ok) Error("Fitting","Fitting failed");
    const ROOT::Fit::FitResult & result = fitter.Result();
   // double chi2 = fitter.Result().Chi2();
    double ndf = 2*3-nPar;
    chi2_ndf = fitter.Result().MinFcnValue()/ndf;  // Minimum value of your function
    
    int status = fitter.Result().Status();
  //  if (status>0 ) {printf("Fit Failed!!!!\n");}
   // if (status>0 || chi2_ndf>10. ) return;
   // cout <<"\033[1;31m Fit Result Chi2:\033[0m"<<chi2_ndf<<endl;
   // result.Print(std::cout);
   
   // Get the covariance matrix
 //  TMatrixDSym covMatrix(5);
  // result.GetCovarianceMatrix(covMatrix); // Matrix for the parameter errors
  // covMatrix.Print();
   
   const double * parFit = result.GetParams();
   const double *FitErr = result.GetErrors();
    
   for (int i = 0; i < nPar; ++i) parFitErr[i] = FitErr[i];
   vertex.SetXYZ(parFit[0], parFit[1], parFit[2]);
   s1 = parFit[3]; s2 = parFit[4];
}

