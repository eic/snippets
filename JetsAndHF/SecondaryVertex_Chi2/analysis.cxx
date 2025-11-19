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

using namespace std;

StPhysicalHelix* gHelix1 = nullptr;
StPhysicalHelix* gHelix2 = nullptr;

const double gPionMass = 0.13957;
const double gKaonMass = 0.493677;

const double twoPi = 2.*3.1415927;
const double eMass = 0.000511;

const double bField = -1.7; // Tesla

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
	
  if(argc!=3 && argc!=1) return 0;

  TString listname;
  TString outname;

  if(argc==1)
    {
      listname  = "/lustrehome/skumar/eic/my_analysis/D0_Sample_Chi2/test.list";
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
  rcTrkCov = new TTreeReaderArray<std::array<float, 21>>{treereader, "CentralCKFTrackParameters.covariance.covariance[21]"};

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

  
    // Create a ROOT file to store the Ntuple
   TFile *file_signal = new TFile("SignalD0.root", "RECREATE");
   TTree *tree_sig = new TTree("treeMLSig", "treeMLSig"); 

  // Define variables to store in the Ntuple
  float d0_pi_sig, d0_k_sig, d0xy_pi_sig, d0xy_k_sig, sum_d0xy_sig, dca_12_sig, dca_D0_sig, decay_length_sig;
  float costheta_sig, costhetaxy_sig, pt_D0_sig, y_D0_sig, mass_D0_sig, sigma_vtx_sig, mult_sig, signif_d0xy_pi_sig, signif_d0xy_k_sig;
  
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
  
  TFile *file_bkg = new TFile("BkgD0.root", "RECREATE");
  TTree *tree_bkg = new TTree("treeMLBkg", "treeMLBkg"); 
  
  // Define variables to store in the Ntuple
  float d0_pi_bkg, d0_k_bkg, d0xy_pi_bkg, d0xy_k_bkg, sum_d0xy_bkg, dca_12_bkg, dca_D0_bkg, decay_length_bkg;
  float costheta_bkg, costhetaxy_bkg, pt_D0_bkg, y_D0_bkg, mass_D0_bkg, sigma_vtx_bkg, mult_bkg, signif_d0xy_pi_bkg, signif_d0xy_k_bkg;
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
  
  int nevents = 0;
  int mult_charged = 0;
  while(treereader.Next())
    {
      if(nevents%1000==0) printf("\n[i] New event %d\n",nevents);
      //if(nevents==20) break;

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
	  TVector3 dcaToVtx = getDcaToVtx(pi_index[i], vertex_rc);
	   std::array<float, 21>& cov_pion = rcTrkCov->At(pi_index[i]);
     int q_pion = rcCharge[pi_index[i]];
	  for(unsigned int j=0; j<k_index.size(); j++)
	    {
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

				// Signed Δv_xy
			  TVector3 dV = (decayVertex-MCVertex_Kaon*0.1)*10; dV.SetZ(0.);   // remove z-component
			  TVector3 pTD0V = parent.Vect(); pTD0V.SetZ(0.);
              double dvxy_signed = dV.Dot(pTD0V.Unit()); 
		      hRes_SVxy_Helixfit->Fill(parent.Pt(), parent.Rapidity(), dvxy_signed);
		      hRes_SVxy_Helixana->Fill(parent.Pt(), parent.Rapidity(), dvxy_signed);
		      
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
			    tree_sig->Fill();
		    }
		  else
		    {
          hEventStat->Fill(6.5);   // Bkg
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
		}
	    }
	}
	      
      nevents++;
    }

  file_signal->cd();  
  tree_sig->Write();
  file_signal->Close();  
  
  file_bkg->cd();  
  tree_bkg->Write();
  file_bkg->Close();

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


