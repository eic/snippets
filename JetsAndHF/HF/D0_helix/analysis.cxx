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

const double twoPi = 2.*3.1415927;
const double eMass = 0.000511;

const double bField = -1.7; // Tesla

TVector3 getDcaToVtx(TVector3 pos, TVector3 mom, int charge, TVector3 vtx);

TLorentzVector getPairParent(TVector3 pos1, TVector3 mom1, int charge1, double mass1,
			     TVector3 pos2, TVector3 mom2, int charge2, double mass2,
			     TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &decayLength, float &V0DcaToVtx);

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

  // ===== read in root files =====
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

  // ===== define histograms =====
  TH1F *hEventStat = new TH1F("hEventStat", "Event statistics", 5, 0, 5);
  hEventStat->GetXaxis()->SetBinLabel(1, "MC events");
  hEventStat->GetXaxis()->SetBinLabel(2, "D0");
  hEventStat->GetXaxis()->SetBinLabel(3, "D0 -> pi+K");
  hEventStat->GetXaxis()->SetBinLabel(4, "Reco D0");

  TH1F *hMcVtxX = new TH1F("hMcVtxX", "x position of MC vertex;x (mm)", 100, -5.05, 4.95);
  TH1F *hMcVtxY = new TH1F("hMcVtxY", "y position of MC vertex;y (mm)", 500, -5.01, 4.99);
  TH1F *hMcVtxZ = new TH1F("hMcVtxZ", "z position of MC vertex;z (mm)", 400, -200, 200);

  TH2F *hMCD0PtRap = new TH2F("hMCD0PtRap", "MC D^{0};y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);
  
  TH3F *h3InvMass = new TH3F("h3InvMass", "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 1.6, 2.0);


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
  TTreeReaderArray<double> mcMomPx = {treereader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomPy = {treereader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomPz = {treereader, "MCParticles.momentum.z"};
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

  // needed for accessing number of hits used for track reconstruction
  TTreeReaderArray<unsigned int> partAssocTrk_begin = {treereader, "ReconstructedChargedParticles.tracks_begin"};
  TTreeReaderArray<int> partAssocTrk_index = {treereader, "_ReconstructedChargedParticles_tracks.index"};
  TTreeReaderArray<int> trkAssocTraj_index = {treereader, "_CentralCKFTracks_trajectory.index"};
  TTreeReaderArray<unsigned int> nMeasurements = {treereader, "CentralCKFTrajectories.nMeasurements"};

  TTreeReaderArray<float> rcTrkLoca = {treereader, "CentralCKFTrackParameters.loc.a"};
  TTreeReaderArray<float> rcTrkLocb = {treereader, "CentralCKFTrackParameters.loc.b"};
  TTreeReaderArray<float> rcTrkqOverP = {treereader, "CentralCKFTrackParameters.qOverP"};
  TTreeReaderArray<float> rcTrkTheta = {treereader, "CentralCKFTrackParameters.theta"};
  TTreeReaderArray<float> rcTrkPhi = {treereader, "CentralCKFTrackParameters.phi"};

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
  
  int nevents = 0;
  while(treereader.Next())
    {
      if(nevents%1000==0) printf("\n[i] New event %d\n",nevents);

      // ===== find MC primary vertex
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

      // ===== get RC primary vertex
      TVector3 vertex_rc(-999., -999., -999.);
      if(prim_vtx_index.GetSize()>0)
	{
	  int rc_vtx_index = prim_vtx_index[0];
	  vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
	}

      // ===== map MC and RC particles
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

      // ===== look for D0 -> pi+K
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
	  hEventStat->Fill(2.5);

	  // D0 kinematics
	  TLorentzVector mc_mom_vec;
	  mc_mom_vec.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
	  double mcRap = mc_mom_vec.Rapidity();
	  double mcPt = mc_mom_vec.Pt();
	  hMCD0PtRap->Fill(mcRap, mcPt);
	}


      // ===== Get reconstructed pions and kaons
      const int pid_mode = 1; // 0 - truth; 1 - realistic
      vector<unsigned int> pi_index;
      vector<unsigned int> k_index;
      pi_index.clear();
      k_index.clear();
      for(unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{

	  // cut on nhits

	  // this assumes 1-to-1-to-1 correspondance between particle-track_trajectory
	  int nhits = nMeasurements[rc_index];
	  
	  // the following implemenation is more robust
	  // int begin_index = partAssocTrk_begin[rc_index]; // take the first track; currently there is 1-to-1 correspondance between track and particle
	  // int track_index = partAssocTrk_index[begin_index];
	  // int traj_index = trkAssocTraj_index[track_index]; // 1-to-1 correspondance between track and trajectory
	  // int nhits = nMeasurements[traj_index];
	  // printf("rc_index = %d, begin_index = %d, track_index = %d, traj_index = %d, nMeasurements = %d = %d\n", rc_index, begin_index, track_index, traj_index, nhits, nMeasurements[rc_index]);

	  if(nhits<4) continue;
	  
	  if(pid_mode==0)
	    {
	      // Use truth PID information
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
	      // Use reconstructed PID
	      if(fabs(rcPdg[rc_index]) == 211) pi_index.push_back(rc_index);
	      if(fabs(rcPdg[rc_index]) == 321) k_index.push_back(rc_index);
	    }
	}

      // ===== pair pion and kaon
      for(unsigned int i=0; i<pi_index.size(); i++)
	{
	  // convert from local to global coordinate
	  // Helix code uses centimeter as the default unit, while ePIC software uses millimeter
	  TVector3 pos1(rcTrkLoca[pi_index[i]] * sin(rcTrkPhi[pi_index[i]]) * -1 * millimeter, rcTrkLoca[pi_index[i]] * cos(rcTrkPhi[pi_index[i]]) * millimeter, rcTrkLocb[pi_index[i]] * millimeter);
	  TVector3 mom1(rcMomPx[pi_index[i]], rcMomPy[pi_index[i]], rcMomPz[pi_index[i]]);
	  int charge1 = rcCharge[pi_index[i]];
	  
	  TVector3 dcaToVtx = getDcaToVtx(pos1, mom1, charge1, vertex_rc);
	  for(unsigned int j=0; j<k_index.size(); j++)
	    {
	      TVector3 pos2(rcTrkLoca[k_index[j]] * sin(rcTrkPhi[k_index[j]]) * -1 * millimeter, rcTrkLoca[k_index[j]] * cos(rcTrkPhi[k_index[j]]) * millimeter, rcTrkLocb[k_index[j]] * millimeter);
	      TVector3 mom2(rcMomPx[k_index[j]], rcMomPy[k_index[j]], rcMomPz[k_index[j]]);
	      int charge2 = rcCharge[k_index[j]];
	  
	      TVector3 dcaToVtx2 = getDcaToVtx(pos2, mom2, charge2, vertex_rc);
	      if(charge1*charge2<0)
		{
		  // -- only look at unlike-sign pi+k pair

		  float dcaDaughters, cosTheta, decayLength, V0DcaToVtx;
		  TLorentzVector parent = getPairParent(pos1, mom1, charge1, gPionMass,
							pos2, mom2, charge2, gKaonMass,
							vertex_rc,
							dcaDaughters, cosTheta, decayLength, V0DcaToVtx);
		  h3InvMass->Fill(parent.Pt(), parent.Rapidity(), parent.M());
		}
	    }
	}	      
      nevents++;
    }

  TFile *outfile = new TFile(outname.Data(), "recreate");

  hEventStat->Write();
  hMcVtxX->Write();
  hMcVtxY->Write();
  hMcVtxZ->Write();
  hMCD0PtRap->Write();
  h3InvMass->Write();
  outfile->Close();

}

//======================================
TVector3 getDcaToVtx(TVector3 pos, TVector3 mom, int charge, TVector3 vtx)
{
  StPhysicalHelix pHelix(mom, pos, bField * tesla, charge);

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);
  
  pHelix.moveOrigin(pHelix.pathLength(vtx_tmp));
  TVector3 dcaToVtx = pHelix.origin() - vtx_tmp;

  dcaToVtx.SetXYZ(dcaToVtx.x()/millimeter, dcaToVtx.y()/millimeter, dcaToVtx.z()/millimeter);
  
  return dcaToVtx;
}

//======================================
TLorentzVector getPairParent(TVector3 pos1, TVector3 mom1, int charge1, double mass1,
			     TVector3 pos2, TVector3 mom2, int charge2, double mass2,
			     TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &decayLength, float &V0DcaToVtx)
{
  // -- get helix  
  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, charge1);
  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla, charge2);

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);

  // -- find two-track DCA point 
  pair<double, double> const ss = p1Helix.pathLengths(p2Helix);
  TVector3 const p1AtDcaToP2 = p1Helix.at(ss.first);
  TVector3 const p2AtDcaToP1 = p2Helix.at(ss.second);
  
  // -- calculate DCA of particle1 to particle2 at their DCA
  dcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag()/millimeter;
	
  // -- calculate Lorentz vector of particle1-particle2 pair
  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * tesla);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * tesla);
  
  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+mass1*mass2));
  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+mass2*mass2));
  
  TLorentzVector parent = p1FourMom + p2FourMom;
	
  // -- calculate decay vertex (secondary or tertiary)
  TVector3 decayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
	
  // -- calculate pointing angle and decay length with respect to primary vertex
  TVector3 vtxToV0 = decayVertex - vtx_tmp;
  float pointingAngle = vtxToV0.Angle(parent.Vect());
  cosTheta = std::cos(pointingAngle);
  decayLength = vtxToV0.Mag()/millimeter;

  // -- calculate V0 DCA to primary vertex
  V0DcaToVtx = decayLength * std::sin(pointingAngle);

  return parent;
}

