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
#include "TH1F.h"
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

using namespace std;

//==================
void analysis(TString listname="", TString outname="")
{
  if(listname=="" || outname=="")
    {
      printf("[e] Missing input list name or output file name.\n");
      return;
    }

  printf("[i] Input list: %s\n", listname.Data());
  printf("[i] Output file: %s\n", outname.Data());

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

  TH1F *hEventStat = new TH1F("hEventStat", "Event statistics", 8, 0, 8);
  hEventStat->GetXaxis()->SetBinLabel(1, "All events");

  TH1F *hMcMult = new TH1F("hMcMult", "MC multiplicity (|#eta| < 3.5);N_{MC}", 50, 0, 50);
  TH1F *hMcVtxX = new TH1F("hMcVtxX", "x position of MC vertex;x (mm)", 100, -5.05, 4.95);
  TH1F *hMcVtxY = new TH1F("hMcVtxY", "y position of MC vertex;y (mm)", 500, -5.01, 4.99);
  TH1F *hMcVtxZ = new TH1F("hMcVtxZ", "z position of MC vertex;z (mm)", 400, -200, 200);
 
  TH1F *hNRecoJet = new TH1F("hNRecoJet", "Number of reconstructed jets with |#eta| < 2.5 and p_{T} > 2;N", 10, 0, 10);
  TH1F *hRecoJetPt = new TH1F("hRecoJetPt", "Reconstructed jets p_{T};p_{T} (GeV/c)", 50, 0, 50);

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
 
  // reconstructed jets
  TTreeReaderArray<int> recoJetType = {treereader, "ReconstructedChargedJets.type"};
  TTreeReaderArray<float> recoJetE = {treereader, "ReconstructedChargedJets.energy"};
  TTreeReaderArray<float> recoJetMomPx = {treereader, "ReconstructedChargedJets.momentum.x"};
  TTreeReaderArray<float> recoJetMomPy = {treereader, "ReconstructedChargedJets.momentum.y"};
  TTreeReaderArray<float> recoJetMomPz = {treereader, "ReconstructedChargedJets.momentum.z"};
  TTreeReaderArray<float> recoJetMass = {treereader, "ReconstructedChargedJets.mass"};

  int nevents = 0;
  while(treereader.Next())
    {
      nevents++;
      if(nevents%1000==0) printf("\n[i] New event %d\n",nevents);
      hEventStat->Fill(0.5);


      //===== find MC primary vertex
      //===== find 4-momentum of incoming electron and p/A
      int nMCPart = mcPartMass.GetSize();
      TVector3 vertex_mc(-999., -999., -999.);
      for(int imc=0; imc<nMCPart; imc++)
	{
	  int status = mcPartGenStatus[imc];
	  int pdg = mcPartPdg[imc];

	  if(status == 4)
	    {
	      double mcE = sqrt(pow(mcMomPx[imc],2)+pow(mcMomPy[imc],2)+pow(mcMomPz[imc],2)+pow(mcPartMass[imc],2));
	      if(pdg == 11)
		{
		  vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
		}
	    }
	}
      hMcVtxX->Fill(vertex_mc.x());
      hMcVtxY->Fill(vertex_mc.y());
      hMcVtxZ->Fill(vertex_mc.z());

      //===== Loop over MC primary particles
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
      hMcMult->Fill(nMcPart);

      //===== select di-jet events
      int nRecoJets = recoJetType.GetSize();
      int nGoodJet = 0;
      for (int i=0; i<nRecoJets; i++)
	{
	  // select jets within acceptance
	  TVector3 jetMom(recoJetMomPx[i],recoJetMomPy[i],recoJetMomPz[i]);
	  if(TMath::Abs(jetMom.PseudoRapidity()) > 2.5) continue;
	  hRecoJetPt->Fill(jetMom.Perp());
	  if(jetMom.Perp()>2) nGoodJet++;
	}
      hNRecoJet->Fill(nGoodJet);
    }

  TFile *outfile = new TFile(outname.Data(), "recreate");

  hEventStat->Write();
  hMcMult->Write();
  hMcVtxX->Write();
  hMcVtxY->Write();
  hMcVtxZ->Write();

  hNRecoJet->Write();
  hRecoJetPt->Write();
  
  outfile->Close();

  return;
}
