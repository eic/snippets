#include <fstream>
#include <map>
#include "TROOT.h"
#include "TClass.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"
#include <string>

using namespace std;

void find_ScatElectron(TString listname = "file.list")
{
  TChain *chain = new TChain("events");
  TTree* meta = NULL;

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
	      if(nfiles==0) meta  = (TTree*)ftmp->Get("podio_metadata");
	      nfiles++;
	    }
	}
    }
  inputstream->close();
  printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());

  TTreeReader treereader(chain);

  // scattered electron
  TTreeReaderArray<int> mc_ScatElec = {treereader, "MCScatteredElectrons_objIdx.index"};
  
  // Reconstructed particles
  TTreeReaderArray<float> rcMomPx = {treereader, "ReconstructedParticles.momentum.x"};
  TTreeReaderArray<float> rcMomPy = {treereader, "ReconstructedParticles.momentum.y"};
  TTreeReaderArray<float> rcMomPz = {treereader, "ReconstructedParticles.momentum.z"};
  TTreeReaderArray<float> rcCharge = {treereader, "ReconstructedParticles.charge"};
  TTreeReaderArray<int> rcPdg = {treereader, "ReconstructedParticles.PDG"};

  TTreeReaderArray<int> rc_index = {treereader, "_ReconstructedParticleAssociations_rec.index"};
  TTreeReaderArray<int> mc_index = {treereader, "_ReconstructedParticleAssociations_sim.index"}; 

  TTreeReaderArray<unsigned int> clusters_begin = {treereader, "ReconstructedParticles.clusters_begin"};
  TTreeReaderArray<unsigned int> clusters_end = {treereader, "ReconstructedParticles.clusters_end"};
  TTreeReaderArray<int> clusters_index = {treereader, "_ReconstructedParticles_clusters.index"};
  TTreeReaderArray<unsigned int> clusters_collectionID = {treereader, "_ReconstructedParticles_clusters.collectionID"};

  TTreeReaderArray<unsigned int> tracks_begin = {treereader, "ReconstructedParticles.tracks_begin"};
  TTreeReaderArray<unsigned int> tracks_end = {treereader, "ReconstructedParticles.tracks_end"};
  TTreeReaderArray<int> tracks_index = {treereader, "_ReconstructedParticles_tracks.index"};
  TTreeReaderArray<unsigned int> tracks_collectionID = {treereader, "_ReconstructedParticles_tracks.collectionID"};

  // cluster collections
  map<string, int> nameToIndex;
  const int nClusterCollections = 3;
  string cluster_names[nClusterCollections] = {"EcalBarrelScFiClusters", "EcalEndcapNClusters", "EcalEndcapPClusters"};
  TTreeReaderArray<float>* clusterE[nClusterCollections], *clusterX[nClusterCollections], *clusterY[nClusterCollections], *clusterZ[nClusterCollections];
  for(int i=0; i<nClusterCollections; i++)
    {
      nameToIndex[cluster_names[i]] = i;
      clusterE[i] = new TTreeReaderArray<float>(treereader, Form("%s.energy", cluster_names[i].c_str()));
      clusterX[i] = new TTreeReaderArray<float>(treereader, Form("%s.position.x", cluster_names[i].c_str()));
      clusterY[i] = new TTreeReaderArray<float>(treereader, Form("%s.position.y", cluster_names[i].c_str()));
      clusterZ[i] = new TTreeReaderArray<float>(treereader, Form("%s.position.z", cluster_names[i].c_str()));
    }

  // Get the relevant leafs from podio_metadata tree
  TLeaf* idLeaf = meta->FindLeaf("events___CollectionTypeInfo.collectionID");
  TLeafC* nameLeafBase = (TLeafC*)meta->FindLeaf("events___CollectionTypeInfo.name");
  meta->GetEntry(0);
  Int_t nCollections = idLeaf->GetLen();
  const char* tmpFile = "podio_metadata.txt";
  {
    TRedirectOutputGuard guard(tmpFile, "w");
    nameLeafBase->PrintValue(nCollections);
  }
  ifstream infile(tmpFile);
  string out;
  getline(infile, out);
  infile.close();

  // parse the out string to save into a map for lookup
  string prefix = "events___CollectionTypeInfo.name =";
  if (out.rfind(prefix, 0) == 0)
    {
      out = out.substr(prefix.size());
    }
  vector<string> names;
  stringstream ss(out);
  string item;

  // split by comma
  while (getline(ss, item, ','))
    {
      // trim leading space
      if (!item.empty() && item[0] == ' ')
	item.erase(0, 1);
      names.push_back(item);
    }

  if (names.size() != nCollections)
    {
      cerr << "Mismatch: names =" << names.size() << " vs IDs =" << nCollections << endl;
      return;
    }
  
  // build map
  unordered_map<UInt_t, string> idToName;
  for (size_t i = 0; i < names.size(); ++i)
    {
      UInt_t cid = static_cast<UInt_t>(idLeaf->GetValue(i));
      idToName[cid] = names[i];
    }

  // Cuts
  const double mIsoR = 0.4;
  const double mEoP_min = 0.8;
  const double mEoP_max = 1.2;
  const double mIsoE = 0.9;
  
  // find scattered electron
  int nevents = 0;
  while(treereader.Next())
    {
      //if(nevents>11) return;
      if(nevents%100==0) printf("\n[i] New event %d\n",nevents);
      printf("\n[i] New event %d\n",nevents);
      nevents++;

      // loop over reconstructed particles
      int scat_elec_index = -1;
      double scat_elec_pt = -1;
      for(int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{
	  //cout << rcMomPx[rc_index] << "  " << chain->FindLeaf("ReconstructedParticles.momentum.x")->GetValue(rc_index) << endl;
	  int nclusters = clusters_end[rc_index] - clusters_begin[rc_index];
	  int ntracks = tracks_end[rc_index] - tracks_begin[rc_index];
	  
	  // at least one cluster and one track
	  if(nclusters == 0 || ntracks == 0) continue;

	  // Negative charge
	  if(rcCharge[rc_index] >= 0) continue;

	  // Calculate associate cluster energy
	  double rcpart_sum_cluster_E = 0;
	  double rcpart_isolation_E = 0;
	  double lead_eta = -999;
	  double lead_phi = -999;
	  double lead_E = 0;
	  for(int ic = clusters_begin[rc_index]; ic < clusters_end[rc_index]; ic++)
	    {
	      int cluster_index = clusters_index[ic];
	      int collectionID = clusters_collectionID[ic];
	      string collection_name = idToName[collectionID];
	      if(!nameToIndex.count(collection_name))
		{
		  cout << "[e] Error: unknown collection name: " << collection_name.c_str() << endl;
		  cout << "[e] Please update the map" << endl;
		  return;
		}
	      int collection_index = nameToIndex[collection_name];
	      double energy = clusterE[collection_index]->At(cluster_index);
	      rcpart_sum_cluster_E += energy;
	      if(energy>lead_E)
		{
		  lead_E = energy;
		  TVector3 pos(clusterX[collection_index]->At(cluster_index),
			       clusterY[collection_index]->At(cluster_index),
			       clusterZ[collection_index]->At(cluster_index));
		  lead_eta = pos.Eta();
		  lead_phi = pos.Phi();
		}
	    }

	  // loop again to get isolation energy
	  for(int j=0; j<rcMomPx.GetSize(); j++)
	    {
	      for(int ic = clusters_begin[j]; ic < clusters_end[j]; ic++)
		{
		  int cluster_index = clusters_index[ic];
		  int collectionID = clusters_collectionID[ic];
		  string collection_name = idToName[collectionID];
		  if(!nameToIndex.count(collection_name))
		    {
		      cout << "[e] Error: unknown collection name: " << collection_name.c_str() << endl;
		      cout << "[e] Please update the map" << endl;
		      return;
		    }
		  int collection_index = nameToIndex[collection_name];
		  TVector3 pos(clusterX[collection_index]->At(cluster_index),
			       clusterY[collection_index]->At(cluster_index),
			       clusterZ[collection_index]->At(cluster_index));

		  double d_eta = pos.Eta() - lead_eta;
		  double d_phi = pos.Phi() - lead_phi;
		  if (d_phi > TMath::Pi()) d_phi-=2*TMath::Pi();
		  if (d_phi < -TMath::Pi()) d_phi+=2*TMath::Pi();

		  double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));
		  // Check if the cluster is within the isolation cone
		  if (dR < mIsoR)
		    {
		      rcpart_isolation_E += clusterE[collection_index]->At(cluster_index);
		    }
		}
	    }

	  //printf("[i] rcpart_sum_cluster_E = %f, rcpart_isolation_E = %f\n", rcpart_sum_cluster_E,  rcpart_isolation_E);

	  // Apply scattered electron ID cuts
	  TVector3 part_mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
	  double recon_EoP = rcpart_sum_cluster_E / part_mom.Mag();
	  double recon_isoE = rcpart_sum_cluster_E / rcpart_isolation_E;
	  if(recon_EoP < mEoP_min || recon_EoP > mEoP_max) continue;
	  if(recon_isoE < mIsoE) continue;

	  // in case of multiple candidate, select the one with largest pt
	  if(part_mom.Pt()>scat_elec_pt)
	    {
	      scat_elec_pt = part_mom.Pt();
	      scat_elec_index = rc_index;
	    } 
	}

      // Sanity check
      if(scat_elec_index>=0)
	{
	  // find corresponding MC particle
	  int scat_elec_mc_index = -1;
	  for(unsigned int i=0; i<rc_index.GetSize(); i++)
	    {
	      if(rc_index[i] == scat_elec_index)
		{
		  scat_elec_mc_index = mc_index[i];
		  break;
		}
	    }

	  // true MC electron index
	  int scat_elec_mc_truth = mc_ScatElec[0];
	  
	  printf("[i] Found a scattered electron candidate: index = %d, pt = %f, mc index = %d, truth = %d\n", scat_elec_index, scat_elec_pt, scat_elec_mc_index, scat_elec_mc_truth);
	}
    }
}
