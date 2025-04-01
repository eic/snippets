// Data model headers
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTFrameReader.h"

// ROOT headers
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// Analysis headers
#include "InclusiveSkim.h"
#include "ElectronID.cc"

void InclusiveSkim() {

	double Ee = 10.;
	double Eh = 100.;

	vector<std::string> inFiles = {"pythia8NCDIS_10x100_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0001.eicrecon.tree.edm4eic.root"};

	auto reader = podio::ROOTFrameReader();
	reader.openFiles(inFiles);

	ElectronID* eFinder = new ElectronID(Ee, Eh);

	TString outFileName = Form("inclusive_skim_%.0fx%.0fGeV.root", Ee, Eh);
	CreateOutputTree(outFileName);


	for(size_t ev = 0; ev < reader.getEntries("events"); ev++) {

		if(ev%1000==0) cout << "Event " << ev << "/" << reader.getEntries("events") << endl;

		ResetVariables();

		const auto event = podio::Frame(reader.readNextEntry("events"));
		eFinder->SetEvent(&event);

		edm4hep::MCParticleCollection e_mc = eFinder->GetMCElectron();

		// Skip events without scattered MC electron
		if(e_mc.size() == 0) continue;

		// Find scattered electrons
		auto e_candidates = eFinder->FindScatteredElectron();
		edm4eic::ReconstructedParticle e_rec;	
		// If there are multiple candidates, 
		// select one with highest pT
		if(e_candidates.size() > 0) {
			positive_eID = 1;			
			e_rec = eFinder->SelectHighestPT(e_candidates);
		}

		// Get momentum vector elements from MC electron
		mc_p = edm4hep::utils::magnitude(e_mc[0].getMomentum());
		mc_eta = edm4hep::utils::eta(e_mc[0].getMomentum());
		mc_phi = edm4hep::utils::angleAzimuthal(e_mc[0].getMomentum());

		// Calculate kinematic variables using MC electron
		TLorentzVector kprime;
		kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, Me);
		CalculateElectronKinematics(Ee, Eh, kprime, mc_xB, mc_Q2, mc_W, mc_y, mc_nu);

		// If electron was identified, get vector elements 
		// and kinematic variables from reconstructed particle
		if(positive_eID) {
			track_p = edm4hep::utils::magnitude(e_rec.getMomentum());
			track_eta = edm4hep::utils::eta(e_rec.getMomentum());
			track_phi = edm4hep::utils::angleAzimuthal(e_rec.getMomentum());

			kprime.SetXYZM(e_rec.getMomentum().x, e_rec.getMomentum().y, e_rec.getMomentum().z, Me);
			CalculateElectronKinematics(Ee, Eh, kprime, e_track_xB, e_track_Q2, e_track_W, e_track_y, e_track_nu);

			// Calculate kinematic variables with calorimeter energy
			TVector3 e3v = kprime.Vect();
			double track_clust_E = eFinder->GetCalorimeterEnergy(e_rec);
			e3v.SetMag(track_clust_E);
			kprime.SetVectM(e3v, Me);
			CalculateElectronKinematics(Ee, Eh, kprime, e_clust_xB, e_clust_Q2, e_clust_W, e_clust_y, e_clust_nu);
		}	


		outTree->Fill();

	}

	outFile->cd();
	outTree->Write(outTree->GetName(), 2);

}


void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T", "Reconstructed and generated information from EICRecon");

	outTree->Branch("positive_eID",		&positive_eID);

	outTree->Branch("mc_p",			&mc_p);
	outTree->Branch("mc_eta",		&mc_eta);
	outTree->Branch("mc_phi",		&mc_phi);

	outTree->Branch("track_p",		&track_p);	
	outTree->Branch("track_eta",		&track_eta);	
	outTree->Branch("track_phi",		&track_phi);

	outTree->Branch("mc_xB",		&mc_xB);
	outTree->Branch("mc_Q2",		&mc_Q2);
	outTree->Branch("mc_W",			&mc_W);
	outTree->Branch("mc_y",			&mc_y);
	outTree->Branch("mc_nu",		&mc_nu);

	outTree->Branch("e_track_xB",		&e_track_xB);
	outTree->Branch("e_track_Q2",		&e_track_Q2);
	outTree->Branch("e_track_W",		&e_track_W);
	outTree->Branch("e_track_y",		&e_track_y);
	outTree->Branch("e_track_nu",		&e_track_nu);

	outTree->Branch("e_clust_xB",		&e_clust_xB);
	outTree->Branch("e_clust_Q2",		&e_clust_Q2);
	outTree->Branch("e_clust_W",		&e_clust_W);
	outTree->Branch("e_clust_y",		&e_clust_y);
	outTree->Branch("e_clust_nu",		&e_clust_nu);

}

void ResetVariables() {

	positive_eID = 0;

	mc_p = -999;
	mc_eta = -999;
	mc_phi = -999;
	
	track_p = -999;
	track_eta = -999;
	track_phi = -999;

	mc_xB = -999;
	mc_Q2 = -999;
	mc_W = -999;
	mc_y = -999;
	mc_nu = -999;
	      
	e_track_xB = -999;
	e_track_Q2 = -999;
	e_track_W = -999;
	e_track_y = -999;
	e_track_nu = -999;

	e_clust_xB = -999;
	e_clust_Q2 = -999;
	e_clust_W = -999;
	e_clust_y = -999;
	e_clust_nu = -999;


}

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, float& xB, float& Q2, float& W, float& y, float& nu) {

		TLorentzVector ki; ki.SetXYZM(0., 0., -fEe, Me);
		TLorentzVector P = GetHadronBeam(fEh);
		TLorentzVector q = ki - kf;
		Q2 = -(q.Dot(q));
		nu = (q.Dot(P))/Mp;
		xB = Q2/(2.*Mp*nu);
		y  = (q.Dot(P))/(ki.Dot(P));
		W  = sqrt(Mp*Mp + (2.*Mp*nu) - Q2);		
}


TLorentzVector GetHadronBeam(double fEh) {
 
	TLorentzVector hadron_beam;
	hadron_beam.SetX(fEh*sin(crossing_angle));
	hadron_beam.SetY(0.);
	hadron_beam.SetZ(fEh*cos(crossing_angle));
	hadron_beam.SetE(std::hypot(fEh, Mp));
	return hadron_beam;

}
