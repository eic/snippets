#include "ElectronID.hh"

#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"

#include <iostream>

ElectronID::ElectronID() {

	mEe = 10.;
	mEh = 100.;
	std::cout << "!!! ElectronID: You have not specified beam energies...defaulting to 10x100 GeV !!!" << std::endl;

	mEoP_min = 0.9;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;


}

ElectronID::ElectronID(double Ee, double Eh) {

	mEe = Ee;
	mEh = Eh;

	mEoP_min = 0.9;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;

}

ElectronID::~ElectronID() {
}


void ElectronID::SetEvent(const podio::Frame* event) {

	mEvent = event;

}


edm4eic::ReconstructedParticleCollection ElectronID::FindScatteredElectron() {

	// Get all the edm4eic objects needed for electron ID
	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
	
	// Create collection for storing scattered electron candidates
	// (subset collection of ReconstructedParticleCollection)
	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
	scatteredElectronCandidates->setSubsetCollection();

	// Loop over all ReconstructedParticles for this event
	for (const auto& reconPart : rcparts) {

		// Require negative particle
		if(reconPart.getCharge() >= 0) continue;

		// Require at least one track and one cluster
		if(reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

		// Calculate rcpart_ member variables for this event
		CalculateParticleValues(reconPart, rcparts);

		// Calculate E/p and isolation fraction for this event
		// Note that the rcpart_ variables are set in CalculateParticleValues
		double recon_EoP = rcpart_sum_cluster_E / edm4hep::utils::magnitude(reconPart.getMomentum());
		double recon_isoE = rcpart_sum_cluster_E / rcpart_isolation_E;
		
		// Apply scattered electron ID cuts
		if(recon_EoP < mEoP_min || recon_EoP > mEoP_max) continue;
		if(recon_isoE < mIsoE) continue;

		// If particle passes cuts, add to output collection
		scatteredElectronCandidates.push_back(reconPart);

	}	

	return scatteredElectronCandidates;

}

edm4hep::MCParticleCollection ElectronID::GetMCElectron() {

	edm4hep::MCParticleCollection meMC;
	meMC->setSubsetCollection();

	auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

	std::vector<edm4hep::MCParticle> mc_electrons;
	
	for(const auto& mcp : mcparts) {
		if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 1) {
			mc_electrons.push_back(mcp);
		}
	}

	// For now, just take first electron
	if(mc_electrons.size() > 0) {
		meMC.push_back(mc_electrons[0]);
	}

	return meMC;
	
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

	edm4hep::MCParticleCollection meMC = GetMCElectron();
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon->setSubsetCollection();

	auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

	for(const auto& assoc : RecoMC) {
		if(assoc.getSim() == meMC[0]) {
			meRecon.push_back(assoc.getRec());
			break;
		}
	}

	return meRecon;

}
	


void ElectronID::CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_sum_cluster_E = 0.;
	rcpart_lead_cluster_E = 0.;
	rcpart_isolation_E = 0.;
	rcpart_deltaH = 0.;

	const edm4eic::Cluster* lead_cluster = nullptr;

	for (const auto& cluster : rcp.getClusters()) {
		rcpart_sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > rcpart_lead_cluster_E) {
			lead_cluster = &cluster;
			rcpart_lead_cluster_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster) return;

	const auto& lead_pos = lead_cluster->getPosition();
	double lead_eta = edm4hep::utils::eta(lead_pos);
	double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);

	for (const auto& other_rcp : rcparts) {
		if (&other_rcp == &rcp) continue;  // Skip the same particle

		for (const auto& other_cluster : other_rcp.getClusters()) {

			const auto& other_pos = other_cluster.getPosition();
			double other_eta = edm4hep::utils::eta(other_pos);
			double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

			double d_eta = other_eta - lead_eta;
			double d_phi = other_phi - lead_phi;

			// Adjust d_phi to be in the range (-pi, pi)
			if (d_phi > M_PI) d_phi-=2*M_PI;
			if (d_phi < -M_PI) d_phi+=2*M_PI;

			double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));

			// Check if the cluster is within the isolation cone
			if (dR < mIsoR) {
				rcpart_isolation_E += other_cluster.getEnergy();
			}
		}
	}

	return;
}

edm4eic::ReconstructedParticle ElectronID::SelectHighestPT(const edm4eic::ReconstructedParticleCollection& ecandidates) {

	edm4eic::ReconstructedParticle erec;
	double max_pT = 0.;
	
	for(auto ecand : ecandidates) {
		double e_pT = edm4hep::utils::magnitudeTransverse(ecand.getMomentum());
		if(e_pT > max_pT) {
			erec = ecand;
			max_pT = e_pT;
		}
	}

	return erec;

}

double ElectronID::GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp) {

	double sum_cluster_E = 0.;
	for (const auto& cluster : rcp.getClusters()) {
		sum_cluster_E += cluster.getEnergy();
	}
	return sum_cluster_E;

}


