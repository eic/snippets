#include "ElectronID.hh"

#include "edm4hep/utils/vector_utils.h"
#include "edm4eic/ClusterCollection.h"

#include <iostream>

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

ElectronID::ElectronID() {

	mEe = 10.;
	mEh = 100.;
	std::cout << "!!! ElectronID: You have not specified beam energies...defaulting to 10x100 GeV !!!" << std::endl;

	mEoP_min = 0.8;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.7;
	mIsoE = 0.9;

	minTrackPoints = 3;

	mEcone = 0.3;
	mEoEH_min = 0.85;

	mPID_veto = 0.62;

	boost = LorentzRotation(); // Initialize to identity
}

ElectronID::ElectronID(double Ee, double Eh) {

	mEe = Ee;
	mEh = Eh;

	mEoP_min = 0.8;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	// mIsoR = 2; // tested for theta > 150
	mIsoR = 0.7; 
	mIsoE = 0.9;

	minTrackPoints = 3;

	mEcone = 0.3;
	mEoEH_min = 0.85;

	mPID_veto = 0.62;

	boost = LorentzRotation(); // Initialize to identity
}

ElectronID::~ElectronID() {
}

void ElectronID::SetEvent(const podio::Frame* event) {

	// std::cout << "** Setting event in ElectronID... " << std::endl;
	mEvent = event;
	eScatIndex = -1;
	hfs_dpt.clear();
	hfs_dpz.clear();
	hfs_de.clear();
	hfs_theta.clear();
	e_det.clear();
	jet_e_det.clear();
	pi_det.clear();
	else_det.clear();
	det_val.clear();
	// std::cout << "** Event set in ElectronID. " << std::endl;

	meMCValid = false;
	meMC = edm4hep::MCParticleCollection();
	meMC.setSubsetCollection();

	return;
}

edm4hep::MCParticle ElectronID::GetMC(edm4eic::ReconstructedParticle e_rec) {

	// std::cout << "Getting MC particle associated with reconstructed particle..." << std::endl;

	const auto& RecoMC = static_cast<const edm4eic::MCRecoParticleAssociationCollection&>(*(mEvent->get("ReconstructedParticleAssociations")));
	for(const auto& assoc : RecoMC) {
		if(assoc.getRec() == e_rec)
			return assoc.getSim();
	}

	return edm4hep::MCParticle();
}

int ElectronID::Check_eID(edm4eic::ReconstructedParticle e_rec) {

	const auto& meMC = GetMCElectron();
	if ( meMC.size() == 0 )
		return 86; // No MC electron found

	const auto& RecoMC = static_cast<const edm4eic::MCRecoParticleAssociationCollection&>(*(mEvent->get("ReconstructedParticleAssociations")));
	for(const auto& assoc : RecoMC) {
		if(assoc.getRec() == e_rec)
		{
		 	if (assoc.getSim() == meMC[0]) 
				return 0;
			else
				return assoc.getSim().getPDG();
		}
	}

	return 86;
}

void ElectronID::CheckClusters() {

	const auto& EcalEndcapNClusters = static_cast<const edm4hep::ClusterCollection&>(*(mEvent->get("EcalEndcapNClusters")));
	const auto& EcalBarrelScFiClusters = static_cast<const edm4hep::ClusterCollection&>(*(mEvent->get("EcalBarrelScFiClusters")));
	const auto& EcalEndcapPClusters = static_cast<const edm4hep::ClusterCollection&>(*(mEvent->get("EcalEndcapPClusters")));

	std::cout << " Number of clusters in EcalEndcapN: " << EcalEndcapNClusters.size() << std::endl;
	std::cout << " Number of clusters in EcalBarrelScFi: " << EcalBarrelScFiClusters.size() << std::endl;
	std::cout << " Number of clusters in EcalEndcapP: " << EcalEndcapPClusters.size() << std::endl;

	return;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindHadronicFinalState(int object_id) {

	// edm4eic::HadronicFinalStateCollection meRecon;
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon.setSubsetCollection();

	// auto& rcparts = mEvent->get<edm4eic::HadronicFinalStateCollection>("HadronicFinalState");
	const auto& rcparts = static_cast<const edm4eic::ReconstructedParticleCollection&>(*(mEvent->get("ReconstructedParticles")));

	for(const auto& mcp : rcparts) {
		if ( mcp.getObjectID().index != object_id )
			meRecon.push_back(mcp);
	}

	return meRecon;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindScatteredElectron() {

	// std::cout << "\nFinding scattered electron candidates..." << std::endl;

	// Get all the edm4eic objects needed for electron ID
	const auto& rcparts = static_cast<const edm4eic::ReconstructedParticleCollection&>(*(mEvent->get("ReconstructedParticles")));
	
	// Create collection for storing scattered electron candidates
	// (subset collection of ReconstructedParticleCollection)
	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
	scatteredElectronCandidates.setSubsetCollection();

	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates_noEoP;
	scatteredElectronCandidates_noEoP.setSubsetCollection();

	// Loop over all ReconstructedParticles for this event
	for (const auto& reconPart : rcparts) {

		// Require at least one track and one cluster
		if(reconPart.getClusters().size() > 0 && reconPart.getTracks().size() > 0) 
		{
			// Require negative particle
			if(reconPart.getCharge() >= 0) 
				continue;

			int n_track_points = reconPart.getTracks()[0].measurements_size();

			// Calculate rcpart_ member variables for this event
			CalculateParticleValues(reconPart, rcparts);

			double trackTheta = edm4hep::utils::anglePolar(reconPart.getMomentum())*(180./M_PI);
			double clusterTheta = GetMomentumVectorFromCluster(reconPart, 0).Theta()*(180./M_PI);

			// Calculate isolation fraction for this event
			double recon_isoE = rcpart_sum_cluster_E / rcpart_isolation_E;

			// Calculate E/p for this event
			double recon_EoP = rcpart_sum_cluster_E / edm4hep::utils::magnitude(reconPart.getMomentum());

			// Calculate E/E+H for this event
			double recon_EoEH = rcpart_sum_cluster_E / (rcpart_sum_cluster_E + rcpart_sum_cluster_H);

			// pID detector look up
			int recon_pID = reconPart.getPDG();
			double pIDgoodness = reconPart.getGoodnessOfPID();

			double Le = 0.0, Lpi = 0.0, Lk = 0.0, Lp = 0.0;

			for (const auto& pid : reconPart.getParticleIDs()) 
			{
				const int apdg = std::abs(pid.getPDG());
				const double L = pid.getLikelihood();
				if (apdg == 11) Le = std::max(Le, L);
				else if (apdg == 211) Lpi = std::max(Lpi, L);
				else if (apdg == 321) Lk = std::max(Lk, L);
				else if (apdg == 2212) Lp = std::max(Lp, L);
			}

			double Lhad = std::max({Lpi, Lk, Lp});
			double Lhadsum = Lpi + Lk + Lp;
			double eFrac = (Le + Lhadsum > 0.0) ? Le / (Le + Lhadsum) : 0.0;
			double hFrac = (Le + Lhad > 0.0) ? Lhad / (Le + Lhad) : 0.0;

			// find real PDG of the particle
			int mc_PDG = Check_eID(reconPart);

			det_val.push_back({mc_PDG, n_track_points, recon_EoP, recon_isoE, recon_EoEH, recon_pID, eFrac, hFrac});

			if ( n_track_points < minTrackPoints ) 
				continue;

			if(recon_isoE < mIsoE) 
				continue;

			if (hFrac > mPID_veto) 
				continue;

			// if (recon_EoEH < mEoEH_min) 
			// 	continue;

			if ( recon_EoP > mEoP_min && recon_EoP < mEoP_max )
			{
				scatteredElectronCandidates.push_back(reconPart);
				continue;
			}
			else if ( (trackTheta > 158 && trackTheta < 162) || (clusterTheta > 22 && clusterTheta < 33) ) // at small angles, it works better to use cluster to estimate theta here
			{
				scatteredElectronCandidates_noEoP.push_back(reconPart);
				continue;
			}
		}
	}	

	// If EoP is found use that, otherwise loosen cuts 
	if (scatteredElectronCandidates.size() > 0)
		return scatteredElectronCandidates;
	else
		return scatteredElectronCandidates_noEoP;
}

edm4hep::MCParticleCollection ElectronID::GetMCHadronicFinalState() {

	edm4hep::MCParticleCollection mhMC;
	mhMC.setSubsetCollection();

	const auto& mcparts = static_cast<const edm4hep::MCParticleCollection&>(*(mEvent->get("MCParticles")));

	const auto& meMC = GetMCElectron();

	bool found_scattered_e = false; 
	for(const auto& mcp : mcparts) {
		if (mcp.getGeneratorStatus() == 1)
		{
			if ( meMC.size() == 0 )
				mhMC.push_back(mcp);
			else if (mcp.getObjectID().index != meMC[0].getObjectID().index ) 
				mhMC.push_back(mcp);	
		}
	}

	return mhMC;
}

const edm4hep::MCParticleCollection& ElectronID::GetMCElectron() {
// updated version: search for status 1 electrons with the shortest generation path to a beam electron. Should also work even if there are multiple beam electrons.

    if (meMCValid)
        return meMC;

    meMC = edm4hep::MCParticleCollection();
    meMC.setSubsetCollection();

    const auto& mcparts = static_cast<const edm4hep::MCParticleCollection&>(*(mEvent->get("MCParticles")));

    std::vector<int> beam_electron_indices;
    for (const auto& p : mcparts) {
        if (p.getPDG() == 11 && p.getGeneratorStatus() == 4)
            beam_electron_indices.push_back(p.getObjectID().index);
    }

    if (beam_electron_indices.empty()) {
        meMCValid = true;
        return meMC;
    }

    std::unordered_map<int, std::pair<int, std::vector<edm4hep::MCParticle>>> per_beam;
    for (int idx : beam_electron_indices)
        per_beam[idx] = {INT_MAX, {}};

    for (const auto& cand : mcparts) {
        if (cand.getPDG() != 11 || cand.getGeneratorStatus() != 1)
            continue;

        std::queue<std::pair<edm4hep::MCParticle, int>> frontier;
        std::unordered_set<int> visited;

        visited.insert(cand.getObjectID().index);
        frontier.emplace(cand, 0);

        while (!frontier.empty()) {
            auto [cur, depth] = frontier.front();
            frontier.pop();

            const int cur_idx = cur.getObjectID().index;

            if (cur.getPDG() == 11 && cur.getGeneratorStatus() == 4
                && per_beam.count(cur_idx)) {

                auto& [min_gen, cands] = per_beam[cur_idx];
                if (depth < min_gen) {
                    min_gen = depth;
                    cands.clear();
                    cands.push_back(cand);
                } else if (depth == min_gen) {
                    cands.push_back(cand);
                }
                break;
            }

            for (auto it = cur.parents_begin(); it != cur.parents_end(); ++it) {
                const auto parent = *it;
                const int parent_idx = parent.getObjectID().index;
                if (parent_idx >= 0 && visited.insert(parent_idx).second)
                    frontier.emplace(parent, depth + 1);
            }
        }
    }

    std::vector<int> sorted_beam = beam_electron_indices;
    std::sort(sorted_beam.begin(), sorted_beam.end());

    for (int beam_idx : sorted_beam) {
        const auto& [min_gen, cands] = per_beam[beam_idx];
        for (const auto& c : cands)
            meMC.push_back(c);
    }

    meMCValid = true;
    return meMC;
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

	// cout << "New process " << endl;

	const auto& meMC = GetMCElectron();
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon.setSubsetCollection();

	if ( meMC.size() == 0 )
		return meRecon; // No MC electron found

	const auto& RecoMC = static_cast<const edm4eic::MCRecoParticleAssociationCollection&>(*(mEvent->get("ReconstructedParticleAssociations")));

	for(const auto& assoc : RecoMC) 
	{
		auto e_candidat = assoc.getSim();

		if(assoc.getSim() == meMC[0]) {
			meRecon.push_back(assoc.getRec());
			// break;
		}
	}

	const auto& rcparts = static_cast<const edm4eic::ReconstructedParticleCollection&>(*(mEvent->get("ReconstructedParticles")));
	edm4hep::Vector3f lead_pos(meMC[0].getEndpoint().x, meMC[0].getEndpoint().y, meMC[0].getEndpoint().z);
	CheckSurroundingClusters(lead_pos, rcparts);

	// std::cout << meRecon.size() << " reconstructed particles associated with the MC electron." << std::endl;

	return meRecon;
}
	
void ElectronID::CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_sum_cluster_E = 0.;
	rcpart_lead_cluster_E = 0.;
	rcpart_isolation_E = 0.;
	rcpart_sum_cluster_H = 0.;

	const edm4eic::Cluster* lead_cluster = nullptr;

	for (const auto& cluster : rcp.getClusters()) {
		rcpart_sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > rcpart_lead_cluster_E) {
			lead_cluster = &cluster;
			rcpart_lead_cluster_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster) return;

	CheckSurroundingClusters(lead_cluster->getPosition(), rcparts);

	return;
}

void ElectronID::CheckSurroundingClusters(const edm4hep::Vector3f& lead_pos,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_isolation_E = 0.;
	rcpart_n_clusters = 0;

	double lead_eta = edm4hep::utils::eta(lead_pos);
	double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);

	for (const auto& other_rcp : rcparts) {
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
				rcpart_n_clusters ++;
			}
		}
	}

	// HCal
	const auto& hcal_barrel = static_cast<const edm4eic::ClusterCollection&>(*(mEvent->get("HcalBarrelClusters")));
	const auto& hcal_f_endcap = static_cast<const edm4eic::ClusterCollection&>(*(mEvent->get("HcalEndcapNClusters")));
	const auto& hcal_b_endcap = static_cast<const edm4eic::ClusterCollection&>(*(mEvent->get("LFHCALClusters")));
	
	for ( const auto& hcal_cluster : hcal_barrel) {
		const auto& other_pos = hcal_cluster.getPosition();
		double other_eta = edm4hep::utils::eta(other_pos);
		double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

		double d_eta = other_eta - lead_eta;
		double d_phi = other_phi - lead_phi;

		// Adjust d_phi to be in the range (-pi, pi)
		if (d_phi > M_PI) d_phi-=2*M_PI;
		if (d_phi < -M_PI) d_phi+=2*M_PI;

		double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));

		// Check if the cluster is within the isolation cone
		if (dR < mEcone) {
			rcpart_sum_cluster_H += hcal_cluster.getEnergy();
			rcpart_n_clusters ++;
		}
	}

	for ( const auto& hcal_cluster : hcal_f_endcap) {
		const auto& other_pos = hcal_cluster.getPosition();
		double other_eta = edm4hep::utils::eta(other_pos);
		double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

		double d_eta = other_eta - lead_eta;
		double d_phi = other_phi - lead_phi;

		// Adjust d_phi to be in the range (-pi, pi)
		if (d_phi > M_PI) d_phi-=2*M_PI;
		if (d_phi < -M_PI) d_phi+=2*M_PI;

		double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));

		// Check if the cluster is within the isolation cone
		if (dR < mEcone) {
			rcpart_sum_cluster_H += hcal_cluster.getEnergy();
			rcpart_n_clusters ++;
		}
	}

	for ( const auto& hcal_cluster : hcal_b_endcap) {
		const auto& other_pos = hcal_cluster.getPosition();
		double other_eta = edm4hep::utils::eta(other_pos);
		double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

		double d_eta = other_eta - lead_eta;
		double d_phi = other_phi - lead_phi;

		// Adjust d_phi to be in the range (-pi, pi)
		if (d_phi > M_PI) d_phi-=2*M_PI;
		if (d_phi < -M_PI) d_phi+=2*M_PI;

		double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));

		// Check if the cluster is within the isolation cone
		if (dR < mEcone) {
			rcpart_sum_cluster_H += hcal_cluster.getEnergy();
			rcpart_n_clusters ++;
		}
	}

	return;
}

void ElectronID::GetEminusPzSum(double &TrackEminusPzSum, double &CalEminusPzSum) {

	const auto& rcparts = static_cast<const edm4eic::ReconstructedParticleCollection&>(*(mEvent->get("ReconstructedParticles")));

	for (const auto& reconPart : rcparts) {

		// Require at least one track and one cluster
		// if(reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

		if(reconPart.getTracks().size() > 0 )
		{
			int n_track_points = reconPart.getTracks()[0].measurements_size();
			if ( n_track_points < minTrackPoints ) continue;

			PxPyPzEVector vT(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, reconPart.getEnergy());
			vT = boost(vT);
			TrackEminusPzSum += (vT.E() - vT.Pz());

			if ( reconPart.getClusters().size() > 0 )
			{
				PxPyPzEVector vC(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, GetCalorimeterEnergy(reconPart));
				// PxPyPzEVector vC = GetMomentumVectorFromCluster(reconPart, reconPart.getMass());
				// if ( reconPart.getTracks().size() > 0 && n_track_points >= minTrackPoints )
					//  vC.SetPxPyPzE(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, GetCalorimeterEnergy(reconPart));
				vC = boost(vC);
				CalEminusPzSum += (vC.E() - vC.Pz());
			}
		}
		else if (reconPart.getClusters().size() > 0 )
		{
			PxPyPzEVector vC = GetMomentumVectorFromCluster(reconPart, 0);
			vC = boost(vC);
			CalEminusPzSum += (vC.E() - vC.E()*std::cos(vC.Theta()));
		}
	}

	// std::cout << " recon E - Pz sum: " << reconEminusPzSum << std::endl;

	return;
}

edm4eic::ReconstructedParticle ElectronID::SelectHighestPT(const edm4eic::ReconstructedParticleCollection& ecandidates) {

	edm4eic::ReconstructedParticle erec;
	double max_pT = 0.;
	
	for(const auto& ecand : ecandidates) {
		double e_pT = ecand.getTracks().size() > 0 ? edm4hep::utils::magnitudeTransverse(ecand.getMomentum()) : GetMomentumVectorFromCluster(ecand, 0.000511).Pt();
		if(e_pT > max_pT) {
			erec = ecand;
			max_pT = e_pT;
		}
	}

	return erec;
}

double ElectronID::GetClusterCone(const edm4eic::ReconstructedParticle& rcp, double frac=1)
{
	const edm4eic::Cluster* lead_cluster = nullptr;
	double sum_cluster_E = 0.;
    double lead_cluster_E = 0.;
	
	for (const auto& cluster : rcp.getClusters()) {
		sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > lead_cluster_E) {
			lead_cluster = &cluster;
			lead_cluster_E = cluster.getEnergy();
		}
	}

	if (!lead_cluster)
        return -1.;

	std::vector<std::pair<double, double>> dr_de;

	for (const auto& cluster : rcp.getClusters())
	{
		if (lead_cluster == &cluster)
			continue;

		double d_eta = edm4hep::utils::eta(cluster.getPosition()) - edm4hep::utils::eta(lead_cluster->getPosition());
		double d_phi = edm4hep::utils::angleAzimuthal(cluster.getPosition()) - edm4hep::utils::angleAzimuthal(lead_cluster->getPosition());

		if (d_phi > M_PI) d_phi -= 2 * M_PI;
		if (d_phi < -M_PI) d_phi += 2 * M_PI;

		double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));
		dr_de.emplace_back(dR, cluster.getEnergy());
	}

	std::sort(dr_de.begin(), dr_de.end(), [](const auto& left, const auto& right) {return left.first < right.first;});

	return dr_de.back().first; // Return the dR of the farthest cluster 
}

double ElectronID::GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp) {

	double sum_cluster_E = 0.;
	for (const auto& cluster : rcp.getClusters()) {
		sum_cluster_E += cluster.getEnergy();
	}
	return sum_cluster_E;

}

double ElectronID::GetClusterTheta(const edm4eic::ReconstructedParticle& rcp) {

	const edm4eic::Cluster* lead_cluster = nullptr;
	double lead_E = 0.;

	for (const auto& cluster : rcp.getClusters()) {
		if(cluster.getEnergy() > lead_E) {
			lead_cluster = &cluster;
			lead_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster) return -999.;

	return lead_cluster->getIntrinsicTheta();
}

PxPyPzEVector ElectronID::GetMomentumVectorFromCluster(const edm4eic::ReconstructedParticle& rcp, double mass) {

	const edm4eic::Cluster* lead_cluster = nullptr;
	double sum_cluster_E = 0.;
	double lead_E = 0.;

	for (const auto& cluster : rcp.getClusters()) {
		sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > lead_E) {
			lead_cluster = &cluster;
			lead_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster)
		cout << "Warning: No clusters found for this reconstructed particle!" << endl;

	if(!lead_cluster) return PxPyPzEVector(0, 0, 0, 0);

	double p = std::sqrt(std::pow(sum_cluster_E, 2) - std::pow(mass, 2));
	double px = p * (lead_cluster->getPosition().x / edm4hep::utils::magnitude(lead_cluster->getPosition()));
	double py = p * (lead_cluster->getPosition().y / edm4hep::utils::magnitude(lead_cluster->getPosition()));
	double pz = p * (lead_cluster->getPosition().z / edm4hep::utils::magnitude(lead_cluster->getPosition()));
	
	return PxPyPzEVector(px, py, pz, sum_cluster_E);
}
