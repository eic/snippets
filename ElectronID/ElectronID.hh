#ifndef ELECTRONID_HH
#define ELECTRONID_HH

#include "podio/Frame.h"

#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"
#include "edm4eic/ClusterCollection.h"

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <algorithm>

class ElectronID{

public:

	ElectronID();
	ElectronID(double Ee, double Eh);
       	~ElectronID();

	inline void SetBeamEnergy(double Ee, double Eh) {mEe = Ee; mEh = Eh;}
	inline void SetEoPMin(double eopmin) {mEoP_min = eopmin;}	
	inline void SetDeltaHMin(double deltahmin) {mDeltaH_min = deltahmin;}	
	inline void SetIsolation(double isor, double isoe) {mIsoR = isor; mIsoE = isoe;}	
	inline void SetMinTrackPoints(int minPoints) { minTrackPoints = minPoints; }

	void SetEvent(const podio::Frame* event); 
	void SetBoost(LorentzRotation fboost) { boost = fboost; }

	int Check_eID(edm4eic::ReconstructedParticle e_rec);
	edm4hep::MCParticle GetMC(edm4eic::ReconstructedParticle e_rec);
	edm4eic::ReconstructedParticleCollection FindHadronicFinalState(int object_id);
	edm4eic::ReconstructedParticleCollection FindScatteredElectron();	
	edm4eic::ReconstructedParticleCollection GetTruthReconElectron();	
	edm4hep::MCParticleCollection GetMCElectron();	
	edm4hep::MCParticleCollection GetMCHadronicFinalState();
	edm4eic::ReconstructedParticle SelectHighestPT(const edm4eic::ReconstructedParticleCollection& rcparts);
	double GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp);
	double GetClusterCone(const edm4eic::ReconstructedParticle& rcp, double frac=1);
	PxPyPzEVector GetMomentumVectorFromCluster(const edm4eic::ReconstructedParticle& rcp, double mass);
	double GetClusterTheta(const edm4eic::ReconstructedParticle& rcp);
	void GetEminusPzSum(double &TrackEminusPzSum, double &CalEminusPzSum);
	void CheckClusters();

	double get_mEoP_min() const { return mEoP_min; }
	double get_mEoP_max() const { return mEoP_max; }
	double get_mDeltaH_min() const { return mDeltaH_min; }
	double get_mDeltaH_max() const { return mDeltaH_max; }
	double get_mIsoR() const { return mIsoR; }
	double get_mIsoE() const { return mIsoE; }
	int GetMinTrackPoints() const { return minTrackPoints; }

	// for HFS QA
	vector<double> hfs_dpt;
	vector<double> hfs_dpz;
	vector<double> hfs_de;
	vector<double> hfs_theta;

	struct DetValues {
		int parType; // 0 for dis electron, rest follow pdg code
		int nTrackPoints;
		double recon_EoP;
		double recon_isoE;
	};
	vector<DetValues> det_val;
	vector<DetValues> e_det;
	vector<DetValues> jet_e_det;
	vector<DetValues> pi_det;
	vector<DetValues> else_det;

	double rcpart_n_clusters;

private:

	const podio::Frame* mEvent;

	double mEe;
	double mEh;
	LorentzRotation boost;

	double mEoP_min;
	double mEoP_max;
	double mDeltaH_min;
	double mDeltaH_max;
	double mIsoR;
	double mIsoE;
	int minTrackPoints = 3;
	
	void CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts);
	void CheckSurroundingClusters(const edm4hep::Vector3f& lead_pos,
		const edm4eic::ReconstructedParticleCollection& rcparts);

	double rcpart_sum_cluster_E;
	double rcpart_lead_cluster_E;
	double rcpart_isolation_E;

	int eScatIndex;
};

#endif
