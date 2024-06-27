// ----------------------------------------------------------------------------
// 'CheckClusterSplittingProcessor.h'
// Derek Anderson
// 04.11.2024
//
// EICrecon plugin to collect some metrics on cluster splitting
// across several calorimeters.
// ----------------------------------------------------------------------------

#ifndef CHECKCLUSTERSPLITTINGPROCESSOR_H
#define CHECKCLUSTERSPLITTINGPROCESSOR_H

// c++ utilities
#include <map>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <cstdlib>
#include <utility>
#include <algorithm>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TDirectory.h>
// edm4hep types
#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
// edm4eic types
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/TrackPoint.h>
#include <edm4eic/TrackSegment.h>
#include <edm4eic/TrackSegmentCollection.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/ReconstructedParticleCollection.h>
// dd4hep utilities
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
// jana libraries
#include <JANA/JEventProcessorSequentialRoot.h>
// eicrecon services
#include "services/rootfile/RootFile_service.h"
#include "services/geometry/dd4hep/DD4hep_service.h"

// make common namespaces implicit
using namespace std;



// convenience types ----------------------------------------------------------

// binning types
typedef tuple<int64_t, float, float> Binning;
typedef pair<string, Binning> AxisDef;

// histogram types
typedef tuple<string, AxisDef> HistDef1D;
typedef tuple<string, AxisDef, AxisDef> HistDef2D;

// definition vectors
typedef vector<AxisDef> VecAxisDef;
typedef vector<HistDef1D> VecHistDef1D;
typedef vector<HistDef2D> VecHistDef2D;

// other types
typedef pair<float, float> FPair;
typedef vector<edm4eic::TrackPoint> TrkVec;
typedef vector<vector<vector<TH1D*>>> TH1Vec;
typedef vector<vector<vector<TH2D*>>> TH2Vec;



// histogram utilities --------------------------------------------------------

namespace Hist {

  // indices/accessors --------------------------------------------------------

  // calorimeters
  enum Calo {NHCal, NECal, BHCal, BECal, PHCal, PECal};

  // histogram axes
  enum EVar {EvtNCl, EvtPE, EvtLE, EvtTE, EvtLT, EvtLP, EvtTP, EvtNA90, EvtNB90, EvtNR};
  enum TVar {TrkRX, TrkRY, TrkRZ, TrkH, TrkF, TrkP, TrkDR, TrkDRS, TrkEP, TrkDep};
  enum CVar {CalNH, CalNP, CalRX, CalRY, CalRZ, CalH, CalF, CalE, CalEps, CalRho, CalEP, CalSig, CalN90, CalDR90};


  // event histogram accessors
  enum EType {EvtAll};
  enum E1D   {EvtNClust, EvtParEne, EvtLeadEne, EvtTotEne, EvtLeadDivTot, EvtLeadDivPar, EvtTotDivPar, EvtNClAb90, EvtNClBe90, EvtNClRec};

  // projection histogram accessors
  enum TType {TrkAll, TrkMatch};
  enum T1D   {TrkPosX, TrkPosY, TrkPosZ, TrkEta, TrkPhi, TrkMom, TrkDeltaR, TrkDeltaRScale, TrkEOverP, TrkDeposit};
  enum T2D   {TrkPosYvsX, TrkMomVsEta, TrkEtaVsPhi, TrkEPvsDR, TrkEPvsDRSig};

  // cluster histogram accessors
  enum CType {CalAll, CalAbove90, CalBelow90};
  enum C1D   {CalNHit, CalNProj, CalPosX, CalPosY, CalPosZ, CalEta, CalPhi, CalEne, CalFraction, CalPurity, CalEOverP, CalSigma, CalNAdd90, CalDRAdd90};
  enum C2D   {CalPosYvsX, CalEneVsEta, CalEneVsPhi, CalEtaVsPhi, CalPurVsEne, CalN90VsEP, CalDR90VsEP};

  // histogram content definitions ------------------------------------------

  // event histogram content
  struct EvtContent {

    // data
    int64_t nClust = 0;
    int64_t nCl90  = 0;
    int64_t nClB90 = 0;
    int64_t nClRec = 0;
    double  ePar   = 0.;
    double  eLead  = 0.;
    double  eTot   = 0.;
    double  lOverT = 0.;
    double  lOverP = 0.;
    double  tOverP = 0.;

  };  // end EvtContent

  // projection histogram content
  struct TrkContent {

    // data
    double rx      = 0.;
    double ry      = 0.;
    double rz      = 0.;
    double eta     = 0.;
    double phi     = 0.;
    double p       = 0.;
    double dr      = 0.;
    double drSig   = 0.;
    double eOverP  = 0.;
    double deposit = 0.;

    // get info from projection
    void GetInfo(edm4eic::TrackPoint& point) {
      rx  = point.position.x;
      ry  = point.position.y;
      rz  = point.position.z;
      eta = edm4hep::utils::eta( point.position );
      phi = atan2( point.position.y, point.position.x );
      p   = edm4hep::utils::magnitude( point.momentum );
      return;
    }  // end 'GetInfo(edm4eic::TrackPoint&)'
  };  // end TrkContent

  // cluster histogram content
  struct CalContent {

    // data
    int64_t nHit    = 0;
    int64_t nProj   = 0;
    int64_t nAdd90  = 0;
    double  rx      = 0.;
    double  ry      = 0.;
    double  rz      = 0.;
    double  eta     = 0.;
    double  phi     = 0.;
    double  ene     = 0.;
    double  frac    = 0.;
    double  purity  = 0.;
    double  eOverP  = 0.;
    double  sigma   = 0.;
    double  drAdd90 = 0.;

    // get info from cluster
    void GetInfo(edm4eic::Cluster& cluster) {
      nHit = cluster.getHits().size();
      rx   = cluster.getPosition().x;
      ry   = cluster.getPosition().y;
      rz   = cluster.getPosition().z;
      eta  = edm4hep::utils::eta( cluster.getPosition() );
      phi  = atan2( cluster.getPosition().y, cluster.getPosition().x );
      ene  = cluster.getEnergy();
      return;
    }  // end 'GetInfo(edm4eic::Cluster&)'
  };  // end CalContent

}  // end Hist namespace



// CheckClusterSplittingProcessor definition ----------------------------------

class CheckClusterSplittingProcessor : public JEventProcessorSequentialRoot {

  // struct to hold user options
  struct Config {
    float  drStep;
    float  drMax;
    string parCollect;
    string projCollect;
    vector<pair<string, string>> vecCalos;
    vector<double> vecAvgEP;
    vector<double> vecSigEP;
  } m_config = {
    .drStep      = 0.05,
    .drMax       = 3.0,
    .parCollect  = "GeneratedParticles",
    .projCollect = "CalorimeterTrackProjections",
    .vecCalos    = {
      make_pair("HcalEndcapNClusters", "HcalEndcapN"),
      make_pair("EcalEndcapNClusters", "EcalEndcapN"),
      make_pair("HcalBarrelClusters",  "HcalBarrel"),
      make_pair("EcalBarrelClusters",  "EcalBarrelImaging"),
      make_pair("LFHCALClusters",      "LFHCAL"),
      make_pair("EcalEndcapPClusters", "EcalEndcapP")
    },
    .vecAvgEP = {
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0
    },
    .vecSigEP = {
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0
    }
  };  // end Config

  public:

    // ctor
    CheckClusterSplittingProcessor() { SetTypeName(NAME_OF_THIS); }

    // public methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

  private:

    // private methods
    void  BuildEvtHistograms();
    void  BuildTrkHistograms();
    void  BuildCalHistograms();
    void  FillEvtHistograms(const int calo, const int type, const Hist::EvtContent& content);
    void  FillTrkHistograms(const int calo, const int type, const Hist::TrkContent& content);
    void  FillCalHistograms(const int calo, const int type, const Hist::CalContent& content);
    void  SaveHistograms();
    void  TryToRecover(const double eStart, const double ePar, bool& didRecover, float& drRecover, int32_t& nRecover);
    void  GetClustersAndDr(const edm4eic::ClusterCollection& others, const Hist::CalContent& reference, const int idRef);
    void  GetProjections(const edm4eic::TrackSegmentCollection& projections, const int calo);
    int   GetCaloID(const int iCalo);
    FPair GetClusterWidths(const edm4eic::Cluster& cluster);
    FPair GetEtaRange(const edm4eic::Cluster& cluster);
    FPair GetPhiRange(const edm4eic::Cluster& cluster);

    // output histograms
    TH1Vec m_vecEvtH1D;
    TH1Vec m_vecTrkH1D;
    TH1Vec m_vecCalH1D;
    TH2Vec m_vecTrkH2D;
    TH2Vec m_vecCalH2D;

    // for matching tracks to clusters
    TrkVec m_vecTrkProj;

    // for adding clusters during split recovery
    map<int, bool>   m_mapClustToChecked;
    map<int, double> m_mapClustToEne;
    map<int, double> m_mapClustToDr;

};  // end CheckClusterSplittingProcessor

#endif

// end ------------------------------------------------------------------------
