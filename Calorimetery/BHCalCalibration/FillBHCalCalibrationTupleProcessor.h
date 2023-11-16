// ----------------------------------------------------------------------------
// 'FillBHCalCalibrationTupleProcessor.h'
// Derek Anderson
// 03.02.2023
//
// A simple JANA plugin to compare the
// reconstructed hit and cluster energy
// in the HCal to simulated particles.
// ----------------------------------------------------------------------------

// C includes
#include <cmath>
// ROOT includes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>  // FIXME update to XYZvectors
#include <TProfile.h>
// JANA includes
#include <JANA/JEventProcessorSequentialRoot.h>
#include <JANA/JEvent.h>
// EDM includes
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/ProtoCluster.h>
#include <edm4eic/Cluster.h>



// JCalibrateWithImagingProcessor definition ----------------------------------

class FillBHCalCalibrationTupleProcessor : public JEventProcessorSequentialRoot {

  // global constants
  enum CONST {
    NCalibVars  = 51,
    NSciFiLayer = 12,
    NImageLayer = 6,
    NRange      = 2,
    NComp       = 3
  };

  private:

    // Data objects we will need from JANA
    PrefetchT<edm4eic::ReconstructedParticle> genParticles       = {this, "GeneratedParticles"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits       = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters      = {this, "HcalBarrelClusters"};
    PrefetchT<edm4eic::Cluster>               bhcalTruthClusters = {this, "HcalBarrelTruthClusters"};
    PrefetchT<edm4eic::CalorimeterHit>        scifiRecHits       = {this, "EcalBarrelScFiRecHits"};
    PrefetchT<edm4eic::CalorimeterHit>        imageRecHits       = {this, "EcalBarrelImagingRecHits"};
    PrefetchT<edm4eic::Cluster>               bemcClusters       = {this, "EcalBarrelImagingMergedClusters"};
    PrefetchT<edm4eic::Cluster>               scifiClusters      = {this, "EcalBarrelScFiClusters"};
    PrefetchT<edm4eic::Cluster>               imageClusters      = {this, "EcalBarrelImagingClusters"};

    // particle histograms
    TH1D *hParChrg                   = nullptr;
    TH1D *hParMass                   = nullptr;
    TH1D *hParEta                    = nullptr;
    TH1D *hParPhi                    = nullptr;
    TH1D *hParEne                    = nullptr;
    TH1D *hParMom                    = nullptr;
    TH1D *hParMomX                   = nullptr;
    TH1D *hParMomY                   = nullptr;
    TH1D *hParMomZ                   = nullptr;
    TH2D *hParEtaVsPhi               = nullptr;
    // bhcal reconstructed hit histograms
    TH1D *hHCalRecHitEta             = nullptr;
    TH1D *hHCalRecHitPhi             = nullptr;
    TH1D *hHCalRecHitEne             = nullptr;
    TH1D *hHCalRecHitPosZ            = nullptr;
    TH1D *hHCalRecHitParDiff         = nullptr;
    TH2D *hHCalRecHitPosYvsX         = nullptr;
    TH2D *hHCalRecHitEtaVsPhi        = nullptr;
    TH2D *hHCalRecHitVsParEne        = nullptr;
    // bhcal cluster hit histograms
    TH1D *hHCalClustHitEta           = nullptr;
    TH1D *hHCalClustHitPhi           = nullptr;
    TH1D *hHCalClustHitEne           = nullptr;
    TH1D *hHCalClustHitPosZ          = nullptr;
    TH1D *hHCalClustHitParDiff       = nullptr;
    TH2D *hHCalClustHitPosYvsX       = nullptr;
    TH2D *hHCalClustHitEtaVsPhi      = nullptr;
    TH2D *hHCalClustHitVsParEne      = nullptr;
    // bhcal reconstructed cluster histograms
    TH1D *hHCalClustEta              = nullptr;
    TH1D *hHCalClustPhi              = nullptr;
    TH1D *hHCalClustEne              = nullptr;
    TH1D *hHCalClustPosZ             = nullptr;
    TH1I *hHCalClustNumHit           = nullptr;
    TH1D *hHCalClustParDiff          = nullptr;
    TH2D *hHCalClustPosYvsX          = nullptr;
    TH2D *hHCalClustEtaVsPhi         = nullptr;
    TH2D *hHCalClustVsParEne         = nullptr;
    // bhcal truth cluster hit histograms
    TH1D *hHCalTruClustHitEta        = nullptr;
    TH1D *hHCalTruClustHitPhi        = nullptr;
    TH1D *hHCalTruClustHitEne        = nullptr;
    TH1D *hHCalTruClustHitPosZ       = nullptr;
    TH1D *hHCalTruClustHitParDiff    = nullptr;
    TH2D *hHCalTruClustHitPosYvsX    = nullptr;
    TH2D *hHCalTruClustHitEtaVsPhi   = nullptr;
    TH2D *hHCalTruClustHitVsParEne   = nullptr;
    // bhcal truth cluster histograms
    TH1D *hHCalTruClustEta           = nullptr;
    TH1D *hHCalTruClustPhi           = nullptr;
    TH1D *hHCalTruClustEne           = nullptr;
    TH1D *hHCalTruClustPosZ          = nullptr;
    TH1I *hHCalTruClustNumHit        = nullptr;
    TH1D *hHCalTruClustParDiff       = nullptr;
    TH2D *hHCalTruClustPosYvsX       = nullptr;
    TH2D *hHCalTruClustEtaVsPhi      = nullptr;
    TH2D *hHCalTruClustVsParEne      = nullptr;
    // bhcal general event-wise histograms
    TH1I *hEvtHCalNumPar             = nullptr;
    // bhcal hit event-wise histograms
    TH1I *hEvtHCalNumHit             = nullptr;
    TH1D *hEvtHCalSumHitEne          = nullptr;
    TH1D *hEvtHCalSumHitDiff         = nullptr;
    TH2D *hEvtHCalSumHitVsPar        = nullptr;
    // bhcal cluster event-wise histograms
    TH1I *hEvtHCalNumClust           = nullptr;
    TH1D *hEvtHCalSumClustEne        = nullptr;
    TH1D *hEvtHCalSumClustDiff       = nullptr;
    TH2I *hEvtHCalNumClustVsHit      = nullptr;
    TH2D *hEvtHCalSumClustVsPar      = nullptr;
    // bhcal lead cluster event-wise histograms
    TH1I *hEvtHCalLeadClustNumHit    = nullptr;
    TH1D *hEvtHCalLeadClustEne       = nullptr;
    TH1D *hEvtHCalLeadClustDiff      = nullptr;
    TH2D *hEvtHCalLeadClustVsPar     = nullptr;
    // bhcal truth cluster event-wise histograms
    TH1I *hEvtHCalNumTruClust        = nullptr;
    TH1D *hEvtHCalSumTruClustEne     = nullptr;
    TH1D *hEvtHCalSumTruClustDiff    = nullptr;
    TH2I *hEvtHCalNumTruClustVsClust = nullptr;
    TH2D *hEvtHCalSumTruClustVsPar   = nullptr;
    // bhcal truth lead cluster event-wise histograms
    TH1I *hEvtHCalLeadTruClustNumHit = nullptr;
    TH1D *hEvtHCalLeadTruClustEne    = nullptr;
    TH1D *hEvtHCalLeadTruClustDiff   = nullptr;
    TH2D *hEvtHCalLeadTruClustVsPar  = nullptr;

    // scifi reconstructed hit histograms
    TH1I *hSciFiRecHitNLayer         = nullptr;
    TH1D *hSciFiRecHitEta            = nullptr;
    TH1D *hSciFiRecHitPhi            = nullptr;
    TH1D *hSciFiRecHitEne            = nullptr;
    TH1D *hSciFiRecHitPosZ           = nullptr;
    TH1D *hSciFiRecHitParDiff        = nullptr;
    TH2D *hSciFiRecHitPosYvsX        = nullptr;
    TH2D *hSciFiRecHitEtaVsPhi       = nullptr;
    TH2D *hSciFiRecHitVsParEne       = nullptr;
    TH2D *hSciFiRecHitEneVsNLayer    = nullptr;
    // image reconstructed hit histograms
    TH1I *hImageRecHitNLayer         = nullptr;
    TH1D *hImageRecHitEta            = nullptr;
    TH1D *hImageRecHitPhi            = nullptr;
    TH1D *hImageRecHitEne            = nullptr;
    TH1D *hImageRecHitPosZ           = nullptr;
    TH1D *hImageRecHitParDiff        = nullptr;
    TH2D *hImageRecHitPosYvsX        = nullptr;
    TH2D *hImageRecHitEtaVsPhi       = nullptr;
    TH2D *hImageRecHitVsParEne       = nullptr;
    TH2D *hImageRecHitEneVsNLayer    = nullptr;
    // bemc reconstructed cluster histograms
    TH1D *hECalClustEta              = nullptr;
    TH1D *hECalClustPhi              = nullptr;
    TH1D *hECalClustEne              = nullptr;
    TH1D *hECalClustPosZ             = nullptr;
    TH1I *hECalClustNumHit           = nullptr;
    TH1D *hECalClustParDiff          = nullptr;
    TH2D *hECalClustPosYvsX          = nullptr;
    TH2D *hECalClustEtaVsPhi         = nullptr;
    TH2D *hECalClustVsParEne         = nullptr;
    // scifi hit event-wise histograms
    TH1D *hEvtSciFiSumEne            = nullptr;
    TH2D *hEvtSciFiSumEneVsNLayer    = nullptr;
    TH2D *hEvtSciFiVsHCalHitSumEne   = nullptr;
    // image hit event-wise histograms
    TH1D *hEvtImageSumEne            = nullptr;
    TH2D *hEvtImageSumEneVsNLayer    = nullptr;
    TH2D *hEvtImageVsHCalHitSumEne   = nullptr;
    // bemc cluster event-wise histograms
    TH1I *hEvtECalNumClust           = nullptr;
    TH1D *hEvtECalSumClustEne        = nullptr;
    TH1D *hEvtECalSumClustDiff       = nullptr;
    TH2D *hEvtECalSumClustVsPar      = nullptr;
    TH2D *hEvtECalVsHCalSumClustEne  = nullptr;
    // bemc lead cluster event-wise histograms
    TH1I *hEvtECalLeadClustNumHit    = nullptr;
    TH1D *hEvtECalLeadClustEne       = nullptr;
    TH1D *hEvtECalLeadClustDiff      = nullptr;
    TH2D *hEvtECalLeadClustVsPar     = nullptr;
    TH2D *hEvtECalVsHCalLeadClustEne = nullptr;

    // ntuple for calibration
    Float_t  varsForCalibration[CONST::NCalibVars];
    TNtuple *ntForCalibration;

  public:

    // ctor
    FillBHCalCalibrationTupleProcessor() { SetTypeName(NAME_OF_THIS); }

    // inherited methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

};  // end FillBHCalCalibrationTupleProcessor definition

// end ------------------------------------------------------------------------
