// ----------------------------------------------------------------------------
// 'FillBHCalCalibrationTupleProcessor.cc'
// Derek Anderson
// 03.02.2023
//
// A simple JANA plugin to compare the
// reconstructed hit and cluster energy
// in the HCal to simulated particles.
// ----------------------------------------------------------------------------

// C includes
#include <vector>
#include <string>
#include <boost/math/special_functions/sign.hpp>
// eicrecon includes 
#include <services/rootfile/RootFile_service.h>
// user includes
#include "FillBHCalCalibrationTupleProcessor.h"

// The following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
  InitJANAPlugin(app);
  app -> Add(new FillBHCalCalibrationTupleProcessor);
  }
}



//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void FillBHCalCalibrationTupleProcessor::InitWithGlobalRootLock(){

  // create directory in output file
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc     -> GetHistFile();
  rootfile -> mkdir("FillBHCalCalibrationTuple") -> cd();

  // initialize bhcal histograms
  const unsigned long nNumBin(200);
  const unsigned long nChrgBin(6);
  const unsigned long nMassBin(1000);
  const unsigned long nPhiBin(60);
  const unsigned long nEtaBin(40);
  const unsigned long nEneBin(200);
  const unsigned long nMomBin(200);
  const unsigned long nPosTrBin(800);
  const unsigned long nPosLoBin(30);
  const unsigned long nDiffBin(200);
  const unsigned long rNumBin[CONST::NRange]   = {0,      200};
  const double        rChrgBin[CONST::NRange]  = {-3.,    3.};
  const double        rMassBin[CONST::NRange]  = {0.,     5.};
  const double        rPhiBin[CONST::NRange]   = {-3.15,  3.15};
  const double        rEtaBin[CONST::NRange]   = {-2.,    2.};
  const double        rEneBin[CONST::NRange]   = {0.,     100.};
  const double        rMomBin[CONST::NRange]   = {-50.,   50.};
  const double        rPosTrBin[CONST::NRange] = {-4000., 4000.};
  const double        rPosLoBin[CONST::NRange] = {-3000., 3000.};
  const double        rDiffBin[CONST::NRange]  = {-5.,    5.};
  // particle histograms
  hParChrg                   = new TH1D("hParChrg",     "Gen. Particles", nChrgBin, rChrgBin[0], rChrgBin[1]);
  hParMass                   = new TH1D("hParMass",     "Gen. Particles", nMassBin, rMassBin[0], rMassBin[1]);
  hParPhi                    = new TH1D("hParPhi",      "Gen. Particles", nPhiBin,  rPhiBin[0],  rPhiBin[1]);
  hParEta                    = new TH1D("hParEta",      "Gen. Particles", nEtaBin,  rEtaBin[0],  rEtaBin[1]);
  hParEne                    = new TH1D("hParEne",      "Gen. Particles", nEneBin,  rEneBin[0],  rEneBin[1]);
  hParMom                    = new TH1D("hParMom",      "Gen. Particles", nEneBin,  rEneBin[0],  rEneBin[1]);
  hParMomX                   = new TH1D("hParMomX",     "Gen. Particles", nMomBin,  rMomBin[0],  rMomBin[1]);
  hParMomY                   = new TH1D("hParMomY",     "Gen. Particles", nMomBin,  rMomBin[0],  rMomBin[1]);
  hParMomZ                   = new TH1D("hParMomZ",     "Gen. Particles", nMomBin,  rMomBin[0],  rMomBin[1]);
  hParEtaVsPhi               = new TH2D("hParEtaVsPhi", "Gen. Particles", nPhiBin,  rPhiBin[0],  rPhiBin[1], nEtaBin, rEtaBin[0], rEtaBin[1]);
  // reco. bhcal hit histograms
  hHCalRecHitPhi             = new TH1D("hHCalRecHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalRecHitEta             = new TH1D("hHCalRecHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalRecHitEne             = new TH1D("hHCalRecHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalRecHitPosZ            = new TH1D("hHCalRecHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalRecHitParDiff         = new TH1D("hHCalRecHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalRecHitPosYvsX         = new TH2D("hHCalRecHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalRecHitEtaVsPhi        = new TH2D("hHCalRecHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalRecHitVsParEne        = new TH2D("hHCalRecHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // bhcal cluster hit histograms
  hHCalClustHitPhi           = new TH1D("hHCalClustHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalClustHitEta           = new TH1D("hHCalClustHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustHitEne           = new TH1D("hHCalClustHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustHitPosZ          = new TH1D("hHCalClustHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustHitParDiff       = new TH1D("hHCalClustHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustHitPosYvsX       = new TH2D("hHCalClustHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustHitEtaVsPhi      = new TH2D("hHCalClustHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustHitVsParEne      = new TH2D("hHCalClustHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. bhcal cluster histograms
  hHCalClustPhi              = new TH1D("hHCalClustPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalClustEta              = new TH1D("hHCalClustEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustEne              = new TH1D("hHCalClustEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustPosZ             = new TH1D("hHCalClustPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustNumHit           = new TH1I("hHCalClustNumHit",   "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalClustParDiff          = new TH1D("hHCalClustParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustPosYvsX          = new TH2D("hHCalClustPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustEtaVsPhi         = new TH2D("hHCalClustEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustVsParEne         = new TH2D("hHCalClustVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // bhcal cluster hit histograms
  hHCalTruClustHitPhi        = new TH1D("hHCalTruClustHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalTruClustHitEta        = new TH1D("hHCalTruClustHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustHitEne        = new TH1D("hHCalTruClustHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustHitPosZ       = new TH1D("hHCalTruClustHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustHitParDiff    = new TH1D("hHCalTruClustHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustHitPosYvsX    = new TH2D("hHCalTruClustHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustHitEtaVsPhi   = new TH2D("hHCalTruClustHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustHitVsParEne   = new TH2D("hHCalTruClustHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth bhcal cluster histograms
  hHCalTruClustPhi           = new TH1D("hHCalTruClustPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalTruClustEta           = new TH1D("hHCalTruClustEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustEne           = new TH1D("hHCalTruClustEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustPosZ          = new TH1D("hHCalTruClustPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustNumHit        = new TH1I("hHCalTruClustNumHit",   "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalTruClustParDiff       = new TH1D("hHCalTruClustParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustPosYvsX       = new TH2D("hHCalTruClustPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustEtaVsPhi      = new TH2D("hHCalTruClustEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustVsParEne      = new TH2D("hHCalTruClustVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // bhcal general event-wise histograms
  hEvtHCalNumPar             = new TH1I("hEvtHCalNumPar", "Barrel HCal", nNumBin, rNumBin[0], rNumBin[1]);
  // bhcal hit event-wise histograms
  hEvtHCalNumHit             = new TH1I("hEvtHCalNumHit",      "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtHCalSumHitEne          = new TH1D("hEvtHCalSumHitEne",   "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtHCalSumHitDiff         = new TH1D("hEvtHCalSumHitDiff",  "Barrel HCal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtHCalSumHitVsPar        = new TH2D("hEvtHCalSumHitVsPar", "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // bhcal cluster event-wise histograms
  hEvtHCalNumClust           = new TH1I("hEvtHCalNumClust",      "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtHCalSumClustEne        = new TH1D("hEvtHCalSumClustEne",   "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtHCalSumClustDiff       = new TH1D("hEvtHCalSumClustDiff",  "Barrel HCal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtHCalNumClustVsHit      = new TH2I("hEvtHCalNumClustVsHit", "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1], nNumBin, rNumBin[0], rNumBin[1]);
  hEvtHCalSumClustVsPar      = new TH2D("hEvtHCalSumClustVsPar", "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // bhcal lead cluster event-wise histograms
  hEvtHCalLeadClustNumHit    = new TH1I("hEvtHCalLeadClustNumHit", "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtHCalLeadClustEne       = new TH1D("hEvtHCalLeadClustEne",    "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtHCalLeadClustDiff      = new TH1D("hEvtHCalLeadClustDiff",   "Barrel HCal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtHCalLeadClustVsPar     = new TH2D("hEvtHCalLeadClustVsPar",  "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // bhcal truth cluster event-wise histograms
  hEvtHCalNumTruClust        = new TH1I("hEvtHCalNumTruClust",        "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtHCalSumTruClustEne     = new TH1D("hEvtHCalSumTruClustEne",     "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtHCalSumTruClustDiff    = new TH1D("hEvtHCalSumTruClustDiff",    "Barrel HCal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtHCalNumTruClustVsClust = new TH2I("hEvtHCalNumTruClustVsClust", "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1], nNumBin, rNumBin[0], rNumBin[1]);
  hEvtHCalSumTruClustVsPar   = new TH2D("hEvtHCalSumTruClustVsPar",   "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // bhcal truth lead cluster event-wise histograms
  hEvtHCalLeadTruClustNumHit = new TH1I("hEvtHCalLeadTruClustNumHit", "Barrel HCal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtHCalLeadTruClustEne    = new TH1D("hEvtHCalLeadTruClustEne",    "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtHCalLeadTruClustDiff   = new TH1D("hEvtHCalLeadTruClustDiff",   "Barrel HCal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtHCalLeadTruClustVsPar  = new TH2D("hEvtHCalLeadTruClustVsPar",  "Barrel HCal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // bhcal particle errors
  hParChrg                   -> Sumw2();
  hParMass                   -> Sumw2();
  hParPhi                    -> Sumw2();
  hParEta                    -> Sumw2();
  hParEne                    -> Sumw2();
  hParMom                    -> Sumw2();
  hParMomX                   -> Sumw2();
  hParMomY                   -> Sumw2();
  hParMomZ                   -> Sumw2();
  hParEtaVsPhi               -> Sumw2();
  // bhcal reco. hit errors
  hHCalRecHitPhi             -> Sumw2();
  hHCalRecHitEta             -> Sumw2();
  hHCalRecHitEne             -> Sumw2();
  hHCalRecHitPosZ            -> Sumw2();
  hHCalRecHitParDiff         -> Sumw2();
  hHCalRecHitPosYvsX         -> Sumw2();
  hHCalRecHitEtaVsPhi        -> Sumw2();
  hHCalRecHitVsParEne        -> Sumw2();
  // bhcal reco. cluster hit errors
  hHCalClustHitPhi           -> Sumw2();
  hHCalClustHitEta           -> Sumw2();
  hHCalClustHitEne           -> Sumw2();
  hHCalClustHitPosZ          -> Sumw2();
  hHCalClustHitParDiff       -> Sumw2();
  hHCalClustHitPosYvsX       -> Sumw2();
  hHCalClustHitEtaVsPhi      -> Sumw2();
  hHCalClustHitVsParEne      -> Sumw2();
  // bhcal reco. cluster errors
  hHCalClustPhi              -> Sumw2();
  hHCalClustEta              -> Sumw2();
  hHCalClustEne              -> Sumw2();
  hHCalClustPosZ             -> Sumw2();
  hHCalClustNumHit           -> Sumw2();
  hHCalClustParDiff          -> Sumw2();
  hHCalClustPosYvsX          -> Sumw2();
  hHCalClustEtaVsPhi         -> Sumw2();
  hHCalClustVsParEne         -> Sumw2();
  // bhcal truth cluster hit errors
  hHCalTruClustHitPhi        -> Sumw2();
  hHCalTruClustHitEta        -> Sumw2();
  hHCalTruClustHitEne        -> Sumw2();
  hHCalTruClustHitPosZ       -> Sumw2();
  hHCalTruClustHitParDiff    -> Sumw2();
  hHCalTruClustHitPosYvsX    -> Sumw2();
  hHCalTruClustHitEtaVsPhi   -> Sumw2();
  hHCalTruClustHitVsParEne   -> Sumw2();
  // bhcal truth cluster errors
  hHCalTruClustPhi           -> Sumw2();
  hHCalTruClustEta           -> Sumw2();
  hHCalTruClustEne           -> Sumw2();
  hHCalTruClustPosZ          -> Sumw2();
  hHCalTruClustNumHit        -> Sumw2();
  hHCalTruClustParDiff       -> Sumw2();
  hHCalTruClustPosYvsX       -> Sumw2();
  hHCalTruClustEtaVsPhi      -> Sumw2();
  hHCalTruClustVsParEne      -> Sumw2();
  // bhcal general event-wise errors
  hEvtHCalNumPar             -> Sumw2();
  // bhcal hit event-wise errors
  hEvtHCalNumHit             -> Sumw2();
  hEvtHCalSumHitEne          -> Sumw2();
  hEvtHCalSumHitDiff         -> Sumw2();
  hEvtHCalSumHitVsPar        -> Sumw2();
  // bhcal cluster event-wise errors
  hEvtHCalNumClust           -> Sumw2();
  hEvtHCalSumClustEne        -> Sumw2();
  hEvtHCalSumClustDiff       -> Sumw2();
  hEvtHCalNumClustVsHit      -> Sumw2();
  hEvtHCalSumClustVsPar      -> Sumw2();
  // bhcal lead cluster event-wise errors
  hEvtHCalLeadClustNumHit    -> Sumw2();
  hEvtHCalLeadClustEne       -> Sumw2();
  hEvtHCalLeadClustDiff      -> Sumw2();
  hEvtHCalLeadClustVsPar     -> Sumw2();
  // bhcal truth cluster event-wise errors
  hEvtHCalNumTruClust        -> Sumw2();
  hEvtHCalSumTruClustEne     -> Sumw2();
  hEvtHCalSumTruClustDiff    -> Sumw2();
  hEvtHCalNumTruClustVsClust -> Sumw2();
  // bhcal truth lead cluster event-wise errors
  hEvtHCalLeadTruClustNumHit -> Sumw2();
  hEvtHCalLeadTruClustEne    -> Sumw2();
  hEvtHCalLeadTruClustDiff   -> Sumw2();
  hEvtHCalLeadTruClustVsPar  -> Sumw2();

  // initialize becal histograms
  const unsigned long nLayerBin(50);
  const unsigned long rLayerBin[CONST::NRange] = {0, 50};
  // scifi reconstructed hit histograms
  hSciFiRecHitNLayer         = new TH1I("hSciFiRecHitNLayer",      "Barrel SciFi", nLayerBin, rLayerBin[0], rLayerBin[1]);
  hSciFiRecHitPhi            = new TH1D("hSciFiRecHitPhi",         "Barrel SciFi", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hSciFiRecHitEta            = new TH1D("hSciFiRecHitEta",         "Barrel SciFi", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hSciFiRecHitEne            = new TH1D("hSciFiRecHitEne",         "Barrel SciFi", nEneBin,   rEneBin[0],   rEneBin[1]);
  hSciFiRecHitPosZ           = new TH1D("hSciFiRecHitPosZ",        "Barrel SciFi", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hSciFiRecHitParDiff        = new TH1D("hSciFiRecHitParDiff",     "Barrel SciFi", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hSciFiRecHitPosYvsX        = new TH2D("hSciFiRecHitPosYvsX",     "Barrel SciFi", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hSciFiRecHitEtaVsPhi       = new TH2D("hSciFiRecHitEtaVsPhi",    "Barrel SciFi", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hSciFiRecHitVsParEne       = new TH2D("hSciFiRecHitVsParEne",    "Barrel SciFi", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hSciFiRecHitEneVsNLayer    = new TH2D("hSciFiRecHitEneVsNLayer", "Barrel SciFi", nLayerBin, rLayerBin[0], rLayerBin[1], nEneBin,   rEneBin[0],   rEneBin[1]);
  // imaging reconstructed hit histograms
  hImageRecHitNLayer         = new TH1I("hImageRecHitNLayer",      "Barrel Image", nLayerBin, rLayerBin[0], rLayerBin[1]);
  hImageRecHitPhi            = new TH1D("hImageRecHitPhi",         "Barrel Image", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hImageRecHitEta            = new TH1D("hImageRecHitEta",         "Barrel Image", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hImageRecHitEne            = new TH1D("hImageRecHitEne",         "Barrel Image", nEneBin,   rEneBin[0],   rEneBin[1]);
  hImageRecHitPosZ           = new TH1D("hImageRecHitPosZ",        "Barrel Image", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hImageRecHitParDiff        = new TH1D("hImageRecHitParDiff",     "Barrel Image", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hImageRecHitPosYvsX        = new TH2D("hImageRecHitPosYvsX",     "Barrel Image", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hImageRecHitEtaVsPhi       = new TH2D("hImageRecHitEtaVsPhi",    "Barrel Image", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hImageRecHitVsParEne       = new TH2D("hImageRecHitVsParEne",    "Barrel Image", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hImageRecHitEneVsNLayer    = new TH2D("hImageRecHitEneVsNLayer", "Barrel Image", nLayerBin, rLayerBin[0], rLayerBin[1], nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. bemc cluster histograms
  hECalClustPhi              = new TH1D("hECalClustPhi",      "Barrel ECal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hECalClustEta              = new TH1D("hECalClustEta",      "Barrel ECal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hECalClustEne              = new TH1D("hECalClustEne",      "Barrel ECal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalClustPosZ             = new TH1D("hECalClustPosZ",     "Barrel ECal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hECalClustNumHit           = new TH1I("hECalClustNumHit",   "Barrel ECal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hECalClustParDiff          = new TH1D("hECalClustParDiff",  "Barrel ECal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalClustPosYvsX          = new TH2D("hECalClustPosYvsX",  "Barrel ECal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hECalClustEtaVsPhi         = new TH2D("hECalClustEtaVsPhi", "Barrel ECal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hECalClustVsParEne         = new TH2D("hECalClustVsParEne", "Barrel ECal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // scifi hit event-wise histogram
  hEvtSciFiSumEne            = new TH1D("hEvtSciFiSumEne",          "Barrel SciFi", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtSciFiSumEneVsNLayer    = new TH2D("hEvtSciFiSumEneVsNLayer",  "Barrel SciFi", nLayerBin, rLayerBin[0], rLayerBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  hEvtSciFiVsHCalHitSumEne   = new TH2D("hEvtSciFiVsHCalHitSumEne", "Barrel SciFi", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin, rEneBin[0], rEneBin[1]);
  // imaging hit event-wise histogram
  hEvtImageSumEne            = new TH1D("hEvtImageSumEne",          "Barrel Image", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtImageSumEneVsNLayer    = new TH2D("hEvtImageSumEneVsNLayer",  "Barrel Image", nLayerBin, rLayerBin[0], rLayerBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  hEvtImageVsHCalHitSumEne   = new TH2D("hEvtImageVsHCalHitSumEne", "Barrel Image", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin, rEneBin[0], rEneBin[1]);
  // bemc cluster event-wise histograms
  hEvtECalNumClust           = new TH1I("hEvtECalNumClust",          "Barrel ECal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtECalSumClustEne        = new TH1D("hEvtECalSumClustEne",       "Barrel ECal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtECalSumClustDiff       = new TH1D("hEvtECalSumClustDiff",      "Barrel ECal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtECalSumClustVsPar      = new TH2D("hEvtECalSumClustVsPar",     "Barrel ECal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  hEvtECalVsHCalSumClustEne  = new TH2D("hEvtECalVsHCalSumClustEne", "Barrel ECal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // bemc lead cluster event-wise histograms
  hEvtECalLeadClustNumHit    = new TH1I("hEvtECalLeadClustNumHit",    "Barrel ECal", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtECalLeadClustEne       = new TH1D("hEvtECalLeadClustEne",       "Barrel ECal", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtECalLeadClustDiff      = new TH1D("hEvtECalLeadClustDiff",      "Barrel ECal", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtECalLeadClustVsPar     = new TH2D("hEvtECalLeadClustVsPar",     "Barrel ECal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  hEvtECalVsHCalLeadClustEne = new TH2D("hEvtECalVsHCalLeadClustEne", "Barrel ECal", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // scifi reconstructed hit errors
  hSciFiRecHitNLayer         -> Sumw2();
  hSciFiRecHitPhi            -> Sumw2();
  hSciFiRecHitEta            -> Sumw2();
  hSciFiRecHitEne            -> Sumw2();
  hSciFiRecHitPosZ           -> Sumw2();
  hSciFiRecHitParDiff        -> Sumw2();
  hSciFiRecHitPosYvsX        -> Sumw2();
  hSciFiRecHitEtaVsPhi       -> Sumw2();
  hSciFiRecHitVsParEne       -> Sumw2();
  hSciFiRecHitEneVsNLayer    -> Sumw2();
  // imaging reconstructed hit errors
  hImageRecHitNLayer         -> Sumw2();
  hImageRecHitPhi            -> Sumw2();
  hImageRecHitEta            -> Sumw2();
  hImageRecHitEne            -> Sumw2();
  hImageRecHitPosZ           -> Sumw2();
  hImageRecHitParDiff        -> Sumw2();
  hImageRecHitPosYvsX        -> Sumw2();
  hImageRecHitEtaVsPhi       -> Sumw2();
  hImageRecHitVsParEne       -> Sumw2();
  hImageRecHitEneVsNLayer    -> Sumw2();
  // bemc reconstructed cluster errors
  hECalClustEta              -> Sumw2();
  hECalClustPhi              -> Sumw2();
  hECalClustEne              -> Sumw2();
  hECalClustPosZ             -> Sumw2();
  hECalClustNumHit           -> Sumw2();
  hECalClustParDiff          -> Sumw2();
  hECalClustPosYvsX          -> Sumw2();
  hECalClustEtaVsPhi         -> Sumw2();
  hECalClustVsParEne         -> Sumw2();
  // scifi hit event-wise histogram
  hEvtSciFiSumEne            -> Sumw2();
  hEvtSciFiSumEneVsNLayer    -> Sumw2();
  hEvtSciFiVsHCalHitSumEne   -> Sumw2();
  // imaging hit event-wise histogram
  hEvtImageSumEne            -> Sumw2();
  hEvtImageSumEneVsNLayer    -> Sumw2();
  hEvtImageVsHCalHitSumEne   -> Sumw2();
  // bemc cluster event-wise errors
  hEvtECalNumClust           -> Sumw2();
  hEvtECalSumClustEne        -> Sumw2();
  hEvtECalSumClustDiff       -> Sumw2();
  hEvtECalSumClustVsPar      -> Sumw2();
  hEvtECalVsHCalLeadClustEne -> Sumw2();
  // bemc lead cluster event-wise errors
  hEvtECalLeadClustNumHit    -> Sumw2();
  hEvtECalLeadClustEne       -> Sumw2();
  hEvtECalLeadClustDiff      -> Sumw2();
  hEvtECalLeadClustVsPar     -> Sumw2();
  hEvtECalVsHCalLeadClustEne -> Sumw2();

  // calibration variables
  const std::vector<std::string> vecCalibVars = {
    "ePar",
    "fracParVsLeadBHCal",
    "fracParVsLeadBEMC",
    "fracParVsSumBHCal",
    "fracParVsSumBEMC",
    "fracLeadBHCalVsBEMC",
    "fracSumBHCalVsBEMC",
    "eLeadBHCal",
    "eLeadBEMC",
    "eSumBHCal",
    "eSumBEMC",
    "diffLeadBHCal",
    "diffLeadBEMC",
    "diffSumBHCal",
    "diffSumBEMC",
    "nHitsLeadBHCal",
    "nHitsLeadBEMC",
    "nClustBHCal",
    "nClustBEMC",
    "hLeadBHCal",
    "hLeadBEMC",
    "fLeadBHCal",
    "fLeadBEMC",
    "eLeadImage",
    "eSumImage",
    "eLeadSciFi",
    "eSumSciFi",
    "nClustImage",
    "nClustSciFi",
    "hLeadImage",
    "hLeadSciFi",
    "fLeadImage",
    "fLeadSciFi",
    "eSumSciFiLayer1",
    "eSumSciFiLayer2",
    "eSumSciFiLayer3",
    "eSumSciFiLayer4",
    "eSumSciFiLayer5",
    "eSumSciFiLayer6",
    "eSumSciFiLayer7",
    "eSumSciFiLayer8",
    "eSumSciFiLayer9",
    "eSumSciFiLayer10",
    "eSumSciFiLayer11",
    "eSumSciFiLayer12",
    "eSumImageLayer1",
    "eSumImageLayer2",
    "eSumImageLayer3",
    "eSumImageLayer4",
    "eSumImageLayer5",
    "eSumImageLayer6"
  };

  // convert to ntuple argument
  std::string argCalibVars("");
  for (size_t iCalibVar = 0; iCalibVar < vecCalibVars.size(); iCalibVar++) {
    argCalibVars.append(vecCalibVars[iCalibVar]);
    if ((iCalibVar + 1) != vecCalibVars.size()) {
      argCalibVars.append(":");
    }
  }

  // ntuple for calibration
  ntForCalibration = new TNtuple("ntForCalibration", "For Calibration", argCalibVars.c_str());
  return;

}  // end 'InitWithGlobalRootLock()'




//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void FillBHCalCalibrationTupleProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // clear array for ntuple
  for (size_t iCalibVar = 0; iCalibVar < NCalibVars; iCalibVar++) {
    varsForCalibration[iCalibVar] = 0.;
  }

  // hit and cluster sums
  double eHCalHitSum(0.);
  double eHCalClustSum(0.);
  double eECalClustSum(0.);
  double eImageClustSum(0.);
  double eSciFiClustSum(0.);
  double eTruHCalClustSum(0.);

  // sum bhcal hit energy
  for (auto bhCalHit : bhcalRecHits()) {
    eHCalHitSum += bhCalHit -> getEnergy();
  }  // end 1st bhcal hit loop

  // if hit sum is 0, skip event
  const bool isHCalHitSumNonzero = (eHCalHitSum > 0.);
  if (!isHCalHitSumNonzero) {
    return;
  }

  // MC particle properties
  float  cMcPar(0.);
  double mMcPar(0.);
  double fMcPar(0.);
  double hMcPar(0.);
  double eMcPar(0.);
  double pTotMcPar(0.);
  double pMcPar[NComp] = {0., 0., 0.};

  // particle loop
  unsigned long nPar(0);
  for (auto par : genParticles()) {

    // grab particle properties
    const auto typePar = par -> getType();
    const auto cPar    = par -> getCharge();
    const auto mPar    = par -> getMass();
    const auto ePar    = par -> getEnergy();
    const auto pParX   = par -> getMomentum().x;
    const auto pParY   = par -> getMomentum().y;
    const auto pParZ   = par -> getMomentum().z;
    const auto pPar    = std::sqrt((pParX * pParX) + (pParY * pParY) + (pParZ * pParZ));
    const auto fPar    = std::atan(pParY / pParX);
    const auto hPar    = std::atanh(pParZ / pPar);

    // select MC particle
    const bool isMcParticle = (typePar == 1);
    if (isMcParticle) {
      cMcPar    = cPar;
      mMcPar    = mPar;
      fMcPar    = fPar;
      hMcPar    = hPar;
      eMcPar    = ePar;
      pMcPar[0] = pParX;
      pMcPar[1] = pParY;
      pMcPar[2] = pParZ;
      pTotMcPar = pPar;
    }
    ++nPar;
  }  // end particle loop

  // fill particle histograms
  hParChrg     -> Fill(cMcPar);
  hParMass     -> Fill(mMcPar);
  hParPhi      -> Fill(fMcPar);
  hParEta      -> Fill(hMcPar);
  hParEne      -> Fill(eMcPar);
  hParMom      -> Fill(pTotMcPar);
  hParMomX     -> Fill(pMcPar[0]);
  hParMomY     -> Fill(pMcPar[1]);
  hParMomZ     -> Fill(pMcPar[2]);
  hParEtaVsPhi -> Fill(fMcPar, hMcPar);

  // reco. bhcal hit loop
  unsigned long nHCalHit(0);
  for (auto bhCalHit : bhcalRecHits()) {

    // grab hit properties
    const auto rHCalHitX   = bhCalHit -> getPosition().x;
    const auto rHCalHitY   = bhCalHit -> getPosition().y;
    const auto rHCalHitZ   = bhCalHit -> getPosition().z;
    const auto eHCalHit    = bhCalHit -> getEnergy();
    const auto rHCalHitS   = std::sqrt((rHCalHitX * rHCalHitX) + (rHCalHitY * rHCalHitY));
    const auto rHCalHitR   = std::sqrt((rHCalHitS * rHCalHitS) + (rHCalHitZ * rHCalHitZ));
    const auto fHCalHit    = boost::math::sign(rHCalHitY) * acos(rHCalHitX / rHCalHitS);
    const auto tHCalHit    = std::acos(rHCalHitZ / rHCalHitR);
    const auto hHCalHit    = (-1.) * std::log(std::atan(tHCalHit / 2.));
    const auto diffHCalHit = (eHCalHit - eMcPar) / eMcPar;

    // fill hit histograms and increment sums/counters
    hHCalRecHitPhi      -> Fill(fHCalHit);
    hHCalRecHitEta      -> Fill(hHCalHit);
    hHCalRecHitEne      -> Fill(eHCalHit);
    hHCalRecHitPosZ     -> Fill(rHCalHitZ);
    hHCalRecHitParDiff  -> Fill(diffHCalHit);
    hHCalRecHitPosYvsX  -> Fill(rHCalHitX, rHCalHitY);
    hHCalRecHitEtaVsPhi -> Fill(fHCalHit, hHCalHit);
    hHCalRecHitVsParEne -> Fill(eMcPar, eHCalHit);
    ++nHCalHit;
  }  // end 2nd bhcal hit loop

  // for highest energy bhcal clusters
  int    iLeadHCalClust(-1);
  int    iLeadTruHCalClust(-1);
  int    nHitLeadHCalClust(-1);
  int    nHitLeadTruHCalClust(-1);
  double hLeadHCalClust(-999.);
  double fLeadHCalClust(-999.);
  double eLeadHCalClust(-999.);
  double eLeadTruHCalClust(-999.);
  double diffLeadHCalClust(-999.);
  double diffLeadTruHCalClust(-999.);

  // get protoclusters
  auto bhCalProtoClusters = event -> Get<edm4eic::ProtoCluster>("HcalBarrelIslandProtoClusters");

  // reco. bhcal cluster loop
  unsigned long iHCalClust(0);
  unsigned long nHCalProto(0);
  unsigned long nHCalClust(0);
  for (auto bhCalClust : bhcalClusters()) {

    // grab cluster properties
    const auto rHCalClustX   = bhCalClust -> getPosition().x;
    const auto rHCalClustY   = bhCalClust -> getPosition().y;
    const auto rHCalClustZ   = bhCalClust -> getPosition().z;
    const auto eHCalClust    = bhCalClust -> getEnergy();
    const auto nHitHCalClust = bhCalClust -> getNhits();
    const auto tHCalClust    = bhCalClust -> getIntrinsicTheta();
    const auto diffHCalClust = (eHCalClust - eMcPar) / eMcPar;

    // calculate cluster eta, phi
    // FIXME update to XYZvectors
    const TVector3 vecPosition(rHCalClustX, rHCalClustY, rHCalClustZ);
    const auto     hHCalClust = vecPosition.Eta();
    const auto     fHCalClust = vecPosition.Phi();
    
    // loop over protoclusters
    unsigned long iHCalProto(0);
    unsigned long nProtoHits(0);
    for (auto bhCalProto : bhCalProtoClusters) {

      // check if proto index is same as cluster index
      // FIXME: this might not be the correct way to match reco and proto clusters
      const bool isSameIndex = (iHCalProto == iHCalClust);
      if (!isSameIndex) continue;

      // loop over hits
      nProtoHits = bhCalProto -> hits_size();
      for (uint32_t iProtoHit = 0; iProtoHit < nProtoHits; iProtoHit++) {

        // get hit
        const auto bhCalProtoHit = bhCalProto -> getHits(iProtoHit);

        // grab hit properties
        const auto rHCalProtoHitX   = bhCalProtoHit.getPosition().x;
        const auto rHCalProtoHitY   = bhCalProtoHit.getPosition().y;
        const auto rHCalProtoHitZ   = bhCalProtoHit.getPosition().z;
        const auto eHCalProtoHit    = bhCalProtoHit.getEnergy();
        const auto rHCalProtoHitS   = std::sqrt((rHCalProtoHitX * rHCalProtoHitX) + (rHCalProtoHitY * rHCalProtoHitY));
        const auto rHCalProtoHitR   = std::sqrt((rHCalProtoHitS * rHCalProtoHitS) + (rHCalProtoHitZ * rHCalProtoHitZ));
        const auto fHCalProtoHit    = boost::math::sign(rHCalProtoHitY) * acos(rHCalProtoHitX / rHCalProtoHitS);
        const auto tHCalProtoHit    = std::acos(rHCalProtoHitZ / rHCalProtoHitR);
        const auto hHCalProtoHit    = (-1.) * std::log(std::atan(tHCalProtoHit / 2.));
        const auto diffHCalProtoHit = (eHCalProtoHit - eMcPar) / eMcPar;

        // fill hit histograms and increment sums/counters
        hHCalClustHitPhi      -> Fill(fHCalProtoHit);
        hHCalClustHitEta      -> Fill(hHCalProtoHit);
        hHCalClustHitEne      -> Fill(eHCalProtoHit);
        hHCalClustHitPosZ     -> Fill(rHCalProtoHitZ);
        hHCalClustHitParDiff  -> Fill(diffHCalProtoHit);
        hHCalClustHitPosYvsX  -> Fill(rHCalProtoHitX, rHCalProtoHitY);
        hHCalClustHitEtaVsPhi -> Fill(fHCalProtoHit, hHCalProtoHit);
        hHCalClustHitVsParEne -> Fill(eMcPar, eHCalProtoHit);
      }
      ++nHCalProto;
    }  // end protocluster loop

    // fill cluster histograms and increment counters
    hHCalClustPhi      -> Fill(fHCalClust);
    hHCalClustEta      -> Fill(hHCalClust);
    hHCalClustEne      -> Fill(eHCalClust);
    hHCalClustPosZ     -> Fill(rHCalClustZ);
    hHCalClustNumHit   -> Fill(nHitHCalClust);
    hHCalClustParDiff  -> Fill(diffHCalClust);
    hHCalClustPosYvsX  -> Fill(rHCalClustX, rHCalClustY);
    hHCalClustEtaVsPhi -> Fill(fHCalClust, hHCalClust);
    hHCalClustVsParEne -> Fill(eMcPar, eHCalClust);
    eHCalClustSum += eHCalClust;
    ++nHCalClust;
    ++iHCalClust;

    // select leading cluster
    const bool isBiggerEne = (eHCalClust > eLeadHCalClust);
    if (isBiggerEne) {
      iLeadHCalClust    = iHCalClust;
      nHitLeadHCalClust = nHitHCalClust;
      hLeadHCalClust    = hHCalClust;
      fLeadHCalClust    = fHCalClust;
      eLeadHCalClust    = eHCalClust;
      diffLeadHCalClust = diffHCalClust;
    }
  }  // end reco. bhcal cluster loop

  // get truth protoclusters
  auto bhCalTruProtoClusters = event -> Get<edm4eic::ProtoCluster>("HcalBarrelTruthProtoClusters");

  // true bhcal cluster loop
  unsigned long iTruHCalClust(0);
  unsigned long nTruHCalProto(0);
  unsigned long nTruHCalClust(0);
  for (auto truthHCalClust : bhcalTruthClusters()) {

    // grab cluster properties
    const auto rTruHCalClustX   = truthHCalClust -> getPosition().x;
    const auto rTruHCalClustY   = truthHCalClust -> getPosition().y;
    const auto rTruHCalClustZ   = truthHCalClust -> getPosition().z;
    const auto eTruHCalClust    = truthHCalClust -> getEnergy();
    const auto nHitTruHCalClust = truthHCalClust -> getNhits();
    const auto fTruHCalClust    = truthHCalClust -> getIntrinsicPhi();
    const auto tTruHCalClust    = truthHCalClust -> getIntrinsicTheta();
    const auto hTruHCalClust    = (-1.) * std::log(std::atan(tTruHCalClust / 2.));
    const auto diffTruHCalClust = (eTruHCalClust - eMcPar) / eMcPar;
    
    // loop over protoclusters
    unsigned long iTruHCalProto(0);
    unsigned long nTruProtoHits(0);
    for (auto bhCalTruProto : bhCalTruProtoClusters) {

      // check if truth proto index is same as truth cluster index
      // FIXME: this might not be the correct way to match reco and proto clusters
      const bool isSameIndex = (iTruHCalProto == iTruHCalClust);
      if (!isSameIndex) continue;

      // loop over hits
      nTruProtoHits = bhCalTruProto -> hits_size();
      for (uint32_t iTruProtoHit = 0; iTruProtoHit < nTruProtoHits; iTruProtoHit++) {

        // get hit
        const auto bhCalTruProtoHit = bhCalTruProto -> getHits(iTruProtoHit);

        // grab hit properties
        const auto rTruHCalProtoHitX   = bhCalTruProtoHit.getPosition().x;
        const auto rTruHCalProtoHitY   = bhCalTruProtoHit.getPosition().y;
        const auto rTruHCalProtoHitZ   = bhCalTruProtoHit.getPosition().z;
        const auto eTruHCalProtoHit    = bhCalTruProtoHit.getEnergy();
        const auto rTruHCalProtoHitS   = std::sqrt((rTruHCalProtoHitX * rTruHCalProtoHitX) + (rTruHCalProtoHitY * rTruHCalProtoHitY));
        const auto rTruHCalProtoHitR   = std::sqrt((rTruHCalProtoHitS * rTruHCalProtoHitS) + (rTruHCalProtoHitZ * rTruHCalProtoHitZ));
        const auto fTruHCalProtoHit    = boost::math::sign(rTruHCalProtoHitY) * acos(rTruHCalProtoHitX / rTruHCalProtoHitS);
        const auto tTruHCalProtoHit    = std::acos(rTruHCalProtoHitZ / rTruHCalProtoHitR);
        const auto hTruHCalProtoHit    = (-1.) * std::log(std::atan(tTruHCalProtoHit / 2.));
        const auto diffTruHCalProtoHit = (eTruHCalProtoHit - eMcPar) / eMcPar;

        // fill hit histograms and increment sums/counters
        hHCalTruClustHitPhi      -> Fill(fTruHCalProtoHit);
        hHCalTruClustHitEta      -> Fill(hTruHCalProtoHit);
        hHCalTruClustHitEne      -> Fill(eTruHCalProtoHit);
        hHCalTruClustHitPosZ     -> Fill(rTruHCalProtoHitZ);
        hHCalTruClustHitParDiff  -> Fill(diffTruHCalProtoHit);
        hHCalTruClustHitPosYvsX  -> Fill(rTruHCalProtoHitX, rTruHCalProtoHitY);
        hHCalTruClustHitEtaVsPhi -> Fill(fTruHCalProtoHit, hTruHCalProtoHit);
        hHCalTruClustHitVsParEne -> Fill(eMcPar, eTruHCalProtoHit);
      }
      ++nTruHCalProto;
    }  // end protocluster loop

    // fill cluster histograms and increment counters
    hHCalTruClustPhi      -> Fill(fTruHCalClust);
    hHCalTruClustEta      -> Fill(hTruHCalClust);
    hHCalTruClustEne      -> Fill(eTruHCalClust);
    hHCalTruClustPosZ     -> Fill(rTruHCalClustZ);
    hHCalTruClustNumHit   -> Fill(nHitTruHCalClust);
    hHCalTruClustParDiff  -> Fill(diffTruHCalClust);
    hHCalTruClustPosYvsX  -> Fill(rTruHCalClustX, rTruHCalClustY);
    hHCalTruClustEtaVsPhi -> Fill(fTruHCalClust, hTruHCalClust);
    hHCalTruClustVsParEne -> Fill(eMcPar, eTruHCalClust);
    eTruHCalClustSum += eTruHCalClust;
    ++nTruHCalClust;

    // select leading cluster
    const bool isBiggerEne = (eTruHCalClust > eLeadTruHCalClust);
    if (isBiggerEne) {
      iLeadTruHCalClust    = iTruHCalClust;
      nHitLeadTruHCalClust = nHitTruHCalClust;
      eLeadTruHCalClust    = eTruHCalClust;
      diffLeadTruHCalClust = diffTruHCalClust;
    }
    ++iTruHCalClust;
  }  // end true bhcal cluster loop

  // for scifi/image hit sums
  double eSciFiHitSum(0.);
  double eImageHitSum(0.);

  double eSciFiHitSumVsNLayer[CONST::NSciFiLayer];
  double eImageHitSumVsNLayer[CONST::NImageLayer];
  for (size_t iSciFi = 0; iSciFi < CONST::NSciFiLayer; iSciFi++) {
    eSciFiHitSumVsNLayer[iSciFi] = 0.;
  }
  for (size_t iImage = 0; iImage < CONST::NImageLayer; iImage++) {
    eImageHitSumVsNLayer[iImage] = 0.;
  }

  // reco. scifi hit loop
  unsigned long nSciFiHit(0);
  for (auto scifiHit : scifiRecHits()) {

    // grab hit properties
    const auto nLayerSciFi  = scifiHit -> getLayer();
    const auto iLayerSciFi  = nLayerSciFi - 1;
    const auto rSciFiHitX   = scifiHit -> getPosition().x;
    const auto rSciFiHitY   = scifiHit -> getPosition().y;
    const auto rSciFiHitZ   = scifiHit -> getPosition().z;
    const auto eSciFiHit    = scifiHit -> getEnergy();
    const auto rSciFiHitS   = std::sqrt((rSciFiHitX * rSciFiHitX) + (rSciFiHitY * rSciFiHitY));
    const auto rSciFiHitR   = std::sqrt((rSciFiHitS * rSciFiHitS) + (rSciFiHitZ * rSciFiHitZ));
    const auto fSciFiHit    = boost::math::sign(rSciFiHitY) * acos(rSciFiHitX / rSciFiHitS);
    const auto tSciFiHit    = std::acos(rSciFiHitZ / rSciFiHitR);
    const auto hSciFiHit    = (-1.) * std::log(std::atan(tSciFiHit / 2.));
    const auto diffSciFiHit = (eSciFiHit - eMcPar) / eMcPar;

    // fill hit histograms
    hSciFiRecHitNLayer      -> Fill(nLayerSciFi);
    hSciFiRecHitPhi         -> Fill(fSciFiHit);
    hSciFiRecHitEta         -> Fill(hSciFiHit);
    hSciFiRecHitEne         -> Fill(eSciFiHit);
    hSciFiRecHitPosZ        -> Fill(rSciFiHitZ);
    hSciFiRecHitParDiff     -> Fill(diffSciFiHit);
    hSciFiRecHitPosYvsX     -> Fill(rSciFiHitX,  rSciFiHitY);
    hSciFiRecHitEtaVsPhi    -> Fill(fSciFiHit,   hSciFiHit);
    hSciFiRecHitVsParEne    -> Fill(eMcPar,      eSciFiHit);
    hSciFiRecHitEneVsNLayer -> Fill(nLayerSciFi, eSciFiHit);

    // increment sums/counters
    eSciFiHitSumVsNLayer[iLayerSciFi] += eSciFiHit;
    eSciFiHitSum                      += eSciFiHit;
    ++nSciFiHit;
  }  // end scifi hit loop

  // reco. image hit loop
  unsigned long nImageHit(0);
  for (auto imageHit : imageRecHits()) {

    // grab hit properties
    const auto nLayerImage  = imageHit -> getLayer();
    const auto iLayerImage  = nLayerImage - 1;
    const auto rImageHitX   = imageHit -> getPosition().x;
    const auto rImageHitY   = imageHit -> getPosition().y;
    const auto rImageHitZ   = imageHit -> getPosition().z;
    const auto eImageHit    = imageHit -> getEnergy();
    const auto rImageHitS   = std::sqrt((rImageHitX * rImageHitX) + (rImageHitY * rImageHitY));
    const auto rImageHitR   = std::sqrt((rImageHitS * rImageHitS) + (rImageHitZ * rImageHitZ));
    const auto fImageHit    = boost::math::sign(rImageHitY) * acos(rImageHitX / rImageHitS);
    const auto tImageHit    = std::acos(rImageHitZ / rImageHitR);
    const auto hImageHit    = (-1.) * std::log(std::atan(tImageHit / 2.));
    const auto diffImageHit = (eImageHit - eMcPar) / eMcPar;

    // fill hit histograms
    hImageRecHitNLayer      -> Fill(nLayerImage);
    hImageRecHitPhi         -> Fill(fImageHit);
    hImageRecHitEta         -> Fill(hImageHit);
    hImageRecHitEne         -> Fill(eImageHit);
    hImageRecHitPosZ        -> Fill(rImageHitZ);
    hImageRecHitParDiff     -> Fill(diffImageHit);
    hImageRecHitPosYvsX     -> Fill(rImageHitX,  rImageHitY);
    hImageRecHitEtaVsPhi    -> Fill(fImageHit,   hImageHit);
    hImageRecHitVsParEne    -> Fill(eMcPar,      eImageHit);
    hImageRecHitEneVsNLayer -> Fill(nLayerImage, eImageHit);

    // increment sums/counters
    eImageHitSumVsNLayer[iLayerImage] += eImageHit;
    eImageHitSum                      += eImageHit;
    ++nImageHit;
  }  // end scifi hit loop

  // for highest energy bemc clusters
  int    iLeadECalClust(-1);
  int    iLeadSciFiClust(-1);
  int    iLeadImageClust(-1);
  int    nHitLeadECalClust(-1);
  int    nHitLeadImageClust(-1);
  int    nHitLeadSciFiClust(-1);
  double hLeadECalClust(-999.);
  double hLeadImageClust(-999.);
  double hLeadSciFiClust(-999.);
  double fLeadECalClust(-999.);
  double fLeadImageClust(-999.);
  double fLeadSciFiClust(-999.);
  double eLeadECalClust(-999.);
  double eLeadImageClust(0.);
  double eLeadSciFiClust(0.);
  double diffLeadECalClust(-999.);
  double diffLeadImageClust(-999.);
  double diffLeadSciFiClust(-999.);

  // reco. bemc cluster loop
  unsigned long iECalClust(0);
  unsigned long iImageClust(0);
  unsigned long iSciFiClust(0);
  unsigned long nECalClust(0);
  unsigned long nImageClust(0);
  unsigned long nSciFiClust(0);
  for (auto bemcClust : bemcClusters()) {

    // grab cluster properties
    const auto rECalClustX   = bemcClust -> getPosition().x;
    const auto rECalClustY   = bemcClust -> getPosition().y;
    const auto rECalClustZ   = bemcClust -> getPosition().z;
    const auto eECalClust    = bemcClust -> getEnergy();
    const auto nHitECalClust = bemcClust -> getNhits();
    const auto tECalClust    = bemcClust -> getIntrinsicTheta();
    const auto diffECalClust = (eECalClust - eMcPar) / eMcPar;

    // calculate cluster eta, phi
    // FIXME update to XYZvectors
    const TVector3 vecPosition(rECalClustX, rECalClustY, rECalClustZ);
    const auto     hECalClust = vecPosition.Eta();
    const auto     fECalClust = vecPosition.Phi();

    // fill cluster histograms and increment counters
    hECalClustPhi      -> Fill(fECalClust);
    hECalClustEta      -> Fill(hECalClust);
    hECalClustEne      -> Fill(eECalClust);
    hECalClustPosZ     -> Fill(rECalClustZ);
    hECalClustNumHit   -> Fill(nHitECalClust);
    hECalClustParDiff  -> Fill(diffECalClust);
    hECalClustPosYvsX  -> Fill(rECalClustX, rECalClustY);
    hECalClustEtaVsPhi -> Fill(fECalClust, hECalClust);
    hECalClustVsParEne -> Fill(eMcPar, eECalClust);
    eECalClustSum += eECalClust;
    ++nECalClust;
    ++iECalClust;

    // select leading cluster
    const bool isBiggerEne = (eECalClust > eLeadECalClust);
    if (isBiggerEne) {
      iLeadECalClust    = iECalClust;
      nHitLeadECalClust = nHitECalClust;
      hLeadECalClust    = hECalClust;
      fLeadECalClust    = fECalClust;
      eLeadECalClust    = eECalClust;
      diffLeadECalClust = diffECalClust;
    }
  }  // end reco. bemc cluster loop

  // loop over scifi clusters
  for (auto scifiClust : scifiClusters()) {

    // grab cluster properties
    const auto rSciFiClustX   = scifiClust -> getPosition().x;
    const auto rSciFiClustY   = scifiClust -> getPosition().y;
    const auto rSciFiClustZ   = scifiClust -> getPosition().z;
    const auto eSciFiClust    = scifiClust -> getEnergy();
    const auto nHitSciFiClust = scifiClust -> getNhits();
    const auto tSciFiClust    = scifiClust -> getIntrinsicTheta();
    const auto diffSciFiClust = (eSciFiClust - eMcPar) / eMcPar;

    // calculate cluster eta, phi
    // FIXME update to XYZvectors
    const TVector3 vecPosition(rSciFiClustX, rSciFiClustY, rSciFiClustZ);
    const auto     hSciFiClust = vecPosition.Eta();
    const auto     fSciFiClust = vecPosition.Phi();

    // increment counters
    eSciFiClustSum += eSciFiClust;
    ++nSciFiClust;
    ++iSciFiClust;

    // select leading cluster
    const bool isBiggerEne = (eSciFiClust > eLeadSciFiClust);
    if (isBiggerEne) {
      iLeadSciFiClust    = iSciFiClust;
      nHitLeadSciFiClust = nHitSciFiClust;
      hLeadSciFiClust    = hSciFiClust;
      fLeadSciFiClust    = fSciFiClust;
      eLeadSciFiClust    = eSciFiClust;
      diffLeadSciFiClust = diffSciFiClust;
    }
  }  // end scifi cluster loop

  // loop over imaging clusters
  for (auto imageClust : imageClusters()) {

    // grab cluster properties
    const auto rImageClustX   = imageClust -> getPosition().x;
    const auto rImageClustY   = imageClust -> getPosition().y;
    const auto rImageClustZ   = imageClust -> getPosition().z;
    const auto eImageClust    = imageClust -> getEnergy();
    const auto nHitImageClust = imageClust -> getNhits();
    const auto tImageClust    = imageClust -> getIntrinsicTheta();
    const auto diffImageClust = (eImageClust - eMcPar) / eMcPar;

    // calculate cluster eta, phi
    // FIXME update to XYZvectors
    const TVector3 vecPosition(rImageClustX, rImageClustY, rImageClustZ);
    const auto     hImageClust = vecPosition.Eta();
    const auto     fImageClust = vecPosition.Phi();

    // increment counters
    eImageClustSum += eImageClust;
    ++nImageClust;
    ++iImageClust;

    // select leading cluster
    const bool isBiggerEne = (eImageClust > eLeadImageClust);
    if (isBiggerEne) {
      iLeadImageClust    = iImageClust;
      nHitLeadImageClust = nHitImageClust;
      hLeadImageClust    = hImageClust;
      fLeadImageClust    = fImageClust;
      eLeadImageClust    = eImageClust;
      diffLeadImageClust = diffImageClust;
    }
  }  // end imaging cluster loop

  // do event-wise calculations
  const auto fracParVsLeadHCal   = eLeadHCalClust / eMcPar;
  const auto fracParVsLeadECal   = eLeadECalClust / eMcPar;
  const auto fracParVsSumHCal    = eHCalClustSum / eMcPar;
  const auto fracParVsSumECal    = eECalClustSum / eMcPar;
  const auto fracLeadHCalVsECal  = eLeadECalClust / (eLeadHCalClust + eLeadECalClust);
  const auto fracSumHCalVsECal   = eECalClustSum / (eHCalClustSum + eECalClustSum);
  const auto diffHCalHitSum      = (eHCalHitSum - eMcPar) / eMcPar;
  const auto diffHCalClustSum    = (eHCalClustSum - eMcPar) / eMcPar;
  const auto diffECalClustSum    = (eECalClustSum - eMcPar) / eMcPar;
  const auto diffTruHCalClustSum = (eTruHCalClustSum - eMcPar) / eMcPar;

  // fill general event-wise bhcal histograms
  hEvtHCalNumPar             -> Fill(nPar);
  // fill hit event-wise bhcal histograms
  hEvtHCalNumHit             -> Fill(nHCalHit);
  hEvtHCalSumHitEne          -> Fill(eHCalHitSum);
  hEvtHCalSumHitDiff         -> Fill(diffHCalHitSum);
  hEvtHCalSumHitVsPar        -> Fill(eMcPar, eHCalHitSum);
  // fill cluster event-wise bhcal histograms
  hEvtHCalNumClust           -> Fill(nHCalClust);
  hEvtHCalSumClustEne        -> Fill(eHCalClustSum);
  hEvtHCalSumClustDiff       -> Fill(diffHCalClustSum);
  hEvtHCalNumClustVsHit      -> Fill(nHCalHit, nHCalClust);
  hEvtHCalSumClustVsPar      -> Fill(eMcPar,   eHCalClustSum);
  // fill lead cluster event-wise bhcal histograms
  hEvtHCalLeadClustNumHit    -> Fill(nHitLeadHCalClust);
  hEvtHCalLeadClustEne       -> Fill(eLeadHCalClust);
  hEvtHCalLeadClustDiff      -> Fill(diffLeadHCalClust);
  hEvtHCalLeadClustVsPar     -> Fill(eMcPar, eLeadHCalClust);
  // fill truth cluster event-wise bhcal histograms
  hEvtHCalNumTruClust        -> Fill(nTruHCalClust);
  hEvtHCalSumTruClustEne     -> Fill(eTruHCalClustSum);
  hEvtHCalSumTruClustDiff    -> Fill(diffTruHCalClustSum);
  hEvtHCalNumTruClustVsClust -> Fill(nHCalClust, nTruHCalClust);
  hEvtHCalSumTruClustVsPar   -> Fill(eMcPar,     eTruHCalClustSum);
  // fill lead truth cluster event-wise bhcal histograms
  hEvtHCalLeadTruClustNumHit -> Fill(nHitLeadTruHCalClust);
  hEvtHCalLeadTruClustEne    -> Fill(eLeadTruHCalClust);
  hEvtHCalLeadTruClustDiff   -> Fill(diffLeadTruHCalClust);
  hEvtHCalLeadTruClustVsPar  -> Fill(eMcPar, eLeadTruHCalClust);

  // fill hit event-wise scifi histograms
  hEvtSciFiSumEne          -> Fill(eSciFiHitSum);
  hEvtSciFiVsHCalHitSumEne -> Fill(eHCalHitSum, eSciFiHitSum);
  for (size_t iSciFi = 0; iSciFi < CONST::NSciFiLayer; iSciFi++) {
    hEvtSciFiSumEneVsNLayer -> Fill(iSciFi + 1, eSciFiHitSumVsNLayer[iSciFi]);
  }

  // fill hit event-wise image histograms
  hEvtImageSumEne          -> Fill(eImageHitSum);
  hEvtImageVsHCalHitSumEne -> Fill(eHCalHitSum, eImageHitSum);
  for (size_t iImage = 0; iImage < CONST::NImageLayer; iImage++) {
    hEvtImageSumEneVsNLayer -> Fill(iImage + 1, eImageHitSumVsNLayer[iImage]);
  }

  // fill cluster event-wise bhcal histograms
  hEvtECalNumClust          -> Fill(nECalClust);
  hEvtECalSumClustEne       -> Fill(eECalClustSum);
  hEvtECalSumClustDiff      -> Fill(diffECalClustSum);
  hEvtECalSumClustVsPar     -> Fill(eMcPar,        eECalClustSum);
  hEvtECalVsHCalSumClustEne -> Fill(eHCalClustSum, eECalClustSum);
  // fill lead cluster event-wise bhcal histograms
  hEvtECalLeadClustNumHit    -> Fill(nHitLeadECalClust);
  hEvtECalLeadClustEne       -> Fill(eLeadECalClust);
  hEvtECalLeadClustDiff      -> Fill(diffLeadECalClust);
  hEvtECalLeadClustVsPar     -> Fill(eMcPar,         eLeadECalClust);
  hEvtECalVsHCalLeadClustEne -> Fill(eLeadHCalClust, eLeadECalClust);

  // set variables for calibration tuple
  varsForCalibration[0]  = (Float_t) eMcPar;
  varsForCalibration[1]  = (Float_t) fracParVsLeadHCal;
  varsForCalibration[2]  = (Float_t) fracParVsLeadECal;
  varsForCalibration[3]  = (Float_t) fracParVsSumHCal;
  varsForCalibration[4]  = (Float_t) fracParVsSumECal;
  varsForCalibration[5]  = (Float_t) fracLeadHCalVsECal;
  varsForCalibration[6]  = (Float_t) fracSumHCalVsECal;
  varsForCalibration[7]  = (Float_t) eLeadHCalClust;
  varsForCalibration[8]  = (Float_t) eLeadECalClust;
  varsForCalibration[9]  = (Float_t) eHCalClustSum;
  varsForCalibration[10] = (Float_t) eECalClustSum;
  varsForCalibration[11] = (Float_t) diffLeadHCalClust;
  varsForCalibration[12] = (Float_t) diffLeadECalClust;
  varsForCalibration[13] = (Float_t) diffHCalClustSum;
  varsForCalibration[14] = (Float_t) diffECalClustSum;
  varsForCalibration[15] = (Float_t) nHitLeadHCalClust;
  varsForCalibration[16] = (Float_t) nHitLeadECalClust;
  varsForCalibration[17] = (Float_t) nHCalClust;
  varsForCalibration[18] = (Float_t) nECalClust;
  varsForCalibration[19] = (Float_t) hLeadHCalClust;
  varsForCalibration[20] = (Float_t) hLeadECalClust;
  varsForCalibration[21] = (Float_t) fLeadHCalClust;
  varsForCalibration[22] = (Float_t) fLeadECalClust;
  varsForCalibration[23] = (Float_t) eLeadImageClust;
  varsForCalibration[24] = (Float_t) eImageClustSum;
  varsForCalibration[25] = (Float_t) eLeadSciFiClust;
  varsForCalibration[26] = (Float_t) eSciFiClustSum;
  varsForCalibration[27] = (Float_t) nSciFiClust;
  varsForCalibration[28] = (Float_t) nImageClust;
  varsForCalibration[29] = (Float_t) hLeadImageClust;
  varsForCalibration[30] = (Float_t) hLeadSciFiClust;
  varsForCalibration[31] = (Float_t) fLeadImageClust;
  varsForCalibration[32] = (Float_t) fLeadSciFiClust;
  varsForCalibration[33] = (Float_t) eSciFiHitSumVsNLayer[0];
  varsForCalibration[34] = (Float_t) eSciFiHitSumVsNLayer[1];
  varsForCalibration[35] = (Float_t) eSciFiHitSumVsNLayer[2];
  varsForCalibration[36] = (Float_t) eSciFiHitSumVsNLayer[3];
  varsForCalibration[37] = (Float_t) eSciFiHitSumVsNLayer[4];
  varsForCalibration[38] = (Float_t) eSciFiHitSumVsNLayer[5];
  varsForCalibration[39] = (Float_t) eSciFiHitSumVsNLayer[6];
  varsForCalibration[40] = (Float_t) eSciFiHitSumVsNLayer[7];
  varsForCalibration[41] = (Float_t) eSciFiHitSumVsNLayer[8];
  varsForCalibration[42] = (Float_t) eSciFiHitSumVsNLayer[9];
  varsForCalibration[43] = (Float_t) eSciFiHitSumVsNLayer[10];
  varsForCalibration[44] = (Float_t) eSciFiHitSumVsNLayer[11];
  varsForCalibration[45] = (Float_t) eImageHitSumVsNLayer[0];
  varsForCalibration[46] = (Float_t) eImageHitSumVsNLayer[1];
  varsForCalibration[47] = (Float_t) eImageHitSumVsNLayer[2];
  varsForCalibration[48] = (Float_t) eImageHitSumVsNLayer[3];
  varsForCalibration[49] = (Float_t) eImageHitSumVsNLayer[4];
  varsForCalibration[50] = (Float_t) eImageHitSumVsNLayer[5];

  // fill tuple
  ntForCalibration -> Fill(varsForCalibration);
  return;

}  // end 'ProcessSequential(std::shared_ptr<JEvent>&)'



//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void FillBHCalCalibrationTupleProcessor::FinishWithGlobalRootLock() {

  // generic axis titles
  const TString sCount("counts");

  // particle axis titles
  const TString sMass("m_{par} [GeV/c^{2}]");
  const TString sCharge("charge");
  const TString sPhiPar("#varphi_{par}");
  const TString sEtaPar("#eta_{Par}");
  const TString sEnePar("E_{par} [GeV]");
  const TString sMomPar("p_{par} [GeV/c]");
  const TString sMomParX("p_{x, par} [GeV/c]");
  const TString sMomParY("p_{y, par} [GeV/c]");
  const TString sMomParZ("p_{z, par} [GeV/c]");
  const TString sNumParEvt("N_{par} per event");

  // hit axis titles
  const TString sPosHitX("x_{hit} [mm]");
  const TString sPosHitY("y_{hit} [mm]");
  const TString sPosHitZ("z_{hit} [mm]");
  const TString sPhiHit("#varphi_{hit}");
  const TString sEtaHit("#eta_{hit}");
  const TString sEneHit("e_{hit} [GeV]");
  const TString sEneHitSum("E^{sum}_{hit} = #Sigmae_{hit} [GeV]");
  const TString sEneHitDiff("#Deltae_{hit} / e_{hit} = (e_{hit} - E_{par}) / e_{hit} [GeV]");
  const TString sEneHitSumDiff("#DeltaE^{sum}_{hit} / E^{sum}_{hit} = (E^{sum}_{hit} - E_{par}) / E^{sum}_{hit} [GeV]");
  const TString sNumHitEvt("N_{hit} per event");

  // reco. cluster axis titles
  const TString sPosClustX("x_{clust} [mm]");
  const TString sPosClustY("y_{clust} [mm]");
  const TString sPosClustZ("z_{clust} [mm]");
  const TString sEneClust("e_{clust} [GeV]");
  const TString sPhiClust("#varphi_{clust}");
  const TString sEtaClust("#eta_{clust}");
  const TString sEneClustSum("E^{sum}_{clust} = #Sigmae_{clust} [GeV]");
  const TString sEneClustDiff("#Deltae_{clust} / e_{clust} = (e_{clust} - E_{par}) / e_{clust} [GeV]");
  const TString sEneClustLead("E^{lead}_{clust} [GeV]");
  const TString sEneClustSumDiff("#DeltaE^{sum}_{clust} / E^{sum}_{clust} = (E^{sum}_{clust} - E_{par}) / E^{sum}_{clust} [GeV]");
  const TString sEneClustLeadDiff("#DeltaE^{lead}_{clust} / E^{lead}_{clust} = (E^{lead}_{clust} - E_{par}) / E^{lead}_{clust} [GeV]");
  const TString sNumHitClust("N_{hit} per cluster");
  const TString sNumClustEvt("N_{clust} per event");

  // truth cluster axis titles
  const TString sPosTruClustX("x_{truth clust} [mm]");
  const TString sPosTruClustY("y_{truth clust} [mm]");
  const TString sPosTruClustZ("z_{truth clust} [mm]");
  const TString sPhiTruClust("#varphi^{truth}_{clust}");
  const TString sEtaTruClust("#eta^{truth}_{clust}");
  const TString sEneTruClust("e^{truth}_{clust} [GeV]");
  const TString sEneTruClustDiff("#Deltae^{truth}_{clust} / e^{truth}_{clust} / (e^{truth}_{clust} - E_{par}) / e^{truth}_{clust} [GeV]");
  const TString sEneTruClustSum("E^{sum/truth}_{clust} = #Sigmae^{truth}_{clust} [GeV]");
  const TString sEneTruClustLead("E^{lead/truth}_{clust} [GeV]");
  const TString sEneTruClustSumDiff("#DeltaE^{sum/truth}_{clust} / E^{sum/truth}_{clust} = (E^{sum/truth}_{clust} - E_{par}) / E^{sum/truth}_{clust} [GeV]");
  const TString sEneTruClustLeadDiff("#DeltaE^{lead/truth}_{clust} / E^{lead/truth}_{clust} = (E^{lead/truth} _{clust} - E_{par}) / E^{lead/truth}_{clust} [GeV]");
  const TString sNumHitTruClust("N_{hit} per truth cluster");
  const TString sNumTruClustEvt("N_{truth clust} per event");

  // set particle axis titles
  hParChrg                   -> GetXaxis() -> SetTitle(sCharge.Data());
  hParChrg                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParMass                   -> GetXaxis() -> SetTitle(sMass.Data());
  hParMass                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParPhi                    -> GetXaxis() -> SetTitle(sPhiPar.Data());
  hParPhi                    -> GetYaxis() -> SetTitle(sCount.Data());
  hParEta                    -> GetXaxis() -> SetTitle(sEtaPar.Data());
  hParEta                    -> GetYaxis() -> SetTitle(sCount.Data());
  hParEne                    -> GetXaxis() -> SetTitle(sEnePar.Data());
  hParEne                    -> GetYaxis() -> SetTitle(sCount.Data());
  hParMom                    -> GetXaxis() -> SetTitle(sMomPar.Data());
  hParMom                    -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomX                   -> GetXaxis() -> SetTitle(sMomParX.Data());
  hParMomX                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomY                   -> GetXaxis() -> SetTitle(sMomParY.Data());
  hParMomY                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomZ                   -> GetXaxis() -> SetTitle(sMomParZ.Data());
  hParMomZ                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParEtaVsPhi               -> GetXaxis() -> SetTitle(sPhiPar.Data());
  hParEtaVsPhi               -> GetYaxis() -> SetTitle(sEtaPar.Data());
  hParEtaVsPhi               -> GetZaxis() -> SetTitle(sCount.Data());
  // set reco. hit bhcal axis titles
  hHCalRecHitPhi             -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalRecHitPhi             -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitEta             -> GetXaxis() -> SetTitle(sEtaHit.Data());
  hHCalRecHitEta             -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitEne             -> GetXaxis() -> SetTitle(sEneHit.Data());
  hHCalRecHitEne             -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitPosZ            -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hHCalRecHitPosZ            -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitParDiff         -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hHCalRecHitParDiff         -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitPosYvsX         -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hHCalRecHitPosYvsX         -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hHCalRecHitPosYvsX         -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalRecHitEtaVsPhi        -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalRecHitEtaVsPhi        -> GetYaxis() -> SetTitle(sEtaHit.Data());
  hHCalRecHitEtaVsPhi        -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalRecHitVsParEne        -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalRecHitVsParEne        -> GetYaxis() -> SetTitle(sEneHit.Data());
  hHCalRecHitVsParEne        -> GetZaxis() -> SetTitle(sCount.Data());
  // set cluster hit bhcal axis titles
  hHCalClustHitPhi           -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalClustHitPhi           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitEta           -> GetXaxis() -> SetTitle(sEtaHit.Data());
  hHCalClustHitEta           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitEne           -> GetXaxis() -> SetTitle(sEneHit.Data());
  hHCalClustHitEne           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitPosZ          -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hHCalClustHitPosZ          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitParDiff       -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hHCalClustHitParDiff       -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitPosYvsX       -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hHCalClustHitPosYvsX       -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hHCalClustHitPosYvsX       -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustHitEtaVsPhi      -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalClustHitEtaVsPhi      -> GetYaxis() -> SetTitle(sEtaHit.Data());
  hHCalClustHitEtaVsPhi      -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustHitVsParEne      -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalClustHitVsParEne      -> GetYaxis() -> SetTitle(sEneHit.Data());
  hHCalClustHitVsParEne      -> GetZaxis() -> SetTitle(sCount.Data());
  // set reco. cluster bhcal axis titles
  hHCalClustPhi              -> GetXaxis() -> SetTitle(sPhiClust.Data());
  hHCalClustPhi              -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustEta              -> GetXaxis() -> SetTitle(sEtaClust.Data());
  hHCalClustEta              -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustEne              -> GetXaxis() -> SetTitle(sEneClust.Data());
  hHCalClustEne              -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustPosZ             -> GetXaxis() -> SetTitle(sPosClustZ.Data());
  hHCalClustPosZ             -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustNumHit           -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hHCalClustNumHit           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustParDiff          -> GetXaxis() -> SetTitle(sEneClustDiff.Data());
  hHCalClustParDiff          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustPosYvsX          -> GetXaxis() -> SetTitle(sPosClustX.Data());
  hHCalClustPosYvsX          -> GetYaxis() -> SetTitle(sPosClustY.Data());
  hHCalClustPosYvsX          -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustEtaVsPhi         -> GetXaxis() -> SetTitle(sPhiClust.Data());
  hHCalClustEtaVsPhi         -> GetYaxis() -> SetTitle(sEtaClust.Data());
  hHCalClustEtaVsPhi         -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustVsParEne         -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalClustVsParEne         -> GetYaxis() -> SetTitle(sEneClust.Data());
  hHCalClustVsParEne         -> GetZaxis() -> SetTitle(sCount.Data());
  // set truth cluster bhcal axis titles
  hHCalTruClustPhi           -> GetXaxis() -> SetTitle(sPhiTruClust.Data());
  hHCalTruClustPhi           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustEta           -> GetXaxis() -> SetTitle(sEtaTruClust.Data());
  hHCalTruClustEta           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustEne           -> GetXaxis() -> SetTitle(sEneTruClust.Data());
  hHCalTruClustEne           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustPosZ          -> GetXaxis() -> SetTitle(sPosTruClustZ.Data());
  hHCalTruClustPosZ          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustNumHit        -> GetXaxis() -> SetTitle(sNumHitTruClust.Data());
  hHCalTruClustNumHit        -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustParDiff       -> GetXaxis() -> SetTitle(sEneTruClustDiff.Data());
  hHCalTruClustParDiff       -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustPosYvsX       -> GetXaxis() -> SetTitle(sPosTruClustX.Data());
  hHCalTruClustPosYvsX       -> GetYaxis() -> SetTitle(sPosTruClustY.Data());
  hHCalTruClustPosYvsX       -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalTruClustEtaVsPhi      -> GetXaxis() -> SetTitle(sPhiTruClust.Data());
  hHCalTruClustEtaVsPhi      -> GetYaxis() -> SetTitle(sEtaTruClust.Data());
  hHCalTruClustEtaVsPhi      -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalTruClustVsParEne      -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalTruClustVsParEne      -> GetYaxis() -> SetTitle(sEneTruClust.Data());
  hHCalTruClustVsParEne      -> GetZaxis() -> SetTitle(sCount.Data());
  // set general event-wise bhcal axis titles
  hEvtHCalNumPar             -> GetXaxis() -> SetTitle(sNumParEvt.Data());
  hEvtHCalNumPar             -> GetYaxis() -> SetTitle(sCount.Data());
  // set hit event-wise bhcal axis titles
  hEvtHCalNumHit             -> GetXaxis() -> SetTitle(sNumHitEvt.Data());
  hEvtHCalNumHit             -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumHitEne          -> GetXaxis() -> SetTitle(sEneHitSum.Data());
  hEvtHCalSumHitEne          -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumHitDiff         -> GetXaxis() -> SetTitle(sEneHitSumDiff.Data());
  hEvtHCalSumHitDiff         -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumHitVsPar        -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalSumHitVsPar        -> GetYaxis() -> SetTitle(sEneHitSum.Data());
  hEvtHCalSumHitVsPar        -> GetZaxis() -> SetTitle(sCount.Data());
  // set cluster event-wise bhcal axis titles
  hEvtHCalNumClust           -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtHCalNumClust           -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumClustEne        -> GetXaxis() -> SetTitle(sEneClustSum.Data());
  hEvtHCalSumClustEne        -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumClustDiff       -> GetXaxis() -> SetTitle(sEneClustSumDiff.Data());
  hEvtHCalSumClustDiff       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalNumClustVsHit      -> GetXaxis() -> SetTitle(sNumHitEvt.Data());
  hEvtHCalNumClustVsHit      -> GetYaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtHCalNumClustVsHit      -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumClustVsPar      -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalSumClustVsPar      -> GetYaxis() -> SetTitle(sEneClustSum.Data());
  hEvtHCalSumClustVsPar      -> GetZaxis() -> SetTitle(sCount.Data());
  // set lead cluster event-wise bhcal axis titles
  hEvtHCalLeadClustNumHit    -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hEvtHCalLeadClustNumHit    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadClustEne       -> GetXaxis() -> SetTitle(sEneClustLead.Data());
  hEvtHCalLeadClustEne       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadClustDiff      -> GetXaxis() -> SetTitle(sEneClustLeadDiff.Data());
  hEvtHCalLeadClustDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadClustVsPar     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalLeadClustVsPar     -> GetYaxis() -> SetTitle(sEneClustLead.Data());
  hEvtHCalLeadClustVsPar     -> GetZaxis() -> SetTitle(sCount.Data());
  // set truth cluster event-wise bhcal axis titles
  hEvtHCalNumTruClust        -> GetXaxis() -> SetTitle(sNumTruClustEvt.Data());
  hEvtHCalNumTruClust        -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumTruClustEne     -> GetXaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtHCalSumTruClustEne     -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumTruClustDiff    -> GetXaxis() -> SetTitle(sEneTruClustSumDiff.Data());
  hEvtHCalSumTruClustDiff    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumTruClustVsPar   -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalSumTruClustVsPar   -> GetYaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtHCalSumTruClustVsPar   -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtHCalNumTruClustVsClust -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtHCalNumTruClustVsClust -> GetYaxis() -> SetTitle(sNumTruClustEvt.Data());
  hEvtHCalNumTruClustVsClust -> GetZaxis() -> SetTitle(sCount.Data());
  // set truth lead cluster event-wise bhcal axis titles
  hEvtHCalLeadTruClustNumHit -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hEvtHCalLeadTruClustNumHit -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadTruClustEne    -> GetXaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtHCalLeadTruClustEne    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadTruClustDiff   -> GetXaxis() -> SetTitle(sEneTruClustLeadDiff.Data());
  hEvtHCalLeadTruClustDiff   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadTruClustVsPar  -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalLeadTruClustVsPar  -> GetYaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtHCalLeadTruClustVsPar  -> GetZaxis() -> SetTitle(sCount.Data());

  // TODO add hit histogram axis titles
  // set reco. cluster bemc axis titles
  hECalClustPhi           -> GetXaxis() -> SetTitle(sPhiClust.Data());
  hECalClustPhi           -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustEta           -> GetXaxis() -> SetTitle(sEtaClust.Data());
  hECalClustEta           -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustEne           -> GetXaxis() -> SetTitle(sEneClust.Data());
  hECalClustEne           -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustPosZ          -> GetXaxis() -> SetTitle(sPosClustZ.Data());
  hECalClustPosZ          -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustNumHit        -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hECalClustNumHit        -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustParDiff       -> GetXaxis() -> SetTitle(sEneClustDiff.Data());
  hECalClustParDiff       -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustPosYvsX       -> GetXaxis() -> SetTitle(sPosClustX.Data());
  hECalClustPosYvsX       -> GetYaxis() -> SetTitle(sPosClustY.Data());
  hECalClustPosYvsX       -> GetZaxis() -> SetTitle(sCount.Data());
  hECalClustEtaVsPhi      -> GetXaxis() -> SetTitle(sPhiClust.Data());
  hECalClustEtaVsPhi      -> GetYaxis() -> SetTitle(sEtaClust.Data());
  hECalClustEtaVsPhi      -> GetZaxis() -> SetTitle(sCount.Data());
  hECalClustVsParEne      -> GetXaxis() -> SetTitle(sEnePar.Data());
  hECalClustVsParEne      -> GetYaxis() -> SetTitle(sEneClust.Data());
  hECalClustVsParEne      -> GetZaxis() -> SetTitle(sCount.Data());
  // set cluster event-wise bemc axis titles
  hEvtECalNumClust        -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtECalNumClust        -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumClustEne     -> GetXaxis() -> SetTitle(sEneClustSum.Data());
  hEvtECalSumClustEne     -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumClustDiff    -> GetXaxis() -> SetTitle(sEneClustSumDiff.Data());
  hEvtECalSumClustDiff    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumClustVsPar   -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalSumClustVsPar   -> GetYaxis() -> SetTitle(sEneClustSum.Data());
  hEvtECalSumClustVsPar   -> GetZaxis() -> SetTitle(sCount.Data());
  // set lead cluster event-wise bemc axis titles
  hEvtECalLeadClustNumHit -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hEvtECalLeadClustNumHit -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadClustEne    -> GetXaxis() -> SetTitle(sEneClustLead.Data());
  hEvtECalLeadClustEne    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadClustDiff   -> GetXaxis() -> SetTitle(sEneClustLeadDiff.Data());
  hEvtECalLeadClustDiff   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadClustVsPar  -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalLeadClustVsPar  -> GetYaxis() -> SetTitle(sEneClustLead.Data());
  hEvtECalLeadClustVsPar  -> GetZaxis() -> SetTitle(sCount.Data());
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
