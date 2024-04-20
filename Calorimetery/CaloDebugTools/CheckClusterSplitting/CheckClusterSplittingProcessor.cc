// ----------------------------------------------------------------------------
// 'CheckClusterSplittingProcessor.cc'
// Derek Anderson
// 04.11.2024
//
// EICrecon plugin to collect some metrics on cluster splitting
// across several calorimeters.
// ----------------------------------------------------------------------------

#define CHECKCLUSTERSPLITTING_CC

// plugin definition
#include "CheckClusterSplittingProcessor.h"

// the following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app -> Add( new CheckClusterSplittingProcessor );
  }
}

// make common namespaces implicit
using namespace std;



// public methods -------------------------------------------------------------

void CheckClusterSplittingProcessor::InitWithGlobalRootLock() {

  // create output directory
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc -> GetHistFile();
  rootfile -> mkdir("CheckClusterSplitting") -> cd();

  // turn on errors
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  // create output histograms
  BuildEvtHistograms();
  BuildTrkHistograms();
  BuildCalHistograms();
  return;

}  // end 'InitWithGlobalRootLock()'



void CheckClusterSplittingProcessor::ProcessSequential(const shared_ptr<const JEvent>& event) {

  // grab particles and track projections
  const auto& particles   = *(event -> GetCollection<edm4eic::ReconstructedParticle>(m_config.parCollect.data()));
  const auto& projections = *(event -> GetCollection<edm4eic::TrackSegment>(m_config.projCollect.data()));

  // get particle energy
  double ePar = 0.;
  for (auto particle : particles) {
    if (particle.getType() == 1) {
      ePar = particle.getEnergy();
      break;
    }
  }  // end particle loop

  // loop over calorimeters
  for (int iCalo = 0; iCalo < m_config.vecCalos.size(); iCalo++) {

    // grab collection
    const auto& clusters = *(event -> GetCollection<edm4eic::Cluster>(m_config.vecCalos[iCalo].first.data()));
    const int   idCalo   = GetCaloID(iCalo);

    // grab relevant track projections
    GetProjections(projections, idCalo);
    if (m_vecTrkProj.size() <= 0) continue;

    // determine lead cluster energy and cluster sum
    double eSum  = 0.;
    double eLead = 0.;
    for (auto cluster : clusters) {
      if (cluster.getEnergy() > eLead) {
        eLead = cluster.getEnergy();
      }
      eSum += cluster.getEnergy();
    }  // end cluster loop
    if (eSum <= 0.) continue;

    // run calculations
    for (auto cluster : clusters) {

      // create struct to hold histogram content
      Hist::CalContent cContent;
      cContent.GetInfo(cluster);

      // skip hit-less clusters
      if (cluster.getHits().size() <= 0) continue;

      // calculate spatial extent of cluster
      FPair widths   = GetClusterWidths(cluster);
      FPair etaRange = GetEtaRange(cluster);
      FPair phiRange = GetPhiRange(cluster);

      // find closest projection
      bool    foundMatch = false;
      double  drMin      = 999.;
      double  drSigMin   = 999.;
      int32_t nProj      = 0;
      int32_t iMatch     = -1;
      int32_t iProject   = 0;
      for (auto project : m_vecTrkProj) {

        // get projection eta/phi
        const double projPhi = atan2(project.position.y, project.position.x);
        const double projEta = edm4hep::utils::eta(project.position);

        // check if in cluster eta/phi extent
        const bool isInEtaRange = ((projEta >= etaRange.first) && (projEta < etaRange.second));
        const bool isInPhiRange = ((projPhi >= phiRange.first) && (projPhi < phiRange.second));
        if (!isInEtaRange || !isInPhiRange) {
          continue;
        } else {
          ++nProj;
        }

        // grab projection info
        Hist::TrkContent tContent;
        tContent.GetInfo(project);

        // get distances
        const double dEta = (projEta - cContent.eta);
        const double dPhi = (projPhi - cContent.phi);

        // grab remaining info
        tContent.dr     = hypot(dEta, dPhi);
        tContent.drSig  = hypot(dEta/widths.first, dPhi/widths.second);
        tContent.eOverP = tContent.p / cluster.getEnergy();

        // fill all projection histograms
        FillTrkHistograms(iCalo, Hist::TType::TrkAll, tContent);

        // if projection is closer to cluster, set current match
        if (tContent.dr <= drMin) {
          drMin      = tContent.dr;
          drSigMin   = tContent.drSig;
          iMatch     = iProject;
          foundMatch = true;
        }
        ++iProject;

      }  // end projection loop

      // grab info of matched cluster
      Hist::TrkContent mContent;
      if (foundMatch) {
        mContent.GetInfo( m_vecTrkProj.at(iMatch) );
        mContent.dr     = drMin;
        mContent.drSig  = drSigMin;
        mContent.eOverP = mContent.p / cluster.getEnergy();
        FillTrkHistograms(iCalo, Hist::TType::TrkMatch, mContent);
      }

      // grab cluster remaining information
      cContent.nProj  = nProj;
      cContent.frac   = cluster.getEnergy() / eSum;
      cContent.purity = cluster.getEnergy() / ePar;
      cContent.eOverP = foundMatch ? mContent.eOverP : -1.;

      // fill cluster histograms
      FillCalHistograms(iCalo, Hist::CType::CalAll, cContent);

    }  // end cluster loop

    // set event level-information
    Hist::EvtContent eContent;
    eContent.nClust = clusters.size();
    eContent.ePar   = ePar;
    eContent.eLead  = eLead;
    eContent.eTot   = eSum;
    eContent.lOverT = eLead / eSum;
    eContent.lOverP = eLead / ePar;
    eContent.tOverP = eSum / ePar;

    // fill event-level histograms
    FillEvtHistograms(iCalo, Hist::EType::EvtAll, eContent);

  }  // end calorimeter loop
  return;

}  // end 'ProcessSequential(shared_ptr<JEvent>&)'



void CheckClusterSplittingProcessor::FinishWithGlobalRootLock() {

  /* nothing to do */
  return;

}  // end 'FinishWithGlobalRootLock()'



// private methods ------------------------------------------------------------

void CheckClusterSplittingProcessor::BuildEvtHistograms() {

  // histogram binning/labales
  VecAxisDef vecEvtAxes = {
    make_pair("N_{clust}",                make_tuple(100, 0., 100.)),
    make_pair("E_{par} [GeV]",            make_tuple(100, 0., 50.)),
    make_pair("E_{lead} [GeV]",           make_tuple(100, 0., 50.)),
    make_pair("#SigmaE_{clust} [GeV]",    make_tuple(100, 0. ,50.)),
    make_pair("E_{lead}/#SigmaE_{clust}", make_tuple(250, 0., 5.)),
    make_pair("E_{lead}/E_{par}",         make_tuple(250, 0., 5.)),
    make_pair("#SigmaE_{clust}/E_{par}",  make_tuple(250, 0., 5.))
  };

  // type labels
  vector<string> vecEvtTypes = {
    "All"
  };

  // histogram definitions
  VecHistDef1D vecEvtDef1D = {
    make_tuple("EvtNClust",         vecEvtAxes[Hist::EVar::EvtNCl]),
    make_tuple("EvtParEne",         vecEvtAxes[Hist::EVar::EvtPE]),
    make_tuple("EvtLeadEne",        vecEvtAxes[Hist::EVar::EvtLE]),
    make_tuple("EvtTotalEne",       vecEvtAxes[Hist::EVar::EvtTE]),
    make_tuple("EvtLeadTotEneFrac", vecEvtAxes[Hist::EVar::EvtLT]),
    make_tuple("EvtLeadParEneFrac", vecEvtAxes[Hist::EVar::EvtLP]),
    make_tuple("EvtTotParEneFrac",  vecEvtAxes[Hist::EVar::EvtTP]),
  };

  // make histograms
  m_vecEvtH1D.resize( m_config.vecCalos.size() );
  for (size_t iCalo = 0; iCalo < m_config.vecCalos.size(); iCalo++) {

    m_vecEvtH1D[iCalo].resize( vecEvtTypes.size() );
    for (size_t iEType = 0; iEType < vecEvtTypes.size(); iEType++) {

      // make 1d histograms
      for (auto evtDef1D : vecEvtDef1D) {

        // make name
        string name("h");
        name += m_config.vecCalos[iCalo].second;
        name += vecEvtTypes[iEType];
        name += get<0>(evtDef1D);

        // make title
        string title(";");
        title += get<1>(evtDef1D).first;
        title += ";counts";

        // create histogram
        m_vecEvtH1D[iCalo][iEType].push_back(
          new TH1D(
            name.data(),
            title.data(),
            get<0>( get<1>(evtDef1D).second ),
            get<1>( get<1>(evtDef1D).second ),
            get<2>( get<1>(evtDef1D).second )
          )
        );
      }  // end 1d hist def loop
    }  // end type loop
  }  // end calo loop
  return;

}  // end 'BuildEvtHistograms()'



void CheckClusterSplittingProcessor::BuildTrkHistograms() {

  // histogram binning/labales
  VecAxisDef vecTrkAxes = {
    make_pair("r_{x} [mm]",      make_tuple(3000, -3000., 3000.)),
    make_pair("r_{y} [mm]",      make_tuple(3000, -3000., 3000.)),
    make_pair("r_{z} [mm]",      make_tuple(3000, -3000., 3000.)),
    make_pair("#eta",            make_tuple(100, -5., 5.)),
    make_pair("#varphi",         make_tuple(180, -3.15, 3.15)),
    make_pair("p [GeV/c]",       make_tuple(100, 0., 50.)),
    make_pair("#Deltar",         make_tuple(500, 0., 5.)),
    make_pair("#Delta#tilde{r}", make_tuple(500, 0., 5.)),
    make_pair("E/p",             make_tuple(250, 0., 5.))
  };

  // type labels
  vector<string> vecTrkTypes = {
    "All",
    "Match"
  };

  // histogram definitions
  VecHistDef1D vecTrkDef1D = {
    make_tuple("TrkPosX",      vecTrkAxes[Hist::TVar::TrkRX]),
    make_tuple("TrkPosY",      vecTrkAxes[Hist::TVar::TrkRY]),
    make_tuple("TrkPosZ",      vecTrkAxes[Hist::TVar::TrkRZ]),
    make_tuple("TrkEta",       vecTrkAxes[Hist::TVar::TrkH]),
    make_tuple("TrkPhi",       vecTrkAxes[Hist::TVar::TrkF]),
    make_tuple("TrkMom",       vecTrkAxes[Hist::TVar::TrkP]),
    make_tuple("TrkDeltaR",    vecTrkAxes[Hist::TVar::TrkDR]),
    make_tuple("TrkDeltaRSig", vecTrkAxes[Hist::TVar::TrkDRS]),
    make_tuple("TrkEoverP",    vecTrkAxes[Hist::TVar::TrkEP])
  };
  VecHistDef2D vecTrkDef2D = {
    make_tuple("TrkPosYvsX",   vecTrkAxes[Hist::TVar::TrkRX],  vecTrkAxes[Hist::TVar::TrkRY]),
    make_tuple("TrkMomVsEta",  vecTrkAxes[Hist::TVar::TrkH],   vecTrkAxes[Hist::TVar::TrkP]),
    make_tuple("TrkEtaVsPhi",  vecTrkAxes[Hist::TVar::TrkF],   vecTrkAxes[Hist::TVar::TrkH]),
    make_tuple("TrkEPvsDR",    vecTrkAxes[Hist::TVar::TrkDR],  vecTrkAxes[Hist::TVar::TrkEP]),
    make_tuple("TrkEPvsDRSig", vecTrkAxes[Hist::TVar::TrkDRS], vecTrkAxes[Hist::TVar::TrkEP])
  };

  // make histograms
  m_vecTrkH1D.resize( m_config.vecCalos.size() );
  m_vecTrkH2D.resize( m_config.vecCalos.size() );
  for (size_t iCalo = 0; iCalo < m_config.vecCalos.size(); iCalo++) {

    m_vecTrkH1D[iCalo].resize( vecTrkTypes.size() );
    m_vecTrkH2D[iCalo].resize( vecTrkTypes.size() );
    for (size_t iTType = 0; iTType < vecTrkTypes.size(); iTType++) {

      // make 1d histograms
      for (auto trkDef1D : vecTrkDef1D) {

        // make name
        string name("h");
        name += m_config.vecCalos[iCalo].second;
        name += vecTrkTypes[iTType];
        name += get<0>(trkDef1D);

        // make title
        string title(";");
        title += get<1>(trkDef1D).first;
        title += ";counts";

        // create histogram
        m_vecTrkH1D[iCalo][iTType].push_back(
          new TH1D(
            name.data(),
            title.data(),
            get<0>( get<1>(trkDef1D).second ),
            get<1>( get<1>(trkDef1D).second ),
            get<2>( get<1>(trkDef1D).second )
          )
        );
      }  // end 1d hist def loop

      // make 2d histograms
      for (auto trkDef2D : vecTrkDef2D) {

        // make name
        string name("h");
        name += m_config.vecCalos[iCalo].second;
        name += vecTrkTypes[iTType];
        name += get<0>(trkDef2D);

        // make title
        string title(";");
        title += get<1>(trkDef2D).first;
        title += ";";
        title += get<2>(trkDef2D).first;
        title += ";counts";

        // create histogram
        m_vecTrkH2D[iCalo][iTType].push_back(
          new TH2D(
            name.data(),
            title.data(),
            get<0>( get<1>(trkDef2D).second ),
            get<1>( get<1>(trkDef2D).second ),
            get<2>( get<1>(trkDef2D).second ),
            get<0>( get<2>(trkDef2D).second ),
            get<1>( get<2>(trkDef2D).second ),
            get<2>( get<2>(trkDef2D).second )
          )
        );
      }  // end 2d hist def loop
    }  // end type loop
  }  // end calo loop
  return;

}  // end 'BuildTrkHistograms()'



void CheckClusterSplittingProcessor::BuildCalHistograms() {

  // histogram binning/labales
  VecAxisDef vecCalAxes = {
    make_pair("N_{hit}",                   make_tuple(50, 0., 50.)),
    make_pair("N_{proj}",                  make_tuple(20, 0., 20.)),
    make_pair("r_{x} [mm]",                make_tuple(3000, -3000., 3000.)),
    make_pair("r_{y} [mm]",                make_tuple(3000, -3000., 3000.)),
    make_pair("r_{z} [mm]",                make_tuple(3000, -3000., 3000.)),
    make_pair("#eta",                      make_tuple(100, -5., 5.)),
    make_pair("#varphi",                   make_tuple(180, -3.15, 3.15)),
    make_pair("E [GeV]",                   make_tuple(100, 0., 50.)),
    make_pair("E_{clust}/#SigmaE_{clust}", make_tuple(250, 0., 5.)),
    make_pair("E_{clust}/E_{par}",         make_tuple(250, 0., 5.)),
    make_pair("E/p",                       make_tuple(250, 0., 5.))
  };

  // type labels
  vector<string> vecCalTypes = {
    "All"
  };

  // histogram definitions
  VecHistDef1D vecCalDef1D = {
    make_tuple("ClustNHit",    vecCalAxes[Hist::CVar::CalNH]),
    make_tuple("ClustNProj",   vecCalAxes[Hist::CVar::CalNP]),
    make_tuple("ClustPosX",    vecCalAxes[Hist::CVar::CalRX]),
    make_tuple("ClustPosY",    vecCalAxes[Hist::CVar::CalRY]),
    make_tuple("ClustPosZ",    vecCalAxes[Hist::CVar::CalRZ]),
    make_tuple("ClustEta",     vecCalAxes[Hist::CVar::CalH]),
    make_tuple("ClustPhi",     vecCalAxes[Hist::CVar::CalF]),
    make_tuple("ClustEne",     vecCalAxes[Hist::CVar::CalE]),
    make_tuple("ClustEpsilon", vecCalAxes[Hist::CVar::CalEps]),
    make_tuple("ClustPurity",  vecCalAxes[Hist::CVar::CalRho]),
    make_tuple("ClustEoverP",  vecCalAxes[Hist::CVar::CalEP])
  };
  VecHistDef2D vecCalDef2D = {
    make_tuple("ClustPosYvsX",  vecCalAxes[Hist::CVar::CalRX], vecCalAxes[Hist::CVar::CalRY]),
    make_tuple("ClustEneVsEta", vecCalAxes[Hist::CVar::CalH],  vecCalAxes[Hist::CVar::CalE]),
    make_tuple("ClustEneVsPhi", vecCalAxes[Hist::CVar::CalF],  vecCalAxes[Hist::CVar::CalF]),
    make_tuple("ClustEtaVsPhi", vecCalAxes[Hist::CVar::CalF],  vecCalAxes[Hist::CVar::CalH]),
    make_tuple("ClustPurVsEne", vecCalAxes[Hist::CVar::CalE],  vecCalAxes[Hist::CVar::CalRho])
  };

  // make histograms
  m_vecCalH1D.resize( m_config.vecCalos.size() );
  m_vecCalH2D.resize( m_config.vecCalos.size() );
  for (size_t iCalo = 0; iCalo < m_config.vecCalos.size(); iCalo++) {

    m_vecCalH1D[iCalo].resize( vecCalTypes.size() );
    m_vecCalH2D[iCalo].resize( vecCalTypes.size() );
    for (size_t iCType = 0; iCType < vecCalTypes.size(); iCType++) {

      // make 1d histograms
      for (auto calDef1D : vecCalDef1D) {

        // make name
        string name("h");
        name += m_config.vecCalos[iCalo].second;
        name += vecCalTypes[iCType];
        name += get<0>(calDef1D);

        // make title
        string title(";");
        title += get<1>(calDef1D).first;
        title += ";counts";

        // create histogram
        m_vecCalH1D[iCalo][iCType].push_back(
          new TH1D(
            name.data(),
            title.data(),
            get<0>( get<1>(calDef1D).second ),
            get<1>( get<1>(calDef1D).second ),
            get<2>( get<1>(calDef1D).second )
          )
        );
      }  // end 1d hist def loop

      // make 2d histograms
      for (auto calDef2D : vecCalDef2D) {

        // make name
        string name("h");
        name += m_config.vecCalos[iCalo].second;
        name += vecCalTypes[iCType];
        name += get<0>(calDef2D);

        // make title
        string title(";");
        title += get<1>(calDef2D).first;
        title += ";";
        title += get<2>(calDef2D).first;
        title += ";counts";

        // create histogram
        m_vecCalH2D[iCalo][iCType].push_back(
          new TH2D(
            name.data(),
            title.data(),
            get<0>( get<1>(calDef2D).second ),
            get<1>( get<1>(calDef2D).second ),
            get<2>( get<1>(calDef2D).second ),
            get<0>( get<2>(calDef2D).second ),
            get<1>( get<2>(calDef2D).second ),
            get<2>( get<2>(calDef2D).second )
          )
        );
      }  // end 2d hist def loop
    }  // end type label loop
  }  // end calo label loop
  return;

}  // end 'BuildHistograms()'



void CheckClusterSplittingProcessor::FillEvtHistograms(const int calo, const int type, const Hist::EvtContent& content) {

  // fill 1d hisgorams
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtNClust)     -> Fill(content.nClust);
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtParEne)     -> Fill(content.ePar);
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtLeadEne)    -> Fill(content.eLead);
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtTotEne)     -> Fill(content.eTot);
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtLeadDivTot) -> Fill(content.lOverT);
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtLeadDivPar) -> Fill(content.lOverP);
  m_vecEvtH1D.at(calo).at(type).at(Hist::E1D::EvtTotDivPar)  -> Fill(content.tOverP);
  return;

}  // end 'FillEvtHistograms(int, int, Hist::EvtContent&)'



void CheckClusterSplittingProcessor::FillTrkHistograms(const int calo, const int type, const Hist::TrkContent& content) {

  // fill 1d histograms
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkPosX)        -> Fill(content.rx);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkPosY)        -> Fill(content.ry);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkPosZ)        -> Fill(content.rz);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkEta)         -> Fill(content.eta);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkPhi)         -> Fill(content.phi);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkMom)         -> Fill(content.p);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkDeltaR)      -> Fill(content.dr);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkDeltaRScale) -> Fill(content.drSig);
  m_vecTrkH1D.at(calo).at(type).at(Hist::T1D::TrkEOverP)      -> Fill(content.eOverP);

  // fill 2d histograms
  m_vecTrkH2D.at(calo).at(type).at(Hist::T2D::TrkPosYvsX)   -> Fill(content.rx,    content.rx);
  m_vecTrkH2D.at(calo).at(type).at(Hist::T2D::TrkMomVsEta)  -> Fill(content.eta,   content.p);
  m_vecTrkH2D.at(calo).at(type).at(Hist::T2D::TrkEtaVsPhi)  -> Fill(content.phi,   content.eta);
  m_vecTrkH2D.at(calo).at(type).at(Hist::T2D::TrkEPvsDR)    -> Fill(content.dr,    content.eOverP);
  m_vecTrkH2D.at(calo).at(type).at(Hist::T2D::TrkEPvsDRSig) -> Fill(content.drSig, content.eOverP);
  return;

}  // end 'FillTrkHisotgrams(int, int, Hist::TrkContent&)'



void CheckClusterSplittingProcessor::FillCalHistograms(const int calo, const int type, const Hist::CalContent& content) {

  // fill 1d histograms
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalNHit)     -> Fill(content.nHit);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalNProj)    -> Fill(content.nProj);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalPosX)     -> Fill(content.rx);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalPosY)     -> Fill(content.ry);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalPosZ)     -> Fill(content.rz);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalEta)      -> Fill(content.eta);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalPhi)      -> Fill(content.phi);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalEne)      -> Fill(content.ene);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalFraction) -> Fill(content.frac);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalPurity)   -> Fill(content.purity);
  m_vecCalH1D.at(calo).at(type).at(Hist::C1D::CalEOverP)   -> Fill(content.eOverP);

  // fill 2d histograms
  m_vecCalH2D.at(calo).at(type).at(Hist::C2D::CalPosYvsX)   -> Fill(content.rx,  content.ry);
  m_vecCalH2D.at(calo).at(type).at(Hist::C2D::CalEneVsEta)  -> Fill(content.eta, content.ene);
  m_vecCalH2D.at(calo).at(type).at(Hist::C2D::CalEneVsPhi)  -> Fill(content.phi, content.ene);
  m_vecCalH2D.at(calo).at(type).at(Hist::C2D::CalEtaVsPhi)  -> Fill(content.phi, content.eta);
  m_vecCalH2D.at(calo).at(type).at(Hist::C2D::CalPurVsEne)  -> Fill(content.ene, content.purity);
  return;

}  // end 'FillCalHistograms(int, int, Hist::CalContent&)'



void CheckClusterSplittingProcessor::SaveHistograms() {

  // make directory and save histograms
  for (size_t iCalo = 0; iCalo < m_config.vecCalos.size(); iCalo++) {

    // save event histograms
    for (size_t iEType = 0; iEType < m_vecEvtH1D[iCalo].size(); iEType++) {
      for (auto evt1D : m_vecEvtH1D[iCalo][iEType]) {
        evt1D -> Write();
      }
    }

    // save track histograms
    for (size_t iTType = 0; iTType < m_vecTrkH1D[iCalo].size(); iTType++) {
      for (auto trk1D : m_vecTrkH1D[iCalo][iTType]) {
        trk1D -> Write();
      }
      for (auto trk2D : m_vecTrkH2D[iCalo][iTType]) {
        trk2D -> Write();
      }
    }

    // save cluster histograms
    for (size_t iCType = 0; iCType < m_vecTrkH1D[iCalo].size(); iCType++) {
      for (auto cal1D : m_vecCalH1D[iCalo][iCType]) {
        cal1D -> Write();
      }
      for (auto cal2D : m_vecCalH2D[iCalo][iCType]) {
        cal2D -> Write();
      }
    }
  }  // end calo loop
  return;

}  // end 'SaveHistograms()'


void CheckClusterSplittingProcessor::GetProjections(const edm4eic::TrackSegmentCollection& projections, const int calo) {

  // make sure vector is empty
  m_vecTrkProj.clear();

  // collect projections
  for (auto project : projections) {
    for (auto point : project.getPoints()) {
      const bool isInSystem  = (point.system  == calo);
      const bool isAtFace    = (point.surface == 1);
      if (isInSystem && isAtFace) {
        m_vecTrkProj.push_back(point);
      }
    }
  }
  return;

}  // end 'GetProjections(edm4eic::TrackSegmentCollection&)'


int CheckClusterSplittingProcessor::GetCaloID(const int iCalo) {

  // get detector element
  auto detector = GetApplication() -> GetService<DD4hep_service>() -> detector();
  auto element  = detector -> detector(m_config.vecCalos.at(iCalo).second);
  return element.id();

}  // end 'GetCaloID(int)'



FPair CheckClusterSplittingProcessor::GetClusterWidths(const edm4eic::Cluster& cluster) {

  // get eta, phi for cluster
  const float hClust = edm4hep::utils::eta( cluster.getPosition() );
  const float fClust = atan2( cluster.getPosition().y, cluster.getPosition().x );

  float hWidth2 = 0.;
  float fWidth2 = 0.;
  for (auto hit : cluster.getHits()) {

    // get eta, phi for hit
    const float hHit = edm4hep::utils::eta( hit.getPosition() );
    const float fHit = atan2( hit.getPosition().y, hit.getPosition().x );

    // increment sum of squares
    hWidth2 += pow(hHit - hClust, 2.);
    fWidth2 += pow(fHit - fClust, 2.);
  }

  // calculate widths
  const float hWidth = sqrt(hWidth2 / cluster.getHits().size());
  const float fWidth = sqrt(fWidth2 / cluster.getHits().size());

  // if nHit = 1, return some min size
  //   FIXME this should be the cell dimensions
  FPair widths  = make_pair(hWidth, fWidth);
  FPair minimum = make_pair(0.05,   0.05);
  return (cluster.getHits().size() > 1) ? widths : minimum;

}  // end 'GetClusterWidths(edm4eic::Cluster&)'



FPair CheckClusterSplittingProcessor::GetEtaRange(const edm4eic::Cluster& cluster) {

  // get hits
  auto hits = cluster.getHits();

  // get hit with min eta
  auto minEtaHit = min_element(
    hits.begin(),
    hits.end(),
    [](const auto minEtaHit1, const auto minEtaHit2) {
      return ( edm4hep::utils::eta(minEtaHit1.getPosition()) < edm4hep::utils::eta(minEtaHit2.getPosition()) );
    }
  );

  // get hit with max eta
  auto maxEtaHit = max_element(
    hits.begin(),
    hits.end(),
    [](const auto maxEtaHit1, const auto maxEtaHit2) {
      return ( edm4hep::utils::eta(maxEtaHit1.getPosition()) < edm4hep::utils::eta(maxEtaHit2.getPosition()) );
    }
  );
  return make_pair(
    edm4hep::utils::eta(minEtaHit -> getPosition()),
    edm4hep::utils::eta(maxEtaHit -> getPosition())
  );

}  // 'GetEtaRange(edm4eic::Cluster& cluster)'



FPair CheckClusterSplittingProcessor::GetPhiRange(const edm4eic::Cluster& cluster) {

  // get hits
  auto hits = cluster.getHits();

  // get hit with min phi
  auto minPhiHit = min_element(
    hits.begin(),
    hits.end(),
    [](const auto minPhiHit1, const auto minPhiHit2) {
      return ( atan2(minPhiHit1.getPosition().y, minPhiHit1.getPosition().x) < atan2(minPhiHit2.getPosition().y, minPhiHit2.getPosition().x) );
    }
  );

  // get hit with max phi
  auto maxPhiHit = max_element(
    hits.begin(),
    hits.end(),
    [](const auto maxPhiHit1, const auto maxPhiHit2) {
      return ( atan2(maxPhiHit1.getPosition().y, maxPhiHit1.getPosition().x) < atan2(maxPhiHit2.getPosition().y, maxPhiHit2.getPosition().x) );
    }
  );
  return make_pair(
    atan2(minPhiHit -> getPosition().y, minPhiHit -> getPosition().x),
    atan2(maxPhiHit -> getPosition().y, maxPhiHit -> getPosition().x)
  );

}  // 'GetPhiRange(edm4eic::Cluster& cluster)'

// end ------------------------------------------------------------------------
