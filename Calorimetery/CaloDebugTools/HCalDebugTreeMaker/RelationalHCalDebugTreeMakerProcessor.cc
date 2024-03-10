// ----------------------------------------------------------------------------
// 'RelationalDebugHCalTreeMaker.cc'
// Derek Anderson
// 03.09.2024
//
// A JANA plugin to construct a tree of HCal clusters
// and related cells and particles.
// ----------------------------------------------------------------------------

#define RELATIONALHCALDEBUGTREEMAKERPROCESSOR_CC

// plugin definition
#include "RelationalHCalDebugTreeMakerProcessor.h"

// The following just makes this a JANA plugin
extern "C" {
    void InitPlugin(JApplication *app) {
        InitJANAPlugin(app);
        app->Add(new RelationalHCalDebugTreeMakerProcessor);
    }
}



// public methods -------------------------------------------------------------

void RelationalHCalDebugTreeMakerProcessor::InitWithGlobalRootLock(){

  //  grab output file
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc     -> GetHistFile();

  // create directory
  m_pluginDir = rootfile -> mkdir(m_config.sPlugin.data());
  m_pluginDir -> cd();

  // reset tree variables
  ResetVariables();

  // initialize event index
  m_evtIndex = 0;

  // initialize output trees
  InitializeDecoder();
  InitializeTree();
  return;

}  // end 'InitWithGlobalLock()'



void RelationalHCalDebugTreeMakerProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // reset tree variables
  ResetVariables();

  // loop over clusters
  int64_t nClust = 0;
  for (auto cluster : m_clusters()) {

    // loop over associated mc particles
    int64_t nAssocToClust = 0;
    for (auto assoc : m_assocs()) {

      // consider only current cluster
      const bool isSameClust = (assoc -> getRecID() == nClust);
      if (!isSameClust) continue;

      // grab associated particle
      auto mcPar = assoc -> getSim();

      FillAssocVariables( mcPar, nClust );
      ++nAssocToClust;

    }  // end association loop

    // grab consituent hits
    const auto clustHits = cluster -> getHits();

    // associate each cell, contributing particle with corresponding cluster
    int64_t nCellsInClust    = 0;
    int64_t nContribToClust = 0;
    for (auto clustHit : clustHits) {

      // get cell ID
      const uint64_t cellID = clustHit.getCellID();

      // get contributing particles
      for (auto simHit : m_simHits()) {

        // identify sim hit based on cell ID
        const bool isSameCell = (cellID == simHit -> getCellID());
        if (!isSameCell) continue;

        // loop over hit contributions
        const auto contribs = simHit -> getContributions();
        for (auto contrib : contribs) {

          // grab contributing particle
          auto mcPar = contrib.getParticle();

          FillContribVariables( mcPar, nClust, cellID );
          ++nContribToClust;

        }  // end contribution loop
      }  // end sim hit loop

      FillCellVariables( clustHit, nClust );
      ++nCellsInClust;

    }  // end cell loop


    FillClusterVariables( cluster, nClust, nCellsInClust, nAssocToClust, nContribToClust );
    ++nClust;

  }  // end cluster loop

  // set event variables
  m_nClust = nClust;

  // fill output tree
  m_outTree -> Fill();

  // increment event index and exit routine
  ++m_evtIndex;
  return;

}  // end 'ProcessSequential(std::shared_ptr<JEvent>&)'



void RelationalHCalDebugTreeMakerProcessor::FinishWithGlobalRootLock() {

  // clean up variables
  ResetVariables();
  return;

}  // end 'FinishWithGlobalRootLock()'



// internal methods -----------------------------------------------------------

void RelationalHCalDebugTreeMakerProcessor::InitializeDecoder() {

  // grab detector
  auto detector = GetApplication() -> GetService<DD4hep_service>() -> detector();

  // make sure readout is available
  dd4hep::IDDescriptor descriptor;
  try {
    descriptor  = detector -> readout("HcalBarrelHits").idSpec();
  } catch (const std::runtime_error &err) {
    throw std::runtime_error("PANIC: readout class is not in output!");
  }

  // grab decoder and test
  std::vector<short> indices(m_config.sFields.size());
  try {
    m_decoder = descriptor.decoder();
    for (size_t iField = 0; iField < m_config.sFields.size(); iField++) {
      indices[iField] = m_decoder -> index(m_config.sFields[iField]);
    }
  } catch (const std::runtime_error &err) {
    throw std::runtime_error("PANIC: something went wrong grabbing the decoder!");
  }
  return;

}  // end 'InitializeDecoder()'



void RelationalHCalDebugTreeMakerProcessor::InitializeTree() {

  // switch to output directory
  m_pluginDir -> cd();

  // instantiate output tree
  m_outTree = new TTree("RelationalDebugTree", "Tree with Cluster-Cell/Particle Relations ");
  m_outTree -> SetDirectory(m_pluginDir);

  // set branches
  m_outTree -> Branch("EvtIndex",          &m_evtIndex, "EvtIndex/I");
  m_outTree -> Branch("EvtNClust",         &m_nClust,   "EvtNClust/I");
  m_outTree -> Branch("ClustIndex",        &m_clustIndex);
  m_outTree -> Branch("ClustNCell",        &m_clustNCells);
  m_outTree -> Branch("ClustNAssocPar",    &m_clustNAssoc);
  m_outTree -> Branch("ClustNContribPar",  &m_clustNContrib);
  m_outTree -> Branch("ClustEne",          &m_clustEne);
  m_outTree -> Branch("ClustEta",          &m_clustEta);
  m_outTree -> Branch("ClustPhi",          &m_clustPhi);
  m_outTree -> Branch("ClustPosX",         &m_clustRX);
  m_outTree -> Branch("ClustPosY",         &m_clustRY);
  m_outTree -> Branch("ClustPosZ",         &m_clustRZ);
  m_outTree -> Branch("ClustTime",         &m_clustTime);
  m_outTree -> Branch("CelllID",           &m_cellID);
  m_outTree -> Branch("CellClustIndex",    &m_cellClustIndex);
  m_outTree -> Branch("CellEne",           &m_cellEne);
  m_outTree -> Branch("CellPosX",          &m_cellRX);
  m_outTree -> Branch("CellPosY",          &m_cellRY);
  m_outTree -> Branch("CellPosZ",          &m_cellRZ);
  m_outTree -> Branch("CellTime",          &m_cellTime);
  m_outTree -> Branch("CellIndexA",        &m_cellIndexA);
  m_outTree -> Branch("CellIndexB",        &m_cellIndexB);
  m_outTree -> Branch("CellIndexC",        &m_cellIndexC);
  m_outTree -> Branch("AssocClustIndex",   &m_assocClustIndex);
  m_outTree -> Branch("AssocGenStat",      &m_assocGenStat);
  m_outTree -> Branch("AssocSimStat",      &m_assocSimStat);
  m_outTree -> Branch("AssocPDG",          &m_assocPDG);
  m_outTree -> Branch("AssocEne",          &m_assocEne);
  m_outTree -> Branch("AssocPhi",          &m_assocPhi);
  m_outTree -> Branch("AssocEta",          &m_assocEta);
  m_outTree -> Branch("AssocMass",         &m_assocMass);
  m_outTree -> Branch("AssocStartVX",      &m_assocStartVX);
  m_outTree -> Branch("AssocStartVY",      &m_assocStartVY);
  m_outTree -> Branch("AssocStartVZ",      &m_assocStartVZ);
  m_outTree -> Branch("AssocStopVX",       &m_assocStopVX);
  m_outTree -> Branch("AssocStopVY",       &m_assocStopVY);
  m_outTree -> Branch("AssocStopVZ",       &m_assocStopVZ);
  m_outTree -> Branch("AssocTime",         &m_assocTime);
  m_outTree -> Branch("ContribClustIndex", &m_contribClustIndex);
  m_outTree -> Branch("ContribCellID",     &m_contribCellID);
  m_outTree -> Branch("ContribGenStat",    &m_contribGenStat);
  m_outTree -> Branch("ContribSimStat",    &m_contribSimStat);
  m_outTree -> Branch("ContribPDG",        &m_contribPDG);
  m_outTree -> Branch("ContribE",          &m_contribEne);
  m_outTree -> Branch("ContribPhi",        &m_contribPhi);
  m_outTree -> Branch("ContribEta",        &m_contribEta);
  m_outTree -> Branch("ContribMass",       &m_contribMass);
  m_outTree -> Branch("ContribStartVX",    &m_contribStartVX);
  m_outTree -> Branch("ContribStartVY",    &m_contribStartVY);
  m_outTree -> Branch("ContribStartVZ",    &m_contribStartVZ);
  m_outTree -> Branch("ContribStopVX",     &m_contribStopVX);
  m_outTree -> Branch("ContribStopVY",     &m_contribStopVY);
  m_outTree -> Branch("ContribStopVZ",     &m_contribStopVZ);
  m_outTree -> Branch("ContribTime",       &m_contribTime);
  return;

} // end 'InitializeTrees()'



void RelationalHCalDebugTreeMakerProcessor::ResetVariables() {

  // reset event variables
  m_nClust = 0;

  // reset cluster variables
  m_clustNCells.clear();
  m_clustNContrib.clear();
  m_clustNAssoc.clear();
  m_clustEne.clear();
  m_clustEta.clear();
  m_clustPhi.clear();
  m_clustRX.clear();
  m_clustRY.clear();
  m_clustRZ.clear();
  m_clustTime.clear();

  // reset cells in cluster variables
  m_cellID.clear();
  m_cellClustIndex.clear();
  m_cellEne.clear();
  m_cellRX.clear();
  m_cellRY.clear();
  m_cellRZ.clear();
  m_cellTime.clear();
  m_cellIndexA.clear();
  m_cellIndexB.clear();
  m_cellIndexC.clear();

  // reset mc particle associated to cluster variables
  m_assocClustIndex.clear();
  m_assocGenStat.clear();
  m_assocSimStat.clear();
  m_assocPDG.clear();
  m_assocEne.clear();
  m_assocPhi.clear();
  m_assocEta.clear();
  m_assocMass.clear();
  m_assocStartVX.clear();
  m_assocStartVY.clear();
  m_assocStartVZ.clear();
  m_assocStopVX.clear();
  m_assocStopVY.clear();
  m_assocStopVZ.clear();
  m_assocTime.clear();

  // reset mc particle contributing to cluster variables
  m_contribClustIndex.clear();
  m_contribCellID.clear();
  m_contribGenStat.clear();
  m_contribSimStat.clear();
  m_contribPDG.clear();
  m_contribEne.clear();
  m_contribPhi.clear();
  m_contribEta.clear();
  m_contribMass.clear();
  m_contribStartVX.clear();
  m_contribStartVY.clear();
  m_contribStartVZ.clear();
  m_contribStopVX.clear();
  m_contribStopVY.clear();
  m_contribStopVZ.clear();
  m_contribTime.clear();
  return;

}  // end 'ResetVariables()'



void RelationalHCalDebugTreeMakerProcessor::FillClusterVariables(
  const edm4eic::Cluster* clust,
  const int64_t iClust,
  const int64_t nCells,
  const int64_t nAssoc,
  const int64_t nContrib
) {

  // calculate eta
  const double theta = clust -> getIntrinsicTheta();
  const double eta   = GetEta(theta);

  m_clustIndex.push_back( iClust );
  m_clustNCells.push_back( nCells );
  m_clustNAssoc.push_back( nAssoc );
  m_clustNContrib.push_back( nContrib );
  m_clustEne.push_back( clust -> getEnergy() );
  m_clustEta.push_back( eta );
  m_clustPhi.push_back( clust -> getIntrinsicPhi() );
  m_clustRX.push_back( clust -> getPosition().x );
  m_clustRY.push_back( clust -> getPosition().y );
  m_clustRZ.push_back( clust -> getPosition().z );
  m_clustTime.push_back( clust -> getTime() );
  return;

}  // end 'FillClusterVariables(edm4eic::Cluster*, int64_t)'



void RelationalHCalDebugTreeMakerProcessor::FillCellVariables(const edm4eic::CalorimeterHit cell, const int64_t iClust) {

  // create vector to hold indices
  std::vector<short> indices(m_const.nFieldMax, -1);

  // get hit indices
  const int64_t cellID = cell.getCellID();
  GetCellIndices(cellID, indices);

  // calculate time
  const double time = GetTime(cell.getTime());

  // set output variables
  m_cellID.push_back( cellID );
  m_cellClustIndex.push_back( iClust );
  m_cellEne.push_back( cell.getEnergy() );
  m_cellRX.push_back( cell.getPosition().x );
  m_cellRY.push_back( cell.getPosition().y );
  m_cellRZ.push_back( cell.getPosition().z );
  m_cellTime.push_back( time );
  m_cellIndexA.push_back( indices[0] );
  m_cellIndexB.push_back( indices[1] );
  m_cellIndexC.push_back( indices[2] );
  return;

}  // end 'FillCellVariables(edm4eic::CalorimeterHit, int64_t)' 



void RelationalHCalDebugTreeMakerProcessor::FillAssocVariables(const edm4hep::MCParticle assoc, const int64_t iClust) {

  // grab particle kinematics
  ROOT::Math::PxPyPzM4D<float> pMC(
    assoc.getMomentum().x,
    assoc.getMomentum().y,
    assoc.getMomentum().z,
    assoc.getMass()
  );

  // set output variables
  m_assocClustIndex.push_back( iClust );
  m_assocGenStat.push_back( assoc.getGeneratorStatus() );
  m_assocSimStat.push_back( assoc.getSimulatorStatus() );
  m_assocPDG.push_back( assoc.getPDG() );
  m_assocEne.push_back( pMC.E() );
  m_assocPhi.push_back( pMC.Phi() );
  m_assocEta.push_back( pMC.Eta() );
  m_assocMass.push_back( assoc.getMass() );
  m_assocStartVX.push_back( assoc.getVertex().x );
  m_assocStartVY.push_back( assoc.getVertex().y );
  m_assocStartVZ.push_back( assoc.getVertex().z );
  m_assocStopVX.push_back( assoc.getEndpoint().x );
  m_assocStopVY.push_back( assoc.getEndpoint().y );
  m_assocStopVZ.push_back( assoc.getEndpoint().z );
  m_assocTime.push_back( assoc.getTime() );
  return;

}  // end 'FillAssocVariables(edm4hep::MCParticle, int64_t)'



void RelationalHCalDebugTreeMakerProcessor::FillContribVariables(
  const edm4hep::MCParticle mc,
  const int64_t iClust,
  const int64_t cellID
) {

  // grab particle kinematics
  ROOT::Math::PxPyPzM4D<float> pMC(
    mc.getMomentum().x,
    mc.getMomentum().y,
    mc.getMomentum().z,
    mc.getMass()
  );

  // set output variables
  m_contribClustIndex.push_back( iClust );
  m_contribCellID.push_back( cellID );
  m_contribGenStat.push_back( mc.getGeneratorStatus() );
  m_contribSimStat.push_back( mc.getSimulatorStatus() );
  m_contribPDG.push_back( mc.getPDG() );
  m_contribEne.push_back( pMC.E() );
  m_contribPhi.push_back( pMC.Phi() );
  m_contribEta.push_back( pMC.Eta() );
  m_contribMass.push_back( mc.getMass() );
  m_contribStartVX.push_back( mc.getVertex().x );
  m_contribStartVY.push_back( mc.getVertex().y );
  m_contribStartVZ.push_back( mc.getVertex().z );
  m_contribStopVX.push_back( mc.getEndpoint().x );
  m_contribStopVY.push_back( mc.getEndpoint().y );
  m_contribStopVZ.push_back( mc.getEndpoint().z );
  m_contribTime.push_back( mc.getTime() );
  return;

}  // end 'FillContribVariables(edm4hep::MCParticle, int64_t, int64_t)'



void RelationalHCalDebugTreeMakerProcessor::GetCellIndices(const int64_t cellID, std::vector<short>& indices) {

  for (size_t iField = 0; iField < m_config.sFields.size(); iField++) {
    if (iField < m_const.nFieldMax) {
      indices.at(iField) = m_decoder -> get(cellID, m_config.sFields[iField].data());
    }
  }
  return;

}  // end 'GetCellIndices(int64_t, vector<short>)'



double RelationalHCalDebugTreeMakerProcessor::GetTime(const double tIn) {

  // set max time
  const double maxTime = std::numeric_limits<double>::max();

  // if above, return max
  double tOut = tIn;
  if (tOut > maxTime) {
    tOut = maxTime;
  }
  return tOut;

}  // end 'GetTime(double)'



double RelationalHCalDebugTreeMakerProcessor::GetEta(const double theta) {

  const double expo = std::tan(theta / 2.);
  const double eta  = -1. * std::log(expo);
  return eta;

}  // end 'GetEta(double)'

// end ------------------------------------------------------------------------

