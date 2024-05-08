// ----------------------------------------------------------------------------
// 'FlatDebugHCalTreeMaker.cc'
// Derek Anderson
// 03.09.2024
//
// A JANA plugin to construct a flat tree of HCal
// cells and clusters, generated particles, and
// all MC particles.
// ----------------------------------------------------------------------------

#define FLATHCALDEBUGTREEMAKERPROCESSOR_CC

// plugin definition
#include "FlatHCalDebugTreeMakerProcessor.h"

// The following just makes this a JANA plugin
extern "C" {
    void InitPlugin(JApplication *app) {
        InitJANAPlugin(app);
        app->Add(new FlatHCalDebugTreeMakerProcessor);
    }
}



// public methods -------------------------------------------------------------

void FlatHCalDebugTreeMakerProcessor::InitWithGlobalRootLock(){

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



void FlatHCalDebugTreeMakerProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // reset tree variables
  ResetVariables();

  // loop over generated particles
  size_t nGenPar = 0;
  for (auto gen : m_genPars()) {

    // only accept truth particles
    const int  typeGen = gen -> getType();
    const bool isTruth = (typeGen == 1);
    if (!isTruth) continue;

    FillGenVariables( gen );
    ++nGenPar;

  }  // end generated particle loop

  // loop over mc particles
  size_t nMCPar = 0;
  for (auto mc : m_mcPars()) {

    FillMCVariables( mc );
    ++nMCPar;

  }  // end mc particle loop

  // loop over bhcal hits
  size_t nCell = 0;
  for (auto cell : m_recHits()) {

    FillCellVariables( cell );
    ++nCell;

  }  // end bhcal hit loop

  // loop over clusters
  int64_t nClust = 0;
  for (auto cluster : m_clusters()) {

    FillClusterVariables( cluster, nClust );
    ++nClust;

  }  // end bhcal cluster loop

  // set event variables
  m_nGen   = nGenPar;
  m_nMC    = nMCPar;
  m_nCells = nCell;
  m_nClust = nClust;

  // fill output tree
  m_outTree -> Fill();

  // increment event index and exit routine
  ++m_evtIndex;
  return;

}  // end 'ProcessSequential(std::shared_ptr<JEvent>&)'



void FlatHCalDebugTreeMakerProcessor::FinishWithGlobalRootLock() {

  // clean up variables
  ResetVariables();
  return;

}  // end 'FinishWithGlobalRootLock()'



// internal methods -----------------------------------------------------------

void FlatHCalDebugTreeMakerProcessor::InitializeDecoder() {

  // grab detector
  auto detector = GetApplication() -> GetService<DD4hep_service>() -> detector();

  // make sure readout is available
  dd4hep::IDDescriptor descriptor;
  try {
    descriptor  = detector -> readout(m_config.sSimHits.data()).idSpec();
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



void FlatHCalDebugTreeMakerProcessor::InitializeTree() {

  // switch to output directory
  m_pluginDir -> cd();

  // instantiate output tree
  m_outTree = new TTree("FlatDebugTree", "Flat Tree with Particles, Cells, and Clusters");
  m_outTree -> SetDirectory(m_pluginDir);

  // set branches
  m_outTree -> Branch("EvtIndex",     &m_evtIndex, "EvtIndex/I");
  m_outTree -> Branch("EvtNGenPar",   &m_nGen,     "EvtNGenPar/I");
  m_outTree -> Branch("EvtNMCPar",    &m_nMC,      "EvtNMCPar/I");
  m_outTree -> Branch("EvtNCell",     &m_nCells,   "EvtNCell/I");
  m_outTree -> Branch("EvtNClust",    &m_nClust,   "EvtNClust/I");
  m_outTree -> Branch("GenType",      &m_genType);
  m_outTree -> Branch("GenPDG",       &m_genPDG);
  m_outTree -> Branch("GenE",         &m_genEne);
  m_outTree -> Branch("GenPhi",       &m_genPhi);
  m_outTree -> Branch("GenEta",       &m_genEta);
  m_outTree -> Branch("GenMass",      &m_genMass);
  m_outTree -> Branch("MCGenStat",    &m_mcGenStat);
  m_outTree -> Branch("MCSimStat",    &m_mcSimStat);
  m_outTree -> Branch("MCPDG",        &m_mcPDG);
  m_outTree -> Branch("MCE",          &m_mcEne);
  m_outTree -> Branch("MCPhi",        &m_mcPhi);
  m_outTree -> Branch("MCEta",        &m_mcEta);
  m_outTree -> Branch("MCMass",       &m_mcMass);
  m_outTree -> Branch("MCStartVx",    &m_mcStartVX);
  m_outTree -> Branch("MCStartVy",    &m_mcStartVY);
  m_outTree -> Branch("MCStartVz",    &m_mcStartVZ);
  m_outTree -> Branch("MCStopVx",     &m_mcStopVX);
  m_outTree -> Branch("MCStopVy",     &m_mcStopVY);
  m_outTree -> Branch("MCStopVz",     &m_mcStopVZ);
  m_outTree -> Branch("MCTime",       &m_mcTime);
  m_outTree -> Branch("CellEne",      &m_cellEne);
  m_outTree -> Branch("CellPosX",     &m_cellRX);
  m_outTree -> Branch("CellPosY",     &m_cellRY);
  m_outTree -> Branch("CellPosZ",     &m_cellRZ);
  m_outTree -> Branch("CellTime",     &m_cellTime);
  m_outTree -> Branch("CellID",       &m_cellID);
  m_outTree -> Branch("CellIndexA",   &m_cellIndexA);
  m_outTree -> Branch("CellIndexB",   &m_cellIndexB);
  m_outTree -> Branch("CellIndexC",   &m_cellIndexC);
  m_outTree -> Branch("CellIndexE",   &m_cellIndexE);
  m_outTree -> Branch("CellIndexF",   &m_cellIndexF);
  m_outTree -> Branch("CellIndexG",   &m_cellIndexG);
  m_outTree -> Branch("ClustIndex",   &m_clustIndex);
  m_outTree -> Branch("ClustNCells",  &m_clustNCells);
  m_outTree -> Branch("ClustEne",     &m_clustEne);
  m_outTree -> Branch("ClustEta",     &m_clustEta);
  m_outTree -> Branch("ClustPhi",     &m_clustPhi);
  m_outTree -> Branch("ClustPosX",    &m_clustRX);
  m_outTree -> Branch("ClustPosY",    &m_clustRY);
  m_outTree -> Branch("ClustPosZ",    &m_clustRZ);
  m_outTree -> Branch("ClustTime",    &m_clustTime);
  return;

} // end 'InitializeTrees()'



void FlatHCalDebugTreeMakerProcessor::ResetVariables() {

  // reset event variables
  m_nGen   = 0;
  m_nMC    = 0;
  m_nCells = 0;
  m_nClust = 0;

  // reset gen particle variables
  m_genType.clear();
  m_genPDG.clear();
  m_genEne.clear();
  m_genPhi.clear();
  m_genEta.clear();
  m_genMass.clear();

  // reset event variables
  m_mcGenStat.clear();
  m_mcSimStat.clear();
  m_mcPDG.clear();
  m_mcEne.clear();
  m_mcPhi.clear();
  m_mcEta.clear();
  m_mcMass.clear();
  m_mcStartVX.clear();
  m_mcStartVY.clear();
  m_mcStartVZ.clear();
  m_mcStopVX.clear();
  m_mcStopVY.clear();
  m_mcStopVZ.clear();
  m_mcTime.clear();

  // reset cell members
  m_cellEne.clear();
  m_cellRX.clear();
  m_cellRY.clear();
  m_cellRZ.clear();
  m_cellTime.clear();
  m_cellIndexA.clear();
  m_cellIndexB.clear();
  m_cellIndexC.clear();
  m_cellIndexE.clear();
  m_cellIndexF.clear();
  m_cellIndexG.clear();

  // reset cluster variables
  m_clustNCells.clear();
  m_clustEne.clear();
  m_clustEta.clear();
  m_clustPhi.clear();
  m_clustRX.clear();
  m_clustRY.clear();
  m_clustRZ.clear();
  m_clustTime.clear();
  return;

}  // end 'ResetVariables()'



void FlatHCalDebugTreeMakerProcessor::FillMCVariables(const edm4hep::MCParticle* mc) {

  // grab particle kinematics
  ROOT::Math::PxPyPzM4D<float> pMC(
    mc -> getMomentum().x,
    mc -> getMomentum().y,
    mc -> getMomentum().z,
    mc -> getMass()
  );

  // set output variables
  m_mcGenStat.push_back( mc -> getGeneratorStatus() );
  m_mcSimStat.push_back( mc -> getSimulatorStatus() );
  m_mcPDG.push_back( mc -> getPDG() );
  m_mcEne.push_back( pMC.E() );
  m_mcPhi.push_back( pMC.Phi() );
  m_mcEta.push_back( pMC.Eta() );
  m_mcMass.push_back( mc -> getMass() );
  m_mcStartVX.push_back( mc -> getVertex().x );
  m_mcStartVY.push_back( mc -> getVertex().y );
  m_mcStartVZ.push_back( mc -> getVertex().z );
  m_mcStopVX.push_back( mc -> getEndpoint().x );
  m_mcStopVY.push_back( mc -> getEndpoint().y );
  m_mcStopVZ.push_back( mc -> getEndpoint().z );
  m_mcTime.push_back( mc -> getTime() );
  return;

}  // end 'FillMCVariables(edm4hep::MCParticle*)'



void FlatHCalDebugTreeMakerProcessor::FillGenVariables(const edm4eic::ReconstructedParticle* gen) {

  // get particle kinematics
  ROOT::Math::PxPyPzE4D<float> pGen(
    gen -> getMomentum().x,
    gen -> getMomentum().y,
    gen -> getMomentum().z,
    gen -> getEnergy()
  );

  // set output variables
  m_genType.push_back( gen -> getType() );
  m_genPDG.push_back( gen -> getPDG() );
  m_genEne.push_back( pGen.E() );
  m_genEta.push_back( pGen.Eta() );
  m_genPhi.push_back( pGen.Phi() );
  m_genMass.push_back( gen -> getMass() );
  return;

}  // end 'FillGenVariables(const edm4eic::ReconstructedParticle*)'



void FlatHCalDebugTreeMakerProcessor::FillCellVariables(const edm4eic::CalorimeterHit* cell) {

  // create vector to hold indices
  std::vector<short> indices(m_const.nFieldMax, -1);

  // get hit indices
  const int64_t cellID = cell -> getCellID();
  GetCellIndices(cellID, indices);

  // calculate time
  const double time = GetTime(cell -> getTime());

  // set output variables
  m_cellEne.push_back( cell -> getEnergy() );
  m_cellRX.push_back( cell -> getPosition().x );
  m_cellRY.push_back( cell -> getPosition().y );
  m_cellRZ.push_back( cell -> getPosition().z );
  m_cellTime.push_back( time );
  m_cellID.push_back( cellID );
  m_cellIndexA.push_back( indices[0] );
  m_cellIndexB.push_back( indices[1] );
  m_cellIndexC.push_back( indices[2] );
  m_cellIndexE.push_back( indices[2] );
  m_cellIndexF.push_back( indices[3] );
  m_cellIndexG.push_back( indices[5] );
  return;

}  // end 'FillCellVariables(edm4eic::CalorimeterHit*)' 



void FlatHCalDebugTreeMakerProcessor::FillClusterVariables(const edm4eic::Cluster* clust, const int64_t iClust) {

  // calculate eta
  const double theta = clust -> getIntrinsicTheta();
  const double eta   = GetEta(theta);

  m_clustIndex.push_back( iClust );
  m_clustNCells.push_back( clust -> getNhits() );
  m_clustEne.push_back( clust -> getEnergy() );
  m_clustEta.push_back( eta );
  m_clustPhi.push_back( clust -> getIntrinsicPhi() );
  m_clustRX.push_back( clust -> getPosition().x );
  m_clustRY.push_back( clust -> getPosition().y );
  m_clustRZ.push_back( clust -> getPosition().z );
  m_clustTime.push_back( clust -> getTime() );
  return;

}  // end 'FillClusterVariables(edm4eic::Cluster*, int64_t)'



void FlatHCalDebugTreeMakerProcessor::GetCellIndices(const int64_t cellID, std::vector<short>& indices) {

  for (size_t iField = 0; iField < m_config.sFields.size(); iField++) {
    if (iField < m_const.nFieldMax) {
      indices.at(iField) = m_decoder -> get(cellID, m_config.sFields[iField].data());
    }
  }
  return;

}  // end 'GetCellIndices(int64_t, vector<short>&)'



double FlatHCalDebugTreeMakerProcessor::GetTime(const double tIn) {

  // set max time
  const double maxTime = std::numeric_limits<double>::max();

  // if above, return max
  double tOut = tIn;
  if (tOut > maxTime) {
    tOut = maxTime;
  }
  return tOut;

}  // end 'GetTime(double)'



double FlatHCalDebugTreeMakerProcessor::GetEta(const double theta) {

  const double expo = std::tan(theta / 2.);
  const double eta  = -1. * std::log(expo);
  return eta;

}  // end 'GetEta(double)'

// end ------------------------------------------------------------------------
