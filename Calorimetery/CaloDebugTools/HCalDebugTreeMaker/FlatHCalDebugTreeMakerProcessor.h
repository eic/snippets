// ----------------------------------------------------------------------------
// 'FlatDebugHCalTreeMaker.h'
// Derek Anderson
// 03.09.2024
//
// A JANA plugin to construct a flat tree of HCal
// cells and clusters, generated particles, and
// all MC particles.
// ----------------------------------------------------------------------------

#ifndef FLATHCALDEBUGTREEMAKERPROCESSOR_H
#define FLATHCALDEBUGTREEMAKERPROCESSOR_H

// c++ utilities
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <cstdlib>
// root libraries
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
// JANA libraries
#include <JANA/JEventProcessor.h>
#include <JANA/JEventProcessorSequentialRoot.h>
// data model types
#include <edm4eic/Cluster.h>
#include <edm4eic/TrackerHit.h>
#include <edm4eic/ProtoCluster.h>
#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/MCRecoClusterParticleAssociation.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHit.h>
// eicrecon services
#include <spdlog/spdlog.h>
#include <services/log/Log_service.h>
#include <services/rootfile/RootFile_service.h>
#include <services/geometry/dd4hep/DD4hep_service.h>
// dd4hep utilities
#include "DD4hep/Objects.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/IDDescriptor.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/CellIDPositionConverter.h"



class FlatHCalDebugTreeMakerProcessor: public JEventProcessorSequentialRoot {

  // struct to hold user options
  struct Config {
    std::string sPlugin;
    std::string sMCPars;
    std::string sGenPars;
    std::string sSimHits;
    std::string sRecHits;
    std::string sClusters;
    std::vector<std::string> sFields;
  } m_config = {
    "FlatHCalDebugTreeMaker",
    "MCParticles",
    "GeneratedParticles",
    "HcalBarrelHits",
    "HcalBarrelRecHits",
    "HcalBarrelClusters",
    {"eta", "phi"}
  };  // end Config

  // struct to hold class-wise constants
  struct Const {
    size_t nFieldMax;
  } m_const = {6};

  public:

    // ctor
    FlatHCalDebugTreeMakerProcessor() { SetTypeName(NAME_OF_THIS); }

    // public methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

  private:

    // internal methods
    void   InitializeDecoder();
    void   InitializeTree();
    void   ResetVariables();
    void   FillMCVariables(const edm4hep::MCParticle* mc);
    void   FillGenVariables(const edm4eic::ReconstructedParticle* gen);
    void   FillCellVariables(const edm4eic::CalorimeterHit* cell);
    void   FillClusterVariables(const edm4eic::Cluster* clust, const int64_t iClust);
    void   GetCellIndices(const int64_t cellID, std::vector<short>& indices);
    double GetTime(const double tIn);
    double GetEta(const double theta);

    // data objects we need from JANA
    PrefetchT<edm4hep::MCParticle>            m_mcPars   = {this, m_config.sMCPars.data()};
    PrefetchT<edm4eic::ReconstructedParticle> m_genPars  = {this, m_config.sGenPars.data()};
    PrefetchT<edm4hep::SimCalorimeterHit>     m_simHits  = {this, m_config.sSimHits.data()};
    PrefetchT<edm4eic::CalorimeterHit>        m_recHits  = {this, m_config.sRecHits.data()};
    PrefetchT<edm4eic::Cluster>               m_clusters = {this, m_config.sClusters.data()};

    // framework members
    dd4hep::BitFieldCoder* m_decoder;

    // i/o members
    TTree*      m_outTree;
    TDirectory* m_pluginDir;

    // event members
    uint64_t m_evtIndex;
    uint64_t m_nGen;
    uint64_t m_nMC;
    uint64_t m_nCells;
    uint64_t m_nClust;

    // generated particle members
    std::vector<int32_t> m_genType;
    std::vector<int32_t> m_genPDG;
    std::vector<float>   m_genEne;
    std::vector<float>   m_genPhi;
    std::vector<float>   m_genEta;
    std::vector<float>   m_genMass;

    // mc particle members
    std::vector<int32_t> m_mcGenStat;
    std::vector<int32_t> m_mcSimStat;
    std::vector<int32_t> m_mcPDG;
    std::vector<float>   m_mcEne;
    std::vector<float>   m_mcPhi;
    std::vector<float>   m_mcEta;
    std::vector<float>   m_mcMass;
    std::vector<float>   m_mcStartVX;
    std::vector<float>   m_mcStartVY;
    std::vector<float>   m_mcStartVZ;
    std::vector<float>   m_mcStopVX;
    std::vector<float>   m_mcStopVY;
    std::vector<float>   m_mcStopVZ;
    std::vector<float>   m_mcTime;

    // cell members
    std::vector<uint64_t> m_cellID;
    std::vector<float>    m_cellEne;
    std::vector<float>    m_cellRX;
    std::vector<float>    m_cellRY;
    std::vector<float>    m_cellRZ;
    std::vector<float>    m_cellTime;
    std::vector<short>    m_cellIndexA;
    std::vector<short>    m_cellIndexB;
    std::vector<short>    m_cellIndexC;
    std::vector<short>    m_cellIndexE;
    std::vector<short>    m_cellIndexF;
    std::vector<short>    m_cellIndexG;

    // cluster members
    std::vector<uint64_t> m_clustIndex;
    std::vector<uint64_t> m_clustNCells;
    std::vector<float>    m_clustEne;
    std::vector<float>    m_clustEta;
    std::vector<float>    m_clustPhi;
    std::vector<float>    m_clustRX;
    std::vector<float>    m_clustRY;
    std::vector<float>    m_clustRZ;
    std::vector<float>    m_clustTime;

};  // end FlatHCalDebugTreeMaker

#endif

// end ------------------------------------------------------------------------
