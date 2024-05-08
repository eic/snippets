// ----------------------------------------------------------------------------
// 'RelationalDebugHCalTreeMaker.h'
// Derek Anderson
// 03.09.2024
//
// A JANA plugin to construct a tree of HCal clusters
// and related cells and particles.
// ----------------------------------------------------------------------------

#ifndef RELATIONALHCALDEBUGTREEMAKERPROCESSOR_H
#define RELATIONALHCALDEBUGTREEMAKERPROCESSOR_H

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



class RelationalHCalDebugTreeMakerProcessor: public JEventProcessorSequentialRoot {

  // struct to hold user options
  struct Config {
    std::string sPlugin;
    std::string sMCPars;
    std::string sGenPars;
    std::string sSimHits;
    std::string sRecHits;
    std::string sClusters;
    std::string sAssocs;
    std::vector<std::string> sFields;
  } m_config = {
    "RelationalHCalDebugTreeMaker",
    "MCParticles",
    "GeneratedParticles",
    "HcalEndcapNHits",
    "HcalEndcapNRecHits",
    "HcalEndcapNClusters",
    "HcalEndcapNClusterAssociations",
    {"layer", "slice"}
  };  // end Config

  // struct to hold class-wise constants
  struct Const {
    size_t nFieldMax;
  } m_const = {6};

  public:

    // ctor
    RelationalHCalDebugTreeMakerProcessor() { SetTypeName(NAME_OF_THIS); }

    // public methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

  private:

    // internal methods
    void   InitializeDecoder();
    void   InitializeTree();
    void   ResetVariables();
    void   FillClusterVariables(const edm4eic::Cluster* clust, const int64_t iClust, const int64_t nCells, const int64_t nAssoc, const int64_t nContrib);
    void   FillCellVariables(const edm4eic::CalorimeterHit cell, const int64_t iClust);
    void   FillAssocVariables(const edm4hep::MCParticle mc, const int64_t iClust);
    void   FillContribVariables(const edm4hep::MCParticle mc, const int64_t iClust, const int64_t cellID);
    void   GetCellIndices(const int64_t cellID, std::vector<short>& indices);
    double GetTime(const double tIn);
    double GetEta(const double theta);

    // data objects we need from JANA
    PrefetchT<edm4hep::MCParticle>                       m_mcPars   = {this, m_config.sMCPars.data()};
    PrefetchT<edm4eic::ReconstructedParticle>            m_genPars  = {this, m_config.sGenPars.data()};
    PrefetchT<edm4hep::SimCalorimeterHit>                m_simHits  = {this, m_config.sSimHits.data()};
    PrefetchT<edm4eic::CalorimeterHit>                   m_recHits  = {this, m_config.sRecHits.data()};
    PrefetchT<edm4eic::Cluster>                          m_clusters = {this, m_config.sClusters.data()};
    PrefetchT<edm4eic::MCRecoClusterParticleAssociation> m_assocs   = {this, m_config.sAssocs.data()};

    // framework members
    dd4hep::BitFieldCoder* m_decoder;

    // i/o members
    TTree*      m_outTree;
    TDirectory* m_pluginDir;

    // event members
    uint64_t m_evtIndex;
    uint64_t m_nClust;

    // cluster members
    std::vector<uint64_t> m_clustIndex;
    std::vector<uint64_t> m_clustNCells;
    std::vector<uint64_t> m_clustNAssoc;
    std::vector<uint64_t> m_clustNContrib;
    std::vector<float>    m_clustEne;
    std::vector<float>    m_clustEta;
    std::vector<float>    m_clustPhi;
    std::vector<float>    m_clustRX;
    std::vector<float>    m_clustRY;
    std::vector<float>    m_clustRZ;
    std::vector<float>    m_clustTime;

    // cells in cluster members
    std::vector<uint64_t> m_cellID;
    std::vector<uint64_t> m_cellClustIndex;
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

    // mc particle associated to cluster members
    std::vector<uint64_t> m_assocClustIndex;
    std::vector<int32_t>  m_assocGenStat;
    std::vector<int32_t>  m_assocSimStat;
    std::vector<int32_t>  m_assocPDG;
    std::vector<float>    m_assocEne;
    std::vector<float>    m_assocPhi;
    std::vector<float>    m_assocEta;
    std::vector<float>    m_assocMass;
    std::vector<float>    m_assocStartVX;
    std::vector<float>    m_assocStartVY;
    std::vector<float>    m_assocStartVZ;
    std::vector<float>    m_assocStopVX;
    std::vector<float>    m_assocStopVY;
    std::vector<float>    m_assocStopVZ;
    std::vector<float>    m_assocTime;

    // mc particle contributing to cluster members
    std::vector<uint64_t> m_contribClustIndex;
    std::vector<uint64_t> m_contribCellID;
    std::vector<int32_t>  m_contribGenStat;
    std::vector<int32_t>  m_contribSimStat;
    std::vector<int32_t>  m_contribPDG;
    std::vector<float>    m_contribEne;
    std::vector<float>    m_contribPhi;
    std::vector<float>    m_contribEta;
    std::vector<float>    m_contribMass;
    std::vector<float>    m_contribStartVX;
    std::vector<float>    m_contribStartVY;
    std::vector<float>    m_contribStartVZ;
    std::vector<float>    m_contribStopVX;
    std::vector<float>    m_contribStopVY;
    std::vector<float>    m_contribStopVZ;
    std::vector<float>    m_contribTime;

};  // end RelationalHCalDebugTreeMaker

#endif

// end ------------------------------------------------------------------------
