/// ===========================================================================
/*! \file    ClusterMergedHitsProcessor.c
 *  \authors Derek Anderson
 *  \date    05.28.2025
 *
 *  A small EICrecon plugin to cluster
 *  merged BHCal hits.
 */
/// ===========================================================================

// plugin definition
#include "ClusterMergedHitsProcessor.h"

// general dependencies
#include <edm4eic/EDM4eicVersion.h>
#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <memory>

// eicrecon-specific components
#include <extensions/jana/JOmniFactoryGeneratorT.h>
#include <factories/calorimetry/CalorimeterClusterRecoCoG_factory.h>
#include <factories/calorimetry/CalorimeterClusterShape_factory.h>
#include <factories/calorimetry/CalorimeterIslandCluster_factory.h>
#include <services/rootfile/RootFile_service.h>



// the following makes this a JANA plugin and adds relevant factories
extern "C" {

  // --------------------------------------------------------------------------
  //! Initialize plugin and add to app
  // --------------------------------------------------------------------------
  /*! This initializes the plugins and wires in
   *  the following factories:
   *    - HcalBarrelMergedHitIslandProtoClusters
   *    - HcalBarrelMergedHitClustersWithoutShapes
   *    - HcalBarrelMergedHitClusters
   */ 
  void InitPlugin(JApplication* app) {

    // supress eicrecon namespace
    using namespace eicrecon;

    // initialize, add to app
    InitJANAPlugin(app);
    app -> Add(new ClusterMergedHitsProcessor);

    // run island clustering on merged hits
    app->Add(
      new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
         "HcalBarrelMergedHitIslandProtoClusters",
        {"HcalBarrelMergedHits"},
        {"HcalBarrelMergedHitIslandProtoClusters"},
        {
          .adjacencyMatrix =
            "("
            "  ( (abs(eta_1 - eta_2) == 1) && (abs(phi_1 - phi_2) == 0) ) ||"
            "  ( (abs(eta_1 - eta_2) == 0) && (abs(phi_1 - phi_2) == 1) ) ||"
            "  ( (abs(eta_1 - eta_2) == 0) && (abs(phi_1 - phi_2) == (320 - 1)) )"
            ") == 1",
          .readout = "HcalBarrelHits",
          .sectorDist = 5.0 * dd4hep::cm,
          .splitCluster = false,
          .minClusterHitEdep = 5.0 * dd4hep::MeV,
          .minClusterCenterEdep = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale = 1.,
        },
        app   // TODO: Remove me once fixed
      )
    );

    // run CoG reconstruction on merged hit clusters
    app->Add(
      new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
         "HcalBarrelMergedHitClustersWithoutShapes",
        {"HcalBarrelMergedHitIslandProtoClusters",
         "HcalBarrelRawHitAssociations"},
        {"HcalBarrelMergedHitClustersWithoutShapes",
         "HcalBarrelMergedHitClusterAssociationsWithoutShapes"},
        {
          .energyWeight = "log",
          .sampFrac = 1.0,
          .logWeightBase = 6.2,
          .enableEtaBounds = false
        },
        app   // TODO: Remove me once fixed
      )
    );

    // run cluster shape calculation on merged hit cluters
    app->Add(
      new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
         "HcalBarrelMergedHitClusters",
        {"HcalBarrelMergedHitClustersWithoutShapes",
         "HcalBarrelMergedHitClusterAssociationsWithoutShapes"},
        {"HcalBarrelMergedHitClusters",
         "HcalBarrelMergedHitClusterAssociations"},
        {
          .energyWeight = "log",
          .logWeightBase = 6.2
        },
        app
      )
    );

  }  // end 'InitPlugin(JApplication*)'

}  // end extern "C"

/// end =======================================================================
