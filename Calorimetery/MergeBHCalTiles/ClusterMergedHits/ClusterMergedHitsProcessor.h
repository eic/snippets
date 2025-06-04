/// ===========================================================================
/*! \file    ClusterMergedHitsProcessor.h
 *  \authors Derek Anderson
 *  \date    05.28.2025
 *
 *  A small EICrecon plugin to cluster
 *  merged BHCal hits.
 */
/// ===========================================================================

#include <JANA/JEventProcessorSequentialRoot.h>
#include <TH2D.h>
#include <TFile.h>



// ============================================================================
//! Cluster Merged BHCal Hits
// ============================================================================
/*! A small EICrecon plugin to cluster
 *  merged BHCal hits.
 */
class ClusterMergedHitsProcessor : public JEventProcessorSequentialRoot {

  public:

    // ctor
    ClusterMergedHitsProcessor() { SetTypeName(NAME_OF_THIS); }

};  // end ClusterMergedHitsProcessor

/// end =======================================================================
