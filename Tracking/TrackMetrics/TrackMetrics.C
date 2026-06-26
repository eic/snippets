// ============================================================
// Build particle -> reco hit associations
//
// Goal:
//   For each MC particle (status==1), save the list of
//   reconstructed tracking hits associated with it.
//
// This will later allow:
//
//   particle -> hits
//   track    -> hits
//
// and therefore:
//   track <-> particle matching via shared hits
// ============================================================

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdint>
#include <iostream>

// ============================================================
// Unique hit key type
// ============================================================
using HitKey = uint64_t;

// ============================================================
// Create globally unique hit key
// ============================================================
inline HitKey makeHitKey(const podio::ObjectID& id) {
    return (uint64_t(id.collectionID) << 32) 
        | uint32_t(id.index);
    }

// ==============================================================
// Create structure for globally unique hit key and quality flag
// ==============================================================
struct TruthHitInfo {
    
    uint64_t recoHitKey;
    int quality;

};

// ==============================================================
// Create structure for Rawhit to Particle with hit quality flag
// ==============================================================
struct RawHitTruthInfo {

    int particleID;
    int quality;

};

// ============================================================
// Function to set plotting style
// ============================================================
void set_style(){
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);
}

// ============================================================
// Main function
// ============================================================
void TrackMetrics(int input_setting = 1){

    // ============================================================
    // Set style
    // ============================================================
    set_style();

    // ============================================================
    // Set maximum number of events allowed to run
    // ============================================================
    int max_events = 1E2;

    // ============================================================
    // Detector config
    // ============================================================
    enum class DetectorType {
        Endcap,
        Barrel
    };

    struct DetectorConfig {
        std::string collection;
        DetectorType type;
    };

    // ============================================================
    // Active detector collections
    // ============================================================
    std::vector<DetectorConfig> detectors = {
        {"SiEndcapTrackerRecHits", DetectorType::Endcap},
        {"SiBarrelVertexRecHits", DetectorType::Barrel},
        {"SiBarrelTrackerRecHits", DetectorType::Barrel},
        {"MPGDBarrelRecHits", DetectorType::Barrel},
        {"OuterMPGDBarrelRecHits", DetectorType::Barrel},
        {"BackwardMPGDEndcapRecHits", DetectorType::Endcap},
        {"ForwardMPGDEndcapRecHits", DetectorType::Endcap},
        {"TOFBarrelRecHits",DetectorType::Barrel},
        {"TOFEndcapRecHits",DetectorType::Endcap}
    };

    // ============================================================
    // Create histograms
    // ============================================================

    // ============================================================
    // Open files and load Reader
    // ============================================================
    std::vector<string> inputfiles;
    podio::ROOTReader r;

    std::string input;
    if (input_setting==0) 
        input = "./input_lists/local_files.txt";
    else if (input_setting==1)
        input = "./input_lists/local_merged_files.txt";
    else if (input_setting==2) 
        input = "./input_lists/list_ForcedBG_10um.txt";
    
    std::ifstream in(input);
    std::string file("");
    while (in >> file){
    	inputfiles.push_back(file);
    }
    r.openFiles(inputfiles);
    in.close();

    // Print run information
    auto nevents = r.getEntries(podio::Category::Event);
    if (nevents > max_events) 
        nevents = max_events;

    cout<<"Running analysis on "<<input<<"!"<<endl;
    cout<<"Analyzing "<< nevents <<" events!"<<endl;

    // ============================================================
    // Event loop
    // ============================================================
    for (unsigned int ievent = 0; ievent < nevents; ievent++) {
        
        if ((ievent % 100) == 0) 
            std::cout << "Processed event "
                      << ievent
                      << std::endl;

        // ========================================================
        // Read event
        // ========================================================
        auto f = 
            podio::Frame( 
                r.readNextEntry(podio::Category::Event) 
            );

        // ========================================================
        // RawHit <-> SimHit associations
        // ========================================================
        auto& raw_hit_assocs =
            f.get<edm4eic::MCRecoTrackerHitAssociationCollection>(
                "CentralTrackingRawHitAssociations"
            );

        // ========================================================
        // Map:
        //
        //   rawHitKey -> RawHitTruthInfo
        //
        // This is built ONCE per event.
        // ========================================================
        std::unordered_map<HitKey, RawHitTruthInfo>
            rawHitTruthMap;

        for (const auto& assoc : raw_hit_assocs) {
            
            auto sim_hit = assoc.getSimHit();
            
            auto mc_particle =
                sim_hit.getParticle();

            int mc_status =
                mc_particle.getGeneratorStatus();

            // ----------------------------------------------------
            // Keep only generator particles
            // ----------------------------------------------------
            if (mc_status != 1)
                continue;

            int particleID =
                mc_particle.getObjectID().index;

            auto raw_hit =
                assoc.getRawHit();

            HitKey rawKey =
                makeHitKey(raw_hit.getObjectID());

            RawHitTruthInfo info;
            
            info.particleID = particleID;
            info.quality   = sim_hit.getQuality();
            
            rawHitTruthMap[rawKey] = info;
        }

        // ========================================================
        // Final structure:
        //
        //   particleID -> set of TruthHitInfo
        //
        // ========================================================
        std::unordered_map<
            int,
            std::vector<TruthHitInfo>
        > particleHits;

        // ========================================================
        // Loop over detector collections
        // ========================================================
        for (const auto& det : detectors) {

            auto& rechit_coll =
                f.get<edm4eic::TrackerHitCollection>(
                    det.collection
                );

            std::cout
                << "Event " << ievent
                << " | Collection: " << det.collection
                << " | rechits: " << rechit_coll.size()
                << std::endl;

            // ====================================================
            // Loop over reconstructed hits
            // ====================================================
            for (const auto& hit : rechit_coll) {

                // ------------------------------------------------
                // Get raw hit associated with reco hit
                // ------------------------------------------------
                auto raw_hit =
                    hit.getRawHit();

                HitKey rawKey =
                    makeHitKey(raw_hit.getObjectID());
                    
                // ------------------------------------------------
                // Find corresponding MC particle
                // ------------------------------------------------
                auto it =
                    rawHitTruthMap.find(rawKey);

                if (it == rawHitTruthMap.end())
                    continue;

                const auto& truth = it->second;
                
                int particleID = truth.particleID;
                int quality = truth.quality;

                // ------------------------------------------------
                // Create reco hit key
                // ------------------------------------------------
                HitKey recoKey =
                    makeHitKey(hit.getObjectID());

                // ------------------------------------------------
                // Save reco hit key and quality
                // ------------------------------------------------
                TruthHitInfo hitInfo;
                
                hitInfo.recoHitKey = recoKey;
                hitInfo.quality    = quality;
                
                particleHits[particleID].push_back(hitInfo);
            }
        }

        // ========================================================
        // Debug printout
        // ========================================================
        std::cout
            << "Event "
            << ievent
            << " has "
            << particleHits.size()
            << " \033[1;34mstatus==1\033[0m"
            << " particles with tracker hits"
            << std::endl;

        // --------------------------------------------------------
        // Print first few particles
        // --------------------------------------------------------
        int nprint = 0;

        for (const auto& [pid, hits] : particleHits) {

            int nQuality0 = 0;

            for (const auto& hit : hits) {

                if (hit.quality == 0)
                nQuality0++;
            }

            std::cout
                << "  Particle "
                << pid
                << " --> "
                << hits.size()
                << " reco hits"
                << " (and "
                << nQuality0
                << " \033[1;32mquality==0\033[0m" 
                << " hits)"
                << std::endl;

            nprint++;

            if (nprint > 10)
                break;
        }

    } // end event loop

    // --------------------------------------------------------
    // Make plots
    // --------------------------------------------------------

} // end main function