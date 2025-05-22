#include "TFile.h"
#include "TH2D.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <cmath>
#include <iostream>

/**
 * ROOT macro for analyzing MC particle distributions
 * Creates pt vs eta histogram from Monte Carlo particle data
 * 
 * @param inFileName  Input ROOT file containing MC data
 * @param outFileName Output ROOT file for histograms
 */
void example_macro(TString inFileName = "input.root", 
                   TString outFileName = "ana.root") {
    
    // Open input file containing the event tree
    TFile* file = new TFile(inFileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inFileName << std::endl;
        return;
    }
    
    // Setup tree reader for MC particle data
    TTreeReader myReader("events", file);
    
    // Define data containers for particle momentum components
    TTreeReaderArray<double> px_mc(myReader, "MCParticles.momentum.x");
    TTreeReaderArray<double> py_mc(myReader, "MCParticles.momentum.y");
    TTreeReaderArray<double> pz_mc(myReader, "MCParticles.momentum.z");
    TTreeReaderArray<int>    status(myReader, "MCParticles.generatorStatus");
    TTreeReaderArray<int>    pdg(myReader, "MCParticles.PDG");
    
    // Create output file and histogram
    TFile* output = new TFile(outFileName, "recreate");
    TH2D* hPtEtaMc = new TH2D("hPtEtaMc", 
                              "MC Particle Distribution;p_{t,gen} (GeV/c);#eta_{gen}", 
                              1500, 0.0, 10.0,    // pt bins: 1500 bins from 0 to 10 GeV/c
                              200, -5.0, 5.0);    // eta bins: 200 bins from -5 to 5
    
    // Event loop
    int nEvents = myReader.GetEntries();
    std::cout << "Total Events: " << nEvents << std::endl;
    
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        myReader.SetEntry(iEvent);
        
        // Progress indicator
        if (iEvent % 1000 == 0) {
            std::cout << "Processing Event: " << iEvent << std::endl;
        }
        
        // Loop over particles in current event
        for (int iParticle = 0; iParticle < status.GetSize(); ++iParticle) {
            
            // Only consider stable particles (status == 1)
            if (status[iParticle] != 1) {
                continue;
            }
            
            // Calculate kinematic variables
            double px = px_mc[iParticle];
            double py = py_mc[iParticle];
            double pz = pz_mc[iParticle];
            
            // Total momentum magnitude
            double p_mc = std::sqrt(px*px + py*py + pz*pz);
            
            // Transverse momentum
            double pt_mc = std::sqrt(px*px + py*py);
            
            // Pseudorapidity (eta)
            double theta_mc = std::acos(pz / p_mc);
            double eta_mc = -std::log(std::tan(theta_mc / 2.0));
            
            // Optional: Azimuthal angle (phi) - currently calculated but not used
            // double phi_mc = std::atan2(py, px);  // Use atan2 for proper quadrant
            
            // Fill histogram
            hPtEtaMc->Fill(pt_mc, eta_mc);
            
        } // End particle loop
        
    } // End event loop
    
    // Write output
    output->cd();
    hPtEtaMc->Write();
    output->Save();
    output->Close();
    
    // Clean up
    file->Close();
    delete file;
    delete output;
    
    std::cout << "Analysis complete. Output saved to: " << outFileName << std::endl;
}
