# Jet Branch Reader Example Macros

A collection of macros which demonstrate how to read the jet branches in the eicrecon output. 

## jetReader_TTreeReader.C

This macro demonstrates how to read the eicrecon output using the TTreeReaderArray functionality in Root. The macro:

- Shows how to access jet quantities from the ReconstructedChargedJets and GeneratedChargedJets branches
- Shows how to access the constituents of jets from the two branches 
- Validates that the energy sum of the jet constituents equal the energy of the parent jet
- Shows how to find the truth particle associated with a particular jet constituent
  - Use to identify jets containing electrons
  - Study PID performance
- Writes several basic histograms to a file

### Prerequisites

- A version of Root that supports TTreeReader (can be local or inside eic-shell)

### Usage

root -l -q jetReader_TTreeReader.C'("/path/to/eicrecon/output/file")'

Note: The specified path can be the xrootd streaming server. For example:

root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.05.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.1950.eicrecon.tree.edm4eic.root

Wildcards are also supported.