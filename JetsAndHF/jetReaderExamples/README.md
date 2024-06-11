# Jet Branch Reader Example Macros

A collection of macros which demonstrate how to read the jet branches in the eicrecon output. 

## jetReader_TTreeReader.C

This macro demonstrates how to read the eicrecon output using the TTreeReaderArray functionality in Root. The macro:

- Shows how to access jet quantities from the ReconstructedChargedJets and GeneratedChargedJets branches
- Shows how to access the constituents of jets from the two branches 
- Validates that the energy sum of the jet constituents equal the energy of the parent jet
- Writes several basic histograms to a file

### Prerequisites

- A version of Root that supports TTreeReader (can be local or inside eic-shell)
- A locally accessible copy of an eicrecon output tree(s)

### Usage

root -l -q jetREader_TTreeReader.C'("/path/to/eicrecon/output/file")'