# Particle Event Simulator and Visualizer

This project is a simulation and visualization tool for generating synthetic particle events in HepMC ASCII format, supporting pion, kaon, and proton species. It uses [ROOT](https://root.cern/) for random generation, histogramming, and fitting.

## Features

- Generate synthetic particle events with randomized momenta and vertices
- Save events in HepMC v3 ASCII format
- Parse and visualize vertex and momentum distributions using ROOT
- Fit distributions with Gaussian and 2D Gaussian functions
- Interactive plotting interface

## Requirements

- Python 3.x
- [ROOT](https://root.cern/install/) (Python bindings must be available, e.g., `import ROOT` should work)

## Getting Started

## Running the Program
```python3 particle_simulator.py```

Options:
1. Pion
2. Kaon
3. Proton

Plot Options:
1. Unfitted Distributions
2. Fitted Distributions
0. Exit


## Output Format
HEMP3 format

```
E 0 0 0 0 0 0
F GenEvent 3 0
F Units 0 0
F MomentumUnit GEV
F LengthUnit MM
F GenParticle 1 1 211 1 0 0 0 0 px py pz E
F GenVertex 1 -1 -1 -1 vx vy vz 0.0
```

## Visualizations
- 2D Vertex Distribution (X vs Y in mm)
- Momentum Magnitude Histogram (|p| in GeV)
- Optional Gaussian/2D Gaussian fits with reduced χ² output


# Using the generated HEPMC files in the npsim:
- Loading the detector geomettry (we are using the default ePIC detector geometry file as an example)

```source /opt/detector/epic-main/bin/thisepic.sh```

- Running the generated detector geometry file "flat_particle_ascii.hepmc" through DD4Hep.

```ddsim --compactFile $DETECTOR_PATH/$DETECTOR_CONFIG.xml --numberOfEvents 10 --inputFiles flat_particle_ascii.hepmc --outputFile envolope_sim_output.root```
