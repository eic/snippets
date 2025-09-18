# Low-momentum tracking efficiency
The codes in this directory generate single-particle simulations at fixed momenta, with the primary goal being to understand the track-finding efficiency for lower-momentum particles. .

Running events through the ePIC simulation
-------------------------------------------
To create a set of simulations for negative muons at fixed momentum values -- 500 MeV/c, 1 GeV/c, 2 GeV/c, 5 GeV/c, 10 GeV/c, 20 GeV/c -- do the following inside the eic-shell container:

```
mkdir output
./run_single_muons_fixed_mom.sh
```

This will create a set of 6 output simulations, each at a fixed momentum, with 1000 events per simulation. The particles will be generated uniformly in eta and phi.

To run a specfic momentum setting only, you can use the relevant command-line argument. For example, to run only the 500 MeV/c setting, do the following inside the eic-shell container:

```
./run_single_muons_fixed_mom.sh 1
```

To run the same setting as above with an adjusted <i> deltaPhiMax </i> seeding parameter, do the following inside the eic-shell container:

```
./run_single_muons_fixed_mom.sh 1 1 0 1000 0.16
```

Running analysis code to extract efficiency
-------------------------------------------
This analysis code extracts the track efficiency as a function of eta and phi for the above events. Results with truth-seeded tracking are shown in blue; results for real-seeded tracking are shown in green. The results are saved to a single pdf file.

```
mkdir plots
root -l -b -q track_eff_fixed_mom.C
```

Comments:

1. For the phi efficiency calculation, only generated particles (events) with |eta| < 3.0 are considered.
2. No track-selection cuts are applied at the analysis level. The ePIC track reconstructed outputs tracks that use at least 3 good measurement hits in the track fit.


