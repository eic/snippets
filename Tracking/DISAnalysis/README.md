# DIS Analysis
The codes in this directory perform analysis on DIS events. They can be run on either locally generated events or the output of the monthly simulation campaigns.

Tracking efficiency in DIS events (DIS_reconstruction.C)
-----------------------------
This analysis code compares the spectra of generated charged particles and reconstructed tracks as a function of pseudorapidity and transverse momentum. It also extracts the efficiency as a function of trasverse momentum.

Hit-based track to MC Particle matching (hit_based_matching.C)
-----------------------------
This analysis code extracts the track purity and resolutions using hit-based track to MC Particle matching for DIS events.

Creating a list of campaign files
-----------------------------
To run the analyses on the output of the simulation campaigns, create a list of campaign files by doing something like the following:

```
mkdir campaign_input
xrdfs dtn-eic.jlab.org ls /volatile/eic/EPIC/RECO/25.05.0/epic_craterlake/DIS/NC/18x275/minQ2=1/ | head -200 | sed 's|^|root://dtn-eic.jlab.org/|g' | tee campaign_input/list_18x275_Q2_1.txt
```

Running the analysis codes
---------------------------
To run the analyses using the simulation campaign files for the 18x275 GeV beam-energy setting with Q2>1 GeV2, do the following

```
mkdir plots
root -l -b -q 'DIS_reconstruction.C(0,1,1)'
root -l -b -q 'hit_based_matching.C(0,1,1)'
```

