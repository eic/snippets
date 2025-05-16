# Hit-based track to MC Particle matching
This analysis code extracts the track purity and resolutions using hit-based track to MC Particle matching for DIS events.

Running the analysis code
--------------------------
For the 18x275 GeV beam-energy setting, with Q2>1 GeV2, if a set of files with names ```Q2_1/eicrecon_*.root``` exist, run the analysis code as

```
mkdir plots
root -l -b -q hit_based_matching.C
```

