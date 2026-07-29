# Extracting Detector Type From Calorimeter (and Tracker) Hits

A short ROOT macro illustrating how to extract the [detector type](https://github.com/AIDASoft/DD4hep/blob/912413d08ca432cc3d84632aa936d83960f6999a/DDCore/include/DD4hep/DetType.h#L40)
from a calorimter (or tracker) hit. This is particularly relevant
for identifying if a calorimter hit/cluster from an _electormagnetic_
or _hadronic_ calorimeter in a generic way.


## Inputs/Outputs

Processes `.edm4hep.root` output from `npsim`. Generates a few
histograms comparing the `(z, r)` and `(eta, r)` distributions
for hits for different types of detectors (trackers, cherenkov,
calorimeters, etc.).


## Dependencies

Needs `PODIO` and `edm4hep`. If you're running in the eic-shell,
then it can be run out-of-the-box with the command below.


## Usage

User options can be set in the `Options` struct defined before
the macro body. And after setting the desired parameters at
the top of the macro, it can be run with:

```bash
root -b -q PlotDetectorTypes.cxx
```
