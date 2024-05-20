# Match Clusters to Track Projections

A short(ish) ROOT macro illustrating two points:

  - How to match calorimeter clusters to track projections using
    distance in eta-phi; and
  - How to use the ROOT-based PODIO Frame Reader.

Reads in EICrecon output (`*.podio.root`). User options can be set in the `Options`
struct defined before the macro body.


## Dependencies

Needs `PODIO`, `edm4hep`, and `edm4eic` to be compiled. If you're running
in the eic-shell, then the commands below will run it out of the box.


## Usage

After setting desired parameters at the top of the macro:

```
root -b -q MatchProjectionsAndClusters.cxx++
```
