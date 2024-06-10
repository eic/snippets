# Read Associations of and Contributions to a Reconstructed Cluster

A short(ish) ROOT macro illustrating three points:

  - How to use the MC Particle-Cluster associations,
  - How to access the contributions (and the particles which
    made the contribution) to a given cluster, and
  - How to use the ROOT-based PODIO Frame Reader.

Reads in EICrecon output (`*.podio.root`, `*.tree.edm4eic.root`). User options can be
set in the `Options` struct defined before the macro body.

**NOTE:** use caution! You need to have the CaloHitContributions saved to the
EICrecon output, i.e. there should be branches labeled `*HitsContributions`
in `events`.


## Dependencies

Needs `PODIO`, `edm4hep`, and `edm4eic` to be compiled. If you're running
in the eic-shell, then the commands below will run it out of the box.


## Usage

After setting desired parameters at the top of the macro:

```
root -b -q MatchProjectionsAndClusters.cxx++
```
