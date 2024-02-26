A simple JANA plugin to prepare trees which consolidate potentially useful information for
debugging hadronic calorimeter software. Downstream of [this repo](https://github.com/ruse-traveler/HCalDebugTreeMaker/tree/main).

Takes an edm4hep file as input, and produces two possible outputs: a "flat" tree and a "relational" tree.

 - **Flat Tree:** a tree which contains select information from all hadronic calorimeter clusters, cells, generated particles, and MC particles (i.e. the truth record).
 - **Relational Tree:** a tree which contains select information from all hadronic clusters and (1) the cells comprising those clusters, (2) the MC particles *associated* to those clusters, and (3) the MC particles *contributing* to those clusters. Belonging relationships between clusters and cells/particles are encoded in indices running parallel to the relevant branches.

(Derived from code by Rederike Bock. Thanks!!)

### JBarrelHCalTreeMaker Usage

To build:

```
# after compiling EICrecon, create plugin by:
eicmkplugin.py JBarrelHCalTreeMaker
cp <path to this repo>/JBarrelHCalTreeMakerProcessor.* $EICrecon_ROOT/JBarrelHCalTreeMaker/

# then compile plugin by:
./eic-build JBarrelHCalTreeMaker -B JBarrelHCalTreeMaker/build

# run EICrecon with macro by:
eicrecon -Pplugins=JBarrelHCalTreeMaker <input edm4hep file>
```

Note that the `eic-build` script is included in the `scripts` directory of this repo for convenience. Options are set by changing the relevant parameters in `JBarrelHCalTreeMaker.h` and then recompiling.
