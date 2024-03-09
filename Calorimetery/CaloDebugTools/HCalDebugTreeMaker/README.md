Two simple JANA plugins to prepare trees which consolidate potentially useful information for
debugging hadronic calorimeter software.

Both takes an edm4hep file as input, and both produce a single TTree as output.

 - `FlatHCalDebugTreeMaker:` produces a tree which contains select information from all hadronic
    calorimeter clusters, cells, generated particles, and MC particles (i.e. the truth record).
 - `RelationalHCalDebugTreeMaker:` a tree which contains select information from all hadronic
    clusters and the following:
     (1) the cells comprising those clusters,
     (2) the MC particles *associated* to those clusters, and
     (3) the MC particles *contributing* to those clusters.

In the case of the latter, belonging relationships between clusters and cells/particles are encoded
in indices running parallel to the relevant branches.

(Initially derived from code by Rederike Bock. Thanks!!)

### {Flat,Relational}HCalDebugTreeMaker Usage

To build either:

```
# after compiling EICrecon, create plugin by:
eicmkplugin.py <tree>HCalDebugTreeMaker
cp <path-to-this-repo>/<tree>HCalDebugTreeMaker.* $EICrecon_ROOT/<tree>HCalDebugTreeMaker/

# then compile plugin by:
./eic-build -B <tree>HCalDebugTreeMaker

# run EICrecon with macro by:
eicrecon -Pplugins=<tree>HCalDebugTreeMaker <input-edm4hep-file>
```

Note that the `eic-build` script is included in the `scripts` directory of this repo for convenience.i
Options are set in the `Config` struct at the top of either plugin header; make changes and then
recompile.
