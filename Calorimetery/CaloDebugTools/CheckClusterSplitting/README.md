# CheckClusterSplitting

An EICrecon plugin to generate histograms of various metrics to check the health
of cluster splitting/merging.

### Usage

To build:

```
# after compiling EICrecon, create plugin by:
eicmkplugin.py CheckClusterSplitting
cp <path-to-this-repo>/CheckClusterSplittingProcessor.* $EICrecon_ROOT/CheckClusterSplitting/

# then compile plugin by:
cmake -S CheckClusterSplitting -B <tree>CheckClusterSplitting/build
cmake --build CheckClusterSplitting/build --target install -- -j8

# run EICrecon with macro by:
eicrecon -Pplugins=CheckClusterSplitting <input-edm4hep-file>
```

Options are set in the `Config` struct at the top of either plugin header; make changes and then
recompile.

 
