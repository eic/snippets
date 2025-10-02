# Using PODIO and Running EICrecon Algorithms in Analysis

A short(ish) ROOT macro illustrating two points:

  - How to use PODIO collections in a standalone ROOT macro,
    especially illustrating how to create new collections and
    save them to an output tree using the `podio::ROOTWriter`
  - How to modify and use EICrecon algorithms in a 
    standalone ROOT Macro.

Reads in EICrecon output (`*.podio.root`, `*.tree.edm4eic.root`,
 etc.).


## To-Do:

There are still a couple to-do items in the example:
  - [ ] Figure out the correct ctor for the jet reconstruction,
        and run it in the event loop.
  - [ ] After fixes in EICrecon, the modified DIS electron
        selection and kinemaitc calculation algorithms
        need to be added/run
  - [ ] And finally, fill the histograms for modified
        algorithms.


## Dependencies

Needs `PODIO`, `edm4hep`, `edm4eic`, and `EICrecon` to be compiled.
If you're running in the eic-shell, then the commands below will
run it out of the box.


## Usage

User options can be set in the `Options` struct defined before the
macro body. And after setting desired parameters at the top of the
macro, it can be run with:

```
source setup.sh
root -b -q RunAlgorithmsOutsideOfReco.cxx
```

**Note:** `setup.sh` only needs to be sourced once!
