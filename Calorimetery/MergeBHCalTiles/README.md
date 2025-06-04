# Merge Tiles in Barrel HCal

An EICrecon plugin and ruby script to merge adjacent tiles in the BHCal and process
the merged tiles up to the `CalorimeterClusterShape` algorithm. These illustrate
two points:

  - The plugin illustrates how to create a custom plugin which wires together
    a new sequence of algorithms that *aren't* run in the default plugins.
  - And the script illustrates how to call this plugin and update the tile-
    merging parameters and adjacency matrix for the BHCal self-consistently.

Reads in either DD4hep output or EICrecon output which contain the `HcalBarrelRecHits`
collection. User options are either set in `RunEICReconWithTileMerging.rb` or from
the commandline (see below).

## Dependencies

Requires all of the usual dependencies for EICrecon (PODIO, EDM4hep, EDM4eic,
DD4hep, etc.)

## Usage

First, the plugin will need to be compiled:

```
/* TODO fill in rest */
```
