# Calculate JES/R With PODIO's Python Interface

A short python script to illustrating how to analyze jets
with PODIO's python interface.

## Inputs/Outputs

Takes EICrecon output as its input and processes the jet
collections `GeneratedChargedJets` and `ReconstructedChargedJets`.
Produces a ROOT file with the histogram used to calculate
the JES/R.

## Dependencies

Needs `PODIO` and `edm4eic` to be run. If you're in the
eic-shell, then it should run out of the box.

## Usage

To run, do:

```
python CalculateJESRWithPODIO.py
```

Input and output files can be specified with the command
line options `-i` and `-o` respectively:

```
python CalculateJESRWithPODIO.py \
    -i <my input file> \
    -o <my ouput file>
```
