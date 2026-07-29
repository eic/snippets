# Extracting Indices From CellIDs

A short ROOT macro illustrating how to extract the value of individual 
fields from a `CellID`.  The `CellID` is a bit field that stores all
of the various indices and IDs, referred to as _fields_ in DD4hep,
that define a unique volume.  The no. of bits in the bit field
allocated to each field is defined in the `readouts` element in each
detector's compact (xml) file.

Really, there are 2 parts to unpacking a CellID. First is to
initialize necessary DD4hep interfaces:

```c++
// 1) entry point into DD4hep interface
static std::unique_ptr<dd4hep::Detector> detector = dd4hep::Detector::make_unique("");
detector->fromCompact("/path/to/my/epic.xml");

// 2) interface for working with readout fields/cell IDs
//   --> argument should be name of readouts element
//       in your xml file
dd4hep::IDDescriptor descriptor = detector->readout("HcalBarrelHits");

// 3) utility for unpacking cell IDs
dd4hep::DDSegmentation::BitFieldCoder* decoder = descriptor.decoder();
```

Then field values can be extracted from a cell ID like so:

```c++ 
value = decoder -> get(hit.getCellID(), "field");
```

Note that while calorimter hits are used as the example, the same
process applies to any object that has a CellID.


## Inputs/Outputs

Processes `.edm4eic.root` output from EICrecon. 


## Dependencies

Needs `PODIO` and `edm4eic`. If you're running in the eic-shell,
then it can be run out-of-the-box with the command below.


## Usage

To run, do:
```
root -b -q ExtractCellIndices.cxx
```
