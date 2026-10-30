This directory contains scripts to create the .gdml files needed for the SVT OB simulation geometry from the CAD files. 

The .STEP, .FCStd, and .gdml files are at https://github.com/eic/epic-data/tree/main/SVTOB_gdml

The .FCStd files are made by loading the .STEP file into FreeCAD and then saving as .FCStd.

The .gdml files can by made using the scripts processL3.sh and processL4.sh. To run these the file freecadcmd_dir must contain the address of the freecadcmd executable. These scripts do the following steps:

Convert the .FCStd file to a directory of .stl files.
Delete the .stl files for components that are not needed.
Remove the parts extruding the ends of the envelope volume.
Merge overlapping carbon-fibre components
Scale the K9 components by 0.97 to avoid overlaps with the carbon fibre
Convert the .stl files to ASCII and then to .gdml
Delete .stl files and leave a folder of .gdml files

