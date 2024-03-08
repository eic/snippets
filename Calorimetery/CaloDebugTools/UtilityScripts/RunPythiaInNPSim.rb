#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunPythiaInNPSim.rb'
# Derek Anderson
# 12.11.2023
#
# Run Pythia8 events (hepmc files) through DDSim.
# -----------------------------------------------------------------------------

# i/o parameters
in_hepmc = "pythia8NCDIS_18x275minQ100xAngleM0025hiDiv1_withBeamEffects_nLines2000K.hepmc"
out_file = "forNewGeoTest.nevt5000.epicmain_onlychcal.py8nc18x275minq100xangleM0025hidiv1.d11m12y2023.edm4hep.root"

# simulation parameters
compact = "$DETECTOR_PATH/epic_hcal_gdml.xml"
numevts = 5000

# run npsim
exec("npsim --compactFile #{compact} --numberOfEvents #{numevts} --inputFiles #{in_hepmc} --outputFile #{out_file}")

# end -------------------------------------------------------------------------
