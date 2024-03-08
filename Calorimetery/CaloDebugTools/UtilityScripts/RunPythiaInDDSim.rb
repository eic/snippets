#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunPythiaInDDSim.rb'
# Derek Anderson
# 06.29.2023
#
# Run Pythia8 events (hepmc files) through DDSim.
# -----------------------------------------------------------------------------

# i/o parameters
in_hepmc = "pythia8NCDIS_18x275_minQ2_1_nLines100K.hepmc"
out_file = "forDynamicRange.py8NCDIS_18x275_minMomTransfer1_100Klines.evt1Kimage.d29m6y2023.edm4hep.root"

# simulation parameters
compact = "$DETECTOR_PATH/$DETECTOR_CONFIG.xml"
numevts = 1000

# run ddsim
exec("ddsim --compactFile #{compact} --numberOfEvents #{numevts} --inputFiles #{in_hepmc} --outputFile #{out_file}")

# end -------------------------------------------------------------------------
