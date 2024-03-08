#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunNPSim.rb'
# Derek Anderson
# 10.30.2023
#
# Run a certain number events in npsim based on
# specified steering and compact files.
# -----------------------------------------------------------------------------

# output file
out_file = "test_ruby_script.edm4hep.root"

# simulation parameters
numevts = 1000
steerer = "../steering/steering.forTowerVsTileCalibCheck_e10th45pim.py"
compact = "$DETECTOR_PATH/epic_imaging.xml"

# run ddsim
exec("npsim --steeringFile #{steerer} --compactFile #{compact} -G -N #{numevts} --outputFile #{out_file}")

# end -------------------------------------------------------------------------
