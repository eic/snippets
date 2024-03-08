#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunDDSim.rb'
# Derek Anderson
# 09.08.2023
#
# Run a certain number events in ddsim based on
# specified steering and compact files.
# -----------------------------------------------------------------------------

# output file
out_file = "test_ruby_script.edm4hep.root"

# simulation parameters
numevts = 100
steerer = "../steering/steering.forTowerVsTileCalibCheck_e10th45pim.py"
compact = "$DETECTOR_PATH/epic_imaging.xml"

# run ddsim
exec("ddsim --steeringFile #{steerer} --compactFile #{compact} -G -N #{numevts} --outputFile #{out_file}")

# end -------------------------------------------------------------------------
