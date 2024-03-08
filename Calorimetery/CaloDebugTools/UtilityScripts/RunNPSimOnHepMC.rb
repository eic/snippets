#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunNPSimOnHepMC.rb'
# Derek Anderson
# 02.22.2024
#
# Run a certain number events from a HepMC file
# in npsim with a specific compact file.
# -----------------------------------------------------------------------------

# output file
out_file = "forBHCalMultiPartSim.e10h11pipXpim.edm4hep.root"

# simulation parameters
numevts = 10000
compact = "$DETECTOR_PATH/epic_bhcal.xml"
hepmc   = "forBHCalMultiPartSim.e10h11pipXpim.hepmc"

# run npsim
exec("npsim -I #{hepmc} -N #{numevts} --compactFile #{compact} --outputFile #{out_file}")

# end -------------------------------------------------------------------------
