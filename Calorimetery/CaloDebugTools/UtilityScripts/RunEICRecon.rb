#!/usr/bin/env ruby
# =============================================================================
# @file   RunEICRecon.rb
# @author Derek Anderson
# @date   05.11.2023
#
# An easy script to interface with EICRecon
# =============================================================================

# i/o parameters
in_ddsim   = "forPFA1.epic26041e10th33pim.edm4hep.root"
out_podio  = "test.eicrecon.root"

# collections to output
out_collect = [
  "HcalBarrelHits",
  "HcalBarrelRawHits",
  "HcalBarrelRawHitLinks",
  "HcalBarrelRecHits",
  "HcalBarrelClusters",
].compact.reject(&:empty?).join(',')

# plugins to run
plugins = [
  "dump_flags",
].compact.reject(&:empty?).join(',')

# additional options
options = [
  "-Peicrecon:LogLevel=debug",
].compact.reject(&:empty?).join(' ')

# set options appropriately
out_opt  = out_collect.empty? ? "" : "-Ppodio:output_collections=#{out_collect}"
plug_opt = plugins.empty? ? "" : "-Pplugins=#{plugins}"
opt_opt  = options.empty? ? "" : "#{options}"

# run EICrecon
exec("eicrecon #{plug_opt} #{out_opt} #{opt_opt} -Ppodio:output_file=#{out_podio} #{in_ddsim}")

# end =========================================================================
