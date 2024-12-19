#!/usr/bin/env ruby
# =============================================================================
# @file   RunEICRecon.rb
# @author Derek Anderson
# @date   05.11.2023
#
# An easy script to interface with EICRecon
# =============================================================================

# i/o parameters
in_ddsim   = "../forTestingJetPodioRelations.edm4hep.root"
out_podio  = "forTestingUserConfigurability_change1V2_withMinJetPtAndMinCstPtChanged_recoJetAlgoSetToGarbage_reco2gen12.podio.root"
out_plugin = "forTestingUserConfigurability_change1V2_withMinJetPtAndMinCstPtChanged_recoJetAlgoSetToGarbage_reco2gen12.plugin.root"

# output collections from EICrecon
out_collect = [
  "GeneratedParticles",
  "GeneratedJets",
  "GeneratedChargedJets",
  "ReconstructedJets",
  "ReconstructedChargedJets"
].compact.reject(&:empty?).join(',')

# plugins to run in EICrecon
plugins = [
  "dump_flags"
].compact.reject(&:empty?).join(',')

# options
options = [
  "-Peicrecon:LogLevel=debug"
].compact.reject(&:empty?).join(' ')

# run EICrecon
exec("eicrecon -Pplugins=#{plugins} -Ppodio:output_collections=#{out_collect} #{options} -Ppodio:output_file=#{out_podio} -Phistsfile=#{out_plugin} #{in_ddsim}")

# end =========================================================================
