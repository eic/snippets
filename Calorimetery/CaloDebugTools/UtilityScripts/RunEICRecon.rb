#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunEICRecon.rb'
# Derek Anderson
# 05.11.2023
#
# An easy script to interface with EICRecon
# -----------------------------------------------------------------------------

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
]

# plugins to run in EICrecon
plugins = [
  "dump_flags"
]

# options
options = [
  "-Preco:ReconstructedJets:minJetPt=2",
  "-Preco:ReconstructedJets:minCstPt=2"
]

# create collection argument
arg_collect = ""
num_collect = out_collect.size
out_collect.each_with_index do |collect, iCollect|
  arg_collect += collect
  arg_collect += "," if iCollect + 1 != num_collect
end

# create plugin argument
arg_plugin = ""
num_plugin = plugins.size
plugins.each_with_index do |plugin, iPlugin|
  arg_plugin += plugin
  arg_plugin += "," if iPlugin + 1 != num_plugin
end

# create options argument
arg_option = ""
num_option = options.size
options.each_with_index do |option, iOption|
  arg_option += option
  arg_option += " " if iOption + 1 != num_option
end

# run EICrecon
exec("eicrecon -Pplugins=#{arg_plugin} -Ppodio:output_collections=#{arg_collect} -Ppodio:output_file=#{out_podio} -Phistsfile=#{out_plugin} #{in_ddsim}")

# end -------------------------------------------------------------------------
