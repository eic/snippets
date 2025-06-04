#!/usr/bin/env ruby
# =============================================================================
# @file   RunEICReconWithTileMerging.rb
# @author Derek Anderson
# @date   12.19.2024
#
# This script generates the appropriate BHCal phi mapping and adjacency
# matrix based on the provided number of tiles to merge into towers.
#
# NOTE: the `ClusterMergedHits` plugin needs to be compiled
# beforehand.
# =============================================================================



# main body of script =========================================================

END {

  # i/o parameters
  in_ddsim   = "output/forBHCalOnlyCheck_rerun.e1pim.file0.d18m5y2025.edm4hep.root"
  out_podio  = "forHoleMystery.bhcalOnly_rerun_nMerge5_e1pim.d27m5y2025.podio.root"
  out_plugin = "forHoleMystery.bhcalOnly_rerun_nMerge5_e1pim.d27m5y2025.plugin.root"

  # output collections from EICrecon
  out_collect = [
    "HcalBarrelRecHits",
    "HcalBarrelMergedHits",
    "HcalBarrelClusters",
    "HcalBarrelMergedHitClusters",
    "HcalBarrelSplitMergeClusters",
    "GeneratedParticles"
  ].compact.reject(&:empty?).join(',')

  # plugins to run in EICrecon
  plugins = [
    "ClusterMergedHits"
  ].compact.reject(&:empty?).join(',')

  # options
  options = [
  ].compact.reject(&:empty?).join(' ')

  # add relevant mapping/matrix
  nmerge = if ARGV.empty? then 5 else ARGV[0] end
  add_map_and_matrix_to_options(nmerge, options)

  # run EICrecon
  exec("eicrecon -Pplugins=#{plugins} -Ppodio:output_collections=#{out_collect} #{options} -Ppodio:output_file=#{out_podio} #{in_ddsim}")

}  # end main body of script



# =============================================================================
# Make phi mapping
# -----------------------------------------------------------------------------
# @brief helper function to generate appropriate phi mapping for the
#   provided number of tiles to merge
#
# @param[in] nmerge number of tiles adjacent in phi to merge
# =============================================================================
def make_phi_mapping(nmerge)

  map = if nmerge.to_i > 1
    "\"phi:phi-(#{nmerge}*((phi/#{nmerge})-floor(phi/#{nmerge})))\""
  else
    "phi:phi"
  end
  return map

end  # end :make_phi_mapping



# =============================================================================
# Make adjacency matrix
# -----------------------------------------------------------------------------
# @brief helper function to generate appropriate adjacency matrix for
#   the provided number of tiles to merge
#
# @param[in] nmerge number of tiles adjacent in phi to merge
# =============================================================================
def make_adjacency_matrix(nmerge)

  # inject number to merge into matrix:
  #   (1) 1st condition: checks for vertically adjacent tiles
  #   (2) 2nd condition: checks for horizontally adjacenct tiles
  #       based on provided number to merge
  #   (3) 3rd condition: checks for tiles adjacent at horizontal
  #       wraparound based on provided number to merge
  # n.b. 320 is the number of tiles per row
  return <<-EOS.gsub(/^[\s\t]*/, '').gsub(/[\s\t]*\n/, ' ').strip
    "(
      ( (abs(eta_1 - eta_2) == 1) && (abs(phi_1 - phi_2) == 0) ) ||
      ( (abs(eta_1 - eta_2) == 0) && (abs(phi_1 - phi_2) == #{nmerge}) ) ||
      ( (abs(eta_1 - eta_2) == 0) && (abs(320 - abs(phi_1 - phi_2)) <= #{nmerge}) )
    ) == 1"
  EOS

end  # end :make_adjacency_matrix



# =============================================================================
# Generate map, matrix and add to options
# -----------------------------------------------------------------------------
# @brief helper function to generate appropriate mapping and adjacency
#   matrix and add to the provided list of EICrecon options
#
# @param[in]  nmerge  number of tiles adjacent in phi to merge
# @param[out] options list of options to append to
# =============================================================================
def add_map_and_matrix_to_options(nmerge, options)

  # generate approriate mapping, add to options
  mapping = make_phi_mapping(nmerge)
  matrix  = make_adjacency_matrix(nmerge)
  options.concat(" -PBHCAL:HcalBarrelMergedHits:fieldTransformations=#{mapping}")
         .concat(" -PClusterMergedHits:HcalBarrelMergedHitIslandProtoClusters:adjacencyMatrix=#{matrix}")

end

# end =========================================================================
