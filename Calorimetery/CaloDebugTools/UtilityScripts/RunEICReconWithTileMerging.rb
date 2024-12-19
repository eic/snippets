#!/usr/bin/env ruby
# =============================================================================
# @file   RunEICReconWithTileMerging.rb
# @author Derek Anderson
# @date   12.19.2024
#
# This script generates the appropriate BHCal phi mapping and adjacency
# matrix based on the provided number of tiles to merge into towers.
# =============================================================================



# main body of script =========================================================

END {

  # i/o parameters
  in_ddsim   = "../input/forBHCalOnlyCheck.e10pim.file0.d30m10y2024.edm4hep.root"
  out_podio  = "forTileMerger.afterReview_change0_removeMatrixMakers.d19m12y2024.podio.root"
  out_plugin = "forTileMerger.afterReview_change0_removeMatrixMakers.d19m12y2024.plugin.root"

  # output collections from EICrecon
  out_collect = [
    "HcalBarrelRecHits",
    "HcalBarrelMergedHits",
    "HcalBarrelClusters",
    "HcalBarrelSplitMergeClusters"
  ].compact.reject(&:empty?).join(',')

  # plugins to run in EICrecon
  plugins = [
    "dump_flags"
  ].compact.reject(&:empty?).join(',')

  # options
  options = [
    "-Pjana:nevents=100",
    "-Peicrecon:LogLevel=trace"
  ].compact.reject(&:empty?).join(' ')

  # add relevant mapping/matrix
  nmerge = if ARGV.empty? then 1 else ARGV[0] end 
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
    '"phi-(#{nmerge}*((phi/#{nmerge})-floor(phi/#{nmerge})))"'
  else
    "phi"
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
  options.concat(" -PBHCAL:HcalBarrelMergedHits:mappings=#{mapping}")
         .concat(" -PBHCAL:HcalBarrelIslandProtoClusters:adjacencyMatrix=#{matrix}")

end

# end =========================================================================
