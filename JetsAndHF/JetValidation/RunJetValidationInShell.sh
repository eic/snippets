#!/bin/bash
# =============================================================================
#! @file    RunJetValidationInShell.sh
#! @authors Derek Anderson (derek.murphy.anderson@protonmail.com)
# -----------------------------------------------------------------------------
#! @brief Script to run JetValidation.C in the eic-shell.
#!
#! @usage If you're already inside eic-shell:
#!     RunJetValidationInShell.sh -p "./plots" \
#!                                -s "my_plots" \
#!                                -l "my_filelist.list" \
#!                                -h "my_hist_file.root" \
#!                                -f 100 \
#!                                -e 1000
#!
#!   Or if you're outside eic-shell:
#!     eic-shell -- RunJetValidationInShell.sh
# =============================================================================

out_path="."
out_suffix="files26071.py8ncdis10x100q100t1000"
file_list="filelists/files26071.py8ncdis10x100q100t1000.list"
hist_file="hists.$out_suffix.root"
num_files=-1
num_events=-1

while getopts "p:s:l:h:f:e:" opt; do
  case $opt in
    p) out_path=$OPTARG ;;
    s) out_suffix=$OPTARG ;;
    l) file_list=$OPTARG ;;
    h) hist_file=$OPTARG ;;
    f) num_files="$OPTARG" ;;
    e) num_events="$OPTARG" ;;
  esac
done

root -b -q "MakeJetValidationHists.C(\"$out_path\", \"$out_suffix\", \"$file_list\", $num_files, $num_events)"
root -b -q "MakeJetValidationPlots.C(\"$out_path\", \"$out_suffix\", \"$hist_file\")"
