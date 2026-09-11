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
#!                                -f 100 \
#!                                -e 1000
#!
#!   Or if you're outside eic-shell:
#!     eic-shell -- RunJetValidationInShell.sh
# =============================================================================

out_path="."
out_suffix="files26071.py8ncdis10x100q100t1000"
file_list="filelists/files26071.py8ncdis10x100q100t1000.list"
num_files=-1
num_events=-1

while getopts "p:s:l:f:e:" opt; do
  case $opt in
    p) out_path=$OPTARG ;;
    s) out_suffix=$OPTARG ;;
    l) file_list=$OPTARG ;;
    f) num_files="$OPTARG" ;;
    e) num_events="$OPTARG" ;;
  esac
done

root -b -q "JetValidation.C(\"$out_path\", \"$out_suffix\", \"$file_list\", $num_files, $num_events)"
