#!/bin/bash
# =============================================================================
#! @file    RunJetValidationInShell.sh
#! @authors Derek Anderson (derek.murphy.anderson@protonmail.com)
# -----------------------------------------------------------------------------
#! @brief Script to run JetValidation.C in the eic-shell.
#!
#! @usage If you're already inside eic-shell:
#!     RunJetValidationInShell.sh -l "my_filelist.list" \
#!                                -n 100 \
#!                                -o "./path/to/my/output"
#!   Or if your outside eic-shell:
#!     eic-shell -- RunJetValidationInShell.sh
# =============================================================================

file_list="filelists/files26071.py8ncdis10x100q100t1000.list"
num_files=10
out_path="."

while getopts "l:u:o:" opt; do
  case $opt in
    l) file_list=$OPTARG ;;
    u) num_files="$OPTARG" ;;
    o) out_path=$OPTARG ;;
  esac
done

root -b -q "JetValidation.C(\"$file_list\", $num_files, \"$out_path\")"
