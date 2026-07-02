#!/bin/bash
# =============================================================================
#! @file    SubmitJetValidationToSlurm.sh
#! @authors Derek Anderson (derek.murphy.anderson@protonmail.com)
# -----------------------------------------------------------------------------
#! @brief Submit this script to run JetValidation.C with slurm. 
#!
#! @usage Update relevant slurm options, eic-shell location,
#!   and shell version below. Then submit with:
#!
#!       sbatch SubmitJetValidationToSlurm.sh
#!
#!   Can also run directly with
#!
#!       ./SubmitJetValidationToSlurm.sh
# =============================================================================
#SBATCH --partition=<my partition>
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --account=eic
#SBATCH --mail-user=<my email>
#SBATCH --mail-type=END,FAIL
#SBATCH --output=<output log name>
#SBATCH --error=<error log name>

eic_shell=$HOME/.bin/eic-shell # NOTE edit me!
version=26.05.0-stable # NOTE edit me!

$eic_shell -v $version -- $PWD/RunJetValidationInShell.sh
