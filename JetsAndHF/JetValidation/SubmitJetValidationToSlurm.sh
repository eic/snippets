#!/bin/bash
# =============================================================================
#! @file    SubmitJetValidationToSlurm.sh
#! @authors Derek Anderson (derek.murphy.anderson@protonmail.com)
# -----------------------------------------------------------------------------
#! @brief Submit this script to run JetValidation.C with slurm. 
#!
#! @usage
#!     sbatch SubmitJetValidationToSlurm.sh
# =============================================================================
#SBATCH --partition=ifarm
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --account=eic
#SBATCH --mail-user=dereka@jlab.org
#SBATCH --mail-type=END,FAIL

eic_shell=$HOME/.bin/eic-shell
version=25.06.0-stable

$eic_shell -v $version -- $PWD/RunJetValidationInShell.sh
