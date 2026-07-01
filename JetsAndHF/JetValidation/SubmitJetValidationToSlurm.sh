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
#SBATCH --output=slurm-26060py8ncdis10x100.out
#SBATCH --error=slurm-26060py8ncdis10x100.err

eic_shell=$HOME/.bin/eic-shell
version=26.05.0-stable

$eic_shell -v $version -- $PWD/RunJetValidationInShell.sh
