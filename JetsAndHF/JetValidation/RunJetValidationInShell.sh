#!/bin/bash
# =============================================================================
#! @file    RunJetValidationInShell.sh
#! @authors Derek Anderson (derek.murphy.anderson@protonmail.com)
# -----------------------------------------------------------------------------
#! @brief Script to run JetValidation.C in the eic-shel.
#!
#! @usage
#!     eic-shell -- RunJetValidationInShell.sh
# =============================================================================

root -b -q JetValidation.C
