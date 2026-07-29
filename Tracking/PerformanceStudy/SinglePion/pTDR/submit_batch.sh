#!/bin/bash
#SBATCH --array=1-42

## SL: uncomment this line to test the script locally before submission.
# SLURM_ARRAY_TASK_ID=5

# Path to the EPIC install that provides setup_local.sh.
setup_dir="parent_path_to_epic"
matmap=""
# matmap="-Pacts:MaterialMap=${setup_dir}/scripts/material_map/material-map.cbor" 
tag="GoldCoating"

nev=10000
nskip=0
mom=("0.5" "1" "2" "5" "10" "15")
mom_index=$((($SLURM_ARRAY_TASK_ID-1) % 6))
momentum="${mom[$mom_index]}"

# Array of theta ranges as strings "min max"
eta_ranges=("-3.5 -3.0" "-3.0 -2.5" "-2.5 -1" "-1 1"  "1 2.5" "2.5 3.0" "3.0 3.5")

eta_index=$((($SLURM_ARRAY_TASK_ID-1) / 6 % 7))
eta="${eta_ranges[$eta_index]}"
eta_min="$(echo $eta | cut -d ' ' -f 1)"
eta_max="$(echo $eta | cut -d ' ' -f 2)"

echo "Running simulation with settings:"
echo "Momentum: $momentum GeV"
echo "Eta Range: $eta_min to $eta_max"
echo "Setup dir: $setup_dir"
echo "Material Map: ${matmap}"

## -----first time using the shifter image, do------
shifter --image=eicweb/eic_xl:26.02.0-stable /opt/local/bin/eic-shell "\
  export eta_min=${eta_min} eta_max=${eta_max} mom=${momentum} matmap=\"${matmap}\"; \
  mkdir -p log; \
  run_sim_nohup.sh \
    ${setup_dir} ${tag} ${nev} ${nskip} ${matmap}"
