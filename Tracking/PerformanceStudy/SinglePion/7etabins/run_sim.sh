#!/bin/bash

DIR=$(cd -- $(dirname -- "${BASH_SOURCE[0]}") &> /dev/null  && pwd)
export LD_LIBRARY_PATH="${DIR}/install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "In directory $PWD!"

# Set ePIC software version
source /opt/detector/epic-main/bin/thisepic.sh

# Read HTCondor ProcID
procid="$1"

# Set Number of events
NEVENTS=1000

# Create Theta/Eta and Momentum arrays
theta_max=(
  3.08121623
  3.04210067
  2.9777899
  2.43656581
  0.70502684
  0.16380276
  0.09949199
)

theta_min=(
  3.04210067
  2.9777899
  2.43656581
  0.70502684
  0.16380276
  0.09949199
  0.06037642
)

eta_low=(
  -3.5
  -3.0
  -2.5
  -1.0
  1.0
  2.5
  3.0
)

eta_hi=(
  -3.0
  -2.5
  -1.0
  1.0
  2.5
  3.0
  3.5
)

momentum=(0.5 1 2 5 10 15)

# Compute indices
Nmom=${#momentum[@]}
Ntheta=${#theta_max[@]}

theta_idx=$(( procid / Nmom ))
mom_idx=$(( procid % Nmom ))

if (( theta_idx < 0 || theta_idx >= Ntheta )); then
  echo "ERROR: ProcID $procid out of range"
  echo "Max ProcID allowed: $(( Ntheta * Nmom - 1 ))"
  exit 1
fi

# Set values for this run
theta_max_val=${theta_max[$theta_idx]}
theta_min_val=${theta_min[$theta_idx]}
eta_low_val=${eta_low[$theta_idx]}
eta_hi_val=${eta_hi[$theta_idx]}
mom_val=${momentum[$mom_idx]}

echo "================================="
echo "ProcID      : $procid"
echo "theta index : $theta_idx"
echo "mom index   : $mom_idx"
echo "theta-max value : $theta_max_val"
echo "theta-min value : $theta_min_val"
echo "eta-low value : $eta_low_val"
echo "eta-hi value : $eta_hi_val"
echo "mom value   : $mom_val"
echo "num events  : $NEVENTS"
echo "================================="

# Run simulation
npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
      --gun.thetaMax $theta_max_val --gun.thetaMin $theta_min_val \
      --gun.momentumMin "${mom_val}*GeV" --gun.momentumMax "${mom_val}*GeV" \
      --gun.particle pi- --numberOfEvents ${NEVENTS} --outputFile output.edm4hep.root

eicrecon -Ppodio:output_file=eicrecon_out.root \
         -Pjana:nevents=${NEVENTS} \
         -Pdd4hep:xml_files=epic_craterlake.xml output.edm4hep.root

# Clean up
mv output.edm4hep.root sim_default_eta_${eta_low_val}_${eta_hi_val}_${mom_val}GeV_${NEVENTS}.edm4hep.root
mv eicrecon_out.root rec_default_eta_${eta_low_val}_${eta_hi_val}_${mom_val}GeV_${NEVENTS}.root

echo "Listing files after running:"
ls
echo ""

echo "Done!"

