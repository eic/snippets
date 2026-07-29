#!/bin/bash

set -e  # exit if any command fails

# ===== SETUP =====
if [ -d "$1" ]; then
    setup_dir="$1"
    shift
else
    echo "ERROR: setup directory '$1' does not exist." >&2
    exit 1
fi
cd "$setup_dir"
source your_script_to_source_and_compile_epic
cd - > /dev/null

# ===== INPUT =====
tag_in=$1
nev=$2
nskip=$3
matmap=$4
detector_xml="${DETECTOR_PATH}/epic_craterlake_tracking_only.xml" ## use tracking_only in dd4hep to save time
detector_xml_full="${DETECTOR_PATH}/epic.xml" ## use full geometry in eicrecon to avoid exccesive amount of error msg from non-tracking detector algorithms.

setup_base="$(basename "$setup_dir")"
tag="${tag_in}${setup_base}_eta_${eta_min}_${eta_max}_mom_${mom}GeV_n${nev}"

# ===== LOG FILE =====
logfile="log/${tag}.log"

# Also copy nohup output to our log file
{
    echo "=== Starting run for $tag ==="
    date
    echo ""

    # ===== STEP 1: npsim =====
    echo "[`date`] Running npsim..."
    npsim --compactFile "$detector_xml" \
          --outputFile "sim_${tag}.root" \
          -N "${nev}" \
          --enableGun --gun.particle pi+ --gun.multiplicity 1 --gun.distribution uniform \
          --gun.etaMin "${eta_min}" --gun.etaMax "${eta_max}" \
          --gun.phiMin "0*degree" --gun.phiMax "360*degree" \
          --gun.momentumMin "${mom}*GeV" --gun.momentumMax "${mom}*GeV"

    echo "[`date`] npsim finished."
    echo ""

    # ===== STEP 2: eicrecon =====
    echo "[`date`] Running eicrecon..."
    eicrecon -Pjana:nevents="${nev}" \
             -Ppodio:output_file="rec_${tag}.root" \
             -Pdd4hep:xml_files="${detector_xml_full}" \
             ${matmap} \
             "sim_${tag}.root"

    echo "[`date`] eicrecon finished."
    echo ""
    date
    echo "=== Run completed for $tag ==="

} > "$logfile" 2>&1
