#!/bin/bash
#Script to prepare 
# shyam kumar; shyam055119@gmail.com
#rucio did list --short "epic:/RECO/26*/*10x250*" 
#Without machine background D0 sample (10x250)
rm test.list
D0sample=true
LcSample=false
DISSample=false
outfile="test.list"

get_files () {
  DID=$1
  RSE=$(rucio list-rules "$DID" | grep -E "BNL-XRD|EIC-XRD" | grep -v "JLAB-TAPE-SE" | awk '{print $5}' | head -n1)

  if [ -z "$RSE" ]; then
    echo "No BNL-XRD/EIC-XRD found for $DID" >&2
    return
  fi

  echo "Using RSE: $RSE for $DID" >&2
  rucio replica list file --pfns --rses "$RSE" "$DID"
}

if $D0sample; then
  echo "Running D0 sample..."
  get_files "epic:/RECO/26.04.1/epic_craterlake/SIDIS/D0_ABCONV/HFsim-PYTHIA/pythia8.312-2.0/ep/10x250/q2_1" > $outfile

elif $LcSample; then
  echo "Running Lc sample..."
  get_files "epic:/RECO/26.04.1/epic_craterlake/SIDIS/Lc_ABCONV/HFsim-PYTHIA/pythia8.312-2.0/ep/10x250/q2_1" > $outfile

elif $DISSample; then
  echo "Running DIS sample..."
  (
    get_files "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_1to10"
    get_files "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_10to100"
    get_files "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_100to1000"
    get_files "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_1000to10000"
  ) | sort -u > $outfile

else
  echo "No sample selected!"
fi
