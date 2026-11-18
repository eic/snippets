#!/bin/bash
#Script to prepare 
# shyam kumar; shyam055119@gmail.com
#rucio did list --short "epic:/RECO/26*/*10x250*" 
#Without machine background:ep (9x130) D0 Sample; and ep (9x130) DIS Sample for ES
rm test.list
D0sample=false
DISSample=true
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
  get_files "epic:/RECO/26.07.1/epic_craterlake/SIDIS/D0_ABCONV/HFsim-PYTHIA/pythia8.316-1.0/ep/9x130/q2_1" > $outfile

elif $DISSample; then
  echo "Running DIS sample..."
  (
    get_files "epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x130/q2_100to1000"
    get_files "epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x130/q2_10to100"
    get_files "epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x130/q2_1to10"
    get_files "epic:/RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/9x130/q2_1000toINF"
  ) | sort -u > $outfile

else
  echo "No sample selected!"
fi
