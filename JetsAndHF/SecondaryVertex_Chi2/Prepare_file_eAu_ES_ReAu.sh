#!/bin/bash
#Script to prepare 
# shyam kumar; shyam055119@gmail.com
#rucio did list --short "epic:/RECO/26*/*10x250*" 
#Without machine background D0 sample for e+p 10x130 and e+Au 10x100
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
  get_files "epic:/RECO/26.04.1/epic_craterlake_without_zdc/SIDIS/D0_ABCONV/HFsim-BeAGLE/BeAGLE1.03.01-2.0/eAu/10x100/q2_1to10000" > $outfile

elif $DISSample; then
  echo "Running DIS sample..."
  (
    get_files "epic:/RECO/26.04.1/epic_craterlake_without_zdc/DIS/BeAGLE1.03.02-2.0/eAu/10x100/q2_1to10"
    get_files "epic:/RECO/26.04.1/epic_craterlake_without_zdc/DIS/BeAGLE1.03.02-2.0/eAu/10x100/q2_10to100"
  ) | sort -u > $outfile

else
  echo "No sample selected!"
fi
