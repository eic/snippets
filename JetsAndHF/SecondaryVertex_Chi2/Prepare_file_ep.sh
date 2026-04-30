#!/bin/bash
#Script to prepare 
# shyam kumar; shyam055119@gmail.com
#rucio did list --short "epic:/RECO/26*/*10x250*" 
#Without machine background D0 sample (10x250)
D0sample=false
LcSample=true
DISSample=false
# Output file
outfile="test.list"

if $D0sample; then
echo "Running D0 sample..."
rucio did content list --short "epic:/RECO/26.04.1/epic_craterlake/SIDIS/D0_ABCONV/HFsim-PYTHIA/pythia8.312-2.0/ep/10x250/q2_1" | sed 's|^epic:/|root://dtn-eic.jlab.org:1094//volatile/eic/EPIC/|' > $outfile
elif $LcSample; then
echo "Running Lc sample..."
rucio did content list --short "epic:/RECO/26.04.1/epic_craterlake/SIDIS/Lc_ABCONV/HFsim-PYTHIA/pythia8.312-2.0/ep/10x250/q2_1" | sed 's|^epic:/|root://dtn-eic.jlab.org:1094//volatile/eic/EPIC/|' > $outfile
elif $DISSample; then
#Without machine background DIS sample (10x250)
(
rucio did content list --short "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_1to10"
rucio did content list --short "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_10to100"
rucio did content list --short "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_100to1000"
rucio did content list --short "epic:/RECO/26.04.1/epic_craterlake/DIS/pythia6.428-1.0/NC/noRad/ep/10x250/q2_1000to10000"
) | sed 's|^epic:/|root://dtn-eic.jlab.org:1094//volatile/eic/EPIC/|' | sort -u > $outfile
else
echo "No sample selected!"
fi
