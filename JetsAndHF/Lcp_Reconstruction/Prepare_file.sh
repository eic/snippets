#!/bin/bash
# Shyam Kumar; INFN Bari, Italy
# put the campaign name as an agrument
xrdfs root://dtn-eic.jlab.org/ ls /volatile/eic/EPIC/RECO/25.04.1/epic_craterlake/SIDIS/D0_ABCONV/pythia8.306-1.1/10x100/q2_1/hiDiv > test.list
sed -i 's|^|root://dtn-eic.jlab.org//|' test.list

