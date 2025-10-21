#!/bin/bash
# example script to generate small sample of background mixed event file locally
# you should always use background files from the simulation campaign data (https://eic.github.io/epic-prod/) when possible 
# Shujie Li, 10.2025

# pre-requists:
#        1.  install SignalBackgroundMerger from https://github.com/eic/HEPMC_Merger
#        2.  recommend to run within eic-shell to use the xrdfs command

## run as: ./run_merger.sh, then follow the options.

# For SR and electron beam gas, please also read README files under JLab ifarm:
#       SynRad: /volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6
#       ESR beam losses: /volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC

## 1000 merged events at 18GeV takes about 2 hours for dd4hep+eicrecon.
## the 5 and 10 GeV mixed events will be 10x slower due to higher SR freq. Use the simulation campaign file if you can. 
nev=100
window=2000 #integration window of 2us --> 2000ns is the default
skip=0 # skip the first N signal events
seed=42
## signal frequency 
#        =0: forced to have one signal per time slice/window, 
#        >0: use the actual signal and background freq in sampling, 
#        <0: (depracated) background with the default minbias (Q2<1 SIDIS) events as signal. 
# Mixing frequencies can be found on ePIC Wiki: https://wiki.bnl.gov/EPIC/index.php?title=Background
# Explanations on mixing scheme at https://github.com/eic/eic.github.io/blob/master/_resources/background_mixed_samples.md
sf=0

if [ "$sf" -gt 0 ]; then
    sig_type="dis"
elif [ "$sf" -eq 0 ]; then
    sig_type="forced"
else
    sig_type="minbias"
fi

###############
## background sources 
## SR, brems, coulomb, touschek, pgas
## generatorStatus ID:
## 2000, 3000, 4000, 5000, 6000
###############
echo "Please select an option:"
echo "1: DIS 18x275"
echo "2: DIS 18x275, no SR"
echo "3: minbias 18x275"
echo "4: DIS 10x275, SR scaled from 18GeV"
echo "5: DIS 5x100, SR scaled from 18GeV"
echo "6: Backgrounds only 18x275"

# Read user input
read -p "Enter option number (1-6): " option
if   [[ "$option" == "1" ]]; then
    ebeam=18
    pbeam=275
    tag=""
    signal="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.hepmc"
    bg_files=(
        # "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/SYNRAD/dataprod_rel_1.0.0/18x275/dataprod_rel_1.0.0_synrad_18x275_run001.hepmc3.tree.root" ## old SR samples with all inner beampipe photons
        "root://dtn-eic.jlab.org//volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6/data/synrad/dataprod_rel_1.0.0/18x275/dataprod_rel_1.0.0_synrad_18x275_run001.preproc_10000repeats.hepmc3.tree.root" ## new SR with outer beampipe photons only, freq=3.3GHz
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electronbrems/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electronbrems_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electroncoulomb/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electroncoulomb_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electrontouschek/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electrontouschek_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/275GeV/pythia8.306-1.0_ProtonBeamGas_275GeV_run001.hepmc3.tree.root"
    )
    freqs=(3300000 18.26 0.86 0.55 22.5) ## kHz
    statuses=(2000 3000 4000 5000 6000)

elif [[ "$option" == "2" ]]; then
    ebeam=18
    pbeam=275
    tag="_noSR_"
    signal="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.hepmc"
    bg_files=(
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electronbrems/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electronbrems_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electroncoulomb/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electroncoulomb_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electrontouschek/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electrontouschek_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/275GeV/pythia8.306-1.0_ProtonBeamGas_275GeV_run001.hepmc3.tree.root"
    )
    freqs=(18.26 0.86 0.55 22.5) ## kHz
    statuses=(3000 4000 5000 6000)

elif [[ "$option" == "3" ]]; then
    ebeam=18
    pbeam=275
    tag="_minbias_"
    signal="pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root"
    bg_files=(
        "dataprod_rel_1.0.0_synrad_18x275_run001.preproc_10000repeats.hepmc3.tree.root" ## new SR with outer beampipe photons only, freq=3.3GHz
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electronbrems/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electronbrems_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electroncoulomb/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electroncoulomb_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electrontouschek/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electrontouschek_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/275GeV/pythia8.306-1.0_ProtonBeamGas_275GeV_run001.hepmc3.tree.root"
    )
    freqs=(3300000 18.26 0.86 0.55 22.5) ## kHz
    statuses=(2000 3000 4000 5000 6000)

elif [[ "$option" == "4" ]]; then
    ebeam=10
    pbeam=275
    tag="_scaled_SR_"
    signal="pythia8NCDIS_10x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.hepmc"

    bg_files=(
    "root://dtn-eic.jlab.org//volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6/data/synrad/dataprod_rel_1.0.0/10x275/dataprod_rel_1.0.0_synrad_10x275_run001.preproc_10000repeats.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electronbrems/dataprod_rel_1.0.1/10x275/dataprod_rel_1.0.1_electronbrems_10x275_50sec.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electroncoulomb/dataprod_rel_1.0.1/10x275/dataprod_rel_1.0.1_electroncoulomb_10x275_50sec.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electrontouschek/dataprod_rel_1.0.1/10x275/dataprod_rel_1.0.1_electrontouschek_10x275_50sec.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/275GeV/pythia8.306-1.0_ProtonBeamGas_275GeV_run001.hepmc3.tree.root"
    )

    freqs=(36608000 172.31 29.56 233.50 32.6) ## kHz
    statuses=(2000 3000 4000 5000 6000)

elif [[ "$option" == "5" ]]; then
    ebeam=5
    pbeam=100
    tag="_scaled_SR_"
    signal="root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/NC/5x100/minQ2=1/pythia8NCDIS_5x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.hepmc3.tree.root"
    bg_files=(
    "root://dtn-eic.jlab.org//volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6/data/synrad/dataprod_rel_1.0.0/5x100/dataprod_rel_1.0.0_synrad_5x100_run001.preproc_10000repeats.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electronbrems/dataprod_rel_1.0.1/5x100/dataprod_rel_1.0.1_electronbrems_5x100_50sec.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electroncoulomb/dataprod_rel_1.0.1/5x100/dataprod_rel_1.0.1_electroncoulomb_5x100_50sec.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electrontouschek/dataprod_rel_1.0.1/5x100/dataprod_rel_1.0.1_electrontouschek_5x100_50sec.hepmc3.tree.root"
    "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/100GeV/pythia8.306-1.0_ProtonBeamGas_100GeV_run001.hepmc3.tree.root"
    )
    statuses=(2000 3000 4000 5000 6000) ## kHz
    freqs=(36608000 328.04 116.57 1112.3 22)

elif [[ "$option" == "6" ]]; then
    ebeam=18
    pbeam=275
    sig_type="bgOnly"
    tag=""
    signal="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.hepmc"
    bg_files=(
        # "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/SYNRAD/dataprod_rel_1.0.0/18x275/dataprod_rel_1.0.0_synrad_18x275_run001.hepmc3.tree.root" ## old SR samples with all inner beampipe photons
        "root://dtn-eic.jlab.org//volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6/data/synrad/dataprod_rel_1.0.0/18x275/dataprod_rel_1.0.0_synrad_18x275_run001.preproc_10000repeats.hepmc3.tree.root" ## new SR with outer beampipe photons only, freq=3.3GHz
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electronbrems/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electronbrems_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electroncoulomb/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electroncoulomb_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC/data/xsuite/electrontouschek/dataprod_rel_1.0.1/18x275/dataprod_rel_1.0.1_electrontouschek_18x275_50sec.hepmc3.tree.root"
        "root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/275GeV/pythia8.306-1.0_ProtonBeamGas_275GeV_run001.hepmc3.tree.root"
    )
    freqs=(3300000 18.26 0.86 0.55 22.5) ## kHz
    statuses=(2000 3000 4000 5000 6000)

else
    echo "Invalid option selected."
    exit 1
fi

## ------------------------------------------------
## assemble command
## ------------------------------------------------
outfile=bgmerged_${sig_type}_${ebeam}x${pbeam}${tag}n${nev}.hepmc3.tree.root
cmd="./SignalBackgroundMerger  -N ${nev} --rngSeed ${seed} -o merged/${outfile}"

## add signal
if [[ "$option" == "6" ]]; then
    cmd+=" -i $signal -sf 0.01 -S $nskip -w $window"
elif [ "$sf" -ge 0 ]; then ## use DIS signal file
    cmd+=" -i $signal -sf $sf -S $nskip -w $window"
else ## use the default from the merger which is the sidis sample
    cmd+=" -i pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root -sf 83 "
fi

## add background
for bg_file in "${bg_files[@]}"; do
    # Extract the filename from the remote path
    filename=$(basename "$bg_file")

    # Check if the file exists locally
    if [ -f "$filename" ]; then
        echo "File '$filename' already exists locally."
    else
        echo "File '$filename' not found locally. Downloading..."
        # Use xrdcp to download the file
        xrdcp "$bg_file" "$filename"

        # Check if the download succeeded
        if [ $? -eq 0 ]; then
            echo "File '$filename' downloaded successfully."
        else
            echo "Failed to download file '$filename'."
            exit
        fi
    fi
done

length=${#bg_files[@]}

# Loop over the arrays using an index
for ((i=0; i<length; i++)); do
    filename=$(basename "${bg_files[i]}")
    cmd+=" -b $filename ${freqs[i]} 0 ${statuses[i]} "
done

# Print the final command
echo "Generated command:"
echo "$cmd"

# Execute the command
mkdir log
mkdir merged

eval "$cmd" >> log/$outfile.log

echo "DONE! See output under log/ and merged/."
echo Please add the following flags to your npsim command:
echo       --physics.alternativeStableStatuses="2001 3001 4001 5001 6001"  --physics.alternativeDecayStatuses="2002 3002 4002 5002 6002"
# ./SignalBackgroundMerger -i DIS_10x275_1.14kHz.hepmc -sf 0 -b proton_beamgas_275GeV_350kHz.hepmc3.tree.root 350 0 2000 -b electron_beamgas_10GeV_3177kHz.hepmc3.tree.root 3177 0 4000 -b synrad_18x275_14000kHz.hepmc.tree.root 14000 0 6000 -o test_new.hepmc -N 100 --rngSeed 42 
