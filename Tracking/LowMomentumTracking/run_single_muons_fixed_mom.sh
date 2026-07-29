#!/usr/bin/bash

source /opt/detector/epic-main/bin/thisepic.sh

#source ../../epic_local/install/bin/thisepic.sh #Use local ePIC version
#source ../../eicrecon_local/install/bin/eicrecon-this.sh #Use local EICRecon version

# Running options
run_set=${1:-0}   # Momentum setting to run
run_npsim=${2:-1} # Run npsim (yes/no)
run_local_map=${3:-0} # 1: Use local material map; 2: Use alternative barrel MPGD digitization
nevents=${4:-1000} # Number of events to run
delta_phi_max=${5:-0.085} # deltaPhiMax seeding parameter

# Setting 1
# Single negative muons at 500 MeV/c from (0,0,0)
if [ "$run_set" -eq 1 ] || [ "$run_set" -eq 0 ]; then
	echo "Running Setting 1!!!"
	echo ""
	
	if [ "$run_npsim" -eq 1 ] ; then
		npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
		--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "0.5*GeV" --gun.momentumMax "0.5*GeV" \
		--numberOfEvents ${nevents} --outputFile output/output_500MeV.edm4hep.root
	fi

	if [ "$run_local_map" -eq 1 ] ; then
		eicrecon -Ppodio:output_file=output/eicrecon_out_500MeV.root \
		-Pjana:nevents=${nevents} \
		-Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
		-Pacts:MaterialMap=material-map.cbor \
		-Pdd4hep:xml_files=epic_craterlake.xml output/output_500MeV.edm4hep.root | tee output/log/out_500MeV.txt
	else
		eicrecon -Ppodio:output_file=output/eicrecon_out_500MeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_500MeV.edm4hep.root | tee output/log/out_500MeV.txt
	fi

fi

# Setting 2
# Single negative muons at 1 GeV/c from (0,0,0)
if [ "$run_set" -eq 2 ] || [ "$run_set" -eq 0 ]; then
	echo "Running Setting 2!!!"
	echo ""

	if [ "$run_npsim" -eq 1 ] ; then
		npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
		--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "1.0*GeV" --gun.momentumMax "1.0*GeV" \
		--numberOfEvents ${nevents} --outputFile output/output_1GeV.edm4hep.root
	fi

	if [ "$run_local_map" -eq 1 ] ; then
                eicrecon -Ppodio:output_file=output/eicrecon_out_1GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pacts:MaterialMap=material-map.cbor \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_1GeV.edm4hep.root | tee output/log/out_1GeV.txt
        else
                eicrecon -Ppodio:output_file=output/eicrecon_out_1GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_1GeV.edm4hep.root | tee output/log/out_1GeV.txt
        fi
	
fi

# Setting 3
# Single negative muons at 2 GeV/c from (0,0,0)
if [ "$run_set" -eq 3 ] || [ "$run_set" -eq 0 ]; then
	echo "Running Setting 3!!!"
	echo ""

	if [ "$run_npsim" -eq 1 ] ; then
		npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
		--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "2.0*GeV" --gun.momentumMax "2.0*GeV" \
		--numberOfEvents ${nevents} --outputFile output/output_2GeV.edm4hep.root
	fi

	if [ "$run_local_map" -eq 1 ] ; then
                eicrecon -Ppodio:output_file=output/eicrecon_out_2GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pacts:MaterialMap=material-map.cbor \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_2GeV.edm4hep.root | tee output/log/out_2GeV.txt
	elif [ "$run_local_map" -eq 2 ] ; then
                eicrecon -Ppodio:output_file=output/eicrecon_out_MpgdDigi_2GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -PMPGD:SiFactoryPattern="0x3" \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_2GeV.edm4hep.root | tee output/log/out_2GeV.txt
        else
                eicrecon -Ppodio:output_file=output/eicrecon_out_2GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_2GeV.edm4hep.root | tee output/log/out_2GeV.txt
        fi

fi

# Setting 4
# Single negative muons at 5 GeV/c from (0,0,0)
if [ "$run_set" -eq 4 ] || [ "$run_set" -eq 0 ]; then
	echo "Running Setting 4!!!"
	echo ""

	if [ "$run_npsim" -eq 1 ] ; then
		npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
		--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "5.0*GeV" --gun.momentumMax "5.0*GeV" \
		--numberOfEvents ${nevents} --outputFile output/output_5GeV.edm4hep.root
	fi

	if [ "$run_local_map" -eq 1 ] ; then
                eicrecon -Ppodio:output_file=output/eicrecon_out_5GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pacts:MaterialMap=material-map.cbor \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_5GeV.edm4hep.root | tee output/log/out_5GeV.txt
        else
                eicrecon -Ppodio:output_file=output/eicrecon_out_5GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_5GeV.edm4hep.root | tee output/log/out_5GeV.txt
        fi

fi

# Setting 5
# Single negative muons at 10 GeV/c from (0,0,0)
if [ "$run_set" -eq 5 ] || [ "$run_set" -eq 0 ]; then
        echo "Running Setting 5!!!"
        echo ""

	if [ "$run_npsim" -eq 1 ] ; then
        	npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
        	--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "10.0*GeV" --gun.momentumMax "10.0*GeV" \
        	--numberOfEvents ${nevents} --outputFile output/output_10GeV.edm4hep.root
	fi

	if [ "$run_local_map" -eq 1 ] ; then
                eicrecon -Ppodio:output_file=output/eicrecon_out_10GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pacts:MaterialMap=material-map.cbor \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_10GeV.edm4hep.root | tee output/log/out_10GeV.txt
        else
                eicrecon -Ppodio:output_file=output/eicrecon_out_10GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_10GeV.edm4hep.root | tee output/log/out_10GeV.txt
        fi

fi

# Setting 6
# Single negative muons at 20 GeV/c from (0,0,0)
if [ "$run_set" -eq 6 ] || [ "$run_set" -eq 0 ]; then
        echo "Running Setting 6!!!"
        echo ""

	if [ "$run_npsim" -eq 1 ] ; then
        	npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
        	--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "20.0*GeV" --gun.momentumMax "20.0*GeV" \
       		--numberOfEvents ${nevents} --outputFile output/output_20GeV.edm4hep.root
	fi

	if [ "$run_local_map" -eq 1 ] ; then
                eicrecon -Ppodio:output_file=output/eicrecon_out_20GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pacts:MaterialMap=material-map.cbor \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_20GeV.edm4hep.root | tee output/log/out_20GeV.txt
        else
                eicrecon -Ppodio:output_file=output/eicrecon_out_20GeV.root \
                -Pjana:nevents=${nevents} \
                -Ptracking:CentralTrackSeedingResults:deltaPhiMax=${delta_phi_max} \
                -Pdd4hep:xml_files=epic_craterlake.xml output/output_20GeV.edm4hep.root | tee output/log/out_20GeV.txt
        fi

fi

echo ""
echo "Done!"
