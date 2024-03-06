#!/bin/bash
# script for material mapping with ACTS python bindings
# please run it in eic-shell, whith DETECTOR_PATH and ACTS_PATH set up 
# run as : ./run_matmap_py.sh 0 --nev 1000 --tag _test
#           use option 0,1,2 step by step to generate the material map (or option -1 to run it at once). option 4 and 5 to produce comparison plots.
# Shujie Li, 03. 2024


echo "=========FOR EPIC Craterlake, 03.2024========="
source detector_setup.sh
ACTS_PATH="" #path to the ACTS source code. see https://github.com/acts-project/acts/tree/main
XML_NAME="epic_craterlake_matmap.xml" #the detector configuration contains all detectors and materials within the volume you want to create the map for
XML_PATH=${PWD}
XML_FILE=${XML_PATH}/${XML_NAME}

nev=1000  #1000
nparticles=1000 #5000 #5000
tag=""
kopt=$1
while [[ $# -gt 1 ]]
do
  key="$2"

  case $key in
    --nev)
      nev=$3
      shift # past value
      shift
      ;;
    --nparticles)
      nparticles=$3
      shift # past value
      shift
      ;;
    --tag)
			tag=$3
			shift
			shift
			;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $2"
      print_the_help
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

recordingFile=geant4_material_tracks${tag}.root
geoFile=geometry-map${tag}.json
matFile=material-map${tag}.json
trackFile=material-map${tag}_tracks.root
propFile=propagation-material${tag}.root

	# ----GEANTINO SCAN------
if  [ "$kopt" == 0 ] || [ "$kopt" -lt 0 ]; then
	# output geant4_material_tracks.root
	# The result of the geantino scan will be a root file containing material tracks. Those contain the direction and production vertex of the geantino, the total material accumulated and all the interaction points in the detector.

python material_recording_ePIC.py -i ${XML_FILE} -n ${nev} -t ${nparticles} -o ${recordingFile}

fi

	#-----MAPPING Configuration-----

if  [ "$kopt" == 1 ]|| [ "$kopt" -lt 0 ]; then


	# map gemoery to geometry-map.json
	python geometry_ePIC.py -i ${XML_FILE} -o ${geoFile}

	# take geometry-map.json and read out config-map.json
	python3 ${ACTS_PATH}/Examples/Scripts/MaterialMapping/writeMapConfig.py ${geoFile} config-map${tag}.json

	# turn on approaches and beampipe surfaces for material mapping
	# you can always manually adjust the mapmaterial flag and binnings in config-map.json
	python3 materialmap_config.py -i config-map${tag}.json -o config-map_new${tag}.json

fi

if  [ "$kopt" == 2 ]|| [ "$kopt" -lt 0 ]; then
	# turn config-map.json into modified geometry-map.json
	python3 ${ACTS_PATH}/Examples/Scripts/MaterialMapping/configureMap.py ${geoFile} config-map_new${tag}.json


	#----MAPPING------------
	# input: geant4_material_tracks.root, geometry-map.json
	# output: material-maps.json or cbor. This is the material map that you want to provide to EICrecon, i.e.  -Pacts:MaterialMap=XXX  .Please --matFile to specify the name and type 
	#         material-maps_tracks.root(recorded steps from geantino, for validation purpose)

	python material_mapping_ePIC.py --xmlFile $XML_FILE --stepFile ${recordingFile} --geoFile ${geoFile} --matFile ${matFile}

fi

## steps below are for validation/debugging
if  [ "$kopt" == 3 ]|| [ "$kopt" -lt 0 ]; then
	#----Prepare validation rootfile--------
	# output propagation-material.root
	python material_validation_ePIC.py --xmlFile $XML_FILE --outputName ${propFile} --matFile ${matFile} -n ${nev} -t ${nparticles}
	# you may see a couple of the error msg below. Should not matter as long as it's not for every event.
	# PropagationA   ERROR     Propagation reached the step count limit of 1000 (did 1000 steps)
	#  PropagationA   ERROR     Step failed with EigenStepperError:3: Step size adjustment exceeds maximum trials

fi

## -------Comparison plots---------
if  [ "$kopt" == 4 ]; then
	rm -rf Validation
	mkdir Validation
	root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map.C'("'$propFile'","'$trackFile'","Validation")'
	root -l -b -q mat_map_local.C'("'$propFile'","'$trackFile'","Validation")'	# 

# fi
# if  [ "$kopt" == 5 ]; then
rm -rf Surfaces
mkdir Surfaces
cd Surfaces
mkdir prop_plot
mkdir map_plot
mkdir ratio_plot
mkdir dist_plot
mkdir 1D_plot
cd ..
root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_ratio.C'("'$propFile'","'$trackFile'",-1,"Surfaces/ratio_plot","Surfaces/prop_plot","Surfaces/map_plot")'
root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_dist.C'("'$trackFile'",-1,"Surfaces/dist_plot")'
root -l -b -q ${ACTS_PATH}/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_1D.C'("'$trackFile'",-1,"Surfaces/1D_plot")'
fi
