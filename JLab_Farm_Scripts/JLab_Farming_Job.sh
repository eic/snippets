#! /bin/bash

### Stephen Kay, University of York
### 22/04/24
### stephen.kay@york.ac.uk
### A script to process some commands within eic-shell. In this example, the simulation is run and then reconstructed.
### Steps to generate files from an event generator can be included in the relevant section.
### This script can be run on its own, or fed into a job on the farm
### VERIFY that it runs interactively for a small sample before submitting a large job!

SimDir="/group/eic/users/${USER}/ePIC/" # Put in the path of your directory here (where your eic-shell is)
echo "Running as ${USER}"
echo "Assuming simulation directory - ${SimDir}"
if [ ! -d $SimDir ]; then   
    echo "!!! WARNING !!!"
    echo "!!! $SimDir - Does not exist - Double check pathing and try again !!!"
    echo "!!! WARNNING !!!"
    exit 1
fi

FileNum=$1 # First arg is the number of files to run
if [[ -z "$1" ]]; then
    echo "I need a number of files to run!"
    echo "Please provide a number of files to run as the first argument"
    exit 2
fi
NumEvents=$2 # Second argument is an output file name
if [[ -z "$2" ]]; then
    echo "I need a number of events to generate per file!"
    echo "Please provide a number of event to generate per file as the second argument"
    exit 3
fi

# Check if an argument was provided for Arg3, if not, set 10
if [[ -z "$3" ]]; then
    Arg3=10
    echo "Arg3 not specified, defaulting to 10"
else
    Arg=$3
fi

# Check if an argument was provided for Arg4, if not, set 10
if [[ -z "$4" ]]; then
    Arg4=10
    echo "Arg4 not specified, defaulting to 10"
else
    Arg4=$4
fi

# Change output path as desired
OutputPath="/volatile/eic/${USER}" # Change as needed, should match JLab_Farming.sh!
export Output_tmp="$OutputPath/im_${FileNum}_${NumEvents}_${Arg3/./p}_${Arg4/./p}"
if [ ! -d "${Output_tmp}" ]; then # Add this in this script too so it can be run interactively
    mkdir $Output_tmp
else
    if [ "$(ls -A $Output_tmp)" ]; then # If directory is NOT empty, prompt a warning
	if [[ "${HOSTNAME}" == *"ifarm"* ]]; then # Only print this warning if running interactively
	    echo "!!!!! Warning, ${Output_tmp} directory exists and is not empty! Files may be overwritten! !!!!!"
	fi
    fi
fi

export EICSHELL=${SimDir}/eic-shell # Point to wherever your eic-shell is!
cd ${SimDir}
# Run EIC shell, generate the events, afterburn them, run the simulation, reconstruct the events
if test -f "${Output_tmp}/ddsimOut_${FileNum}_${NumEvents}.edm4hep.root"; then # If simulatuon file exists, only run the reconstruction
    # Should add additional checks on this, check how old files is (if before 2024, run anyway, if small, re-rerun etc)
    echo "${Output_tmp}/ddsimOut_${FileNum}_${NumEvents}.edm4hep.root already exists, running reconstruction only"
# Init_Env.sh is a script to source various things when in eic-shell
cat <<EOF | $EICSHELL
source Init_Env.sh
cd $Output_tmp
eicrecon -Pjana:nevents=${NumEvents} -Phistsfile=EICReconOut_${FileNum}_${NumEvents}.root ${Output_tmp}/ddsimOut_${FileNum}_${NumEvents}.edm4hep.root 
sleep 5

echo; echo; echo "Reconstruction finished, output file is - ${Output_tmp}/EICReconOut_${FileNum}_${NumEvents}.root"; echo; echo;
EOF
# If simulation file doen't exist, generate files, run simulation and then reconstruction
else
cat <<EOF | $EICSHELL
source Init_Env.sh

# Add any event generation steps here

# Add any input files in as needed
npsim -v 4 --inputFiles PATH_TO_INPUT_FILE --outputFile ${Output_tmp}/ddsimOut_${FileNum}_${NumEvents}.edm4hep.root --compactFile ${SimDir}/epic/epic_ip6.xml -N ${NumEvents}
sleep 5

# Run reconstruction, add any plugins you want to use
echo; echo; echo "Simulation finished, running reconstruction."; echo; echo;
cd $Output_tmp
eicrecon -Pplugins=PLUGINS -Pjana:nevents=${NumEvents} -Phistsfile=EICReconOut_${FileNum}_${NumEvents}.root ${Output_tmp}/ddsimOut_${FileNum}_${NumEvents}.edm4hep.root 
sleep 5

echo; echo; echo "Reconstruction finished, output file is - ${Output_tmp}/EICReconOut_${FileNum}_${NumEvents}.root"; echo; echo;

EOF
fi

exit 0
