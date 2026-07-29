#!/bin/bash

submit=$1

if [ -z $submit ]; then
    echo "[i] Set it to submit = 0 for local test"
    submit=0
fi

##################################################
### Local test ###
##################################################
if [ $submit -eq 0 ]; then
    echo "[i] Running locally"
    root -b -q analysis.C\(\"test.list\",\"test.output.root\"\)
    exit
fi

###############################################################
### Batch production for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    config=DIS_10x100_minQ2_1_25.10.0
    
    echo "[i] Submit batch jobs for sample ${config}"
    pwd=$PWD
	
    odir=$pwd/output
    logdir=$odir/log
    if [ ! -d $odir ]; then
	mkdir -pv $odir
    fi
    rm -rf $odir/*
    mkdir $logdir

    executable=job_run.sh
    cp -v ${executable} $odir/.
    cp -v analysis.C $odir/.
    
    # Initialising Condor File
    condor_file=CondorFile_submit
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${odir}/${executable}" >> ${condor_file}
    echo "GetEnv  =  True" >> ${condor_file}
    echo "Arguments = \$(oodir) \$(inputfile) \$(outputfile) "  >> ${condor_file}

    echo "log = ${logdir}/log_\$(number).log"  >> ${condor_file}
    echo "error = ${logdir}/log_\$(number).err"  >> ${condor_file}
    echo "output = ${logdir}/log_\$(number).out"  >> ${condor_file}

    echo "" >> ${condor_file}
    echo "queue number, oodir, inputfile, outputfile from (" >> ${condor_file}

    # change this line as needed
    files=`ls $pwd/input_files/$config/subList*`
    for file in $files; do
	listNum=`basename ${file} | sed "s/.list//g" | cut -f 2 -d _`
	OutFile=${odir}/output_${listNum}.root
	echo ${listNum}, ${odir}, ${file}, ${OutFile} >> ${condor_file}
    done
    echo ")" >> ${condor_file}
    mv ${condor_file} $odir/.
    cd $odir

    #submit condor jobs
    condor_submit ${condor_file}
    cd $pwd
fi

