This respository hosts scripts and macros for illustrating how to submit batch jobs on RCF through Condor.

Here are the steps
* Log into your EIC account on RCF
* Clone this respository
* Optional: prepare file lists for the simulation campaign you want to analyze
  * Use ```input_files/prep_file_list.sh```
  * Remember to modify the script to fetch the simulation files you want
* Submit jobs
  * This needs to be done outside of eic-shell
  * Replace the path to eic-shell to yours in ```job_run.sh``` 
  * To run a local test, use ```./submit_jobs.sh 0```
  * To submit batch jobs, use ```./submit_jobs.sh 1```
     * It creates a file called "CondorFile_submit", and ```condor_submit CondorFile_submit``` submits the jobs, which execute ```job_run.sh``` to run ```analysis.C```
     * If you run through you own file lists, remember to make changes accordingly in ```submit_jobs.sh```
