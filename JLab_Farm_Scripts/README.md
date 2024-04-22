# JLab Farm Scripts
A small collection of scripts for creating and running jobs on the JLab farm.

Remember, verify that your script runs and works on an interactive node before submitting a large set of jobs. Check that your output is as expected.
*Please* read the comments in the scripts!

## JLab_Farming.sh
This script creates and submits jobs to a new swif2 workflow. As is, it takes 4 arguments.
1. The number of files to run.
2. The number of events to run per file.
3. A blank argument 3, defaults to 10.
4. A blank argument 4, defaults to 10.

The intention with args 3/4 is to demonstrate how to add and utilise further arguments you may need in your script.
Check the paths on lines 10, 58 and 80 carefully. Modify as needed.
Insert the script you want to run on line 103.

## JLab_Farming_Job.sh
This script is the job actually submitted to the farm to run. The arguments are the same as JLab_Farming.sh.
Make sure the paths match between the two files.
As above, check file paths carefully. Test a small sample first!
Note that this script calls Init_Env.sh as the first command after entering eic-shell.

## Init_Env.sh

This simply sources a series of things you will need to run your simulation/reconstruction.
Change as needed for your desired configuration.
