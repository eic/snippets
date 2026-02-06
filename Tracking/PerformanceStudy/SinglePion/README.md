# Scripts to generate preTDR (sec 8.3.3) style plots for single particle tracking performance study
Shujie Li  at lbl dot gov, Oct 2025
* this only works for single particle (pi+ by default) simualtion. Only take first track per event to bypass the track-particle matching.

## standard performance study
Folder __7etabins__: a complete workflow to generate simulation files, calculate tracking efficiency and resolutions bin by bin, and plot them v.s. momentum.

## single file quality inspection
__plot_single_resol.py__: for a given local / simulation campaign file, plot pull distributions, resolutions, efficiency.

__plot_single_measurements.py__: for a given local / simulation campaign file,
* plot the number of simhits/measurements/outliers per detector (the example shows only barrel trackers, endcap plots can be produced in a similar fashion).
* check the chi2, and residual distributions for all good measurements. 

__generate_install_script.py__: helper script to install missing packages to run a given python script. 