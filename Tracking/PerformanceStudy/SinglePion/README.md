# Scripts to generate preTDR (sec 8.3.3) style plots for single particle tracking performance study
Shujie Li  at lbl dot gov, Oct 2025
* this only works for single particle (pi+ by default) simualtion. Only take first track per event to bypass the track-particle matching.
  
Steps:
0. generate rec rootfiles at a range of eta slices (e.g. -3.5, -2.5, -1, 1, 2.5, 3.5) and fixed momentum. Default output file name format:
   `rec_{setting}_eta_{eta_list[ii]:g}_{eta_list[ii+1]:g}_{mom}GeV_{nev}.root`
   One can also adjust this code to use the single particle simulation campaign file. 
1. run run_performance_study.py to pull out bin-by-bin efficiency and dp/p, theta, phi, DCA_r resolutions. By default, results are saved at `eff_out.txt` and `resol_out_whole.txt`, one row = one momentum and eta slice.
2. run run_performnace_plot.py to generate plots.
