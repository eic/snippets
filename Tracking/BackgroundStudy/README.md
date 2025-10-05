# Example scripts for tracking performance study with beam background
* See https://indico.bnl.gov/event/29753/
* contact: Shujie Li, shujieli@lbl.gov

## run_merger.sh: 
example script to generate small sample of background mixed hepmc file locally

you should always use background files from the simulation campaign data (https://eic.github.io/epic-prod/) when possible 

### pre-requists:
   1.  install SignalBackgroundMerger from https://github.com/eic/HEPMC_Merger
   2.  recommend to run within eic-shell to use the xrdfs command

 run as: ./run_merger.sh, then follow the options.
 
* For SR and electron beam gas, please also read README files under JLab ifarm (contact: Andrii Natochii):
    SynRad: /volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6
    ESR beam losses: /volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC

* 1000 merged events at 18GeV takes about 2 hours for dd4hep+eicrecon. The 5 and 10 GeV mixed events will be 10x slower due to higher SR freq. Use the simulation campaign file if you can. 

* You need then run npsim with the generated hepmc file and corresponding flags.


## epic_analysis-ak_example.ipynb
example python script to pull out background particle distributions, track efficiency and purity etc.  It uses uproot (>2.7), seaborn, particle, and a few other python modules. You may need to pip install them. 

Hopefully one can convert desired functions to plain python script, of ROOT/C easily with chatgpt. 

* please run this within eic-shell if you want to use any podio modules:
    1. in eic-shell, type jupyer-lab to start the server (may take a minute)
    2. copy the localhost link to your browser or editor to open the notebook.

* key functions: get_traj_hits, get_part_hits. Most interesting plots are under subsection "efficiency" and "track purity"

* suggestion to general physics analysis: only use tracks with at least 4 hits, e.g.  __CentralCKFTrajectories.nMeaurements>3__ 
    

