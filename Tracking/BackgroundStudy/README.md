# Example scripts for tracking performance study with beam background
* See https://indico.bnl.gov/event/29753/
* contact: Shujie Li, shujieli@lbl.gov
### Note to general physics analysis:
only use tracks with at least 4 hits, e.g.  __CentralCKFTrajectories.nMeaurements>3__ 

## run_merger.sh: 
This is an example script to generate small sample of background mixed hepmc file locally.

Pleae always use background files from the simulation campaign data (https://eic.github.io/epic-prod/) when possible 

### pre-requists:
   1.  install SignalBackgroundMerger from https://github.com/eic/HEPMC_Merger
   2.  recommend to run within eic-shell to use the xrdfs command

 run as: ./run_merger.sh, then follow the options.

* Mixing frequencies can be found on [ePIC background Wiki](https://wiki.bnl.gov/EPIC/index.php?title=Background)
* Explanations on [mixing scheme](https://github.com/eic/eic.github.io/blob/master/_resources/background_mixed_samples.md)

* For SR and electron beam gas, please also read README files under JLab ifarm (contact: Andrii Natochii):
    SynRad: /volatile/eic/andrii/SynradG4_HepMC_Files_SR_on_IP6
    ESR beam losses: /volatile/eic/andrii/Xsuite_HepMC_Files_ESR_BeamLoss_on_ePIC

* 1000 merged events at 18GeV takes about 2 hours for dd4hep+eicrecon. The 5 and 10 GeV mixed events will be 10x slower due to higher SR freq. Use the simulation campaign file if you can. 

* You need then run npsim with the generated hepmc file and corresponding flags.


## epic_analysis-ak_example.ipynb
example python script to pull out background particle distributions, hit rate per detector area, and track efficiency and purity etc.  It uses uproot (>2.7), seaborn, particle, and a few other python modules. You may need to pip install them. 

Please navigate the jupyter notebook with section titles. One should be able to convert relevant code to plain python script, or ROOT/C easily with chatgpt. 

* The ___inspect...___ section does not depend on podio
   1. plot mcparticle generator status (signal or backgrounds)
   2. plot number of detector hit on a given Silicon tracker surface per ms per unit area, see [this talk](https://indico.bnl.gov/event/29602/contributions/112847/attachments/64502/110757/background_rate_09022025.pdf)
* to run the ___hit-based analysis___ section:
   1. in eic-shell, type jupyer-lab to start the server (may take a minute)
   2. copy the localhost link to your browser or editor to open the notebook
   3. select the kernel provided by eic-shell to use podio modules

   key functions: get_traj_hits, get_part_hits. Most interesting plots are under subsection "efficiency" and "track purity"

    

