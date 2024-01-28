## Time imaging reconstruction
The following can be done either within or outside ```eic-shell```.
For time imaging reconstruction we need an input file with the same number of events for two different particle types. Let's consider the case of pi+ and K+. 
- Use ```npsim``` to generate 25000 events separately for each of the two cases ("pi+" and "kaon+" for ```--gun.particle```). Rename the ```sim.edm4hep.root``` file to match the particle type.
- Run ```eicrecon``` with the two plugins ```hpDIRCsimHits``` and ```hpDIRCrawHits``` over each of the simulated file. Rename the ```eicrecon.root``` file to match the particle type.
- Merge the two ```eicrecon``` root files. For example,
  ```bash
  hadd eicrecon_30_theta_mix_pik_50000_events.root EICReconOut_pi+_25000evts_theta_30deg.root EICReconOut_kaon+_25000evts_theta_30deg.root
  ```
- To create PDFs(Probability Density Functions) run the macro ```createPdf_epic.cpp``` on the input file (merged ```eicrecon``` file). This will create a ROOT file with extension ```.pdf.root``` and will contain PDFs for time imaging reconstruction.
  ```
  root -q -b 'createPdf_epic.cpp("eicrecon_30_theta_mix_pik_50000_events.root")'
  ```
- For time imaging reconstruction use the input file and the PDF file as arguments and run the macro ```recoPdf_epic.cpp```. 
```
root -q -b 'recoPdf_epic.cpp("eicrecon_30_theta_mix_pik_50000_events.root","eicrecon_30_theta_mix_pik_50000_events.pdf.root")'
```
- This will create a ROOT file with extension ```.root_r.root```. This file contains a TTree called ```reco``` in which the photon yield (```nph```) and separation power (```sep```) are saved and can be used to evaluate the hpDIRC performance.
- To plot the hit pattern run the macro ```draw_hp.C``` with the input file as the argument.
```
root -q -b 'draw_hp.C("eicrecon_30_theta_mix_pik_50000_events.root")'
```
