# Instructions to use for BNL HTCondor farm:

1. Login to eic node
2. Run this command with different parameters (it will use default if not specified)
```bash
MOMENTUM=3 ./run_submit.sh
```
3. You could modify DEFAULT parameters inside the script.
4. Example in loop:
``` bash
momenta=(1 2 5)
for mom in "${momenta[@]}"; do
    export MOMENTUM=$mom
    ./run.sh
done
```
5. This will set everything for you, including `eic-shell`, `epic` repository and `eicrecon`.
And run in batch mode full chain:
```bash
ddsim SIMFILE
eicrecon RECOFILE
root -l -n yourMacro.C (RECOFILE, ANAFILE)
```
6. There is a possibility to wait for jobs using another script `condor_control.sh` which monitors your batch jobs and continues when everything is finished( e.g. you could merge output root analized files which contain your own trees/histograms)

7. The script is a bit long (~400 lines) yet it does fully automated procedure in single file. 




Hi all, here is some basic info on nHcal/EPIC beginners

* * *
EPIC Analysis Workflow: Simulation → Reconstruction → Analysis
================================================================
 Overview
-----------

Three-step workflow for EPIC analysis:

1.  **Simulate**: Use `ddsim` or `npsim` to generate simulated events.
    
2.  **Reconstruct**: Process simulated data with `EICRecon`.
    
3.  **Analyze**: Run your custom **ROOT macros** for physics analysis.
    
* * *
 Getting Started
------------------

*    [EPIC Get Started](https://eic.github.io/documentation/getstarted.html)
    
*   [EIC Tutorials Overview](https://eic.github.io/documentation/tutorials.html)
    
*   [SSH Guide for STAR/EIC - password-less authorization](https://star-juniors.github.io/software/ssh/)

* * *

Step 1: **Simulation**
-------------------------

*   [eic/epic](https://github.com/eic/epic) — detector geometry, compact files, configurations
###  Tools: `ddsim` / `npsim` with Geant4

*   [Single Particle Simulation Tutorial](https://eic.github.io/tutorial-simulations-using-ddsim-and-geant4/01-single-particle-simulations/index.html)
    
*   [DD4hep Manual](https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManual.pdf)
    

### Detector Description Files for nHcal

From the `epic/` repository:

*   [`backward.xml`](https://github.com/eic/epic/blob/main/compact/hcal/backward.xml) — full geometry description
    
*   [`definitions.xml`](https://github.com/eic/epic/blob/main/compact/definitions.xml#L649-L654) — constants  which are used in `backward.xml`
    
*   [`PolyhedraEndcapCalorimeter2_geo.cpp`](https://github.com/eic/epic/blob/main/src/PolyhedraEndcapCalorimeter2_geo.cpp) — geometry implementation code
* * *

 Step 2: **Reconstruction**
-----------------------------

###  Tool: `EICRecon`

*    [EICRecon GitHub Repository](https://github.com/eic/EICrecon)
    

### Parameter Configuration

*   [`EICRecon/src/detectors/EHCAL/EHCAL.cc`](https://github.com/eic/EICrecon/blob/main/src/detectors/EHCAL/EHCAL.cc) — example of parameter tuning for reconstruction (clustering parameters)    

* * *

Step 3: **Analysis**
-----------------------

### Tool: Custom ROOT Macros

Once you have reconstructed output, analyze it using your own ROOT-based analysis scripts.

* * *

My Workflow at BNL SDCC ssh:
-----------------------

I use SSH connection via [`VS Code`](https://star-juniors.github.io/software/vscode.html) which setups the whole environment as you are working on your own laptop. See the picture.

- [SSH Guide for STAR/EIC - password-less authorization](https://star-juniors.github.io/software/ssh/)
- [Tunnel guide](https://star-juniors.github.io/software/vs-code-tunnel.html)
- And one can you use [`tmux`](https://pragmaticpineapple.com/gentle-guide-to-get-started-with-tmux/) for detaching process from SSH connection to be run in background while you are disconnected from SSH.

![image(1)](https://github.com/user-attachments/assets/d9a119ee-eff0-42b8-84f8-6411ae0c2b64)
