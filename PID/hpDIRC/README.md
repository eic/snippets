# Instructions for ePIC hpDIRC simulation and reconstruction

Please follow the [tutorial](https://eic.github.io/tutorial-setting-up-environment/02-eic-shell/index.html) on setting up the EIC environment first.
The following assumes you're running in ```eic-shell```. 

## Detector geometry

The hpDIRC detector geometry is based on [DD4hep](https://github.com/AIDASoft/DD4hep). The following repository is a fork of [eic/epic](https://github.com/eic/epic) and has a branch that contains a version of hpDIRC geometry and all the material properties required for Cherenkov photon propagation. Get a copy of it using:
```bash
git clone https://github.com/niwgit/epic.git
cd epic
git checkout box_envelope
```
## Compilation
The following commands build and install the geometry to the ```install``` directory (can be replaced by any directory path).
```bash
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build
cmake --install build
```
Load the geometry using 
```bash
source install/setup.sh
```
## Modifying detector geometry parameters
The source code of hpDIRC geometry is located in ```src/DIRC_geo.cpp```. The parameters used in that file are defined in the compact file ```compact/pid/dirc.xml```. The wavelength-dependent properties such as refractive index and absorption length of the materials are defined in ```compact/optical_materials.xml```.

Use the following command to see a list of compact files with detector geometry configurations. The file ```epic_dirc_only.xml``` is used to include only the hpDIRC detector in the simulation.
```bash
ls $DETECTOR_PATH
```
```${DETECTOR_PATH}/compact/pid/dirc.xml``` can be used to modify the parameters defining the hpDIRC geometry. For example, change the parameter ```DIRCBox_count``` to 1 to use only 1 bar box (default is 12 bar boxes).

## Viewing the geometry
```bash
dd_web_display --export $DETECTOR_PATH/epic_dirc_only.xml
```
The above command creates an output ROOT file called ```detector_geometry.root```. Use [jsroot](https://root.cern/js/) to open the ROOT file and view the geometry.

## Simulation
Geant4 simulations with DD4hep are done via ```npsim```. For simulating events with one charged particle per event as the primary particle, we use a particle gun with ``-G`` flag. The number of events is specified with the ```-N``` flag. The paticle gun has options to set the momentum, vertex position and polar and azimuthal angles of the charged particle. 

The following example is for running 100 events with 6 GeV/c pi+ track of polar angle 30 deg and azimuthal angle 355.334 deg. The azimuthal angle was selected such that the charged particle hits the center of one of the middle bars in the bar box when the magnetic field is included. 
```bash
npsim --runType batch --compactFile $DETECTOR_PATH/epic_dirc_only.xml -G -N 100  --gun.particle "pi+" --gun.momentumMin 6*GeV --gun.momentumMax 6*GeV --gun.phiMin 355.334*deg --gun.phiMax 355.334*deg --gun.thetaMin 30*deg --gun.thetaMax 30*deg --gun.distribution uniform --gun.position 0*cm,0*cm,0*cm --part.userParticleHandler='' --outputFile sim.edm4hep.root
```
This will create an output ROOT file called ```sim.edm4hep.root```. The file name can be modified but it should contain the extension ```.edm4hep.root```. This file will contain all the hits on the DIRC sensitive detector (MCP plane) as ```DIRCBarHits```.  These are ```SimTrackerHit``` objects that are defined in the [EDM4hep](https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml#L227) data model. The truth information of each particle simulated is saved in [MCParticle](https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml#L158) objects. The primary particle can be selected using ```MCParticles.generatorStatus==1```.

To see the event display use ```--runType vis``` with a visualization macro file (e.g. ```--macroFile myvis.mac```) and ```--enableGun``` or ```-G``` option. Then in the Geant4 interactive mode run the event with ```/run/beamOn 1```.
![alt text](https://github.com/eic/snippets/blob/main/PID/hpDIRC/pic/event_display_myvis_box_envelope.png)

## Applying quantum efficiency (QE) of sensors
Currently, the ```npsim``` output file contains information on all the optical photons that hit the MCP plane. So the number of hits per event (track) can be high as 800-1000. But in reality, the detected number of photons will be much lower (50-150) due to the photon detection efficiency of the sensors. This efficiency is wavelength-dependent and ideally should be applied during the simulation to avoid tracking optical photons that are not detected. But at the moment these efficiencies are applied with the algorithm [PhotoMultiplierHitDigi](https://github.com/eic/EICrecon/blob/main/src/algorithms/digi/PhotoMultiplierHitDigi.cc) of ```eicrecon```. The configuration parameters of this algorithm for hpDIRC is set using https://github.com/eic/EICrecon/blob/main/src/detectors/DIRC/DIRC.cc and the quantum efficiencies for sensors are included in ```eicrecon``` with the default plugin ```DIRC```. This plugin creates collections ```DIRCRawHits``` and ```DIRCRawHitsAssociations```.

## Using plugins with eicrecon

 First you'll need to clone the repository EICrecon and follow the build instructions. 
```bash
git clone https://github.com/eic/EICrecon.git
cd EICrecon
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build
cmake --install build
```

After that source the EICrecon environment.
```bash
source install/bin/eicrecon-this.sh
```

Instructions on creating a plugin to make custom histograms/trees are given in the tutorial https://eic.github.io/tutorial-jana2/03-end-user-plugin/index.html. As an example,
```bash
eicmkplugin.py  hpDIRCrawHits
```

For the hpDIRC, we use two plugins (https://github.com/eic/snippets/tree/main/PID/hpDIRC/plugins) with ```eicrecon```. 
- ```hpDIRCsimHits``` : This plugin is used to create histograms using simulated hits before applying sensor QE. For example, we can look at the wavelength distribution of Cherenkov photons.
- ```hpDIRCrawHits``` : This plugin uses hits after applying QE. Additionally, it removes the hits that are in gaps between MCPs. The main purpose of this plugin is to save information such as the number of hits per event, time, mcp id and pixel id of each hit in a TTree called ```dirctree```.

The following shows how to build the plugin ```hpDIRCrawHits``` that will be installed in the directory ```EICrecon_MY```. The ```$EICrecon_MY$``` environment variable should be set to that directory so that when ```eicrecon``` runs the plugin will be identified. This assumes the files ```hpDIRCrawHitsProcessor.h```, ```hpDIRCrawHitsProcessor.cc``` and ```CMakelists.txt``` are within the source directory ```hpDIRCrawHits```.
```bash
mkdir EICrecon_MY
export EICrecon_MY=${PWD}/EICrecon_MY
cmake -S hpDIRCrawHits -B hpDIRCrawHits/build
cmake --build hpDIRCrawHits/build --target install
```
To run the plugins with ```eicrecon``` over the simulation file ```sim.edm4hep.root``` which was created earlier using ```npsim``` do:
```bash
eicrecon -Pplugins=hpDIRCsimHits,hpDIRCrawHits -Ppodio:output_collections=DIRCRawHits,DIRCRawHitsAssociations sim.edm4hep.root
```
This will create two output files ```eicrecon.root``` and ```podio_output.root```. The file ```eicrecon.root``` should contain two directories ```hpDIRCsimHits``` and ```hpDIRCrawHits``` corresponding to the plugins we used. The histograms and trees defined in the plugins can be found in these directories. The ```dirctree``` in the ```hpDIRCrawHits``` directory will be used to obtain information required for the time imaging reconstruction.
