# Analysis code for studying track measurement, outlier, and hole hits

Generating single negative muon events
----------------------------------------
In the container, do the following:

```
mkdir output

source /opt/detector/epic-main/bin/thisepic.sh

npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "2.0*GeV" --gun.momentumMax "2.0*GeV" \
--numberOfEvents 10000 --gun.position 0,0,0 --outputFile output/output_2GeV.edm4hep.root

eicrecon -Ppodio:output_file=output/eicrecon_out_2GeV.root -Pjana:nevents=10000 -Pdd4hep:xml_files=epic_craterlake.xml output/output_2GeV.edm4hep.root
```

This will generate 10,000 single negative muon events with 2 GeV/c momentum, eta = [-4,4], and phi = [0,2Pi].

Running the analysis code
---------------------
The analysis code needs to be run in the container, since it uses PODIO classes. If you want to print out information for every event, set the ```print_evt_info``` variable to true. (N.B. There sometimes is an issue when using the CLING interpreter. So, the code should be run in compiled mode.)

```
mkdir plots
ln -sf output/output_2GeV.edm4hep.root hit_matching.input.root
source /opt/detector/epic-main/bin/thisepic.sh
root -l -b -q hit_matching.C+ 
```

If your `compactFile` happens to include so-called 2DStrip MPGDs, you may want to boost the processing of those by supplying a geometry file:

```
ln -sf detector_geometry.2DStrip.root hit_matching.geometry.root
```
where the geometry file can be obtained with:

```
dd_web_display -o detector_geometry.2DStrip.root --export $DETECTOR_PATH/epic_craterlake_tracking_2DStrip.xml
```

