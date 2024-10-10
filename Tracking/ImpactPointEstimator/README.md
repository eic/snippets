# Simple analysis code using ImpactPointEstimator

Generating events from (x,y,z) = (1,0,0)mm
-------------------------------------------
In the container, do the following:

```
mkdir output
mkdir output/log

source /opt/detector/epic-main/bin/thisepic.sh

npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --enableGun --gun.distribution 'eta' \
--gun.thetaMax 3.106 --gun.thetaMin 0.036 --gun.momentumMin "0.5*GeV" --gun.momentumMax "20*GeV" \
--numberOfEvents 10000 --gun.position 1,0,0 --outputFile output/output_1_0_0.edm4hep.root

eicrecon -Ppodio:output_file=output/eicrecon_out_1_0_0.root -Pjana:nevents=10000 -Pdd4hep:xml_files=epic_craterlake.xml output/output_1_0_0.edm4hep.root | tee output/log/out_1_0_0.txt
```

Running the analysis code
---------------------
The analysis code needs to be run in the container.

```
mkdir plots
source /opt/detector/epic-main/bin/thisepic.sh
root -l -b -q pca_global_impactpoint.C
```

