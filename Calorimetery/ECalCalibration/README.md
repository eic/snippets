This directory contains scripts that calculate slope parameter for EcalBarrelScFiCalibration. Here are the steps:

1. Extract timeWalkCorrectionParameters and lightSpeedParameters from Minho's script. This script ONLY extract the slope parameter for EcalBarrelScFiCalibration and nothing else.
2. Start your eic-shell environment.
3. Install and source your own EICrecon with the correct timeWalkCorrectionParameters and lightSpeedParameters put into EcalBarrelScFiCalibration class in src/detectors/BEMC/BEMC.cc. At the time this README is created, some modifications necessary to run this script (e.g. https://github.com/eic/EICrecon/pull/2770, https://github.com/eic/EICrecon/pull/2789 and a new CalorimeterHitDigi class for which pull request is not created yet) have not been merged to main. You may skip this step if you are running it when all those changes are accepted into the main branch.
4. Execute 'snakemake -j <insert-number-of-CPUs>'. It generates 'plot/theta80_100/gamma/zRes0mm/UncalibratedERecoVsETruth_CALOROCCluster_Single_energy10_20000.png' as output.
5. Open that png and read the equation of linear fit on the legend on the subplot on the left. The inverse of the slope parameter on the legend is the slope for EcalBarrelScFiCalibration. Note that you need to take the inverse of what you are reading from the png. 
