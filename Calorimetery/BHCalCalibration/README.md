# BHCalCalibration

Here you'll find code associated with the ML-assissted BHCal Calibration workflow.

### FillBHCalCalibrationTuple Input

Since this is an eicrecon plugin, it runs on `.edm4hep.root` files.  However, it is currently designed
to only work on single particle events.

### FillBHCalCalibrationTuple Usage

After compiling `EICrecon`, create and compile the plugin with:

```
eicmkplugin.py FillBHCalCalibrationTuple
cp <path to this repo>/plugin/FillBHCalCalibrationTupleProcessor.* ./FillBHCalCalibrationTuple/
cmake -S JCalibrateHcal -B FillBHCalCalibrationTuple/build
cmake --build FillBHCalCalibrationTuple/build --target install
```

Then run it with:

```
eicrecon -Pplugins=FillBHCalCalibrationTuple <input edm4hep file>
```
