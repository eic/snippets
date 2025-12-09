# Helix

contact person: marr@bnl.gov

This example utilizes helix swimming method to find two-track DCA
which can be used for reconstructing D0 decay vertices. 

Usage: 
```
eic-shell
make
./analysis
```

If you want to use the code locally, please make sure you have ROOT
installed. In particular, if your ROOT is compiled with a recent c++ verion, you
might need to replace the following two lines in StHelix.h
```
  if (!::finite(mDipAngle    )) 	return   11;
  if (!::finite(mCurvature   )) 	return   12; 
```
with
```
  if (!std::isfinite(mDipAngle    )) 	return   11;
  if (!std::isfinite(mCurvature   )) 	return   12;
```
