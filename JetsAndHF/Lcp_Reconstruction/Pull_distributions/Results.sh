#!/bin/bash
# Script; shyam kumar; INFN Bari, Italy
# shyam.kumar@ba.infn.it; shyam055119@gmail.com
rm *.png
root -b -l -q compare_Pulls.C'("hPullVtxX","Pull (PVx)")'
root -b -l -q compare_Pulls.C'("hPullVtxY","Pull (PVy)")'
root -b -l -q compare_Pulls.C'("hPullVtxZ","Pull (PVz)")'
root -b -l -q compare_Pulls.C'("hPullVtxX_1","Pull (PVx)[Recotrack=1]")'
root -b -l -q compare_Pulls.C'("hPullVtxY_1","Pull (PVy)[Recotrack=1]")'
root -b -l -q compare_Pulls.C'("hPullVtxZ_1","Pull (PVz)[Recotrack=1]")'
root -b -l -q compare_Pulls.C'("hPullVtxX_2","Pull (PVx)[Recotrack>1]")'
root -b -l -q compare_Pulls.C'("hPullVtxY_2","Pull (PVy)[Recotrack>1]")'
root -b -l -q compare_Pulls.C'("hPullVtxZ_2","Pull (PVz)[Recotrack>1]")'
root -b -l -q compare_Pulls.C'("hResVtxX","PVx_{Rec}-PVx_{Mc}(mm) [Recotrack>1]")'
root -b -l -q compare_Pulls.C'("hResVtxY","PVy_{Rec}-PVy_{Mc}(mm) [Recotrack>1]")'
root -b -l -q compare_Pulls.C'("hResVtxZ","PVz_{Rec}-PVz_{Mc}(mm) [Recotrack>1]")'



