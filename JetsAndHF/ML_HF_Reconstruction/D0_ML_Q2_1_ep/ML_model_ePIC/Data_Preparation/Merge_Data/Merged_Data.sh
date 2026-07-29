#!/bin/bash
#script to run the machine learning differential in pT and y
# Shyam Kumar; INFN Bari, Italy
rm *.png *.root
rm -rf Output/ && mkdir Output/
nEvents_D0Sample=$1 # Simulated ~1747M
nTotalEvents=$2 #10fb^-1 (ep)
cd Merge_Data/
y_arr=(-3.0 -1.0 1.0 3.0) 
pt_arr=(0.0 1.0 2.0 5.0 10.0)
for ((i=0; i<${#y_arr[@]}-1; i++)); do
for ((j=0; j<${#pt_arr[@]}-1; j++)); do
rm final_merged_${y_arr[i]}_${y_arr[i+1]}_${pt_arr[j]}_${pt_arr[j+1]}.root
root -b -l -q Sample_signal.C'('${y_arr[i]}','${y_arr[i+1]}','${pt_arr[j]}','${pt_arr[j+1]}','$nEvents_D0Sample','$nTotalEvents')'
root -b -l -q Sample_bkg.C'('${y_arr[i]}','${y_arr[i+1]}','${pt_arr[j]}','${pt_arr[j+1]}','$nTotalEvents')'
hadd final_merged_${y_arr[i]}_${y_arr[i+1]}_${pt_arr[j]}_${pt_arr[j+1]}.root  merged_bkg_y_D0_${y_arr[i]}_${y_arr[i+1]}_pt_D0_${pt_arr[j]}_${pt_arr[j+1]}.root merged_signal_y_D0_${y_arr[i]}_${y_arr[i+1]}_pt_D0_${pt_arr[j]}_${pt_arr[j+1]}.root
done
done
root -b -l -q fitD0withStudentT.C



