#!/bin/bash
#script to run the machine learning differential in pT and y
# Shyam Kumar; INFN Bari, Italy;
# shyam.kumar@ba.infn.it; shyam055119@gmail.com
# Supported by The FAIR Spoke 6 Project, funded by the NextGenerationEU program in Italy
rm *.png
prepare_sample=true  
preselection_cuts='(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100.'
nEvents_D0Sample=1.747e+9 # Simulated ~1747M
nTotalEvents=4.72162e+9 #10fb^-1 (ep)
y_arr_ML=(-1.0 1.0 3.0) 
pt_arr_ML=(0.0 1.0 2.0)

if $prepare_sample; then
find . -type f \( -name "*.root" -o -name "*.png" -o -name "*.txt" \) -delete # delete all root, txt, png files in directory
cd Data_Preparation/Filtered_D0Sample/
root -b -l -q Create_Signal_with_Cuts.C'("'"$preselection_cuts"'")'
root -b -l -q Create_bkg_with_Cuts.C'("'"$preselection_cuts"'")'
cd ../Filtered_DISSample
root -b -l -q Create_Signal_with_Cuts.C'("'"$preselection_cuts"'")'
root -b -l -q Create_bkg_with_Cuts.C'("'"$preselection_cuts"'")'
cd ../../
# Signal and Backgrounds for ML
hadd Data_Preparation/SignalD0.root  Data_Preparation/Filtered_D0Sample/SignalD0.root Data_Preparation/Filtered_DISSample/SignalD0.root 
cp Data_Preparation/Filtered_DISSample/BkgD0.root Data_Preparation/
fi
# Prepare data using resampling and applying proper scaling factor
cd Data_Preparation/Merge_Data/
source Merged_Data.sh "$nEvents_D0Sample" "$nTotalEvents"
cd ../../
# RunML differential in y and pT
for ((i=0; i<${#y_arr_ML[@]}-1; i++)); do
for ((j=0; j<${#pt_arr_ML[@]}-1; j++)); do
rm -rf ML_Output_Optuna_${y_arr_ML[i]}_${y_arr_ML[i+1]}_${pt_arr_ML[j]}_${pt_arr_ML[j+1]} && mkdir ML_Output_Optuna_${y_arr_ML[i]}_${y_arr_ML[i+1]}_${pt_arr_ML[j]}_${pt_arr_ML[j+1]}
python3 machine_learning_Final.py --ymin ${y_arr_ML[i]} --ymax ${y_arr_ML[i+1]} --ptmin ${pt_arr_ML[j]} --ptmax ${pt_arr_ML[j+1]}
root -b -l -q draw_soverb_significance.C'('${y_arr_ML[i]}', '${y_arr_ML[i+1]}','${pt_arr_ML[j]}','${pt_arr_ML[j+1]}')'
root -b -l -q Superimpose_BDT_efficiencies.C'('${y_arr_ML[i]}', '${y_arr_ML[i+1]}','${pt_arr_ML[j]}','${pt_arr_ML[j+1]}')'
done
done

