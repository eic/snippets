#!/usr/bin/env bash
# Convert .FCStd file to directory of .stl files
source FCStd_to_stl.sh L4.FCStd
# Remove .stl files that are not needed (silicon, or components outside envelope)
rm L4/2025*
rm L4/a*
rm L4/A*
rm L4/FIB*
rm L4/LEC*
rm L4/MT*
rm L4/Part*
rm L4/REC*
rm L4/U*
# Remove parts of components which extrude L4envelope in z direction
./cut_all_stl.sh L4 L4cut L4envelope.stl
# Merge overlapping carbon fibre components (offset skins by 10um so new mesh is watertight)
python merge_stl.py L4cut/Carbon_Compsite_Central_Spar_L006.stl L4cut/Carbon_Compsite_Central_Spar_L007.stl L4cut/Carbon_Composite.stl
python translate_stl.py L4cut/Carbon_Composite_Stave_Skin_Top_-_Overlap001.stl L4cut/Carbon_Composite_Stave_Skin_Top_-_Overlap001.stl --dy -0.01
python merge_stl.py L4cut/Carbon_Composite.stl L4cut/Carbon_Composite_Stave_Skin_Top_-_Overlap001.stl L4cut/Carbon_Composite.stl
python translate_stl.py L4cut/Carbon_Composite_Stave_Skin_Bottom_-_Overlap001.stl L4cut/Carbon_Composite_Stave_Skin_Bottom_-_Overlap001.stl --dy 0.01
python merge_stl.py L4cut/Carbon_Composite.stl L4cut/Carbon_Composite_Stave_Skin_Bottom_-_Overlap001.stl L4cut/Carbon_Composite.stl 
rm L4cut/Carbon_Composite_Stave_Skin_Bottom_-_Overlap001.stl  
rm L4cut/Carbon_Compsite_Central_Spar_L007.stl
rm L4cut/Carbon_Composite_Stave_Skin_Top_-_Overlap001.stl     
rm L4cut/Carbon_Compsite_Central_Spar_L006.stl
# Scale K9 components by 0.97 to avoid overlaps with carbon fibre
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant012.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant012.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant013.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant013.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant014.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant014.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant015.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant015.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant016.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant016.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant017.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant017.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant018.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant018.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant019.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant019.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant020.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant020.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant021.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant021.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant022.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant022.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_C-Channel_Variant023.stl L4cut/K9_Split_CrossBrace_C-Channel_Variant023.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_End_Varient002.stl L4cut/K9_Split_CrossBrace_End_Varient002.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_End_Varient003.stl L4cut/K9_Split_CrossBrace_End_Varient003.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_Sensor_4-007.stl L4cut/K9_Split_CrossBrace_Sensor_4-007.stl --scale 0.97
python scale_stl.py L4cut/K9_Split_CrossBrace_Sensor_4-008.stl L4cut/K9_Split_CrossBrace_Sensor_4-008.stl --scale 0.97
python scale_stl.py L4cut/MirrorK9_Split_CrossBrace_End_Varient002.stl L4cut/MirrorK9_Split_CrossBrace_End_Varient002.stl --scale 0.97
python scale_stl.py L4cut/MirrorK9_Split_CrossBrace_End_Varient003.stl L4cut/MirrorK9_Split_CrossBrace_End_Varient003.stl --scale 0.97
# Convert .stl files to ASCII and then to .gdml convertL4.sh contains offset parameters to set origin to centre of stave
./convertL4.sh L4cut
# Delete .stl files and rename directory of .gdml file to L4
rm L4cut/*.stl
rm -rf L4
mv L4cut L4
