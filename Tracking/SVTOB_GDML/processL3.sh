#!/usr/bin/env bash
# Convert .FCStd file to directory of .stl files
source FCStd_to_stl.sh L3.FCStd
# Remove .stl files that are not needed (silicon, or components outside envelope)
rm L3/A*
rm L3/B*
rm L3/CAP*
rm L3/Clip*
rm L3/DF*
rm L3/LEC*
rm L3/MO*
rm L3/MT*
rm L3/Opto*
rm L3/P*
rm L3/R*
rm L3/T*
rm L3/U*
# Remove parts of components which extrude L3envelope in z direction
./cut_all_stl.sh L3 L3cut L3envelope.stl
# Merge overlapping carbon fibre components (offset skins by 10um so new mesh is watertight)
python merge_stl.py L3cut/Carbon_Compsite_Central_Spar_L3.stl L3/Carbon_Compsite_Central_Spar_L004.stl L3cut/Carbon_Composite.stl
python merge_stl.py L3cut/Carbon_Composite.stl L3cut/Carbon_Composite_Stave_Skin_L3_updated_Top.stl L3cut/Carbon_Composite.stl 
python translate_stl.py L3cut/Carbon_Composite_Stave_Skin_L3_updated_Top.stl L3cut/Carbon_Composite_Top.stl --dy -0.01
python merge_stl.py L3cut/Carbon_Composite.stl L3cut/Carbon_Composite_Top.stl L3cut/Carbon_Composite.stl
python translate_stl.py L3cut/Carbon_Composite_Stave_Skin_L3_updated_Bottom.stl L3cut/Carbon_Composite_Bottom.stl --dy 0.01
python merge_stl.py L3cut/Carbon_Composite.stl L3cut/Carbon_Composite_Bottom.stl L3cut/Carbon_Composite.stl
rm L3cut/Carbon_Composite_Bottom.stl                        
rm L3cut/Carbon_Compsite_Central_Spar_L3.stl
rm L3cut/Carbon_Composite_Stave_Skin_L3_updated_Bottom.stl  
rm L3cut/Carbon_Composite_Top.stl
rm L3cut/Carbon_Composite_Stave_Skin_L3_updated_Top.stl     
rm L3cut/Carbon_Compsite_Central_Spar_L004.stl
# Scale K9 components by 0.97 to avoid overlaps with carbon fibre
python scale_stl.py L3cut/K9_Split_CrossBrace_L004.stl L3cut/K9_Split_CrossBrace_L004.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_L005.stl L3cut/K9_Split_CrossBrace_L005.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_L006.stl L3cut/K9_Split_CrossBrace_L006.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_L3_End001.stl L3cut/K9_Split_CrossBrace_L3_End001.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_L3_End.stl L3cut/K9_Split_CrossBrace_L3_End.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_L3.stl L3cut/K9_Split_CrossBrace_L3.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_Sensor_4-077.stl L3cut/K9_Split_CrossBrace_Sensor_4-077.stl --scale 0.97
python scale_stl.py L3cut/K9_Split_CrossBrace_Sensor_4-078.stl L3cut/K9_Split_CrossBrace_Sensor_4-078.stl --scale 0.97
python scale_stl.py L3cut/MirrorK9_Split_CrossBrace1_L3_End.stl L3cut/MirrorK9_Split_CrossBrace1_L3_End.stl --scale 0.97
python scale_stl.py L3cut/MirrorK9_Split_CrossBrace2_L3_End.stl L3cut/MirrorK9_Split_CrossBrace2_L3_End.stl --scale 0.97
# Convert .stl files to ASCII and then to .gdml convertL3.sh contains offset parameters to set origin to centre of stave
./convertL3.sh L3cut
# Delete .stl files and rename directory of .gdml file to L3
rm L3cut/*.stl
rm -rf L3
mv L3cut L3

