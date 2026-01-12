# A Machine Learning Framework for Heavy-Flavor Analysis in the ePIC Experiment

Supported by The FAIR (Future Artificial Intelligence Research) Spoke 6 Project, funded by the NextGenerationEU program in Italy

based on hipe4ml:https://doi.org/10.5281/zenodo.5070131

Shyam Kumar; INFN Bari, Italy
shyam.kumar@ba.infn.it; shyam055119@gmail.com

```
Instructions to Run the code:
1. Install hipe4ml
pip install hipe4ml (Ubuntu) or brew install libomp (Mac)
2. This requires signal and background features
Run D0 reconstruction for D0 Sample and DIS Samples
Put the output root files in the directories: D0_Sample and DIS_Sample 
3. Go to the directory: ML_model_ePIC (I added the root files from ep 10x100, Q^2>1 GeV^2, Campaign: 25.04.1 as an example)
source Run_ML_withData.sh
```

## Tutorial for D0 Reconstruction in the ePIC experiment

If you want to get familiar with hipe4ml, the following tutorials are available:

| Type | Link |
| -------------- | ------------- |
| Binary classifier | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1JtKXPSfRBXTOu4-yUd_h7fOENDyWLCDm?usp=sharing) |
