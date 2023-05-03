### Introduction

This Github repository contains all code necessary for reproducing the analysis for our manuscript about imputation of MFC data. In general, all imputation and pre-processing was performed in R after which Jupyter Notebooks (Python) were used to visualize and interpret the resulting imputed FCS files. The following sections will explain which scripts were used for what purpose.



### Pre-processing and imputation (R)

**Transformation and generation of pseudo-tubes**

Pre-processing and splitting of pre-gated (T-cell) data was performed using the script ```Split_MFC.R```. For the real-life data, the ```Preprocess_MFC_RL.R``` was used. For Mosmann and Nilsonn datasets, the ```Split_Mosmann.R``` and ```Split_Nilsson.R``` files were used.

In this script, pre-gated FCS files are transformed with split into pseudo-tubes. This results in 3 different FCS files with either the suffix _gt (original, transformed data), _ff1 (pseudotube 1) or _ff2 (pseudotube 2).



**Imputation of pseudo-tubes**

Imputation of pseudo-tubes is performed in the script ```Impute_MFC.R```.  For the other datasets, please find the other R scripts which are named after the respective dataset.

Please account for the fact that this script does not impute for Infinicyt and requires these files to be already present beforehand. 

The script ```CytoBackBone_MOD.R``` contains a modified version of the ```merge``` function which is used in ```Impute_MFC.R```.



**Generation of gating files**

Gating was performed on aggregated files of ground truth and imputed data. These were generated using the script ```Create_MFC_Aggregates.R```



**Generation of gating labels**

The script ```Get_MFC_Gating.R```  uses FlowJo workspace files in combination with the aggregated fcs files to generate csv files containing the gating labels.



### Where to find the code for every plot in the manuscript

All analyses were performed in Jupyter notebooks. 



Figure 2: ```Distance_Plots.ipynb```

Figure 3: ```Density_Plots.ipynb```

Figure 4: ```Density_Plots.ipynb```

Figure 5: ```Distance_Plots.ipynb```

Figure 6: ```Backbone_Experiment.ipynb```

Figure 7: ```Labeling_Plots.ipynb```



Supplemental Figure 3: ```Distance_Plots.ipynb```

Supplemental Figure 4: ```Infinicyt_Plots.ipynb```

Supplemental Figure 5:  ```Distance_Plots.ipynb```

Supplemental Figure 6:  ```Distance_Plots.ipynb```

Supplemental Figure 7:  ``Distance_Plots.ipynb``

Supplemental Figure 8: ```FlowSOM_analysis.ipynb```

Supplemental Figure 9: ```Labeling_Plots.ipynb```

Supplemental Figure 10: ```FlowSOM_analysis.ipynb```
