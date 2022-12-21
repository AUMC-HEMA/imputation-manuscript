### Introduction

This Github repository contains all code necessary for reproducing the analysis for our manuscript about imputation of MFC data. In general, all imputation and pre-processing was performed in R after which Jupyter Notebooks (Python) were used to visualize and interpret the resulting imputed FCS files. The following sections will explain which scripts were used for what purpose.



### Pre-processing and imputation (R)

**Transformation and generation of pseudo-tubes**

Pre-processing and splitting of pre-gated (T-cell) data was performed using the script ```Split_MFC.R```. In this script, pre-gated FCS files are transformed with unique co-factors and split into pseudo-tubes. This results in 3 different FCS files with either the suffix _gt (original, transformed data), _ff1 (pseudotube 1) or _ff2 (pseudotube 2).



**Imputation of pseudo-tubes**

Imputation of pseudo-tubes is performed in the script ```Impute_MFC.R```. Please account for the fact that this script does not impute for Infinicyt and requires these files to be already present beforehand. 

The script ```CytoBackBone_MOD.R``` contains a modified version of the ```merge``` function which is used in ```Impute_MFC.R```.

Additionally, the script ```KNN.R``` imputes split files using the FNN library using k=1 for a comparison between KNN and Infinicyt as described in the supplementary figure.



**Generation of gating files**

Gating was performed on aggregated files of ground truth and imputed data. These were generated using the script ```Create_MFC_Aggregates.R```



**Generation of gating labels**

The scripts ```Get_MFC_Gating.R```  and ```Get_SingleMarker_MFC_Gating.R``` usesFlowJo workspace files in combination with the aggregated fcs files to generate csv files containing the gating labels.



### Analysis and interpretation of results (Python)

All analyses were performed in Jupyter notebooks. 

* Figure 3 and 4 were generated in ```Density_Plots.ipynb```
* Figures 2, 5 and Supplemental Figures 5 and 6 were generated in ```Distance_Plots.ipynb```
* Supplemental Figure 4 was generated in ```Infinicyt_Plots.ipynb```
* Figures 6 and 7 and Supplemental Figures 7 and 8 were generated in ```Labeling_Plots.ipynb```

