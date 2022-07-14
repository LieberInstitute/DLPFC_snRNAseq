# Code


### 01_align/

### 02_cellranger_metrics/

`00_metrics_functions.R`  

`01_Tran_et_al_metrics.R`  

`02_cellranger_metrics.R`  
Extract Cell Ranger QC metrics from `metrics_summary.csv` files. Save table in 
`..processed-data/02_cellranger_metrics/tran_metrics.csv`  

`03_combine_metrics.R`  
Compare Cell Ranger metrics with previous studies, create boxplots

### 03_build_sce/
`01_build_basic_sce.R`

`02_get_droplet_scores.R`

`03_droplet_qc.R`

`04_GLM_PCA.R`

`05_harmony_correction.R`

`06_cluster.R`

`07_hierarchical_cluster.R`

`08_cluster_annotation.R`

`09_find_markers.R`
