# Code


### 01_align  
Use cellranger to process and align data

### 02_cellranger_metrics  

`00_metrics_functions.R`  

`01_Tran_et_al_metrics.R`  

`02_cellranger_metrics.R`  
Extract Cell Ranger QC metrics from `metrics_summary.csv` files. Save table in 
`..processed-data/02_cellranger_metrics/tran_metrics.csv`  

`03_combine_metrics.R`  
Compare Cell Ranger metrics with previous studies, create boxplots

### 03_build_sce 
1. `01_build_basic_sce.R` Build sce object with all output from cellRanger in 

2. `02_get_droplet_scores.R` For each sample run `DropletUtils::emptyDrops()` to find empty droplets

3. `03_droplet_qc.R` Read in scores from `emptyDrops()` filer out empty droplets, preform other QC with `scater::isOutlier` (precent mito, total coutns, detected featuers) and drop failing nulci, compute doublet scores for assesment after clustering

4. `04_GLM_PCA.R` Calculate reduced dimensions were calculated with generalized linear models principal component analysis

5. `05_harmony_correction.R` Correct for batch effects with harmony, used sample as group variable after testing 

6. `06_cluster.R` First step of hierarchical clustering, use `igraph::cluster_walktrap()` to find prelim clusters

7. `07_hierarchical_cluster.R` Second step of hierarchical clustering, pseudobulk and hierarchical cluster prelim clusters

8. `08_cluster_annotation.R`  Explore QC metrics and cell type identity of clusters, plot marker genes, assign cell types to clusters
*layer level annotations were added in [spatilDLPFC](https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/analysis/12_spatial_registration_sn)

9. `09_find_markers.R` Find marker gense for HC clusters, not used much down stream


### 04_synapse_upload  
Code relevant to preparing data & metadata tables to upload to [synapse](www.synapse.org/)

### 05_explore_sce  
0. `00_pub_metrics.R` check metrics needed for manuscript

1. `01_reduced_dim_plots.R` Plot UMAP and TSNE plots with cell type annotataions

2. `02_cellType_prop.R` Plot and explore cell type proportions by sample and other breakdowns, make composition bar plots

3. `03_explore_markers.R` Plot marker genes including heatmaps and violin plots

4. `04_3D_UMAP.R` Try a 3D UMAP of the data (Fun!)
 
5. `05_azimuth_validation.R` Applied reference-based mapping tool Azimuth
 
6. `06_explore_azimuth_annotations.R` Check out correspondence of cell clusters from our annotaitons vs. Azimuth

7. `07_convert_SCE_HDF5.R` Export data in portable form for iSEE app and sharing with collaborators

8. `08_pseudobulk_cellType.R` Pseudobulk data and store for various applications

9. `09_SCZ_interaction_expression.R`

### 06_iSEE_apps 
Build and launch [iSEE app ](https://libd.shinyapps.io/spatialDLPFC_snRNA-seq/)

### 07_ccc
Code relevant to cell-to-cell communication analysis
