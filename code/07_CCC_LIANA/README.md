## Cell-cell communication analysis with LIANA

## LIANA
A collection of cell-cell communication methods and databases and provides consensus of CCC methods result

* More details: https://saezlab.github.io/liana/
* The anlaysis was conducted with LIANA v0.1.8

## To run the analysis
You can source start_CCC.R, which will create jobs for each sample(R3)/brain piece(R2) using `ccc_batch_config.sh` that in turn run `CCC_withLIANA`.

After having the results, one can run `Vis_snRNA_CCC.R` to create heatmaps for visualization.

Some experimental code that summarizing CCC analysis are in `Archive/interpret_CCC.R`
