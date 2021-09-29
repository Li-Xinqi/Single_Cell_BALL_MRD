# Single_Cell_BALL_MRD
Codes and trained models for Paper 'Elucidating Minimal Residual Disease of Pediatric B cell Acute Lymphoblastic Leukemia by Single Cell Analysis'.

#### B cell differentiation stage classifier

The B cell differentiation stage classifier  was trained by the one-class logistic regression classifier to distinguish different differentiation stages (HSC/LMPP, CLP, proB, preBI, preBII, immatureB, matureB and activatedB).

#### Non-leukemic/leukemic cell classifier

The Non-leukemic/leukemic cell classifier was trained by a classical binary logistic regression to distinguish non-leukemic and leukemic cells in B-ALL samples.

#### Notes: Codes for the two classifiers are in the folder "Codes", and the trained models are in the folder "Data".

#### R packages and versions

```
R version 3.6.1 (2019-07-05)
Packages and versions:
 [1] magrittr_1.5                Seurat_3.1.5                
 [3] Rcpp_1.0.5                  knitr_1.28                 
 [5] DescTools_0.99.36           pROC_1.16.2                
 [7] survminer_0.4.7             survival_3.1-12            
 [9] gelnet_1.2.1                sampling_2.8               
[11] caret_6.0-86                lattice_0.20-41            
[13] ggalluvial_0.11.3           scales_1.1.1               
[15] VennDiagram_1.6.20          futile.logger_1.4.3        
[17] org.Hs.eg.db_3.10.0         AnnotationDbi_1.48.0       
[19] clusterProfiler_3.14.3      pheatmap_1.0.12            
[21] plyr_1.8.6                  cowplot_1.0.0              
[23] DropletUtils_1.6.1          SingleCellExperiment_1.8.0 
[25] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
[27] BiocParallel_1.20.1         matrixStats_0.56.0         
[29] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
[31] IRanges_2.20.2              S4Vectors_0.24.4           
[33] reshape2_1.4.4              scCancer_2.1.0             
[35] monocle_2.14.0              DDRTree_0.1.5              
[37] irlba_2.3.3                 VGAM_1.1-3                 
[39] Biobase_2.46.0              BiocGenerics_0.32.0        
[41] Matrix_1.2-18               ggpubr_0.3.0               
[43] ggplot2_3.3.1               dplyr_1.0.0                
```
