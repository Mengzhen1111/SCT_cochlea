# Spatial Transcriptomics Analysis Pipeline
**Authors:Mengzhen yan, Chengdu Institute of Biology, Chinese Academy of Sciences, Chengdu 610213, China
**Corresponding author:Jun Li, lijun0309@wchscu.cn
This repository contains all custom scripts and workflows required to reproduce the main figures and supporting analyses of the manuscript:
> **Spatiotemporal transcriptomic analyses reveal molecular gradient patterning during development and the tonotopic organization along the cochlear axis**
##Session Information
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Nebulosa_1.0.1        cowplot_1.2.0         ggforce_0.5.0        
 [4] semla_1.4.2           ComplexHeatmap_2.20.0 colorRamp2_0.1.0     
 [7] RColorBrewer_1.1-3    clustree_0.5.1        ggraph_2.2.2         
[10] readxl_1.4.5          patchwork_1.3.2       png_0.1-8            
[13] config_0.3.2          reshape2_1.4.5        gridExtra_2.3        
[16] ggpubr_0.6.3          spacexr_2.2.1         Matrix_1.7-3         
[19] lubridate_1.9.5       forcats_1.0.1         stringr_1.6.0        
[22] dplyr_1.1.4           purrr_1.0.2           readr_2.2.0          
[25] tidyr_1.3.2           tibble_3.3.0          ggplot2_4.0.3        
[28] tidyverse_2.0.0       Seurat_5.4.0          SeuratObject_5.4.0   
[31] sp_2.2-1             

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.23            splines_4.4.1              
  [3] later_1.4.8                 cellranger_1.1.0           
  [5] polyclip_1.10-7             fastDummies_1.7.5          
  [7] lifecycle_1.0.5             rstatix_0.7.3              
  [9] doParallel_1.0.17           globals_0.19.1             
 [11] lattice_0.22-7              MASS_7.3-65                
 [13] backports_1.5.0             magrittr_2.0.3             
 [15] plotly_4.12.0               yaml_2.3.12                
 [17] httpuv_1.6.17               otel_0.2.0                 
 [19] sctransform_0.4.3           spam_2.11-3                
 [21] sessioninfo_1.2.3           pkgbuild_1.4.8             
 [23] spatstat.sparse_3.1-0       reticulate_1.45.0          
 [25] pbapply_1.7-4               zlibbioc_1.50.0            
 [27] abind_1.4-8                 pkgload_1.5.0              
 [29] GenomicRanges_1.56.2        Rtsne_0.17                 
 [31] BiocGenerics_0.56.0         pracma_2.4.6               
 [33] tweenr_2.0.3                GenomeInfoDbData_1.2.12    
 [35] circlize_0.4.17             IRanges_2.38.1             
 [37] S4Vectors_0.42.1            ggrepel_0.9.8              
 [39] irlba_2.3.7                 listenv_0.10.1             
 [41] spatstat.utils_3.2-2        goftest_1.2-3              
 [43] RSpectra_0.16-2             spatstat.random_3.4-5      
 [45] fitdistrplus_1.2-6          parallelly_1.47.0          
 [47] DelayedArray_0.30.1         codetools_0.2-20           
 [49] shape_1.4.6.1               tidyselect_1.2.1           
 [51] UCSC.utils_1.0.0            farver_2.1.2               
 [53] viridis_0.6.5               stats4_4.4.1               
 [55] matrixStats_1.5.0           spatstat.explore_3.8-0     
 [57] jsonlite_2.0.0              GetoptLong_1.1.0           
 [59] ks_1.15.1                   ellipsis_0.3.2             
 [61] tidygraph_1.3.1             progressr_0.19.0           
 [63] Formula_1.2-5               ggridges_0.5.7             
 [65] survival_3.8-6              iterators_1.0.14           
 [67] dbscan_1.2.4                foreach_1.5.2              
 [69] tools_4.4.1                 ica_1.0-3                  
 [71] Rcpp_1.1.1                  glue_1.8.1                 
 [73] SparseArray_1.4.8           MatrixGenerics_1.16.0      
 [75] usethis_3.2.1               GenomeInfoDb_1.40.1        
 [77] withr_3.0.2                 fastmap_1.2.0              
 [79] shinyjs_2.1.1               digest_0.6.39              
 [81] timechange_0.4.0            R6_2.6.1                   
 [83] mime_0.13                   colorspace_2.1-1           
 [85] scattermore_1.2             tensor_1.5.1               
 [87] dichromat_2.0-0.1           spatstat.data_3.1-9        
 [89] generics_0.1.4              data.table_1.18.2.1        
 [91] S4Arrays_1.4.1              graphlayouts_1.2.3         
 [93] httr_1.4.8                  htmlwidgets_1.6.4          
 [95] uwot_0.2.4                  pkgconfig_2.0.3            
 [97] gtable_0.3.6                lmtest_0.9-40              
 [99] S7_0.2.2                    XVector_0.44.0             
[101] SingleCellExperiment_1.26.0 htmltools_0.5.9            
[103] carData_3.0-6               dotCall64_1.2              
[105] clue_0.3-67                 Biobase_2.64.0             
[107] scales_1.4.0                spatstat.univar_3.1-7      
[109] rstudioapi_0.18.0           tzdb_0.5.0                 
[111] rjson_0.2.23                nlme_3.1-168               
[113] cachem_1.1.0                zoo_1.8-15                 
[115] GlobalOptions_0.1.3         KernSmooth_2.23-26         
[117] vipor_0.4.7                 parallel_4.4.1             
[119] miniUI_0.1.2                ggrastr_1.0.2              
[121] pillar_1.11.1               vctrs_0.7.3                
[123] RANN_2.6.2                  promises_1.5.0             
[125] car_3.1-5                   xtable_1.8-8               
[127] cluster_2.1.8.2             beeswarm_0.4.0             
[129] zeallot_0.2.0               magick_2.9.1               
[131] mvtnorm_1.3-7               cli_3.6.6                  
[133] compiler_4.4.1              rlang_1.2.0                
[135] crayon_1.5.3                future.apply_1.20.2        
[137] ggsignif_0.6.4              mclust_6.1.2               
[139] ggbeeswarm_0.7.3            plyr_1.8.9                 
[141] fs_2.0.1                    stringi_1.8.4              
[143] viridisLite_0.4.3           deldir_2.0-4               
[145] lazyeval_0.2.2              devtools_2.5.0             
[147] spatstat.geom_3.7-3         RcppHNSW_0.6.0             
[149] hms_1.1.4                   future_1.70.0              
[151] shiny_1.13.0                SummarizedExperiment_1.34.0
[153] ROCR_1.0-12                 igraph_2.0.3               
[155] broom_1.0.12                memoise_2.0.1   

## Data Availability
All raw and reference data are publicly available:
The raw single-cell spatial transcriptomic FASTQ files have been deposited at https://ngdc.cncb.ac.cn/gsa under accession CRA019653.
All original histology images (H&E) required to align transcriptomic spots to anatomical positions have been deposited in the Figshare (10.6084/m9.figshare.33025988).
The mouse inner ear single-cell RNA-seq data at E16 and P7 stage from publication can be accessible in the Gene Expression Omnibus data repository under accession code (GSE137299).
The mouse inner ear single-cell RNA-seq data at adult stage from our publication can be accessible in the National Genomics Data Center under the number (CRA015540).

## Raw Data Processing

**Important:** The raw sequencing data for this study were generated using the **BMKMANU S series** platform. 
The initial processing from FASTQ to expression matrix was performed using the official **`BSTMatrix`** software.
This software can be obtained from http://www.bmkmanu.com/portfolio/tools.

The complete processing pipeline, including barcode/UMI identification, STAR alignment, and tissue detection was executed with the following command:

```bash
./BSTMatrix -c config.txt -s 0

##config.txt
FQ1          /path/to/read_1.fq.gz       
FQ2          /path/to/read_2.fq.gz   
FLU          /path/to/flu_info.txt   
HE     /path/to/HE.tif                              
 GenomeVer    
INDEX      /path/to/STAR/index/dir/       
GFF          /path/to/ref/gene/gff3/file               
FEATURE          /path/to/features.tsv  
OUTDIR    /path/to/result/dir/        
PREFIX    outfile-prefix      
BCType                  V2                             
BCThreads              8                              
Sjdboverhang        100                  
STARThreads          8 
