
# Scripts for the manuscipt comparing of spatial reordering across single-cell protocols using datasets on the liver lobule


This repository contains 

1. Code to reproduce the analyses and figures in "Reproducibility across single-cell RNA-seq protocols for spatial ordering analysise". This includes all preprocessing steps and the generation of figures in the main text and supplement. 

    * Full-length dataset is available at GSE116140, processed data is included as a supplement to the manuscript.
    * Download UMI dataset summary NIHMS70855-supplement-Supplementary_Table_3.xlsx and counts NIHMS70855-supplement-Supplementary_Table_1.zip from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321580/
    * Download 10X dataset at: https://figshare.com/articles/Robject_files_for_tissues_processed_by_Seurat/5821263
    * Run all codes in the the Preprocess folder first, then the scripts in the Analysis folder in the indicated order. Directory set-up is a main folder called scSpatialReconstructCompare-Paper with subfolders: DATA, RDATA, PLOTS, and OUT.
  
2. A Shiny application exists here to interactively explore the correlations as shown in Figure 4 and Supp. Figure 4 for any KEGG category.


# Instructions to view data in Shiny app


1. Assuming you have R installed already, you will need to have these packages installed:

```
library(shiny)
library(shinyFiles)
library(scales)
library(viridis)

if (!requireNamespace("clusterProfiler", quietly=TRUE))
install.packages("clusterProfiler")

if (!requireNamespace("org.Mm.eg.db", quietly=TRUE))
install.packages("org.Mm.eg.db")

```

2. To run the Shiny App, type in:

```
runGitHub("rhondabacher/scSpatialReconstructCompare-Paper")
```

### Questions

Any questions may be reported using Github issues: https://github.com/rhondabacher/scSpatialReconstructCompare-Paper/issues

or emailing Rhonda Bacher at rbacher@ufl.edu

