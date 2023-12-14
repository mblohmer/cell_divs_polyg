# Tracking Cell Divisions

This is the github repo for the upcoming paper.

## Introduction
This github repository contains the steps to recreate all analyses presented in the manuscript. Prerendered HTML files with the R code and resulting figures can be found under "rendered-analysis". 
The R scripts to redo the analysis from the raw data is stored at "unrendered-analysis". To run the code, R version >=4.1.2 and a program to run jupyter notebooks (e.g using <a href="https://docs.jupyter.org/en/latest/install/notebook-classic.html">Anaconda</a> or <a href="https://code.visualstudio.com/docs/datascience/jupyter-notebooks">VSCode</a> ) is needed.

## Reproducing the analysis

These R packages are needed to generate the data:

```r
analysis.packages <- c( "data.table","tidyverse","ggtree",
                        "phytools","tidytree", "phytools", 
                        "readxl", "lubridate", "ape", 
                        "rstatix", "treeio", "survival",
                        "survminer", "contsurvplot", "pammtools",
                        "devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

new.analysis.packages <- analysis.packages[!(analysis.packages %in% installed.packages()[,"Package"])]
if(length(new.analysis.packages)>0) BiocManager::install(new.analysis.packages)
```

To graphically display the analysis, these R packages are needed:

```r
plot.packages <- c( "ggpubr", "patchwork", "pheatmap", 
                    "ggrepel", "wesanderson", "nicolash2/ggdendroplot",
                    "ggbeeswarm", "ggbreak", "ggdist", 
                    "grid", "scales", "gridExtra",
		    "zoo", "ggfortify")

new.plot.packages <- plot.packages[!(plot.packages %in% installed.packages()[,"Package"])]
if(length(new.plot.packages)>0) BiocManager::install(new.plot.packages)
```

