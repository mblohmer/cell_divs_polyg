# Quantifying cell divisions along evolutionary lineages in cancer

This is the Github repository  for the paper *doi.xxx*. It contains **A)** all data analyzed in this paper and **B)** the steps to recreate all analyses presented.

![Quantifying cell divisions along evolutionary lineages in cancer](/data/title_plot.png "Title plot")

## Data structure
Raw data needed to reproduce all parts of this paper can be found in the *data* directory. The *data* directory contains both raw polyguanine profiles and metadata.  

Raw polyguanine profiles in the subdirectories *multiple_lung_cancers*, *hnssc*, and *in_vitro_evolution* were newly generated for this paper. Raw polyguanine profiles in the subdirectories *colon_iscs*, *natgen_2020*, and *science_2017* were part of previous publications.  
Each subdirectory, e.g. *multiple_lung_cancers*, contains further subdirectories for each patient. Within each patient's subdirectory, another subdirectory termed *repre_repli_data* contains the raw fluorescent intensity data. Every text file contains data on one genotyped polyguanine repeat. The columns in these files denote the sample name. Each row denotes an allele of a distinct length, consecutively listed from the shortest observed allele to the longest observed allele.  

Newly generated metadata for this paper is contained in the files *Supplementary_tables.xlsx* (clinical data on patients with multifocal lung cancer), *cell_division_numbers.xlsx* (number of cell divisions separating  each sample in the in vitro evolution experiment), and *polyG_chr_positions.tsv* (mapping polyguanine repeats into chromosomes). The rest of the tables are supplementary information from other publications. The data directory *mutation_burden* contains previously published mutation burden data.

## Reproducing the analysis

Pre-rendered HTML files with the R code and resulting figures can be found under *rendered-analysis*. 
The R scripts to redo the analysis, starting from the raw data, is stored at *analysis*. To run the code, R version >=4.1.2 and a program to run Jupyter Notebooks (e.g <a href="https://docs.jupyter.org/en/latest/install/notebook-classic.html">Anaconda</a> or <a href="https://code.visualstudio.com/docs/datascience/jupyter-notebooks">VSCode</a> ) is needed.

These R packages are needed to generate the results:

```r
analysis.packages <- c( "data.table","tidyverse","ggtree",
                        "phytools","tidytree", "phytools", 
                        "readxl", "lubridate", "ape", 
                        "rstatix", "treeio", "survival",
                        "survminer", "contsurvplot", "pammtools",
                        "devtools", "parallel")

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
		            "zoo", "ggfortify", "dendextend"
                    "ggtext")

new.plot.packages <- plot.packages[!(plot.packages %in% installed.packages()[,"Package"])]
if(length(new.plot.packages)>0) BiocManager::install(new.plot.packages)
```
All data plots will be saved in the corresponding subdirectory under *plots*. 
Intermediary results files are in saved under *results*.  