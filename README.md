# acNMF (v1.0.0)
This repository describes the acNMF method, along with supporting code for subsequent analyses.
<br><br>
![Alt Text](images/acNMF_schematic.png)

## acNMF Method
A notebook of the generic method with detailed instructions can be found at rchapple2.github.io/acNMF/.<br><br>

i) The first step of the acNMF method is to perform a random split of counts x cells single cell gene expression matrix to obtain two roughly equal datasets.<br>
ii) Second, [cNMF](https://github.com/dylkot/cNMF/tree/master) is performed over a wide range of ranks on each of these independent data splits. <br>
iii) Graph therory is next used to identify similar latent factors generated in each cNMF analysis. For each rank, jaccard similarity scores across many jaccard lengths are calculated between the two datasets, and if statistical significance is reached they are represented as node pairs in the network.  Community detection algorithms then identify groups/pairs of interconnected nodes and are recorded.<br>
iv) Finally, the results are plotted across all ranks and jaccard lengths.  Inflection point detection of stable curves identify the most suitable rank and jaccard length for that particular dataset.<br>

## acNMF on Simulated Data
The scripts to reproduce the acNMF analysis on simulated data can be found in the Simulated Data folder.  These should be run in the following order. <br><br>
**Simulate.ipynb** Jupyter notebook modified from [cNMF publication](https://github.com/dylkot/cNMF/blob/master/Tutorials/analyze_simulated_example_data.ipynb)<br>
**acNMF_input.R** Convert file into acNMF compatable format<br>
**cNMF_runscript and cnmf_v2.0.py** Modified [cNMF](https://github.com/dylkot/cNMF/tree/master) code that enables faster runtimes on HPC environments<br>
**acNMF_output.R** Calculates Jaccard Similarity and plots results<br>

## acNMF on Neuroblastoma scRNA-seq Datasets
The scripts to reproduce this analysis are found in the acNMF Method folder, and should be run in the following order.<br><br>
**input_processing.R** This code is representative of the processing required for a single dataset to become compatable for acNMF.  This code was performed on each dataset in our analysis separately.<br>
**cnmf_splitrun.py** This code is used to run [cNMF](https://github.com/dylkot/cNMF/tree/master) on each data split generated from the input_processing.R script.<br>
**post_cNMF_analysis.R**  This script calculates the Jaccard index on the cNMF output, plots the results across all ranks, and chooses the most appropriate rank for each dataset. <br>
**cNMF_meta_analysis.R**  This script is used to generate the network graph in which similar gene expression programs from independent datasets are represented as interconnected nodes in the network.<br> 

## Post-acNMF Analyses
The scripts to reproduce these analysis are found in the Post-acNMF Analysis folder. <br><br>
**inferCNV.R** Conducts inferCNV analysis using a pre-defined reference and the subclustering module.<br>
**DESeq2.R** Performs DESeq2 and GSEA on pseudobulked mouse scRNA-seq dataset.<br>
**velocyto.bsub and scvelo.py** Creates loom file and performs RNA velocity analysis, respectively.<br>
**Monocle.R** Calculates pseudotime trajectory plot for mouse dataset.<br>

## Scripts for Shiny App and Accompanying Gene Expression Progam Reports
The scripts for the Shiny App and GEP reports are found in the Reports folder.<br><br>
### GEP Reports
**literatureCuratedGenes.R** Contains all gene sets to which each gene expression program is compared.  Although this file was generated from the neuroblastoma literature, this file can be modified to include any gene set from any domain. <br>
**commonCodeForSummaries.R** This file loads all the data that is ubiquitously required for each GEP report.<br>
**gepSummaryKnitr.R** This file contains all of the Knitr code to generate the HTML report including all interactive figures, statistics, and gene set comparisons.<br><br>
To generate the reports: 
```{r, message=F}
Rscript scriptToKnitAllPrograms.R
```

### Shiny App
**navbarpage_srcdata.R** This file is needed to create the metadata tables and additional figures that are included in the NB_meta_shiny app.
**NB_meta_analysis_Shiny_app.R** Generates the Shiny app.
