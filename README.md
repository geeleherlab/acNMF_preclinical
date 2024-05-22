# acNMF (v1.0.0)
This repository describes the acNMF method, along with supporting code for subsequent analyses.
<br><br>
![Alt Text](images/acNMF_schematic.png)

## acNMF Method
A notebook of the generic method with detailed instructions can be found at rchapple2.github.io/acNMF/.<br><br>

## acNMF on Simulated Data
The scripts to reproduce the acNMF analysis on simulated data can be found in the Simulated Data folder <br><br>
**Simulate.ipynb** Jupyter notebook modified from [cNMF publication](https://github.com/dylkot/cNMF/blob/master/Tutorials/analyze_simulated_example_data.ipynb)<br>
**acNMF_input.R** Convert file into acNMF compatable format<br>
**cNMF_runscript and cnmf_v2.0.py** Modified [cNMF](https://github.com/dylkot/cNMF/tree/master) code that enables faster runtimes on HPC environments<br>
**acNMF_output.R** Calculates Jaccard Similarity and plots results<br>

## acNMF on Neuroblastoma scRNA-seq Datasets
The scripts to reproduce this analysis are found in the acNMF Method folder. <br><br>
**input_processing.R** This code is representative of the processing required for a single dataset to become compatable for acNMF.  This code was performed on each dataset in our analysis separately.<br>
**

## Post-acNMF Analyses
The scripts to reproduce these analysis are found in the Post-acNMF Analysis folder. <br><br>
**inferCNV.R** Conducts inferCNV analysis using a pre-defined reference and the subclustering module.<br>
**DESeq2.R** Performs DESeq2 and GSEA on pseudobulked mouse scRNA-seq dataset.<br>
**velocyto.bsub and scvelo.py** Creates loom file and performs RNA velocity analysis, respectively.<br>
**Monocle.R** Calculates pseudotime trajectory plot for mouse dataset.<br>

## Scripts for Shiny App and Accompanying Gene Expression Progam Reports
The scripts for the Shiny App and GEP reports are found in the Reports folder<br><br>

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
