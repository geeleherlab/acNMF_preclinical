# acNMF (v1.0.0)
This repository describes the acNMF method, along with supporting code for subsequent analyses.
<br><br>
![Alt Text](images/acNMF_schematic.png)

## acNMF on Simulated Data
The scripts to reproduce the acNMF analysis on simulated data can be found in the Simulated Data folder <br><br>
**Simulate.ipynb** Jupyter notebook modified from [cNMF publication](https://github.com/dylkot/cNMF/blob/master/Tutorials/analyze_simulated_example_data.ipynb)
**acNMF_input.R** Convert file into acNMF compatable format
**cNMF_runscript and cnmf_v2.0.py** Modified [cNMF](https://github.com/dylkot/cNMF/tree/master) code that enables faster runtimes on HPC environments
**acNMF_output.R** Calculates Jaccard Similarity and plots results

## Post-acNMF Analyses
### RNA Velocity
### inferCNV
### DESeq2

## Scripts for Shiny App and Accompanying Gene Expression Progam Reports
