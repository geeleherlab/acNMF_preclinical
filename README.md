# Classification-MES_like_cells
- Note that you will need to download essential data files from <link to OSF> in order to run the scripts. 

Here I provide scripts to run InferCNV and classification models on three datasets: GOSH (PD46693), Dong (T200, T214, T69, T230), and Mouse (NB831, NB837, NB839, NB847, NB849, NB853, NB856, NB883). 

The scripts perform on sample basis, results were saved to Output folder.

Order of scripts and output: 
1. InferCNV_\<sample\>.R:
  - reads cnmf RDS files, metadat file, and gene order file as input
  - outputs cellname.RDS, preliminary InferCNV results, and InferCNV results with HMM
2. Classification_\<sample\>.R:
  - reads cellname.RDS, preliminary InferCNV results, and InferCNV results with HMM, cnmf RDS files, gene order file as input
  - saves classification models and balanced accuracy on validation data
  - plots boxplots of GLM scores on test data and validation data 
