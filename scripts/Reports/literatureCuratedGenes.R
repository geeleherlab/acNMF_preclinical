## This file contains manually curated gene sets for neuroblastoma, obtained from the existing literature.
## It cotains both "large" (literatureMarkers_largeGeneSets) and "small" (literatureMarkers_smallGeneSets) gene sets, these are displayed differently in the final report.


# List elements are formatted so that the description is in the first element.

# This list contains information from the following papers:
# Yuan: bioRxiv, 
# Dong: Cancer Cell, 
# Olsen: bioRxiv, 
# Kildisiute: Science Advances, 
# Kamavena: Nature Genetics
# Jansky: Nature Genetics
# Kinker: Nature Genetics

literatureMarkers_smallGeneSets <- list(

"Proliferating Cells"=c("These are canonical marker genes of cycling cells, obtained from Fig. 3 of Hsiao et al. (PMID 32312741)", "CDK1", "UBE2C", "TOP2A", "H4C5", "H4C3", "MKI67"), # These don't mark a sepcific cell cycle phase, but rather can broadly be used to predict cell cycle phase (mostly high in G2/M I think),... https://genome.cshlp.org/content/30/4/611.full


## *Yuan* (Dan Carter) smart-seq neuroblastoma paper, Bioxiv (https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1),
# Yuan cell types (Fig 1),
"Symphathoadrenal cells (Yuan)"=c("Obtained from Fig. 1 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "SOX10", "PHOX2B", "TUBB3", "CHGA"), # Dong put SOX10 in schwannCellPrecursor in Fig 3, and PHOX2B in chromaffinCell.

"Mesenchymal cells (Yuan)"=c("Obtained from Fig. 1 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "COL1A2", "COL6A1", "COL3A1", "VIM"),

"Immune cells (Yuan)"=c("Obtained from Fig. 1 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "PTPRC", "CXCR4"),

# Yuan cell types (Fig 6),
"Cytotoxic T cells (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "PRF1"),

"Other T cells (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "CD4", "CCR7", "SELL", "TCF7", "LEF1", "FOXP3"), # mixutre of T cells, e.g. CD4 = helper

"B cells (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "CD79A", "CD79b"),

"Plasmacytoid dendritic cells (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "LILRA4", "LAMP5", "IRF4"),

"Inflammatory Macrophages (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "LYZ", "CSTA", "VCAN", "FCN1"),

"Noninflammatory Macrophages (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "CD163", "VSIG4", "MS4A7", "MAF"),

"Adrenal Cortex (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "CYP21A2", "CYP11B1"),

"Endothelial Cells (Yuan)"=c("Obtained from Fig. 6 of Yuan et al. https://www.biorxiv.org/content/10.1101/2020.05.15.097469v1", "PECAM1", "VWF", "ENG"),

## *Dong* et al, Cancer Cell.
# Fig 3.(c),
"Schwann cell precursors (Dong)"=c("Obtained from Fig. 3c of Dong et al. (PMID 32946775).", "PLP1", "SOX10", "FOXD3"),
"Chromaffin Cells (Dong)"=c("Obtained from Fig. 3c of Dong et al. (PMID 32946775). Jansky et al (PMID) criticize this choice of markers, stating: the neuronal markers NPY, PRPH, NTRK1 and ISL1 are used to annotate chromaffin cells; these were clearly expressed in neuroblasts in our analyses.", "PHOX2B", "PHOX2A", "ISL1", "NYP", "PRPH", "NTRK1"), # PHOX2B is in "sympathoadrenal_yuan"

"Symphatoblasts (Dong)"=c("Obtained from Fig. 3c of Dong et al. (PMID 32946775). Jansky et al (PMID) criticize this choice of markers, stating:  Sympathoblasts were identified by Dong et al.28 based on the expression of CARTPT and INSM1, which we found expressed in chromaffin cells", "CARTPT", "INSM1"),

"Broad sympathetic marker (Dong)"=c("Obtained from Fig. 3c of Dong et al. (PMID 32946775).", "TH", "CHGB", "DBH", "GATA2", "GATA3", "HAND2"), # DBH is referred to as "adreneergic" in kharchenco paper.

# Fig. 5(A), - the Meta Programs, derived for each tumor separately, then combined, I think... IS THE ENTIRE META-PROGRAM SOMEHOW AVAILABLE?....
"Translation (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "EIF3E", "EIF3L", "EIF3F"), 

"Stress response (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "FOS", "JUNB", "JUN"),

"Chromaffin cell-like (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "NPY", "PRPH", "NTRK1"),

"Chromaffin cell development (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "GATA3", "PHOX2B", "HAND2"),

"Meta5 unclear (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "STMN1", "STMN2"), # not really discribed, genes seem like they may be cell cycle associated (directly involved in the cytoskeleton, I think),

"Cell cycle (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "TOP2A", "MKI67", "PCNA"),

"Meta_7 undefined sympathetic (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "TUBB2A", "CHGA", "DBH"), 

"Meta_8 undefined undifferentiated neuronal (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "SOX11", "CTNNB1"), # Sox11 is very highly and specifically expressed in fetal brain in GTEx. CTNNB1 is very broadly expressed, plays a role in cell adhesion, definitely involved in neural tube formation.

"Meta_9 ATRX (Dong)"=c("Markers obtained in computationally inferred meta program (by NMF applied to their neuroblastoma sc-RNA-seq data) from Fig. 5(a) of Dong et al. (PMID 32946775).", "ZFHX3", "ATRX"),





## *Olsen* (Kharchecko), biorxiv (https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1.full.pdf+html).
# From their main text. (pg. 3),
"Mesenchymal (Olsen)"=c("Markers listed in on page 3 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these markers are stated as being expressed in a group of mesenchymal-like cells (some at least claimed to be neuroblastoma cells) shown in their UMAP representation on Fig. 1D. They reference the two Nature Genetics paper.", "PRRX1", "LEPR", "PDGFRA", "DCN"), 

"Adrenergic (Olsen)"=c("Markers listed in on page 3 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these markers are for adrenergic genes", "TH", "DBH"),

# pg. 4
"Schwann cell precursor (Olsen)"=c("Markers listed in on page 4 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these markers are stated as being characteristic of Schwann cell lineage and Furlan et al. (Science 2017, PMID 28684471) is referenced to support this. These cell types seem to be claimed to form a bridge between the mesenchymal and adrenergic cell types.", "SOX10", "S100B", "ERBB3", "FOXD3", "PLP1"),

"Myelination (Olsen)"=c("Markers listed in on page 4 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these markers are state as being involved in Myelination and are expressed in the Schwann cell precursor cluster of their UMAP (Fig. 1D), these are suggested to represent more mature myelinating Schwann cells. MPZ is Myelin Protien Zero.", "LGI4", "PMP2", "MPZ"),

"Axon guidance (Olsen)"=c("Markers listed in on page 4 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these markers are for axon guidance instead of myelination", "SEMA3B"),

"Nerves and neurofiliment (Olsen)"=c("Markers listed in on page 4 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - NF200 is expressed in neurofiliments and is found close to SCPs in normal tissue, NRG1 is secreted by neurons, to which ERBB3 (on glia) binds.", "NF200", "NRG1"),


# pg. 6
"Normal Schwann cell precurors, but not Neuroblastoma Scps (Olsen)"=c("Based on their neuroblastoma data and stated on page 6 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1", "PHOX2B", "POSTN"),
"Neuroblastoma SCPs but not  normal SCPs (Olsen)"=c("The same as above", "GFRA3", "ITGB8"),


"Autonomic neurons mouse (olsen)"=c("Stated on page 6 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - seems to be referencing previous mouse lineage tracing studies in mouse which sucggested SCPs give rise to autonomic neurons and mesenchymal cells", "PHOX2B", "ASCL1"),
"Mesencchymal populations mouse (Olsen)"=c("Stated on page 6 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - seems to be referencing previous mouse lineage tracing studies in mouse which sucggested SCPs give rise to autonomic neurons and mesenchymal cells", "PRRX1", "FLI1"),
"Mixed phenotype identified at the fate split point mouse (olsen)"=c("Stated on page 6 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - seems to be referencing previous mouse lineage tracing studies in mouse which sucggested SCPs give rise to autonomic neurons and mesenchymal cells. markers are of the split point where SCPs differentiate.", "PRRX1", "PHOX2B", "SOX10"),


"Adrenergic neuroblatoma (Olsen)"=c("Stated on page 6 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - neuroblastoma adrenergic markers derived from their data, which are state as resembling sympathoblasts and lacking chromffin markers such as PNMT", "PRPH", "GAP43"), # they also mention here that these "lacked chromaffin markers such as PNMT"

# pg. 7 (they argue for heterogeneity in the adrenergic compartment, this is likely consistent with what we see),
"Adrenergic proliferating cluster (Olsen)"=c("Stated on pages 6 and 7 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - it is argued that upon reanalysis of the adrenergic neuroblastoma cells, there is further heterogeneity, grouped into proliferating, mature and immature cell populations. These genes were found as markers for proliferating adrenergic cells", "MKI67", "EZH2"),
"Adrenergic mature (olsen)"=c("Stated on pages 6 and 7 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - it is argued that upon reanalysis of the adrenergic neuroblastoma cells, there is further heterogeneity, grouped into proliferating, mature and immature cell populations. These markers are specific to the mature adrenergic cells", "NTRK1", "PRPH", "GAP43"), # the state these are "markers characteristic of differentiation into sympathetic neurons". NTRK1 is a landmark of spontaneous regression.
"Adrenergic immature (Olsen)"=c("Stated on pages 6 and 7 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - it is argued that upon reanalysis of the adrenergic neuroblastoma cells, there is further heterogeneity, grouped into proliferating, mature and immature cell populations. These is a marker listed for the immature adrenergic cells. SOX11 is referenced as important in early phases of pro-adrenergic differentiation (cites PMID 20147379) whereas SOX4 appears later. SOX11 is also associated with poor outcome in neuroblatoma bulk RNA-seq data (Olsen analysis).", "SOX11"),  

# pg. 8 (heterogeneity in the mysenchymal compartment - malignant and non-malignant cells (fibroblasts, pericytes),),
"Mesemchymal validated (Olsen)"=c("Stated on pages 7 and 8 of the main text of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - it is stated that a mesenchymal-like neuroblastoma cell state is validate by combined immunofluorescence staining for PDGFRA and FISH for PPMID gain, a gene on 17q.", "PDGFRA"),





# Fig. 1
# (a),
"Mesenchymal Fig 1D (Olsen)"=c("Selected mesenchymal marker genes shown in Fig. 1D of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these are highly expressed in their mesenchymal cluster on their UMAP.", "COL1A2", "LUM", "DCN", "COL3A1", "COL1A1", "COL6A2", "MGP", "SPARC", "CALD1", "BGN", "PRRX1", "LEPR", "PDGFRA"),

"Schwann cell precursor-like Fig 1D (Olsen)"=c("Similar to above, but for selected SCP marker genes.", "S100B", "PLP1", "SOX10", "ERBB3"),

"Adrenergic Fig 1D (Olsen)"=c("Similar to above, but for selected adrenergic marker genes", "STMN2", "RTN1", "MAP1B", "MLLT11", "UCHL1", "BEX1", "CD24", "DBH", "ELAVL4", "RGS5", "TH", "ISL1", "NRG1"),

# Fig. 2
# (b),
"SCPs mouse unique (Olsen)"=c("Selected list of SCP marker genes shown in Fig. 2B of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - these genes were unique to mouse SCPs as determined by the linearge tracing study Furlan et al. (Science 2017, PMID 28684471).", "SCRG1", "H1FX", "CRABP1", "HNRNPA1P48", "PHOX2B", "POSTN", "ASCL1"),
"SCPs mouse and neuroblastoma (Olsen)"=c("Similar to above, but for genes common to neuroblastoma SCP-like cells (as determined from their UMAP plot) and mouse SCP cells.", "PLP1", "S100B", "CDH19", "GPM6B", "FXYD1", "MPZ", "NRXN1", "DST", "ERBB3", "ABCA8", "CNN3", "CRYAB", "FOXD3-AS1", "SOX10", "FOXD3", "PMP2", "LGI4", "SEMA3B"),
"SCPs neuroblastoma unique (Olsen)"=c("Similar to above, but for genes unique to neuroblastoma SCP-like cells. Note FN1 is often used as a mesenchymal marker gene.", "IFITM3", "GFRA3", "SELENOM", "ITGB8", "ANXA1", "FN1", "TIMP3"),

# Fig. 3
"Bridge region mesencchymal (Olsen)"=c("Marker genes shown in Fig. 3A of Olsen et al. https://www.biorxiv.org/content/10.1101/2020.05.04.077057v1 - Cells were ordered by pseudotime in the bridge region, this first group of cells were from the mesenchymal part of the trajectory.", "BGN", "COL1A2", "PDGFRA", "LUM", "CALD1", "PRRX1", "FOSB1", "MYC", "KLF4", "MAFF", "FGF1", "TNNT2", "NPNT", "PODXL"),
"Bridge region mesenchymal-SCP transition (Olsen)"=c("As above but for cells in the mesenchymal transitioning to SCP region", "IFITM3", "B2M", "S100A10"),
"Bridge region SCP-like (Olsen)"=c("As above but for the SCP-like region", "MPZ", "SOX10", "S100B", "ERBB3", "CDH19", "FXYD1", "ABCA8"),
"Bridge region SCP-adrenergic transition (Olsen)"=c("As above but for cells in the mesenchymal transitioning to SCP-adrenergic region", "SRP14", "HMGB1", "MDK", "PCBP2"),
"Bridge region adrenergic transition (Olsen)"=c("As above but for cells in the mesenchymal transitioning to adrenergic region", "STMN2", "MLLT11", "RGS5", "PHOX2B"),




## Kildisiute Science advances (science advances).
# Literature curated markers (with actual citations!), are available in their Supplementary Table 2.
"Vascular (Kildisiute)"=c("Vascular markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "PLVAP", "KDR", "PTPRB", "PECAM1"), 

"Erythroblast (Kildisiute)"=c("Erythroblast markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "HBG1", "HBG2", "HBB"), 

"Leukocyte (Kildisiute)"=c("Leukocyes markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "PTPRC"), 

"Adrenal cortex (Kildisiute)"=c("Adrenal cortex markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "STAR", "MC2R"), 

"Mesenchymal (Kildisiute)"=c("Mesenchymal markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "TCF21", "PDGFRB"), 

"Schwann cell precursor (Kildisiute)"=c("SCP markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "SOX10", "MPZ", "PLP1", "ERBB3"), 

"Bridge (Kildisiute)"=c("Bridge cell (cites Furlan 2017) markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx.", "DLL3", "TLX2"), 

"Sympathobasts or chromaffin (Kildisiute)"=c("Sympathoblast and chromaffin markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "TH", "DBH", "PHOX2B", "CHGB"), 

"Sympathoblasts (Kildisiute)"=c("Sympathoblasts  markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "GAP43", "BCL2", "NPY"), 

"Chromaffin (Kildisiute)"=c("Chromaffin markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "PNMT"), 

"Neuroblastoma (Kildisiute)"=c("Neuroblastoma markers obtained from Kildisiute et al, Supplmenentary Table 2, references supporting these genes are provided in Supp Table S2 of Kildisiute et al (PMID 33547074) https://www.science.org/doi/suppl/10.1126/sciadv.abd3311/suppl_file/abd3311_tables_s1_to_s12.xlsx", "PHOX2B", "PHOX2A", "MYCN"),



## Kameneva, Nature Genetics, PMID 33833454
# Fig. 1
"Subepicardial mesenchyme (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception.", "COL3A1", "PRRX1", "TWIST1", "TWIST2", "TBX18", "COL1A1", "COL1A2"),

"Kidney Primordium (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception. ", "PAX2", "LYPD1", "LHX1"),

"Adrenal gland cortex (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception. ", "STAR", "NR5A1", "CYP17A1", "CYP11A1", "CYP21A2"),

"Schwann cell precursors (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception. The authors data suggest these SCPs can give rise to symphathoblasts, which can then give rise to chromaffin cells. However, these SCPs can also directly give rise to chromaffin cells (i.e there are two bridges/paths).", "SOX10", "PLP1", "FOXD3"), 

"Symphathoblasts (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception. Several of these markers are also (more weakly) expressed in chromaffin cells. The authors data suggest Symphthoblasts give rise to chromaffin cells.", "STMN2", "PHOX2B", "ISL1", "PRPH", "ELAVL3", "ELAVL4"),

"Chromaffin (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception. These markers are expressed at a lower level in symphathoblasts.", "CHGA", "PNMT", "TH"),

"Endothelium (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception.", "PECAM1", "KDR", "CAVIN2", "FLT1", "EGFL7", "PRCP"),

"Liver Primordium (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception.", "HNF4A", "AHSG", "ITIH1", "ALDOB", "VTN"),

"HSCs and immune cells (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception.", "RGS10", "FCGR1A", "CD163", "AIF1"),

"Erythroid cells (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception.", "HBB", "ALAS2"),

"Melanocytes (Kameneva)"=c("Marker gene were obtained from Fig. 1D of Kameneva et al (PMID 33833454). These genes were used by the authors to annotate each cell type in their human fetal adrenal scRNA-seq data obtained 6, 8, 9, 11, 12 and 14 weeks post conception.", "MITF", "DCT", "PMEL", "TYR"),

# Fig. 2
"Schwann cell precursor Fig2 (Kamenva)"=c("Marker genes were obtained from Fig. 2H of Kameneva et al (PMID 33833454). These genes were used to show cells transitioning from a SCP to a symphathoblast state in their human fetal single-cell RNA-seq data.", "SOX10", "PLP1", "FOXD3", "FABP7", "S100B", "ERBB3", "NGFR", "MBP", "MPZ", "COL2A1", "POSTN", "MOXD1", "GAS7"),

"Symphathoblasts Fig2 (Kamenva)"=c("Marker genes were obtained from Fig. 2H of Kameneva et al (PMID 33833454) for symphathoblasts (differentiating from SCPs)", "STMN2", "HAND2", "ELAVL4", "STMN4", "ISL1", "PRPH", "ELAVL2", "HMX1"),

"Chromaffin cells Fig2 (Kamenva)"=c("Marker genes were obtained from Fig. 2H of Kameneva et al (PMID 33833454) for chromaffin cells (differentiating from symphathoblasts)", "CHGA", "CHGB", "INSM1", "PENK", "PNMT", "SLC3D3"),




## Marker genes from Hanemaaijer et al (PNAS 2021, PMID 28684471) (Hans Clevers lab) 
## This gene list is from their Table S2 and contains an intersection of their identified marker genes with those identified in Furlan et al.
"Medulla (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for genes defining the broad group labelled adrenal medulla, which included 7 subclusters (SCP, Bridge, Committed Progenitor, N Chromaffin, E Chromaffin, Neuroblast) - see UMAP on their Fig1B for cluster assignments.", "PHOX2A", "CYB561", "CHGA", "CHGB", "DISP2", "SLC18A1", "EML5", "HAND2", "DBH", "UCHL1", "TH", "NNAT", "PCSK1N", "MAP1B", "DDC", "GATA3", "PHOX2B"),

"Schwann cell Precursor (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the SCP subcluster, which is part of the Adrenal Medulla cluster", "ERBB3", "FOXD3", "SERPINE2", "ATP1A2", "FABP7", "PLP1", "CRYAB", "SOX10", "TTYH1", "PTPRZ1", "LMO4", "DAGLA", "NGFR", "CST3", "CHL1", "CNP"), 

"Bridge (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Bridge subcluster, which is part of the Adrenal Medulla cluster.", "NFASC", "MIAT", "TLX2", "CKB", "TBX20", "DLL3", "LDHB", "SOX11", "GSE1", "WDR6", "DPYSL3", "RCC2", "ASCL1", "CDKN1C", "HTR3A"),

"Committed Progenitor (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Committed Progenitor subcluster (seems to mean progenitor of chromaffin cells), which is part of the Adrenal Medulla cluster.", "MEG3", "NRK", "AGTR2", "SYTL4", "NDUFA4L2", "DGKK", "C2CD4B", "RIAN", "UNC5C", "NTRK3", "PARM1", "RGS4", "SLC18A2", "ATP6V1B2", "ZDBF2", "FAM155A", "MIRG", "RGS5", "SYT1"),

"Neuroblast (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Neuroblast subcluster, which is part of the Adrenal Medulla cluster.", "TUBB3", "ELAVL3", "NEFL", "NEFM", "ISL1", "RTN1", "ELAVL4", "CCND1", "STMN2", "BASP1", "STMN1", "PRPH", "INA"),

"Pre-E Chromafin (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Pre Epinepherine Chromaffin subcluster, which is part of the Adrenal Medulla cluster.", "ARHGEF28", "FAM163A", "LRP11", "NELFCD", "PCP4", "SLC31A1", "TLE1", "ARNT2", "SLC38A4", "TMX4", "RPS6KA6", "JUN"),

"N Chromafin (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Norepinepherine Chromaffin subcluster, which is part of the Adrenal Medulla cluster.", "GNAS", "PTPRN", "SNAP25", "PCLO", "SCG3", "CACNA2D1", "PPFIA2", "SLCO3A1", "CXCL14", "SYN2", "SLC35D3", "ADCYAP1R1", "NAP1L5", "SCG5", "LGR5", "C1QL1", "CELF4"),

"E Chromaffin (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Epinepherine Chromaffin subcluster, which is part of the Adrenal Medulla cluster.", "DLK1", "RESP18", "RAB3C", "RAB3B", "GATA2", "NPY", "ANK3"),

"Cortex (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for their broad Cortex group, which is composed of Fetal Zone, Definitive Zone and Adr Primordium subclusters", "KCNK3H", "SPE1", "HSPD1"),

"Adrenal Premordium (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Adrenal Premordium  subcluster, which is part of the Cortex cluster.", "TPI1", "NPM1", "TK1", "MIF", "RPS2"),

"Fetal Zone (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Fetal Zone subcluster, which is part of the Cortex cluster.", "GRAMD1B"),

"Definitive Zone (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the Definitive Zone subcluster, which is part of the Cortex cluster.", "CPE", "CALN1", "GREB1", "CHCHD10"),

"Stroma (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the broad Stroma cluster", "COL5A1", "COL3A1", "FN1", "GPC3", "COL6A1"),

"Stromal 1 (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for Stroma subcluster", "HMGA2", "MMP14", "TGFB2", "CHD3", "COL27A1", "IGSF3", "NFIB"),

"Stromal 2 (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for Stroma subcluster", "POSTN", "COL14A1", "GAS1", "FBN2", "FSTL1", "COL5A2", "IGFBP5", "GSN", "ITM2A", "SFRP1", "PENK", "FBN1"),

"Stromal 3 (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for Stroma subcluster", "NOTCH3", "OLFML3", "PIK3R1"),

"Endothelium (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the broad Endothelium cluster", "SPARCL1"),

"Angioblast (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for Angioblast subcluster of Endothelium cluster", "NES", "MAP4K4", "PREX1"),

"Endothelial subclusters 2-6 (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for Various endothelial subclusters of Endothelium clusters", "GJA1", "DUSP6", "SHROOM2", "SPRY1", "CDK1", "IVNS1ABP", "IGFBP3", "TM4SF1", "JAG1", "FBLN2", "MECOM", "HES1", "MCF2L", "FAM102A", "MARCKSL1", "LFNG", "ADAMTS1", "THBD", "RHOB", "CLIC4"),

"Immune (Hanemaaijer)"=c("Marker genes obtained from Supplementary Table SD of Hanemaaijer et al (PMID 33500353). The authors generated single-cell RNA-seq data (sort-seq, 2,229 cells total) from mouse adrenal glads at E13.5, E14.5, E17.5, E18.5, P1 and P5. These were marker genes that matched with a similar dataset generated by Furlan et al (PMID 28684471). This particular set of markers are for the broad Immune cluster", "LCP2", "TMSB4X"),




## Jansky paper (from Frank Westerman lab) states "Sympathoblasts were identified by Dong et al.28 based on the expression of CARTPT and INSM1, which we found expressed in chromaffin cells. Conversely, the neuronal markers NPY, PRPH, NTRK1 and ISL1 are used to annotate chromaffin cells; these were clearly expressed in neuroblasts in our analyses. Expression of many other sympathetic neuronal marker genes in our neuroblasts supports our cell type annotation, including NEFM, STMN2, SYN3 and ALK27,29,30,31."

## Jansky looked at data up to 17 weeks post conception, which is later than the other papers, and is where the neuroblast cell type seems to emerge properly!!!

## Jansky, Nature Genetics... This paper described a neuroblastic cell type, and made no metion of sympathoblasts (although they do have "Connecting progenitor cells" and Early chromaffin cells, which would appear to be in the right place for sympathoblasts.).

"Hepatocytes (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that ALB+ Hepatocytes were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D. The main text also speculated that these hepatocytes were probably derived from neighboring tissue, and hence artifactual, but there is no evidence presented for this.", "ALB", "AFP", "ASGR1"),

"Muscle progenitor (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that PAX7+ Muscle progenitor were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "PAX7", "CDON", "PAX3"),

"Myofibroblasts (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that ACTA2+ Myofibroblasts were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "ACTA2", "TAGLN", "MYH11"),

"Myocytes (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that MYH3+ Myocytes were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "MYH3", "TTN", "MYH8"),

"Adrenal cortical (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that CYP11A1+ Adrenal cortical cells were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "CYP11A1", "CYP11B1", "NR5A1"),

"Endothelial cells (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that PTPRB+ Endothelial cells were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "PTPRB", "FLT1", "EGFL7"),

"Mesenchymal cells (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that COL1A1+ Mesenchymal were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "COL1A1", "COL1A2", "PDGFRB"),

"Immune cells (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that PTPRC+ Immune cells were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "PTPRC", "CD247", "ITGAM"),

"Erythrocytes (Jansky)"=c("Mentioned in the main text (Jansky et al, Nature Genetics (2021)), page 1, that HBA2+ Erythrocytes were identified in their fetal adrenal glands, this is also highlighted in the UMAP plot on their Fig 1B. Additional genes also shown in their Extended data Figure 2D.", "HBA2", "ANK1", "HBG2"),


# Jansky Fig. 1E, adrenal medulalla markers from their data... (note PNMT and NEFM)
"Schwann cell precursors (Jansky Fig1E)"=c("Heatmap in Fig. 1E, marker genes derived from Jansks scRNA-seq data for adrenal medulla cell types.", "PLP1", "MPZ", "SOX10", "CDH19", "ERBB3"),

"Bridge (Jansky Fig1E)"=c("Heatmap in Fig. 1E, marker genes derived from Jansks scRNA-seq data for adrenal medulla cell types.", "ERBB4", "CDH9", "CTTNBP2", "ASCL1"),

"Chromaffin and connecting progenitor cells (Jansky Fig1E)"=c("Heatmap in Fig. 1E, marker genes derived from Jansks scRNA-seq data for adrenal medulla cell types.", "DBH", "TH", "CHGA", "DDC"),

"Late chromaffin cells (Jansky Fig1E)"=c("Heatmap in Fig. 1E, marker genes derived from Jansks scRNA-seq data for adrenal medulla cell types. Note: this late chromaffin cell signature also includes PNMT, which is an enzyme catalyzing methylation of norepinephrine to form epinephrine; this marker is absent from the ealier less differentiated chromaffin cells (which may have been referred to as sympathoblasts in other datasets, I think)", "DBH", "TH", "CHGA", "DDC", "PNMT"),

"Neuroblasts (Jansky Fig1E)"=c("Heatmap in Fig. 1E, marker genes derived from Jansks scRNA-seq data for adrenal medulla cell types.", "NEFM", "GAP43", "STMN2", "ISL1", "ALK", "SYN3", "IL7"),

"Late Neuroblasts (Jansky Fig1E)"=c("Heatmap in Fig. 1E, marker genes derived from Jansks scRNA-seq data for adrenal medulla cell types. Note: the Jansky Late Neuroblast signature is the same as the earlier Neuroblast signature, but NEFM and ALK are now absent (meaning lower expressed at later stages of differentiation)", "GAP43", "STMN2", "ISL1", "SYN3", "IL7"),











## Cell lines pan cancer signle cell paper (Regev\Tirosh). Kinker et al. (Nature genetics, 2020)
# Cellular programs identified in Figure X. I.e. programs shared across all cell lines.
"Proteasomal degradation (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, this is for proteasomal degradataion", "PSMA3", "PSMA4", "PSME2", "PSMB3", "PSMC2"),

"Protein maturation (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers.", "HSPA5", "HSPA8", "OS9", "RPN2", "PDIA3"),

"Stress response (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers.", "SQSTM1", "DDIT3", "DDIT3", "ATF3", "GADD45A", "GADD45B"),

"EpiSen (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, an epithelial senescence associated program.", "S100AB", "S100A9", "SLPI", "SPRR1B", "LCN2", "CLDN4", "AQP3"),

"p53 Dependent Senescence (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers.", "CDKN1A", "TP53I3", "NEAT1", "TP53TG1"),

"EMT III (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, this EMT program was enriched in non-cycling cells.", "LAMA3", "LAMB3", "JUP", "ITGA2", "COL17A1"),

"IFN Response (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, this program contained interferon response genes.", "ISG15", "ISG20", "IF16", "IF144", "OASL", "IFIT1", "IFIT2", "IFIT3"),

"EMT II (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, this EMT program was enriched in HNSCC cell lines", "VIM", "FN1", "AXL", "COL5A1"),

"EMT I melanoma (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, this EMT program was melanoma specific.", "CYR61", "RRAS", "IL8", "TPM1", "CAV1", "CTGF"),

"Pigmentation (Kinker)"=c("These marker genes were obtained in a pan cancer cell lines NMF analysis of scRNA-eq data in Kinker et al (PMID 33128048) - 10 metaprograms were recovered that manifest across cancers, Skin pigmentation genes, identified in melanoma.", "MITF", "PMEL", "MLANA", "DCT"),


## Cancer-associated fibroblast Nature Cancer review (Lavie et. al. 2022)
# Table of marker genes for CAF states

"Myofibroblastic CAF" = c("These marker genes were curated across cancer subtypes in multiple organ systems as reviewed in Lavie et. al. (PMID 35883004) and contain myofibroblastic specific CAF genes", "ACTA2", "TAGLN", "POSTN", "MYL9", "MYH11", "MYLK", "TPM1", "TPM2", "MMP11", "MMP2", "HOPX", "BGN", "DCN", "IGFBP7", "ITGA7", "LUM", "FN1", "ACTG2", "VCAN", "MEF2C", "COL8A1", "COL15A1", "COL10A1", "COL4A1", "COL14A1", "COL13A1", "COL5A1", "COL5A2", "COL63A", "COL1A1", "COL3A1", "COL11A1", "COL1A2", "COL12A1", "RGS5", "IGFBP3", "THY1", "THBS2", "THBS1", "MYL8", "CNN2", "CNN3", "TNC", "TMEM119", "TGFBR1", "TGFBR2", "TGFB1", "TGFB2", "WNT5A", "PGF", "VEGFA"),

"Inflammatory CAF" = c("These marker genes were curated across cancer subtypes in multiple organ systems as reviewed in Lavie et. al. (PMID 35883004) and contain inflammatory specific CAF genes", "IGF1", "PDGFD", "PDPN", "PDGFRA", "PDGFRB", "IGFBP6", "LY6C1", "IL6", "IL33", "IL8", "IL1R1", "IL1", "IL10", "CXCL1", "CXCL12", "CXCL13", "CXCL1", "CXCL9", "CXCL10", "CXCL2", "CXCL14", "CCL2", "CCL7", "CCL8", "C3", "C4B", "C1S1", "C1S2", "C1S", "C1R", "C7", "CFD", "LIF", "C1QA", "CIQB", "C1QC", "CFB", "SERPING1", "SLP1", "SAA1", "HGF", "CCL21", "CFD", "SOD2", "ADH1B", "GPX3", "RGMA", "SCARA5", "SCARA3", "GMCSF"),

"Antigen Presenting CAF" = c("These marker genes were curated across cancer subtypes in multiple organ systems as reviewed in Lavie et. al. (PMID 35883004) and contain antigen presenting specific CAF genes", "CD74", "S100A4", "S100A3", "S100A10", "S100A13", "SLPI", "SPPI", "CD73", "CCL19", "CCL5", "H2-AA", "H2-AB1", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "CD83", "H2-EB1", "KRT8", "KRT18", "FSP1", "KRT19", "MSLN", "UPK3B", "LRNN4"),

## Myeloid Subsets
# Table of marker genes for MDSCs (PMIC 33526920)

"PMN-MDSC" = c("These marker genes were curated for MDSC subtypes as reviewed in Veglia et. al. (PMID 33526920)", "ITGAM", "CD15", "CD84", "GM-CSF", "VEGF", "IL6", "IL1B", "HIF1A", "TNFRSF10B", "SLC27A2", "LOX1", "CD36", "CD244", "ARG1", "STAT1", "STAT3", "STAT6", "IRF1", "S100A8", "S100A9", "ANXA1", "LYZ2", "PTGS2", "ARG2", "TGFB1", "CSF1", "IL4R"),

"M-MDSC" = c("These marker genes were curated for MDSC subtypes as reviewed in Veglia et. al. (PMID 33526920)", "ITGAM", "CD14", "CD84", "M-CSF", "VEGF", "HIF1A", "TNFRSF10B", "CD36", "CXCR1", "CD274", "IL1B", "IL10", "TGFB1", "S100A9", "S100A8", "ARG1", "ARG2", "NOS2", "VEGFA", "TNF", "STAT3", "IL6"),

# Table of marker genes for Macrophages

"M1 Macrophage" = c("These genes were collated from multiple sources", "IL1B", "IL6","IL33", "TNF", "CD14", "CD80", "ITGAX", "CD86", "IL12", "IL23", "IL18", "CCL2", "CCL3", "CCL4", "CCL5", "CXCL8", "CXCL9", "CXCL10", "CXCL11", "CXCL16", "TLR4"),

"M2 Macrophage" = c("These genes were collated from multiple sources", "TGFB1", "CCL18", "IL10", "IL4", "IL13", "CD14", "CD163", "CD204", "CD206", "VEGF", "VEGFA", "ARG1", "CCl17", "CCL18", "CCL22", "CCL24", "IL2RA", "CXCR1", "CXCR2", "CD23")

)



# From Nat Gen paper, Van Grongingen.
literatureMarkers_largeGeneSets <- list(
"VanGroningen Adrenergic Genes"=c("Adrenergic marker genes from Supplementary Table 2 of Van Groningen et al. Nature Genetics 2017. These genes were identified by differential expression analysis of mesenchymal-like and adrenergic-like neuroblastoma cell lines.", "ABCA3", "ABCB1", "ABLIM1", "ACOT7", "ACTL6B", "ACVR1B", "ADAM22", "ADCYAP1R1", "ADGRB3", "ADRBK2", "AGTPBP1", "AHSA1", "AKAP1", "AKAP12", "ALK", "ANK2", "ANKRD46", "ANP32A", "AP1S2", "ARHGEF7", "ARL6IP1", "ASCL1", "ASRGL1", "ATCAY", "ATL1", "ATP6V0E2", "ATP6V1B2", "AUTS2", "BEND4", "BEX1", "BEX2", "BIRC5", "BMP7", "BMPR1B", "C11orf95", "C14orf132", "C3orf14", "C4orf48", "C7orf55", "CACNA1B", "CACNA2D2", "CADM1", "CAMSAP1", "CCDC167", "CCND1", "CCNI", "CCP110", "CCSAP", "CD200", "CDC42EP3", "CDCA5", "CDKN2C", "CDKN3", "CELF2", "CENPU", "CENPV", "CEP44", "CERK", "CETN3", "CHGA", "CHGB", "CHML", "CHRNA3", "CKB", "CLASP2", "CLGN", "CRH", "CRMP1", "CSE1L", "CXADR", "CXCR4", "CYFIP2", "CYGB", "DACH1", "DAPK1", "DBH", "DCX", "DDC", "DDX39A", "DIABLO", "DKK1", "DLK1", "DNAJB1", "DNAJC6", "DNAJC9", "DNER", "DPYSL2", "DPYSL3", "DPYSL5", "DTD1", "DUSP4", "EEF1A2", "EIF1B", "ELAVL2", "ELAVL3", "ELAVL4", "EML4", "EML6", "ENDOG", "ENO2", "EPB41L4A-AS1", "ESRRG", "EVL", "EXOC5", "EYA1", "FABP6", "FAM107B", "FAM155A", "FAM163A", "FAM167A", "FAM169A", "FAM171B", "FAM60A", "FAXC", "FBLL1", "FBXO8", "FEV", "FHOD3", "FIGNL1", "FKBP1B", "FKBP4", "FOXM1", "FOXO3", "FSD1", "FZD3", "GABRB3", "GAL", "GAP43", "GATA2", "GATA3", "GCH1", "GDAP1", "GDAP1L1", "GDI1", "GDPD1", "GGCT", "GGH", "GLCCI1", "GLDC", "GLRX", "GMNN", "GNB1", "GNG4", "GPR22", "GPR27", "GRB10", "GRIA2", "H1FX", "HAND1", "HAND2-AS1", "HES6", "HEY1", "HK2", "HMGA1", "HMP19", "HN1", "HNRNPA0", "HS6ST2", "ICA1", "IGFBPL1", "IGSF3", "INA", "INO80C", "INSM1", "INSM2", "IRS2", "ISL1", "KDM1A", "KIAA1211", "KIDINS220", "KIF15", "KIF1A", "KIF21A", "KIF2A", "KIF5C", "KLC1", "KLF13", "KLF7", "KLHL13", "KLHL23", "KNSTRN", "L1CAM", "LEPROTL1", "LIN28B", "LINC00888", "LMO3", "LOC100507194", "LOC101928409", "LRRTM2", "LSM3", "LSM4", "LYN", "MAGI3", "MANEAL", "MAP1B", "MAP2", "MAP6", "MAPK8", "MAPT", "MARCH11", "MCM2", "MCM6", "MCM7", "MIAT", "MMD", "MRPL48", "MSH6", "MSI2", "MTCL1", "MXI1", "MYBL2", "MYEF2", "MYO5A", "MYRIP", "NANOS1", "NAP1L5", "NAPB", "NARS2", "NBEA", "NCAM1", "NCAN", "NCOA7", "NCS1", "NEFL", "NEFM", "NELFCD", "NELL2", "NET1", "NFIL3", "NGRN", "NMNAT2", "NNAT", "NOL4", "NPTX2", "NPY", "NRCAM", "NRSN1", "NSG1", "NUDT11", "NUF2", "NUSAP1", "OLA1", "OLFM1", "PARP6", "PBK", "PBX3", "PDK1", "PEG3", "PHF21B", "PHOX2A", "PHOX2B", "PHPT1", "PHYHIPL", "PIK3R1", "PKIA", "PLPPR5", "PNMA2", "POLB", "POPDC3", "PPM1E", "PPP1R9A", "PPP2R3C", "PRC1", "PRCD", "PRIM1", "PRPH", "PRSS12", "PRSS3", "PTS", "QDPR", "RAB33A", "RAB6B", "RALGDS", "RANBP1", "RBBP8", "RBMS3", "RBP1", "REC8", "REEP1", "RET", "RFC4", "RGS17", "RGS5", "RIMBP2", "RIMS3", "RNF144A", "RNF150", "RNF157", "RNF165", "RNFT2", "RPS6KA2", "RRM2", "RTN1", "RTN2", "RUFY3", "RUNDC3A", "RUNDC3B", "SATB1", "SBK1", "SCAMP5", "SCG2", "SCG3", "SCN3A", "SEC11C", "SEPT3", "SEPT6", "SERP2", "SETD7", "SHC3", "SHD", "SIX3", "SLC10A4", "SLC35G2", "SLIT1", "SLIT3", "SNAP25", "SNAP91", "SOX11", "ST3GAL6", "STMN2", "STMN4", "STRA6", "STXBP1", "SV2C", "SYNPO2", "SYT1", "SYT4", "TACC2", "TAGLN3", "TBC1D30", "TBPL1", "TCEAL7", "TDG", "TENM4", "TFAP2B", "TH", "THSD7A", "TIAM1", "TMEM108", "TMEM178B", "TMEM97", "TMOD1", "TMOD2", "TMTC4", "TOX2", "TRAP1", "TSPAN13", "TSPAN7", "TTC8", "TUB", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4B", "UBE2C", "UBE2T", "UCP2", "UNC79", "VRK1", "ZNF195", "ZNF22", "ZNF24", "ZNF512", "ZNF536", "ZNF704", "ZNF711", "ZNF738", "ZNF91", "ZWILCH"),

"VanGroningen Mesenchymal Genes"=c("Mesenchymal marker genes from Supplementary Table 2 of Van Groningen et al. Nature Genetics 2017. These genes were identified by differential expression analysis of mesenchymal-like and adrenergic-like neuroblastoma cell lines.", "A2M", "ABRACL", "ACADVL", "ACAP2", "ACTA2", "ACTN1", "ADAM19", "ADAM9", "ADAMTS5", "ADGRE5", "ADGRG6", "AEBP1", "AJUBA", "ALDH1A3", "AMMECR1", "ANTXR1", "ANXA1", "ANXA2", "ANXA5", "ANXA6", "APOE", "APP", "ARHGAP1", "ARHGEF40", "ARL1", "ARL4A", "ARMCX2", "ARPC1B", "ASPH", "ATP10D", "ATP1B1", "ATP2B1", "ATP2B4", "ATP6V0E1", "ATP8B2", "ATXN1", "B2M", "BAG3", "BGN", "BMP5", "BNC2", "BOC", "BTN3A2", "C1orf198", "C1orf54", "C4orf32", "C6orf120", "CALD1", "CALU", "CAPN2", "CAPN6", "CBFB", "CBLB", "CCDC80", "CD164", "CD44", "CD59", "CD63", "CDH11", "CETN2", "CFH", "CFI", "CILP", "CKAP4", "CLIC4", "CMTM3", "CMTM6", "CNN3", "COL11A1", "COL12A1", "COL1A1", "COL27A1", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COPA", "CPED1", "CPS1", "CRABP2", "CREB3L2", "CREG1", "CRELD2", "CRISPLD1", "CRTAP", "CSRP1", "CTDSP2", "CTNNA1", "CTSB", "CTSC", "CTSO", "CXCL12", "CYBRD1", "CYFIP1", "CYP26A1", "CYR61", "DCAF6", "DDOST", "DDR2", "DESI2", "DKK3", "DLC1", "DLX1", "DLX2", "DMD", "DNAJC1", "DNAJC10", "DNAJC3", "DNM3OS", "DPY19L1", "DSE", "DUSP14", "DUSP5", "DUSP6", "EDEM1", "EDNRA", "EFEMP2", "EGFR", "EGR1", "EGR3", "EHD2", "ELAVL1", "ELF1", "ELK3", "ELK4", "EMILIN1", "EMP1", "ENAH", "EPHA3", "EPS8", "ERBIN", "ERLIN1", "ERRFI1", "ETS1", "EVA1A", "EXT1", "EXTL2", "F2R", "F2RL2", "FAM102B", "FAM114A1", "FAM120A", "FAM129A", "FAM3C", "FAM43A", "FAM46A", "FAT1", "FBN1", "FBN2", "FGFR1", "FIBIN", "FILIP1L", "FKBP14", "FLNA", "FLRT2", "FMOD", "FN1", "FNDC3B", "FSTL1", "FUCA2", "FZD1", "FZD2", "FZD7", "GABRR1", "GALNT10", "GAS1", "GAS2", "GDF15", "GJA1", "GNAI1", "GNG12", "GNS", "GORAB", "GPC6", "GPR137B", "GPX8", "GRN", "GSN", "HES1", "HEXB", "HIBADH", "HIPK3", "HIST1H2AC", "HIST1H2BK", "HLA-A", "HLA-B", "HLA-C", "HLA-F", "HLX", "HNMT", "HOMER1", "HS3ST3A1", "HSP90B1", "HSPA5", "HSPB1", "HTRA1", "HYOU1", "ID1", "ID3", "IFI16", "IFITM2", "IFITM3", "IGF2R", "IGFBP5", "IGFBP6", "IL13RA1", "IL6ST", "INSIG1", "IQGAP2", "ITGA10", "ITGA4", "ITGAV", "ITGB1", "ITM2B", "ITM2C", "ITPR1", "ITPRIPL2", "JAK1", "JAM3", "KANK2", "KCNK2", "KCTD12", "KDELC2", "KDELR2", "KDELR3", "KDM5B", "KIAA1462", "KIF13A", "KIRREL", "KLF10", "KLF4", "KLF6", "L3HYPDH", "LAMB1", "LAMC1", "LAMP1", "LAPTM4A", "LASP1", "LATS2", "LEPROT", "LGALS1", "LHFP", "LHX8", "LIFR", "LIPA", "LITAF", "LIX1L", "LMAN1", "LMNA", "LOXL2", "LPP", "LRP10", "LRRC17", "LRRC8C", "LTBP1", "LUZP1", "MAGT1", "MAML2", "MAN2A1", "MANF", "MBD2", "MBNL1", "MBTPS1", "MEOX1", "MEOX2", "MEST", "MGAT2", "MGP", "MGST1", "MICAL2", "MMP2", "MOB1A", "MRC2", "MXRA5", "MYADM", "MYDGF", "MYL12A", "MYL12B", "MYLIP", "NANS", "NBR1", "NEK7", "NES", "NFIA", "NFIC", "NID1", "NID2", "NOTCH2", "NOTCH2NL", "NPC2", "NPTN", "NQO1", "NR3C1", "NRP1", "OGFRL1", "OLFML2A", "OLFML2B", "OLFML3", "OSTC", "P4HA1", "PALLD", "PAPSS2", "PCDH18", "PCOLCE2", "PCSK5", "PDE3A", "PDE7B", "PDGFC", "PDIA3", "PDIA4", "PDIA6", "PDLIM1", "PEA15", "PEAK1", "PHLDA3", "PHLDB2", "PHTF2", "PIAS3", "PLAGL1", "PLEKHA2", "PLEKHH2", "PLK2", "PLOD2", "PLOD3", "PLPP1", "PLS3", "PLSCR1", "PLSCR4", "PLXDC2", "POLR2L", "PON2", "POSTN", "PPIB", "PPIC", "PPT1", "PRCP", "PRDM6", "PRDX4", "PRDX6", "PROM1", "PRRX1", "PTBP1", "PTGER4", "PTGFRN", "PTN", "PTPN14", "PTPRG", "PTPRK", "PTRF", "PXDC1", "PXDN", "PYGL", "QKI", "QSOX1", "RAB13", "RAB29", "RAB31", "RAP1A", "RAP1B", "RBMS1", "RCN1", "RECK", "REST", "RGL1", "RGS10", "RGS3", "RHOC", "RHOJ", "RIN2", "RIT1", "RNFT1", "RNH1", "ROBO1", "ROR1", "RRBP1", "S1PR3", "SASH1", "SCPEP1", "SCRG1", "SDC2", "SDC4", "SDCBP", "SDF4", "SEC14L1", "SEL1L3", "SEMA3C", "SEMA3F", "SEPT10", "SERPINE2", "SERPINH1", "SFT2D1", "SFT2D2", "SGK1", "SH3BGRL", "SHC1", "SHROOM3", "SIX1", "SIX4", "SKIL", "SLC16A4", "SLC30A1", "SLC30A7", "SLC35F5", "SLC38A2", "SLC38A6", "SLC39A14", "SMAD3", "SNAI2", "SNAP23", "SOSTDC1", "SOX9", "SPARC", "SPARCL1", "SPATA20", "SPCS3", "SPRED1", "SPRY1", "SPRY4", "SPRY4-IT1", "SQSTM1", "SRPX", "SSBP4", "SSR1", "SSR3", "STAT1", "STAT3", "STEAP1", "STK38L", "SUCLG2", "SURF4", "SVIL", "SYDE1", "SYNJ2", "SYPL1", "TCF7L2", "TFE3", "TFPI", "TGFB1I1", "TGFBR2", "THBS1", "TIMP1", "TJP1", "TM4SF1", "TM9SF2", "TMBIM4", "TMED9", "TMEFF2", "TMEM263", "TMEM50A", "TMEM87B", "TNC", "TNFRSF12A", "TNFRSF1A", "TNMD", "TNS1", "TOR1AIP1", "TPBG", "TPM1", "TPM2", "TRAM1", "TRAM2", "TRIL", "TRIM5", "TSC22D2", "TSC22D3", "TSPAN4", "TUBB6", "TWSG1", "TXNDC12", "UAP1", "UGDH", "VCL", "VIM", "WIPI1", "WLS", "WNT5A", "WWTR1", "YAP1", "ZCCHC24", "ZFP36L1", "ZNF217"),

## /addedCode below.

"Descartes adrenocortical markers"=c("Top 50 marker genes of adrenocortical cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/adrenocortical/in/adrenal)", "MSMO1", "MC2R", "CYP11A1", "CYP21A2", "DHCR24", "HMGCR", "LDLR", "DHCR7", "CYP17A1", "HMGCS1", "PAPSS2", "SLC2A14", "TM7SF2", "FRMD5", "FREM2", "CYP11B1", "FDXR", "LINC00473", "SCARB1", "STAR", "SH3BP5", "NPC1", "GRAMD1B", "HSPD1", "JAKMIP2", "GSTA4", "PDE10A", "AC009410.1", "POR", "SLC1A2", "DFNB31", "BAIAP2L1", "SULT2A1", "SLC16A9", "HSPE1", "SCAP", "SH3PXD2B", "SGCZ", "ERN1", "INHA", "FDX1", "RP11-1101K5.1", "CLU", "CYB5B", "PEG3", "IGF1R", "FDPS", "APOC1", "DNER", "LINC01059"),

"Descartes chromaffin markers"=c("Top 50 marker genes of chromaffin cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/chromaffin/in/adrenal)", "GAP43", "MARCH11", "PRPH", "NPY", "HMX1", "HS3ST5", "NTRK1", "PTCHD1", "TMEFF2", "IL7", "RP11-543F8.2", "EYA4", "STMN4", "TMEM132C", "RP11-40F8.2", "RPH3A", "TUBB2B", "TUBB2A", "ALK", "KCNB2", "RP11-161D15.3", "PLXNA4", "LINC00693", "FAT3", "RBFOX1", "GREM1", "MLLT11", "ANKFN1", "CNTFR", "ELAVL2", "REEP1", "CCND1", "SLC44A5", "MAB21L1", "RP11-509E10.1", "GAL", "EYA1", "SLC6A2", "RYR2", "MAP1B", "ISL1", "BASP1", "MAB21L2", "TUBA1A", "RP11-420N3.2", "CNKSR2", "RGMB", "STMN2", "SYNPO2", "EPHA6"),

"Descartes Vascular_endothelial markers"=c("Top 50 marker genes of Vascular_endothelial cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/vascular_endothelial/in/adrenal)", "GALNT15", "CALCRL", "FLT4", "DNASE1L3", "SLCO2A1", "FCGR2B", "F8", "TM4SF18", "IRX3", "CXORF36", "CYP26B1", "TMEM88", "KDR", "ROBO4", "PODXL", "NR5A2", "ESM1", "NOTCH4", "CEACAM1", "ELTD1", "SOX18", "AC011526.1", "PLVAP", "ECSCR", "BTNL9", "NPR1", "CDH5", "KANK3", "CRHBP", "PTPRB", "SHANK3", "RP11-768F21.1", "TEK", "CHRM3", "ARHGAP29", "RAMP2", "PPAP2B", "MMRN2", "RASIP1", "APLNR", "EHD3", "C8ORF4", "MYRIP", "TIE1", "EFNB2", "ID1", "SHE", "HYAL2", "CDH13", "CLDN5"),

"Descartes stromal markers"=c("Top 50 marker genes of stromal cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/stromal/in/adrenal)", "SULT1E1", "HHIP", "MGP", "SFRP2", "DCN", "RP11-383H13.1", "ITGA11", "OGN", "DKK2", "RSPO3", "COL3A1", "LUM", "CCDC102B", "SCARA5", "PRRX1", "ABCC9", "PCOLCE", "LAMC3", "LRRC17", "PDGFRA", "ZNF385D", "ADAMTSL3", "LOX", "CLDN11", "ABCA6", "COL6A3", "EDNRA", "PAMR1", "C7", "BICC1", "COL1A1", "MXRA5", "GLI2", "COL12A1", "GAS2", "POSTN", "CD248", "PCDH18", "IGFBP3", "FREM1", "CCDC80", "ACTA2", "CDH11", "ISLR", "PRICKLE1", "FNDC1", "ELN", "COL1A2", "ADAMTS2", "COL27A1"),

"Descartes sympathoblasts markers"=c("Top 50 marker genes of sympathoblasts cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/sympathoblasts/in/adrenal)", "SLC24A2", "PENK", "ARC", "SLC35F3", "ST18", "KSR2", "CHGA", "CNTN3", "CCSER1", "CNTNAP5", "PNMT", "PACRG", "CDH18", "SORCS3", "AGBL4-IT1", "GALNTL6", "AGBL4", "PCSK2", "KCTD16", "DGKK", "LINC00632", "GRM7", "RP11-176N18.2", "C1QL1", "SPOCK3", "KIAA1244", "TBX20", "FAM155A-IT1", "FAM155A", "NTNG1", "LAMA3", "SLC18A1", "CDH12", "INSM1", "GCH1", "ROBO1", "TENM1", "HTATSF1", "MGAT4C", "EML6", "CHGB", "TMEM130", "PTPLAD1", "FGF14-IT1", "FGF14", "SCG2", "UNC80", "PCSK1N", "GRID2", "TIAM1"),

"Descartes erythroblasts markers"=c("Top 50 marker genes of erythroblasts cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/erythroblasts/in/adrenal)", "C17ORF99", "SLC4A1", "AHSP", "GYPA", "HBM", "EPB42", "SPTA1", "GYPB", "ALAS2", "TMEM56", "HBZ", "HBA1", "HBG2", "HBA2", "HBB", "GYPE", "SLC25A21", "HBG1", "SLC25A37", "TMCC2", "SELENBP1", "RHCE", "FECH", "CPOX", "ABCB10", "TRAK2", "XPO7", "TFR2", "ANK1", "TSPAN5", "CAT", "SPTB", "RHD", "HEMGN", "HECTD4", "BLVRB", "RP5-964N17.1", "SPECC1", "SOX6", "GCLC", "EPB41", "DENND4A", "RGS6", "GYPC", "MICAL2", "RHAG", "SNCA", "RAPGEF2", "CR1L", "MARCH3"),

"Descartes myeloid markers"=c("Top 50 marker genes of myeloid cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/myeloid/in/adrenal)", "MS4A7", "MS4A4A", "CD163", "CD14", "VSIG4", "C1QB", "C1QA", "C1QC", "MS4A6A", "CD163L1", "MS4A4E", "SLCO2B1", "FGD2", "ADAP2", "SPP1", "CYBB", "MPEG1", "LGMN", "MSR1", "HCK", "SLC1A3", "CSF1R", "CPVL", "IFNGR1", "HLA-DPA1", "CST3", "MERTK", "HLA-DRB1", "HRH1", "HLA-DRA", "FMN1", "RBPJ", "CTSB", "MARCH1", "WWP1", "FGL2", "CD74", "CTSS", "TGFBI", "ATP8B4", "CTSD", "RGL1", "AXL", "RNASE1", "ITPR2", "SLC9A9", "PTPRE", "CTSC", "SFMBT2", "ABCA1"),

"Descartes Schwann markers"=c("Top 50 marker genes of Schwann cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/schwann/in/adrenal)", "MPZ", "SOX10", "CDH19", "OLFML2A", "PLP1", "SLC35F1", "PTPRZ1", "TRPM3", "ERBB3", "XKR4", "PPP2R2B", "GRIK3", "COL18A1", "SCN7A", "NRXN1", "SFRP1", "ABCA8", "KCTD12", "LRRTM4", "SOX5", "NRXN3", "VIM", "DST", "MDGA2", "NLGN4X", "PMP22", "COL25A1", "IL1RAPL1", "ADAMTS5", "FIGN", "EGFLAM", "ZNF536", "EDNRB", "LAMA4", "LAMC1", "PTN", "MARCKS", "PAG1", "VCAN", "FAM134B", "IL1RAPL2", "SORCS1", "ERBB4", "STARD13", "GAS7", "PLCE1", "LAMB1", "COL5A2", "HMGA2", "GFRA3"),

"Descartes Megakaryocytes markers"=c("Top 50 marker genes of Megakaryocytes cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/megakaryocytes/in/adrenal)", "GP9", "PF4", "PPBP", "GP1BA", "ITGA2B", "ITGB3", "RP11-556I14.2", "P2RX1", "THBS1", "PLEK", "FERMT3", "ZYX", "TUBB1", "MMRN1", "LTBP1", "ARHGAP6", "ACTN1", "LIMS1", "TLN1", "RAB27B", "TGFB1", "FLNA", "VCL", "RAP1B", "PSTPIP2", "SLC24A3", "MED12L", "PRKAR2B", "PDE3A", "ANGPT1", "BIN2", "MYH9", "DOK6", "ACTB", "MCTP1", "CD84", "SLC2A3", "TRPC6", "TMSB4X", "GSN", "TPM4", "INPP4B", "STOM", "FLI1", "CD9", "SPN", "MYLK", "UBASH3B", "STON2", "HIPK2"),


"Descartes Lyphoid markers"=c("Top 50 marker genes of Lyphoid cells in the Decartes fetal adrenal single cell map (https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/cell/lymphoid/in/adrenal)", "SKAP1", "BACH2", "FAM65B", "ARHGAP15", "PTPRC", "ETS1", "ANKRD44", "MCTP2", "LCP1", "HLA-B", "PRKCH", "HLA-C", "IKZF1", "HLA-A", "LINC00299", "FYN", "SCML4", "FOXP1", "MSN", "PDE3B", "SP100", "GNG2", "LEF1", "MBNL1", "SAMD3", "WIPF1", "NCALD", "SORL1", "CELF2", "PITPNC1", "CCND3", "ARID5B", "DOCK10", "B2M", "PLEKHA2", "CD44", "RCSD1", "RAP1GAP2", "ABLIM1", "BCL2", "TOX", "STK39", "TMSB10", "IFI16", "ARHGDIB", "CCL5", "NKG7", "EVL", "ITPKB", "ANKRD44-IT1")

)



# this set of literature marker genes (and learned markers) was obtained from cellTypist (https://www.celltypist.org/encyclopedia), from this specific excel spreadsheet: (https://github.com/Teichlab/celltypist_wiki/blob/main/atlases/Pan_Immune_CellTypist/encyclopedia/encyclopedia_table.xlsx). 
literatureMarkers_celltypistImmuneGeneSets <- list(

# Curated markers
"B cells: B cells (curated markers)"=c("B lymphocytes with diverse cell surface immunoglobulin receptors recognising specific antigenic epitopes and mediating humoral immunity", "CD79A", "MS4A1", "CD19"), 

"B cells: Follicular B cells (curated markers)"=c("resting mature B lymphocytes found in the primary and secondary lymphoid follicles and participating in T cell-dependent immune responses", "CXCR5", "TNFRSF13B", "CD22"), 

"B cells: Germinal center B cells (curated markers)"=c("proliferating mature B cells that undergo somatic hypermutation and class-switch recombination in secondary lymphoid organs", "POU2AF1", "CD40", "SUGCT"), 

"B cells: Memory B cells (curated markers)"=c("long-lived mature B lymphocytes which are formed within germinal centers following primary infection and selected for higher-affinity immunoglobulin", "CR2", "CD27", "MS4A1"), 

"B cells: Naive B cells (curated markers)"=c("mature B lymphocytes which express cell-surface IgM and IgD and have not been exposed to/activated by antigens", "IGHM", "IGHD", "TCL1A"), 

"B cells: Transitional B cells (curated markers)"=c("immature B cell precursors in the bone marrow which connect Pre-B cells with mature naive B cells and are subject to the process of B cell selection", "CD24", "MYO1C", "MS4A1"), 

"B-cell lineage: Large pre-B cells (curated markers)"=c("proliferative B lymphocyte precursors derived from Pro-B cells and expressing membrane μ chains with surrogate light chains in their receptors", "MME", "CD24", "MKI67"), 

"B-cell lineage: Small pre-B cells (curated markers)"=c("non-proliferative B lymphocyte precursors derived from Pro-B cells and expressing membrane μ chains with surrogate light chains in their receptors", "MME", "CD24", "IGLL5"), 

"B-cell lineage: Pre-pro-B cells (curated markers)"=c("the earliest primitive B cell progenitors which express CD34 and SPINK2", "have germline Ig genes and give rise to Pro-B cells", "IL7R", "ZCCHC7", "RAG1"), 

"B-cell lineage: Pro-B cells (curated markers)"=c("early B lymphocyte progenitors undergoing D-J joining on the H chain chromosome and joining of a V segment to the rearranged D-J", "MME", "DNTT", "IGLL1"), 

"Cycling cells: Cycling B cells (curated markers)"=c("proliferating B lymphocytes", "MKI67", "TOP2A", "CD19"), 

"Cycling cells: Cycling DCs (curated markers)"=c("proliferating dendritic cells", "MKI67", "TOP2A", "CLEC10A"), 

"Cycling cells: Cycling gamma-delta T cells (curated markers)"=c("proliferating gamma-delta (γδ) T lymphocytes", "MKI67", "TOP2A", "TRDC"), 

"Cycling cells: Cycling monocytes (curated markers)"=c("proliferating monocytes", "MKI67", "TOP2A", "S100A9"), 

"Cycling cells: Cycling NK cells (curated markers)"=c("proliferating natural killer cells", "MKI67", "TOP2A", "GNLY"), 

"Cycling cells: Cycling T cells (curated markers)"=c("proliferating T lymphocytes", "MKI67", "TOP2A", "CD3D"), 

"DC: DC (curated markers)"=c("dendrite-shaped cells which process and present antigens to naïve T cells and induce adaptive immune responses ", "CD1C", "FCER1A", "CLEC10A"), 

"DC: DC1 (curated markers)"=c("conventional type 1 dendritic cells which constitute a rare dendritic cell population with superior cross-presentation ability", "BATF3", "CADM1", "CLEC9A"), 

"DC: DC2 (curated markers)"=c("conventional type 2 dendritic cells which constitute the major subset of myeloid dendritic cells", "CLEC10A", "FCER1A", "CD1C"), 

"DC: DC3 (curated markers)"=c("a dendritic cell subtype found in bone marrow which shares features with monocytes and conventional type 2 dendritic cells", "CLEC10A", "FCER1A", "S100A9"), 

"DC: Migratory DCs (curated markers)"=c("migratory dendritic cells that transport antigens to the draining lymph nodes during both homeostatic conditions and infections", "EBI3", "CCR7", "CCL19"), 

"DC: Transitional DC (curated markers)"=c("immature dendritic cells which follow dendritic cell precursors and are committed to the mature dendritic cells", "CLEC10A", "KLF4", "AXL"), 

"DC precursor: DC precursor (curated markers)"=c("dendritic cell precursors that differentiate into dendritic cells", "IRF8", "CLEC10A", "PCLAF"), 

"Double-negative thymocytes: Double-negative thymocytes (curated markers)"=c("the early immature thymocytes from the thymus in the double negative (co-receptors CD4- and CD8-) stage", "FXYD2", "HES1", "CD99"), 

"Double-positive thymocytes: Double-positive thymocytes (curated markers)"=c("immature thymocytes capable of binding MHC class I or II in the double positive (co-receptors CD4+ and CD8+) stage", "CD1A", "CD8A", "SMPD3"), 

"Early MK: Early MK (curated markers)"=c("early megakaryocytes which are committed to megakaryocyte lineage cells", "PF4", "CMTM5", "SELP"), 

"Endothelial cells: Endothelial cells (curated markers)"=c("single-layered lining cells that constitute the interior vasculature", "form part of the microenvironment and fulfill different functions", "CLDN5", "PLVAP", "SPARCL1"), 

"Epithelial cells: Epithelial cells (curated markers)"=c("highly specialised cells which are located in the outer layer of stroma and form the cellular components of the epithelium", "KRT15", "KRT17", "CCL19"), 

"Erythrocytes: Erythrocytes (curated markers)"=c("biconcave enucleated red blood cells filled with hemoglobin to transport oxygen and carbon dioxide between the lungs and tissues", "HBM", "HBQ1", "GYPA"), 

"Erythroid: Early erythroid (curated markers)"=c("early erythroid cells which are committed to erythroid lineage cells", "APOC1", "TESPA1", "GATA2"), 

"Erythroid: Late erythroid (curated markers)"=c("late erythroid cells following early and middle erythroid cells", "GYPA", "GYPB", "HBA1"), 

"Erythroid: Mid erythroid (curated markers)"=c("middle erythroid cells in fetal liver", "blood and bone marrow following early erythroid cells and maturing into late erythroid cells", "PRDX2", "KCNH2", "GYPA"), 

"ETP: ETP (curated markers)"=c("early thymic progenitors migrated to the thymus before turning into double-negative thymocytes or other cells", "ACY3", "SPINK2", "CD34"), 

"Fibroblasts: Fibroblasts (curated markers)"=c("the most common cells of connective tissues which synthesize the extracellular matrix and collagen to maintain tissue homeostasis", "COL1A1", "COL1A2", "DCN"), 

"Granulocytes: Granulocytes (curated markers)"=c("cells of innate immune system with specific granules in the cytoplasm including neutrophils", "eosinophils and basophils", "KIT", "CPA3", "SLC45A3"), 

"Granulocytes: Neutrophils (curated markers)"=c("the most abundant type of granulocytes that contains distinctive cytoplasmic granules and forms an essential part of the innate immune system", "LCN2", "ORM1", "MMP8"), 

"HSC/MPP: CMP (curated markers)"=c("highly proliferative common myeloid progenitors which later give rise to granulocyte-monocyte progenitors and megakaryocyte-erythroid progenitors", "MPO", "CTSG", "FLT3"), 

"HSC/MPP: Early lymphoid/T lymphoid (curated markers)"=c("early lymphoid/T lymphocytes with lymphocyte potential in the fetal liver before T cells emerged from the thymus", "ETS1", "IL32", "GATA3"), 

"HSC/MPP: ELP (curated markers)"=c("early lymphoid progenitors which originate from hematopoietic stem cells of the bone marrow and migrate to the thymus", "FLT3", "LTB", "RUNX2"), 

"HSC/MPP: GMP (curated markers)"=c("hematopoietic granulocyte-monocyte progenitors that are committed to the granulocyte and monocyte lineage cells", "ELANE", "MPO", "PRTN3"), 

"HSC/MPP: HSC/MPP (curated markers)"=c("hematopoietic stem cells and multipotent progenitor cells with the potential of differentiating into different blood cells", "CD34", "SPINK2", "CRHBP"), 

"HSC/MPP: Megakaryocyte-erythroid-mast cell progenitor (curated markers)"=c("shared progenitors in the fetal liver which originate from common myeloid progenitors and differentiate into megakaryocytes", "erythroid cells or mast cells", "GATA2", "FCER1A", "GATA1"), 

"HSC/MPP: MEMP (curated markers)"=c("shared progenitors which are derived from common myeloid progenitors and will differentiate into megakaryocytes", "erythroid cells or mast cells", "GATA2", "ITGA2B", "GATA1"), 

"HSC/MPP: Neutrophil-myeloid progenitor (curated markers)"=c("progenitors of neutrophils and myeloid cells which are transitioned from hematopoietic stem cells and multipotent progenitors", "LYZ", "MPO", "SERPINB1"), 

"ILC: ILC (curated markers)"=c("specialised innate immune cells from the lymphoid lineage but without antigen-specific T-cell receptors on the surface", "S100A13", "TLE1", "AREG"), 

"ILC: ILC1 (curated markers)"=c("innate lymphoid cell subpopulation I that is non-cytotoxic and has overlapping phenotypes and functions with natural killer cells", "CXCR3", "CD3D", "IKZF3"), 

"ILC: ILC2 (curated markers)"=c("innate lymphoid cell subpopulation II that promotes type 2 inflammation and is involved in immune response against large extracellular pathogens", "GATA3", "KLRG1", "HPGDS"), 

"ILC: ILC3 (curated markers)"=c("innate lymphoid cell subpopulation III that is required for host defense against specific extracellular bacteria and fungi", "IL4I1", "RORC", "KIT"), 

"ILC: CD16+ NK cells (curated markers)"=c("CD16+ granular lymphocytes that play protective roles against both infectious pathogens and cancer using antibody-dependent cell-mediated cytotoxicity", "GNLY", "FCGR3A", "NKG7"), 

"ILC: CD16- NK cells (curated markers)"=c("CD16- granular lymphocytes that play protective roles against both infectious pathogens and cancer using antibody-dependent cell-mediated cytotoxicity", "GNLY", "CD160", "NKG7"), 

"ILC: NK cells (curated markers)"=c("granular lymphocytes that play protective roles against both infectious pathogens and cancer using antibody-dependent cell-mediated cytotoxicity", "GNLY", "XCL2", "NKG7"), 

"ILC: Transitional NK (curated markers)"=c("immature natural killer cells which originate from natural killer cell precursors and are committed to mature natural killer cells", "CCL5", "GZMK", "FCGR3A"), 

"ILC precursor: ILC precursor (curated markers)"=c("innate lymphoid cell precursors which give rise to innate lymphoid cells and lack the characteristics of their mature progenies", "LST1", "HPN", "SCN1B"), 

"Macrophages: Hofbauer cells (curated markers)"=c("primitive placental resident macrophages with granules and vacuoles found in placenta particularly during early pregnancy", "LYVE1", "SCIN", "F13A1"), 

"Macrophages: Kidney-resident macrophages (curated markers)"=c("long-lived macrophages resident in the kidney under non-inflammatory conditions and maintaining homeostasis and resolving inflammation ", "SEPP1", "CSF2RA", "CD74"), 

"Macrophages: Kupffer cells (curated markers)"=c("resident macrophages in liver under non-inflammatory conditions which line the hepatic sinusoids and are involved in erythrocyte clearance", "CD5L", "VCAM1", "CETP"), 

"Macrophages: Macrophages (curated markers)"=c("specialised mononuclear phagocytic cells which recognise", "engulf and degrade cellular debris and pathogens in various tissues", "C1QA", "CD68", "TREM2"), 

"Mast cells: Mast cells (curated markers)"=c("long-lived innate migrant cells found in most tissues with many large basophilic granules rich in histamine and heparin", "TPSAB1", "TPSB2", "CPA3"), 

"Megakaryocyte precursor: Megakaryocyte precursor (curated markers)"=c("megakaryocyte precursors in the bone marrow that are committed to megakaryocytes", "GATA2", "PRSS57", "CAVIN1"), 

"Megakaryocytes/platelets: Megakaryocytes/platelets (curated markers)"=c("large multinucleated cells with numerous azurophilic granules and with mature blood platelets released from the cytoplasm", "CMTM5", "ITGA2B", "PF4"), 

"MNP: MNP (curated markers)"=c("mononuclear phagocytes including dendritic cells", "monocytes and macrophages that have in common the property of phagocytosis", "HLA-DPA1", "HLA-DQA1", "FGL2"), 

"Mono-mac: Mono-mac (curated markers)"=c("mononuclear phagocytes including monocytes and macrophages that are involved in the regulation of innate and adaptive immunity", "TYROBP", "C1QC", "HMOX1"), 

"Monocyte precursor: Monocyte precursor (curated markers)"=c("monocyte precursors that are committed to the monocytes", "LYZ", "VCAN", "S100A9"), 

"Monocytes: Classical monocytes (curated markers)"=c("CD14+ myeloid mononuclear recirculating leukocytes that are capable of differentiating into macrophages and myeloid lineage dendritic cells", "S100A9", "CD14", "S100A12"), 

"Monocytes: Non-classical monocytes (curated markers)"=c("CD16+ myeloid mononuclear recirculating leukocytes that are capable of differentiating into macrophages and myeloid lineage dendritic cells", "FCGR3A", "C1QA", "CX3CR1"), 

"Monocytes: Monocytes (curated markers)"=c("myeloid mononuclear recirculating leukocytes that are capable of differentiating into macrophages and myeloid lineage dendritic cells", "S100A9", "LYZ", "FCN1"), 

"Myelocytes: Myelocytes (curated markers)"=c("early granulocyte precursors that are derive from promyelocytes and later mature into metamyelocytes", "S100A8", "FTL", "CSTA"), 

"pDC: pDC (curated markers)"=c("rare plasmacytoid dendritic cell subpopulation which serves as the major source of type I interferons when the body is infected by a virus", "IL3RA", "LILRA4", "PLD4"), 

"pDC precursor: pDC precursor (curated markers)"=c("precursors of plasmacytoid dendritic cells which have intermixed lymphoid and myeloid origins and give rise to plasmacytoid dendritic cells", "IL3RA", "CCDC50", "IRF7"), 

"Plasma cells: Plasma cells (curated markers)"=c("B-lymphocyte white blood cells capable of secreting large quantities of immunoglobulins or antibodies", "JCHAIN", "MZB1", "XBP1"), 

"Promyelocytes: Promyelocytes (curated markers)"=c("early granulocyte progenitors in the bone marrow which are derived from the myeloblasts and develop into the myelocytes", "CTSG", "MPO", "ELANE"), 

"T cells: CD8a/a (curated markers)"=c("unconventional T lymphocytes with CD8 alpha/alpha homodimers binding MHC molecules", "ZNF683", "GNG4", "PDCD1"), 

"T cells: CD8a/b(entry) (curated markers)"=c("T lymphocytes with CD8 alpha/beta heterodimers in the late double-positive (entry) stage during T cell development", "TOX2", "SATB1", "CCR9"), 

"T cells: Follicular helper T cells (curated markers)"=c("specialised CD4+ T cell subpopulation found in the periphery within B cell follicles and stimulating the cognate follicular B cells for class-switching", "PDCD1", "ICOS", "CXCR5"), 

"T cells: gamma-delta T cells (curated markers)"=c("unconventional T lymphocyte subpopulation expressing a gamma-delta T cell receptor complex on the surface to recognise antigens", "TRDC", "TRGC1", "CCL5"), 

"T cells: MAIT cells (curated markers)"=c("mucosal-associated invariant T cells which have semi-invariant T-cell receptors and are restricted by the MHC I-like molecule MR1", "KLRB1", "SLC4A10", "TRAV1-2"), 

"T cells: Memory CD4+ cytotoxic T cells (curated markers)"=c("CD4+ memory T cells which have cytotoxic activities by secreting granzymes and perforin and by killing target cells on the basis of MHC class II", "GZMK", "CD4", "IL10"), 

"T cells: NKT cells (curated markers)"=c("distinct T lymphocyte subpopulation that expresses both alpha/beta T-cell receptors and surface receptors of natural killer cells", "NKG7", "GNLY", "CD8A"), 

"T cells: Regulatory T cells (curated markers)"=c("T cell subpopulation which modulates immune responses and regulates other cells through direct cell-cell contact and cytokine release", "CTLA4", "IL2RA", "FOXP3"), 

"T cells: T(agonist) (curated markers)"=c("unconventional T cell subpopulation in the thymus which expresses MIR155HG and shares some signatures with differentiating regulatory T cells", "MIR155HG", "BIRC3", "SMS"), 

"T cells: Tcm/Naive cytotoxic T cells (curated markers)"=c("CD8+ cytotoxic T lymphocytes mainly localized in secondary lymphoid tissues and sustaining the responses by proliferating and producing new effectors", "CD8A", "CCR7", "SELL"), 

"T cells: Tcm/Naive helper T cells (curated markers)"=c("CD4+ helper T lymphocytes mainly localized in secondary lymphoid tissues and sustaining the responses by proliferating and producing new effectors", "CCR7", "SELL", "CD4"), 

"T cells: Tem/Trm cytotoxic T cells (curated markers)"=c("CD8+ cytotoxic T lymphocytes mainly localized in lymphoid and peripheral tissues and presenting an immediate", "but not sustained", "immune defense", "GZMK", "CD8A", "CCL5"), 

"T cells: Tem/Temra cytotoxic T cells (curated markers)"=c("terminally differentiated CD8+ cytotoxic T lymphocytes with effector memory phenotypes that re-express CD45RA ", "CX3CR1", "GZMB", "GNLY"), 

"T cells: Trm cytotoxic T cells (curated markers)"=c("tissue resident CD8+ cytotoxic T lymphocytes mainly localized in epithelial tissues such as gut", "ITGA1", "ITGAE", "CXCR6"), 

"T cells: Tem/Effector helper T cells (curated markers)"=c("CD4+ helper T lymphocytes mainly localized in lymphoid and peripheral tissues and presenting an immediate", "but not sustained", "immune defense", "KLRB1", "AQP3", "ITGB1"), 

"T cells: Tem/Effector helper T cells PD1+ (curated markers)"=c("CD4+ helper T lymphocyte subpopulation in the thymus which features the expression of programmed cell death protein 1 (PD-1)", "PDCD1", "CD4", "CTLA4"), 

"T cells: Treg(diff) (curated markers)"=c("unconventional T lymphocyte subpopulation in the thymus which connects αβ T cells and canonical regulatory T cells and is still differentiating", "CD27", "CCR7", "IKZF4"), 

"T cells: Type 1 helper T cells (curated markers)"=c("CD4+ helper T lymphocyte subpopulation which is capable of producing interferon-gamma and modulating cell-mediated immune responses", "CCL5", "CXCR3", "TBX21"), 

"T cells: Type 17 helper T cells (curated markers)"=c("CD4+ helper T lymphocyte subpopulation which is capable of producing interleukin 17 (IL-17) and mediating protective immunity and autoimmunity", "IL7R", "CCR6", "ZBTB16"),


# Model markers
"B cells: B cells (model markers)"=c("B lymphocytes with diverse cell surface immunoglobulin receptors recognising specific antigenic epitopes and mediating humoral immunity", "HBG2", "HBA2", "HLA-DQA2", "HLA-DRB1", "HLA-DRA", "MALAT1", "CD74", "CD79A", "HBA1", "HSPA1B"), 

"B cells: Follicular B cells (model markers)"=c("resting mature B lymphocytes found in the primary and secondary lymphoid follicles and participating in T cell-dependent immune responses", "MT-RNR2", "IGHA1", "JCHAIN", "GAS5", "RPL41", "CD83", "OR2A25", "C11orf72", "VPREB3", "RP11-705C15.5"), 

"B cells: Germinal center B cells (model markers)"=c("proliferating mature B cells that undergo somatic hypermutation and class-switch recombination in secondary lymphoid organs", "SOST", "SERPINA9", "GLYATL2", "AICDA", "REG1A", "IFIT1B", "LINC01644", "BMP3", "RP11-589C21.6", "IL17A"), 

"B cells: Memory B cells (model markers)"=c("long-lived mature B lymphocytes which are formed within germinal centers following primary infection and selected for higher-affinity immunoglobulin", "RPS17", "RACK1", "CD74", "MIR1244-2", "HLA-DRA", "IGKC", "HES1", "MT-ATP8", "MT-ND3", "IGHA2"), 

"B cells: Naive B cells (model markers)"=c("mature B lymphocytes which express cell-surface IgM and IgD and have not been exposed to/activated by antigens", "TCL1A", "IGHD", "CD74", "IGHM", "YBX3", "RPS17", "IGKC", "RPL41", "RACK1", "IL4R"), 

"B cells: Transitional B cells (model markers)"=c("immature B cell precursors in the bone marrow which connect Pre-B cells with mature naive B cells and are subject to the process of B cell selection", "MIR1244-2", "GNG3", "C11orf72", "LINC01644", "LINC02206", "IFIT1B", "OR2A25", "GLYATL2", "IGHV5-78", "LINC01709"), 

"B-cell lineage: Large pre-B cells (model markers)"=c("proliferative B lymphocyte precursors derived from Pro-B cells and expressing membrane μ chains with surrogate light chains in their receptors", "IGLL1", "VPREB1", "LINC01483", "CMA1", "LINC01644", "RP11-589C21.6", "CTAG2", "ZFHX4-AS1", "IL1RAPL1", "GPIHBP1"), 

"B-cell lineage: Small pre-B cells (model markers)"=c("non-proliferative B lymphocyte precursors derived from Pro-B cells and expressing membrane μ chains with surrogate light chains in their receptors", "CD79B", "OR2A25", "MIR1244-2", "FCGR2C", "IL22RA2", "TLDC2", "AICDA", "CD207", "IGHE", "LINC02206"), 

"B-cell lineage: Pre-pro-B cells (model markers)"=c("the earliest primitive B cell progenitors which express CD34 and SPINK2", "have germline Ig genes and give rise to Pro-B cells", "ZFHX4-AS1", "VPREB1", "FCGR2C", "SIGLEC17P", "RP11-84C10.2", "LINC02202", "C11orf72", "HLA-DPB2", "OR2A25", "LINC01709"), 

"B-cell lineage: Pro-B cells (model markers)"=c("early B lymphocyte progenitors undergoing D-J joining on the H chain chromosome and joining of a V segment to the rearranged D-J", "VPREB1", "DNTT", "XCR1", "LINC01013", "SOCS2-AS1", "RP11-701P16.2", "IGHE", "RP11-84C10.2", "RP11-589C21.6", "OR2A25"), 

"Cycling cells: Cycling B cells (model markers)"=c("proliferating B lymphocytes", "IGHV5-78", "RGS13", "LINC01709", "CD207", "KIAA0087", "C11orf72", "GRIN1", "OR2A25", "LCNL1", "LINC01644"), 

"Cycling cells: Cycling DCs (model markers)"=c("proliferating dendritic cells", "SIGLEC17P", "MMP12", "TCL6", "IGHE", "IGHV5-78", "SNHG26", "AP001059.6", "RP11-701P16.2", "SNX29P1", "DRAIC"), 

"Cycling cells: Cycling gamma-delta T cells (model markers)"=c("proliferating gamma-delta (γδ) T lymphocytes", "OR2A25", "DRAIC", "LCNL1", "C11orf72", "FCGR2C", "SYCP1", "LINC01644", "SIGLEC17P", "LINC02202", "PMCH"), 

"Cycling cells: Cycling monocytes (model markers)"=c("proliferating monocytes", "FCGR2C", "RP11-589C21.6", "AJ271736.10", "IL22RA2", "LINC02202", "OR2A25", "LINC01644", "TCL1B", "LDLRAD2", "IL22"), 

"Cycling cells: Cycling NK cells (model markers)"=c("proliferating natural killer cells", "PCLAF", "IGKC", "KL", "IGHV5-78", "C11orf72", "L3MBTL4-AS1", "LINC01644", "KIAA0087", "GRIN1", "IGHE"), 

"Cycling cells: Cycling T cells (model markers)"=c("proliferating T lymphocytes", "PCLAF", "CD207", "DUSP26", "KIAA0087", "IGHV5-78", "CALB2", "SOST", "HLA-DPB2", "SIGLEC17P", "IL22RA2"), 

"DC: DC (model markers)"=c("dendrite-shaped cells which process and present antigens to naïve T cells and induce adaptive immune responses ", "RGS1", "FCER1A", "RP11-1275H24.3", "LINC01644", "SOCS2-AS1", "OR2A25", "IGHV5-78", "SIGLEC17P", "LCNL1", "SYCP1"), 

"DC: DC1 (model markers)"=c("conventional type 1 dendritic cells which constitute a rare dendritic cell population with superior cross-presentation ability", "CLEC9A", "WFDC21P", "SYCP1", "HLA-DPB2", "FOXH1", "RPS3AP34", "OSTN-AS1", "RP11-701P16.2", "LINC01644", "IL1RAPL1"), 

"DC: DC2 (model markers)"=c("conventional type 2 dendritic cells which constitute the major subset of myeloid dendritic cells", "CST3", "RPL41", "CLEC10A", "JAML", "HLA-DRB5", "RPS17", "HLA-DQB1", "MMP9", "IL22RA2", "AGRP"), 

"DC: DC3 (model markers)"=c("a dendritic cell subtype found in bone marrow which shares features with monocytes and conventional type 2 dendritic cells", "FOXH1", "FCN1", "C11orf72", "BMP3", "GPIHBP1", "LINC01644", "FCER1A", "OGDHL", "LINC02206", "OR2A25"), 

"DC: Migratory DCs (model markers)"=c("migratory dendritic cells that transport antigens to the draining lymph nodes during both homeostatic conditions and infections", "BX255923.3", "LINC01644", "NCCRP1", "IGHV5-78", "ENTHD1", "SIGLEC17P", "WDR49", "OR2A25", "GPIHBP1", "CIB3"), 

"DC: Transitional DC (model markers)"=c("immature dendritic cells which follow dendritic cell precursors and are committed to the mature dendritic cells", "SCT", "FOXH1", "SOST", "SYCP1", "HLA-DPB2", "LINC01709", "TCL1B", "IGHE", "RP11-1275H24.3", "GPIHBP1"), 

"DC precursor: DC precursor (model markers)"=c("dendritic cell precursors that differentiate into dendritic cells", "PCLAF", "RNASE2", "GPIHBP1", "RPS3AP34", "C11orf72", "OGDHL", "IGHE", "SOCS2-AS1", "RP11-589C21.6", "LINC01644"), 

"Double-negative thymocytes: Double-negative thymocytes (model markers)"=c("the early immature thymocytes from the thymus in the double negative (co-receptors CD4- and CD8-) stage", "CD99", "IGLL1", "PTCRA", "FXYD2", "PDLIM1", "MFAP4", "MZB1", "JCHAIN", "CD7", "SELL"), 

"Double-positive thymocytes: Double-positive thymocytes (model markers)"=c("immature thymocytes capable of binding MHC class I or II in the double positive (co-receptors CD4+ and CD8+) stage", "SMPD3", "SH3TC1", "ELOVL4", "CD1B", "CD8B", "RP11-144L1.4", "ARPP21", "CD52", "RAG2", "CHRNA3"), 

"Early MK: Early MK (model markers)"=c("early megakaryocytes which are committed to megakaryocyte lineage cells", "FOXH1", "CCDC175", "OR2A25", "CIB3", "AC098614.1", "C11orf72", "IGHV5-78", "PCP2", "IGHE", "CD207"), 

"Endothelial cells: Endothelial cells (model markers)"=c("single-layered lining cells that constitute the interior vasculature", "form part of the microenvironment and fulfill different functions", "SELE", "RAMP3", "ACKR1", "SOX17", "ADGRL4", "VWF", "JAM2", "MMRN2", "SELP", "RP11-536O18.1"), 

"Epithelial cells: Epithelial cells (model markers)"=c("highly specialised cells which are located in the outer layer of stroma and form the cellular components of the epithelium", "ASCL1", "PSMB11", "FOXN1", "WFDC2", "PAX1", "FOXG1", "COL17A1", "DSP", "CDH3", "TBATA"), 

"Erythrocytes: Erythrocytes (model markers)"=c("biconcave enucleated red blood cells filled with hemoglobin to transport oxygen and carbon dioxide between the lungs and tissues", "RP11-797H7.5", "CTA-363E6.6", "CTA-392E5.1", "RP11-1275H24.3", "LCNL1", "LINC01644", "IGHV5-78", "ZNF812P", "IL1RAPL1", "GPIHBP1"), 

"Erythroid: Early erythroid (model markers)"=c("early erythroid cells which are committed to erythroid lineage cells", "FAM178B", "APOC1", "CD207", "KIAA0087", "PCLAF", "OR2A25", "LINC01644", "PRSS57", "IL17A", "PRG2"), 

"Erythroid: Late erythroid (model markers)"=c("late erythroid cells following early and middle erythroid cells", "IFIT1B", "PMCH", "OR2A25", "LINC01644", "SLC12A3", "DRAIC", "ABCC13", "CTAG2", "FOXH1", "CD207"), 

"Erythroid: Mid erythroid (model markers)"=c("middle erythroid cells in fetal liver", "blood and bone marrow following early erythroid cells and maturing into late erythroid cells", "IGHE", "AC098614.1", "ZNF812P", "LINC01644", "XCR1", "C11orf72", "RP11-589C21.6", "RPS3AP34", "HOXB-AS3", "IL17A"), 

"ETP: ETP (model markers)"=c("early thymic progenitors migrated to the thymus before turning into double-negative thymocytes or other cells", "ACY3", "LINC01644", "WDR49", "FOXH1", "GPIHBP1", "RP11-589C21.6", "ZFHX4-AS1", "KIAA0087", "LINC01709", "PENK"), 

"Fibroblasts: Fibroblasts (model markers)"=c("the most common cells of connective tissues which synthesize the extracellular matrix and collagen to maintain tissue homeostasis", "SMOC2", "SFRP1", "NTRK2", "PRRX1", "EBF2", "OLFML1", "MXRA5", "F10", "ANGPTL1", "PDGFRA"), 

"Granulocytes: Granulocytes (model markers)"=c("cells of innate immune system with specific granules in the cytoplasm including neutrophils", "eosinophils and basophils", "SIGLEC17P", "MS4A2", "IGHV5-78", "ZFHX4-AS1", "RP11-589C21.6", "OSTN-AS1", "MAK", "GPIHBP1", "SLC10A5", "KIAA0087"), 

"Granulocytes: Neutrophils (model markers)"=c("the most abundant type of granulocytes that contains distinctive cytoplasmic granules and forms an essential part of the innate immune system", "LINC01644", "GATA3-AS1", "OR2A25", "TCL6", "IL22RA2", "IGHV5-78", "SIGLEC17P", "IGHE", "LINC01709", "FCGR2C"), 

"HSC/MPP: CMP (model markers)"=c("highly proliferative common myeloid progenitors which later give rise to granulocyte-monocyte progenitors and megakaryocyte-erythroid progenitors", "MPO", "SIGLEC17P", "SORD2P", "AC098614.1", "CD207", "CPA3", "SYCP1", "ZFHX4-AS1", "IL22", "OR2A25"), 

"HSC/MPP: Early lymphoid/T lymphoid (model markers)"=c("early lymphoid/T lymphocytes with lymphocyte potential in the fetal liver before T cells emerged from the thymus", "HBA2", "HBG2", "RACK1", "ATP5F1E", "HBA1", "IGHV5-78", "SIGLEC17P", "SOST", "ST18", "RP11-505E24.2"), 

"HSC/MPP: ELP (model markers)"=c("early lymphoid progenitors which originate from hematopoietic stem cells of the bone marrow and migrate to the thymus", "IGHV5-78", "MIR1-1HG-AS1", "SIGLEC17P", "GPIHBP1", "LINC01644", "REG1A", "OSTN-AS1", "KIAA0087", "FOXH1", "IGHE"), 

"HSC/MPP: GMP (model markers)"=c("hematopoietic granulocyte-monocyte progenitors that are committed to the granulocyte and monocyte lineage cells", "FCGR2C", "MPO", "OR2A25", "LINC02202", "LINC02227", "FOXP1-IT1", "ZNF321P", "LINC01943", "OSTN-AS1", "ELANE"), 

"HSC/MPP: HSC/MPP (model markers)"=c("hematopoietic stem cells and multipotent progenitor cells with the potential of differentiating into different blood cells", "AVP", "SPINK2", "ZNF321P", "LINC01709", "GATA3-AS1", "DTX2P1-UPK3BP1-PMS2P11", "OR2A25", "CLEC9A", "GRIN1", "ERICH1-AS1"), 

"HSC/MPP: Megakaryocyte-erythroid-mast cell progenitor (model markers)"=c("shared progenitors in the fetal liver which originate from common myeloid progenitors and differentiate into megakaryocytes", "erythroid cells or mast cells", "FCER1A", "CKB", "HBD", "IGHE", "DRAIC", "HLA-DPB2", "OR2A25", "KIAA0087", "C11orf72", "ZFHX4-AS1"), 

"HSC/MPP: MEMP (model markers)"=c("shared progenitors which are derived from common myeloid progenitors and will differentiate into megakaryocytes", "erythroid cells or mast cells", "CTAG2", "LINC01644", "FOXH1", "GIHCG", "SIGLEC17P", "RP11-701P16.2", "LINC02206", "IGHE", "OR2A25", "CALB2"), 

"HSC/MPP: Neutrophil-myeloid progenitor (model markers)"=c("progenitors of neutrophils and myeloid cells which are transitioned from hematopoietic stem cells and multipotent progenitors", "MPO", "MS4A3", "LINC01709", "OSTN-AS1", "KIAA0087", "KCNE5", "CUX2", "RP11-589C21.6", "C11orf72", "GPIHBP1"), 

"ILC: ILC (model markers)"=c("specialised innate immune cells from the lymphoid lineage but without antigen-specific T-cell receptors on the surface", "AREG", "IL7R", "IGHV5-78", "XCR1", "KLRB1", "SYCP1", "RP11-589C21.6", "SIGLEC17P", "LINC01644", "GRIN1"), 

"ILC: ILC1 (model markers)"=c("innate lymphoid cell subpopulation I that is non-cytotoxic and has overlapping phenotypes and functions with natural killer cells", "SIGLEC17P", "LINC01013", "LINC01644", "LINC02202", "KIAA0087", "SYCP1", "OSTN-AS1", "REG1A", "IGHE", "IL22"), 

"ILC: ILC2 (model markers)"=c("innate lymphoid cell subpopulation II that promotes type 2 inflammation and is involved in immune response against large extracellular pathogens", "IGHE", "IL17RB", "LINC01709", "MIR144", "FCGR2C", "TRDJ1", "OR2A25", "LINC02227", "SNORA67", "Metazoa_SRP"), 

"ILC: ILC3 (model markers)"=c("innate lymphoid cell subpopulation III that is required for host defense against specific extracellular bacteria and fungi", "GATA3-AS1", "FOXH1", "EIF2S2P4", "CD207", "RP11-677O4.2", "LINC02227", "LINC01644", "RP11-1275H24.3", "SIGLEC17P", "OSTN-AS1"), 

"ILC: CD16+ NK cells (model markers)"=c("CD16+ granular lymphocytes that play protective roles against both infectious pathogens and cancer using antibody-dependent cell-mediated cytotoxicity", "FCER1G", "GNLY", "IGFBP7", "TYROBP", "GZMB", "MYOM2", "PTGDS", "FCGR3A", "PRSS57", "ADAMTS1"), 

"ILC: CD16- NK cells (model markers)"=c("CD16- granular lymphocytes that play protective roles against both infectious pathogens and cancer using antibody-dependent cell-mediated cytotoxicity", "TYROBP", "LDB2", "CCL3", "FCER1G", "NKG7", "GSTP1", "CLIC3", "IRF8", "CXCL3", "KLRB1"), 

"ILC: NK cells (model markers)"=c("granular lymphocytes that play protective roles against both infectious pathogens and cancer using antibody-dependent cell-mediated cytotoxicity", "NKG7", "CST3", "GNLY", "SPINK2", "HMOX1", "MALAT1", "IGFBP4", "HBG2", "STMN1", "KLHL23"), 

"ILC: Transitional NK (model markers)"=c("immature natural killer cells which originate from natural killer cell precursors and are committed to mature natural killer cells", "GNLY", "TYROBP", "KLRB1", "PRDM16", "FOXH1", "AGRP", "SYCP1", "LINC02206", "RP11-423H2.3", "GPIHBP1"), 

"ILC precursor: ILC precursor (model markers)"=c("innate lymphoid cell precursors which give rise to innate lymphoid cells and lack the characteristics of their mature progenies", "RBM24", "GPIHBP1", "IGHE", "KIAA0087", "HPN", "PCDH9-AS1", "HUNK", "OR2A25", "PDZK1", "RP11-423H2.3"), 

"Macrophages: Hofbauer cells (model markers)"=c("primitive placental resident macrophages with granules and vacuoles found in placenta particularly during early pregnancy", "CGA", "PAGE4", "XAGE3", "GAPLINC", "PANX2", "GATA3-AS1", "LINC01644", "LYVE1", "ZFHX4-AS1", "KIAA0087"), 

"Macrophages: Kidney-resident macrophages (model markers)"=c("long-lived macrophages resident in the kidney under non-inflammatory conditions and maintaining homeostasis and resolving inflammation ", "IGHV5-78", "LINC02206", "KIAA0087", "C11orf72", "TCL1B", "LINC01644", "IGHE", "FCGR2C", "KLRF2", "RP11-1275H24.3"), 

"Macrophages: Kupffer cells (model markers)"=c("resident macrophages in liver under non-inflammatory conditions which line the hepatic sinusoids and are involved in erythrocyte clearance", "TIMD4", "CETP", "SNORD3D", "EGR2", "PENK", "LINC01644", "ECSCR", "ARRDC3", "TCHH", "CD5L"), 

"Macrophages: Macrophages (model markers)"=c("specialised mononuclear phagocytic cells which recognise", "engulf and degrade cellular debris and pathogens in various tissues", "GNLY", "APOE", "C1QA", "CD5L", "FTL", "FTH1", "C1QB", "APOC1", "SELENOP", "RNASE1"), 

"Mast cells: Mast cells (model markers)"=c("long-lived innate migrant cells found in most tissues with many large basophilic granules rich in histamine and heparin", "TPSAB1", "CPA3", "RP11-354E11.2", "CMA1", "IGHE", "IL1RL1", "SOST", "CCNA1", "KRT1", "OR2A25"), 

"Megakaryocyte precursor: Megakaryocyte precursor (model markers)"=c("megakaryocyte precursors in the bone marrow that are committed to megakaryocytes", "ADRA2A", "PCDH9-AS1", "FCER1A", "OR2A25", "SLC10A5", "IL22", "IGHV5-78", "C11orf72", "KIAA0087", "RPL7P23"), 

"Megakaryocytes/platelets: Megakaryocytes/platelets (model markers)"=c("large multinucleated cells with numerous azurophilic granules and with mature blood platelets released from the cytoplasm", "RP11-423H2.3", "LY6G6F-LY6G6D", "PF4V1", "RP11-1275H24.3", "SOST", "LINC02227", "OR2A25", "GLYATL2", "IGHV5-78", "SDAD1P1"), 

"MNP: MNP (model markers)"=c("mononuclear phagocytes including dendritic cells", "monocytes and macrophages that have in common the property of phagocytosis", "MTRNR2L10", "IGHV5-78", "LINC01644", "TCL6", "SOCS2-AS1", "GATA3-AS1", "LINC01483", "IL17A", "LINC01943", "RPS3AP34"), 

"Mono-mac: Mono-mac (model markers)"=c("mononuclear phagocytes including monocytes and macrophages that are involved in the regulation of innate and adaptive immunity", "NKG7", "KLRB1", "HLA-DRA", "KRT13", "FTH1", "SPRR3", "C1QC", "SAT1", "LYZ", "S100A2"), 

"Monocyte precursor: Monocyte precursor (model markers)"=c("monocyte precursors that are committed to the monocytes", "RETN", "RPS3AP34", "LYZ", "FCN1", "LINC01644", "IGHE", "ST18", "LCNL1", "IGHV5-78", "CUX2"), 

"Monocytes: Classical monocytes (model markers)"=c("CD14+ myeloid mononuclear recirculating leukocytes that are capable of differentiating into macrophages and myeloid lineage dendritic cells", "S100A9", "LINC02206", "EEF1A1", "CD207", "TCL1B", "CST7", "OR2A25", "SOST", "IL22", "MIR1-1HG-AS1"), 

"Monocytes: Non-classical monocytes (model markers)"=c("CD16+ myeloid mononuclear recirculating leukocytes that are capable of differentiating into macrophages and myeloid lineage dendritic cells", "RHOC", "HES4", "MS4A7", "DNAJA4", "ABI3", "BAG3", "CUX2", "IGHE", "LINC02206", "SIGLEC17P"), 

"Monocytes: Monocytes (model markers)"=c("myeloid mononuclear recirculating leukocytes that are capable of differentiating into macrophages and myeloid lineage dendritic cells", "TYROBP", "FTL", "IGKC", "TRAC", "RPS17", "SAT1", "NEAT1", "LYZ", "HLA-DPB1", "IGLC2"), 

"Myelocytes: Myelocytes (model markers)"=c("early granulocyte precursors that are derive from promyelocytes and later mature into metamyelocytes", "LYZ", "FTH1", "MT-RNR2", "S100A6", "SRGN", "FTL", "TYROBP", "ACTB", "GPIHBP1", "MT-ND3"), 

"pDC: pDC (model markers)"=c("rare plasmacytoid dendritic cell subpopulation which serves as the major source of type I interferons when the body is infected by a virus", "SCT", "FOXH1", "LCNL1", "SYCP1", "OR2A25", "LINC01709", "DRAIC", "TRIM71", "OGDHL", "LINC02227"), 

"pDC precursor: pDC precursor (model markers)"=c("precursors of plasmacytoid dendritic cells which have intermixed lymphoid and myeloid origins and give rise to plasmacytoid dendritic cells", "MYBPH", "DRAIC", "LINC01709", "SIGLEC17P", "GPR15", "REG1A", "AGRP", "IL22", "C11orf72", "SOST"), 

"Plasma cells: Plasma cells (model markers)"=c("B-lymphocyte white blood cells capable of secreting large quantities of immunoglobulins or antibodies", "JCHAIN", "IGHA2", "MZB1", "XBP1", "IGKC", "RP11-354E11.2", "RPL23AP18", "TGFBR3L", "C11orf72", "ST18"), 

"Promyelocytes: Promyelocytes (model markers)"=c("early granulocyte progenitors in the bone marrow which are derived from the myeloblasts and develop into the myelocytes", "PRTN3", "ELANE", "FCGR2C", "AZU1", "XCR1", "LINC01644", "LINC01709", "FOXH1", "CD207", "GPIHBP1"), 

"T cells: CD8a/a (model markers)"=c("unconventional T lymphocytes with CD8 alpha/alpha homodimers binding MHC molecules", "NUCB2", "CD27", "LEF1", "ZNF683", "MALAT1", "PRKCH", "CD8A", "GPIHBP1", "TRGC2", "CTSW"), 

"T cells: CD8a/b(entry) (model markers)"=c("T lymphocytes with CD8 alpha/beta heterodimers in the late double-positive (entry) stage during T cell development", "CD1E", "ITM2A", "MALAT1", "SATB1", "RPL10", "CD44", "FTL", "TMSB10", "RPS27", "MT-ND2"), 

"T cells: Follicular helper T cells (model markers)"=c("specialised CD4+ T cell subpopulation found in the periphery within B cell follicles and stimulating the cognate follicular B cells for class-switching", "KLRB1", "SOCS3", "CH25H", "IGHV5-78", "HOXB5", "LINC01644", "SIGLEC17P", "GPIHBP1", "ZFHX4-AS1", "FOXH1"), 

"T cells: gamma-delta T cells (model markers)"=c("unconventional T lymphocyte subpopulation expressing a gamma-delta T cell receptor complex on the surface to recognise antigens", "KIR2DL4", "KLRC2", "KCNK10", "LINC02227", "SYCP1", "KIAA0087", "GRIN1", "ZFHX4-AS1", "SIGLEC17P", "LINC01644"), 

"T cells: MAIT cells (model markers)"=c("mucosal-associated invariant T cells which have semi-invariant T-cell receptors and are restricted by the MHC I-like molecule MR1", "KLRB1", "IL7R", "GZMK", "LINC01709", "CEBPD", "LINC01871", "KLRG1", "IGHV5-78", "BMP3", "NCR3"), 

"T cells: Memory CD4+ cytotoxic T cells (model markers)"=c("CD4+ memory T cells which have cytotoxic activities by secreting granzymes and perforin and by killing target cells on the basis of MHC class II", "GZMK", "GZMA", "SOST", "IL22RA2", "SIGLEC17P", "IL26", "RP11-84C10.2", "PENK", "PPP3CB-AS1", "MIR144"), 

"T cells: NKT cells (model markers)"=c("distinct T lymphocyte subpopulation that expresses both alpha/beta T-cell receptors and surface receptors of natural killer cells", "NKG7", "GNLY", "IL32", "PRDM16", "ITM2C", "LINC02227", "CD52", "RPS3AP34", "TCL1B", "FOXH1"), 

"T cells: Regulatory T cells (model markers)"=c("T cell subpopulation which modulates immune responses and regulates other cells through direct cell-cell contact and cytokine release", "MS4A6A", "IL32", "TIMD4", "PTPRCAP", "RTKN2", "TIGIT", "MIR4435-2HG", "RGS1", "SRGN", "FOXH1"), 

"T cells: T(agonist) (model markers)"=c("unconventional T cell subpopulation in the thymus which expresses MIR155HG and shares some signatures with differentiating regulatory T cells", "HLA-B", "CD74", "SMS", "COTL1", "PTPRC", "TUBB", "BIRC3", "B2M", "ITM2A", "HLA-A"), 

"T cells: Tcm/Naive cytotoxic T cells (model markers)"=c("CD8+ cytotoxic T lymphocytes mainly localized in secondary lymphoid tissues and sustaining the responses by proliferating and producing new effectors", "CD8B", "CD8A", "LINC02446", "B2M", "FTL", "HLA-C", "RPS2", "GAS5", "HLA-A", "HLA-B"), 

"T cells: Tcm/Naive helper T cells (model markers)"=c("CD4+ helper T lymphocytes mainly localized in secondary lymphoid tissues and sustaining the responses by proliferating and producing new effectors", "IGFBP1", "MT-RNR2", "FCGRT", "TMSB10", "SNHG3", "HBG1", "RPS18", "LTB", "LINC01709", "CTSL"), 

"T cells: Tem/Trm cytotoxic T cells (model markers)"=c("CD8+ cytotoxic T lymphocytes mainly localized in lymphoid and peripheral tissues and presenting an immediate", "but not sustained", "immune defense", "CCL5", "CD8B", "CD8A", "GZMK", "RPS29", "RPS2", "ATP5F1E", "CST7", "IGHM", "NKG7"), 

"T cells: Tem/Temra cytotoxic T cells (model markers)"=c("terminally differentiated CD8+ cytotoxic T lymphocytes with effector memory phenotypes that re-express CD45RA ", "GZMH", "RP11-466H18.1", "FCGR2C", "SIGLEC17P", "PZP", "IL7R", "OR2A25", "MT-ND3", "RP11-1275H24.3", "SNX29P1"), 

"T cells: Trm cytotoxic T cells (model markers)"=c("tissue resident CD8+ cytotoxic T lymphocytes mainly localized in epithelial tissues such as gut", "CCL5", "CD8A", "IGHA1", "CCL4", "RPLP2", "IFIT1B", "CD8B", "IGHV5-78", "IL7R", "LINC01644"), 

"T cells: Tem/Effector helper T cells (model markers)"=c("CD4+ helper T lymphocytes mainly localized in lymphoid and peripheral tissues and presenting an immediate", "but not sustained", "immune defense", "ANXA1", "DNTT", "MIR1244-2", "KLRB1", "S100A4", "DONSON", "AQP3", "LTB", "RP11-589C21.6", "SFTPA1"), 

"T cells: Tem/Effector helper T cells PD1+ (model markers)"=c("CD4+ helper T lymphocyte subpopulation in the thymus which features the expression of programmed cell death protein 1 (PD-1)", "FOXH1", "RP11-589C21.6", "SOST", "CTAG2", "PANX2", "IGHV5-78", "ZFHX4-AS1", "SOCS2-AS1", "KIAA0087", "INPP5J"), 

"T cells: Treg(diff) (model markers)"=c("unconventional T lymphocyte subpopulation in the thymus which connects αβ T cells and canonical regulatory T cells and is still differentiating", "S100A4", "HLA-B", "RPS26", "RGS1", "HLA-A", "CD1E", "TSC22D3", "PTPRC", "FCGR2C", "SIGLEC17P"), 

"T cells: Type 1 helper T cells (model markers)"=c("CD4+ helper T lymphocyte subpopulation which is capable of producing interferon-gamma and modulating cell-mediated immune responses", "CCL5", "KLRB1", "SIGLEC17P", "FCGR2C", "RPL17", "CEBPB", "PRDM16", "ANXA1", "FOXH1", "TCL1B"), 

"T cells: Type 17 helper T cells (model markers)"=c("CD4+ helper T lymphocyte subpopulation which is capable of producing interleukin 17 (IL-17) and mediating protective immunity and autoimmunity", "KLRB1", "ANXA1", "S100A4", "CEBPD", "FCGR2C", "MIR144", "IGHE", "LINC02245", "LINC02206", "LINC02202")

)

