# This file contains code that is useful for some combination of dataset level, community level and GEP level summary reports.


programNames <- c("1. Naive B Cell", 
"2. Inflammatory:cytokine signalling I", 
"3. Normal Adrenal Fibroblasts", 
"4. T Cell Activation", 
"5. Neuroblastoma: Adrenergic I (pre-chromaffin)", 
"6. Hepatocyte-like Cell", 
"7. Endothelial: Vascular", 
"8. Monocyte:M-MDSC", 
"9. Adrenocortical Cells", 
"10. Erythroblasts", 
"11. Plasma Cells", 
"12. GMP", 
"13. Cell Cycle (G2-M)", 
"14. Mature NKT Cells", 
"15. Megakaryocyte", 
"16. Mononuclear Phagocyte", 
"17. Undifferentiated lymphocyte",
"18. Antigen Presentation",
"19. Cancer Associated Fibroblast: Myofibroblastic (IGF1+)", 
"20. Keratin Expression", 
"21. Mast Cells", 
"22. Neuroblastoma: MYCN", 
"23. Megakaryocyte-erythroid precursors", 
"24. M2 Macrophage", 
"25. Cancer Associated Fibroblast: Myofibroblast (POSTN+)", 
"26. Pro-B Cells", 
"27. Endothelial: Lymphatic", 
"28. TNFa signalling via NFkB", 
"29. Endothelial", 
"30. Plasmacytoid Dendritic Cell", 
"31. Myelocyte", 
"32. Regulatory T Cell", 
"33. Monocyte: Inflammatory", 
"34. Natural killer cell (CD16+)", 
"35. Cancer Associated Fibroblast: Myofibroblast (Contractile)", 
"36. Endothelial: Angiogenesis",
"37. Naive T Cells", 
"38. CD8+ T Cells", 
"39. Neuroblastoma: Adrenergic II (pre-neuronal like)",
"40. Fibroblast: Inflammatory", 
"41. Inflammatory:cytokine signalling II")

programDescriptions <- c("1.", 
"2.", 
"3.", 
"4.", 
"5.", 
"6.", 
"7.", 
"8.", 
"9.", 
"10", 
"11.", 
"12.", 
"13.", 
"14.", 
"15.", 
"16.", 
"17.",
"18.",
"19.", 
"20.", 
"21.", 
"22.", 
"23.", 
"24.", 
"25.", 
"26.", 
"27.", 
"28.", 
"29.", 
"30.", 
"31.", 
"32.", 
"33.", 
"34.", 
"35.", 
"36.",
"37.", 
"38.", 
"39.",
"40.", 
"41.")


#+ loadLibraries, echo=FALSE, message = FALSE
library(ggplot2)
library(plotly)
library(gapminder)
library(knitr)
library(htmlwidgets)
library(org.Hs.eg.db)
library(M3C)
library(GSA)
library(RColorBrewer)
library(jaccard)
library(org.Hs.eg.db) ## /addedCode remember to install it if you don't have it already
library(vroom)
library(SeuratDisk)
library(BioCircos)
library(webshot)
library(infercnv)
library(viridisLite)
library(tidyr)


## Load the required data

data.dir.meta <- "/Volumes/clusterhome/rchapple/scRNAseq/cNMF_meta_analysis/jaccard/all_pairwise/unfiltered_data/redo/"
suppdata.dir <- "/Users/rchapple/Documents/scRNASeq/cNMF_meta_analysis/R/annotation/V3.0/human/GOSH/GOSH_data/"
indata.dir <- "/Volumes/clusterhome/rchapple/scRNAseq/human/analyses/cNMF/input/benchmark/Kildisiute_2021/GOSH/h5ad/"
outdata.dir <- "/Volumes/clusterhome/rchapple/scRNAseq/human/analyses/cNMF/output/benchmark/unfiltered_data/Kildisiute_GOSH/h5ad/"

params <- readRDS(paste0(data.dir.meta, "optimized_community_paramters.RDS"))

#comms <- readRDS(paste0(data.dir.meta, "community_geps_filtered.RDS"))
comms <- readRDS(paste0(data.dir.meta, "community_geps.RDS"))
comms <- comms[["Kildisiute_GOSH"]]


commsSorted <- sort(comms[, thisProgNumber])

# Load the raw count data.

split0 <- LoadH5Seurat(paste0(indata.dir, "GOSH.split0.h5Seurat"), assays = "counts")
split0 <- split0@assays$RNA@counts
split1 <- LoadH5Seurat(paste0(indata.dir, "GOSH.split1.h5Seurat"), assays = "counts")
split1 <- split1@assays$RNA@counts

allData <- rbind(t(split0), t(split1))


noDetectableExpression <- colnames(allData)[apply(allData, 2, function(thecol)return(sum(thecol) == 0))] ## /addedCode, want to list undetectable genes....


## /addedCode, TPM block below.
# Create TPM data.
countsMatNb <- t(allData)
geneLensFromBiomart <- read.delim(paste0(suppdata.dir,"mart_export.txt")) # I got this from BioMart.
geneLists <- split(geneLensFromBiomart, geneLensFromBiomart[,"Gene.name"]) # for each gene pull out the maximum CDS length.
theMaxes <- sapply(geneLists, function(a)return(max(a[,"Transcript.length..including.UTRs.and.CDS."])))
hasLengthAndExpression <- names(theMaxes)[names(theMaxes) %in% rownames(countsMatNb)]
countsMatNb_subset <- countsMatNb[hasLengthAndExpression, ]
txLengths <- theMaxes[hasLengthAndExpression]

# function for calculating TPM from count data and lengths. (Credit: Mike Love on Biostars)
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
allDataTpm <- t(tpm3(countsMatNb_subset, txLengths))



# Load the H matrices

split0_H <- vroom(file = paste0(outdata.dir, "cNMF_split0/cNMF_split0.usages.k_", params$Kildisiute_GOSH$Rank, ".dt_0_10.consensus.txt"), show_col_types = FALSE)
split0_H <- tibble::column_to_rownames(split0_H, 1)
split1_H <- vroom(file = paste0(outdata.dir, "cNMF_split1/cNMF_split1.usages.k_", params$Kildisiute_GOSH$Rank, ".dt_0_10.consensus.txt"), show_col_types = FALSE)
split1_H <- tibble::column_to_rownames(split1_H, 1)

# Load the list mapping community number to column of H (i.e. H includes all programs for this split, not just the ones reproduced in split0 & split1).

communitiesOther <- readRDS(paste0(data.dir.meta, "Kildisiute_GOSH_gepindex.RDS"))

# Load the Singler output for cell type annotations and delta values (which are a type of confidence score assocaited with the cell type annotations)
singlerCellLabels <- readRDS(paste0(suppdata.dir, "GOSH.sc.label.scores.RDS"))
singlerDeltaValues <- readRDS(paste0(suppdata.dir, "delta_values.RDS"))


# Load the original cluster annotations from this study.
theMetaData <- read.delim(paste0(suppdata.dir, "GOSH_metadata.csv"), as.is=T, sep = ",")
boxplotCategories <- theMetaData[, c("SampleName", "Annotation")] # Create the catagories for the boxplot (this will vary based on the names provided in theMetaData)

# Create the text string for hovering in the UMAP (this contains the generic text strings, the activity score will be added when the UMAP is created). Having a generic string will mean this can be created for each dataset based on whatever additional information is available (the Activity Score will always be available and added later).
genericCellHoverText <- paste("CellID: ", rownames(theMetaData), 
                              "\nPatient ID: ", theMetaData[, 1], 
                              "\nOriginal Annotation: ", theMetaData[, "Annotation"], 
                              "\nTotal Read Counts: ", theMetaData[, "nCount_RNA"], 
                              "\nGenes Detected: ", theMetaData[, "nFeature_RNA"], 
                              "\nmtGenes: ", round(theMetaData[, "mtGenes"], digits=2), 
                              "\nHeat Shock Genes: ", round(theMetaData[, "hspGenes"], digits=2), 
                              "\nRibosomal Genes: ", round(theMetaData[, "riboGenes"], digits=2), sep="")
names(genericCellHoverText) <- rownames(theMetaData)


GOSH <- LoadH5Seurat("/Volumes/clusterhome/rchapple/scRNAseq/human/analyses/cNMF/input/Kildisiute_2021/unfiltered_data/umap/GOSH.h5Seurat", assays = "counts")
umapOut <- GOSH@reductions$umap@cell.embeddings
