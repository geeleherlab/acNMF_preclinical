library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

mouse <- readRDS(file = 'NB_DBHiCre.RDS')
Idents(mouse) <- mouse@meta.data$seurat_clusters
mouse <- RenameIdents(object = mouse,
                      `0` = "Neuroendocrine",
                      `1` = "Neuroendocrine",
                      `2` = "Neuroendocrine",
                      `3` = "Neuroendocrine",
                      `4` = "Neuroendocrine",
                      `5` = "Neuroendocrine",
                      `6` = "Neuroendocrine",
                      `7` = "Neuroendocrine",
                      `8` = "Neuroendocrine",
                      `9` = "Myeloid",
                      `10` = "Neuroendocrine",
                      `11` = "Neuroendocrine",
                      `12` = "Neuroendocrine",
                      `13` = "Neuroendocrine",
                      `14` = "Neuroendocrine",
                      `15` = "B Cell",
                      `16` = "Neuroendocrine",
                      `17` = "Neuroendocrine",
                      `18` = "Neuroendocrine",
                      `19` = "Mesenchymal",
                      `20` = "T Cell",
                      `21` = "Schwannian Stromal Cell",
                      `22` = "Neuroendocrine",
                      `23` = "Neuroendocrine",
                      `24` = "Neuroendocrine",
                      `25` = "Neuroendocrine"
)

treatment <- factor(mouse@meta.data$drug)
mouse_id <- factor(mouse@meta.data$orig.ident)
cell_type <- factor(as.vector(Idents(mouse)))
metadata <- NULL
metadata$Treatment <- treatment
metadata$Mouse_ID <- mouse_id
metadata$Cell_Type <- cell_type

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = mouse@assays$RNA@counts), colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "mouse_id", "annotation", "treatment")]

# Named vector of celltypes
celltypes <- purrr::set_names(levels(sce$Cell_Type))
ncelltypes <- length(celltypes)
# Named vector of sample names
ids <- purrr::set_names(levels(sce$Mouse_ID))
nids <- length(ids)
#Named vector of annotation
trt <- purrr::set_names(levels(sce$Treatment))
ntrt <- length(trt)

#Sample level metadata
## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$Mouse_ID))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in ids vector
m <- match(ids, sce$Mouse_ID)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>% select(-"Cell_Type")

########################################
# Calculate quality control (QC) metrics
dim(sce)
[1] 20982 68253
sce.qc <- perCellQCMetrics(sce)
sce.qc$is_outlier <- isOutlier(metric = sce.qc$total, nmads = 2, type = "both", log = TRUE)
# Remove outlier cells
sce.filtered <- sce[, sce.qc$is_outlier == FALSE]
dim(sce.filtered)
[1] 20982 65111
## Remove lowly expressed genes which have less than 10 cells with any counts
sce.filtered <- sce.filtered[rowSums(counts(sce.filtered) > 1) >= 10, ]
dim(sce.filtered)
[1] 15776 65111

####################################################
# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce.filtered)[, c("Cell_Type", "Mouse_ID")]


# Aggregate across cluster-sample groups
sce.agg <- aggregate.Matrix(t(counts(sce.filtered)), groupings = groups$Mouse_ID, fun = "sum") 

########################################
#DESeq2
#######################################

# Get sample names for each of the cell type clusters
# prep. data.frame for plotting
get_sample_ids <- function(x){sce.agg.split[[x]] %>% colnames()}
de_samples <- map(1:length(celltypes), get_sample_ids) %>% unlist()

# Get cluster IDs for each of the samples
samples_list <- map(1:length(celltypes), get_sample_ids)
get_celltype_ids <- function(x){rep(names(sce.agg.split)[x], each = length(samples_list[[x]]))}
de_celltype_ids <- map(1:length(celltypes), get_celltype_ids) %>% unlist()

# Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(Cell_Type = de_celltype_ids, Mouse_ID = de_samples)
gg_df <- left_join(gg_df, ei[, c("Mouse_ID", "Treatment")]) 

metadata.new <- gg_df %>% dplyr::select(Cell_Type, Mouse_ID, Treatment) 
metadata.new   

# # Subset the metadata to only the Neuroendocrine cells
celltype_metadata <- metadata.new[which(metadata.new$Cell_Type == "Neuroendocrine"), ]
head(celltype_metadata)


# Assign the rownames of the metadata to be the sample IDs
rownames(celltype_metadata) <- celltype_metadata$Mouse_ID
head(celltype_metadata)

#RHC 230302
celltype_metadata <- celltype_metadata[, -1]

# Subset the counts to only the Neuroendocrine cells
counts <- t(sce.agg)

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(celltype_metadata) == colnames(counts))   


dds <- DESeqDataSetFromMatrix(counts, 
                              colData = celltype_metadata, 
                              design = ~ Treatment)

#################################3
#Visualization
#################################3

#PCA

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, intgroup = "Treatment")

#Heirarchical clustering

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = celltype_metadata[, c("Treatment"), drop=F])

#############################
#Run DESeq2

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Output results of Wald test for contrast for stim vs ctrl
levels(celltype_metadata$Treatment)[2]
[1] "vehicle"
levels(celltype_metadata$Treatment)[1]
[1] "cisplatin"

contrast <- c("Treatment", levels(celltype_metadata$Treatment)[2], levels(celltype_metadata$Treatment)[1])

res <- results(dds, contrast = contrast, alpha = 0.05)
resultsNames(dds)
[1] "Intercept"      "Treatment_vehicle_vs_cisplatin"

res2 <- lfcShrink(dds, contrast =  contrast, coef = 2, type = "ashr")
# > using 'ashr' for LFC shrinkage. If used in published research, please cite:
# >   Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
# >   https://doi.org/10.1093/biostatistics/kxw041

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res2 %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

# Check results output
res_tbl

# Write all results to file
write.csv(res_tbl, paste0("results/alldata/Neuroendocrine_", levels(celltype_metadata$Treatment)[2], "_vs_", levels(celltype_metadata$Treatment)[1], "_all_genes.csv"), quote = FALSE, row.names = FALSE)

# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)

# Check significant genes output
sig_res

# Write significant results to file
write.csv(sig_res, paste0("results/alldata/Neuroedocrine_", levels(celltype_metadata$Treatment)[2], "_vs_", levels(celltype_metadata$Treatment)[1], "_sig_genes.csv"), quote = FALSE, row.names = FALSE)

###################
#Visualization

## ggplot of top genes
normalized_counts <- counts(dds, normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n=20)
top20_sig_norm <- data.frame(normalized_counts) %>% rownames_to_column(var = "gene") %>% dplyr::filter(gene %in% top20_sig_genes)
gathered_top20_sig <- top20_sig_norm %>% gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
gathered_top20_sig <- inner_join(ei[, c("Mouse_ID", "Treatment" )], gathered_top20_sig, by = c("Mouse_ID" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


#Heatmap
# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Set a color palette
#heat_colors <- brewer.pal(6, "YlOrRd")
library(viridis)
hc <- viridis(12)
# Run pheatmap using the metadata data frame for the annotation

ancol <- list(Mouse_ID = c(NB831 = "#2D1160FF", NB837 = "#51127CFF", NB839 = "#721F81FF", NB847 = "#932B80FF", NB849 = "#B63679FF", NB853 = "#D8456CFF", NB855 = "#F1605DFF", NB856 = "#FB8861FF", NB864 ="#FEAF77FF", NB883 = "#FED799FF", NB887 = "#FCFDBFFF"),
              Treatment = c(vehicle = "#EE302E", cisplatin = "#76C376"))
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = hc, 
         annotation_colors = ancol,
         cluster_rows = T, 
         show_rownames = F,
         annotation = celltype_metadata[, c("Mouse_ID", "Treatment")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20, 
         cutree_rows = 2, cutree_cols = 2)  

#Volcano plot
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_tbl %>% mutate(threshold = padj < 0.01 & abs(log2FoldChange) >= 1.5)

## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of treated relative to control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,30)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  


#Enrichr

library(enrichR)

dbs.full <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", 
         "GO_Cellular_Component_2021", 
         "GO_Biological_Process_2021", 
         "MSigDB_Oncogenic_Signatures", 
         "MSigDB_Hallmark_2020", 
         "KEGG_2021_Human", 
         "Descartes_Cell_Types_and_Tissue_2021", 
         "Drug_Perturbations_from_GEO_up",
         "Drug_Perturbations_from_GEO_down",
         "Disease_Perturbations_from_GEO_up",
         "Disease_Perturbations_from_GEO_down",
         "LINCS_L1000_Chem_Pert_Consensus_Sigs", 
         "LINCS_L1000_CRISPR_KO_Consensus_Sigs")


#Repressed
rep.index <- which(sig_res$log2FoldChange[1:500] > 0)
enrichr.drug.repressed <- enrichr(sig_res$gene[rep.index], dbs)
#Induced
ind.index <- which(sig_res$log2FoldChange[1:500] < 0)
enrichr.drug.induced <- enrichr(sig_res$gene[ind.index], dbs)

saveRDS(enrichr.drug.induced, file = 'results/alldata/enrichr.drug.induced.RDS')
saveRDS(enrichr.drug.repressed, file = 'results/alldata/enrichr.drug.repressed.RDS')


##############
#GSEA
#############
library(fgsea)

GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=500) %>%  ## maximum gene set size
                        #nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) 
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}




plot_geneset_clusters = function( gs_results, GO_file, min.sz = 4, main="GSEA clusters"){
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  
  myGO = fgsea::gmtPathways(GO_file)
  df = matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
  rownames(df) = colnames(df) = gs_results$pathway
  
  for ( i in 1:nrow(gs_results)) {
    genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
    for (j in 1:nrow(gs_results)) {
      genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
      ## Jaccards distance  1 - (intersection / union )
      overlap = sum(!is.na(match(genesI, genesJ )))
      jaccards = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = 1-jaccards
    }
  }
  
  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(df)
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
  ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it
  
  gs_results$Cluster = clust
  gs_results = gs_results[gs_results$Cluster != 0, ]
  
  ## select gene sets to label for each clusters
  bests = gs_results %>%  
    group_by( Cluster ) %>% 
    top_n(wt = abs(size), n = 1) %>% 
    .$pathway
  ## determine cluster order for plotting
  clust_ords = gs_results %>% 
    group_by( Cluster ) %>% 
    summarise("Average" = NES ) %>% 
    arrange(desc(Average)) %>% 
    .$Cluster %>% 
    unique
  
  gs_results$Cluster = factor(gs_results$Cluster, levels = clust_ords)
  
  gs_results$Label = ""
  gs_results$Label[gs_results$pathway %in% bests ] = gs_results$pathway[gs_results$pathway %in% bests ]
  gs_results$Label = str_remove(gs_results$Label, "GO_")
  gs_results$Label = tolower(gs_results$Label)
  
  g1 = ggplot(gs_results, aes(x = Cluster, y = NES, label = Label )) +
    geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
    scale_size_continuous(range = c(0.5,5)) +
    geom_text_repel( force = 2, max.overlaps = Inf) +
    ggtitle(main)
  
  return(g1)
}

#Gene list
gene_list <- res_tbl$log2FoldChange
names(gene_list) <- str_to_title(tolower(res_tbl$gene))
gene_list <- sort(gene_list, decreasing = T)

gene_list <- sig_res$log2FoldChange
names(gene_list) <- str_to_title(tolower(sig_res$gene))
gene_list <- sort(gene_list, decreasing = T)

#GSEA MSigDB
#Hallmark (MH)

GO_file <- "GSEA/Neuroendocrine/preranked_list/MSigDB/alldata/MH/mh.all.v2022.1.Mm.symbols.hallmark.gmt"
gsea_res = GSEA(gene_list, GO_file, pval = 0.05)
dim(gsea_res$Results)

plot_geneset_clusters( gs_results = gsea_res$Results[gsea_res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )

#Curated gene sets (M2)
GO_file <- "GSEA/preranked_list/MSigDB/alldata/M2/m2.all.v2022.1.Mm.symbols.gmt"
gsea_res = GSEA(gene_list, GO_file, pval = 0.05)
dim(gsea_res$Results)

plot_geneset_clusters( gs_results = gsea_res$Results[gsea_res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )


#Curate gene sets (M2) plus ADRN MES
admes <- read.delim(file = "GSEA/ADRN_MES.csv", sep = ",", header = F)
colnames(admes) <- c("GENE", "ADMES")
adrn <- admes %>% filter(ADMES == "ADRN")
mes <- admes %>% filter(ADMES == "MES")
adrn <- str_to_title(tolower(adrn$GENE))
mes <- str_to_title(tolower(mes$GENE))
adrn <- c("ADRN", "doi:10.1038/ng.3889", adrn)
mes <- c("MES", "doi:10.1038/ng.3889", mes)
write_delim(as.data.frame(adrn), file = "GSEA/ADRN.txt", delim = '\t')
write_delim(as.data.frame(mes), file = "GSEA/MES.txt", delim = '\t')

GO_file <- "GSEA/MSigDB/alldata/M2/m2.all.v2022.1.Mm.symbols.ADRNMES.gmt"
gsea_res = GSEA(gene_list, GO_file, pval = 0.05)
dim(gsea_res$Results)

plot_geneset_clusters( gs_results = gsea_res$Results[gsea_res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )

#Ontology (M5)
GO_file <- "GSEA/MSigDB/alldata/M5/m5.all.v2022.1.Mm.symbols.gmt"
gsea_res = GSEA(gene_list, GO_file, pval = 0.05)
dim(gsea_res$Results)

plot_geneset_clusters( gs_results = gsea_res$Results[gsea_res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )



#New GSEA - for gsea plots
BiocManager::install("clusterProfiler")
library(org.Mm.eg.db)
library(clusterProfiler)

gse <- gseGO(gene_list, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Mm.eg.db", eps = 1e-300)
saveRDS(gse, file = "GSEA/clusterProfiler_GSEA_output.RDS")
gseaplot(gse, geneSetID = 1)

#Better

library(msigdbr)
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% select(gs_name, gene_symbol)

em2 <- GSEA(geneList = gene_list, TERM2GENE = m_t2g)

library(enrichplot)
gseaplot2(em2, geneSetID = 1, title = em2@result$Description[1])
gseaplot2(em2, geneSetID = 8, title = em2@result$Description[8])
gseaplot2(em2, geneSetID = 5, title = em2@result$Description[5])


#Exporting data for GSEA sample level permutation
eds <- as.matrix(normalized_counts)
rownames(eds) <- str_to_title(tolower(rownames(eds)))
write.table(as.data.frame(eds), file = "GSEA/alldata/eds.txt", sep = "\t")

mesadrn.chip <- as.data.frame(cbind(paste0("Probe_", seq(1:length(rownames(counts)))), rownames(counts), "na"))
colnames(mesadrn.chip) <- c("Probe Set ID", "Gene Symbol", "Gene Title")
write.table(mesadrn.chip, file = "GSEA/alldata/mesadrn.chip", sep = "\t")

mesadrn.gmx <- read.delim(file = "GSEA/Neuroendocrine/MES_ADRN.gmx")
mesadrn.gmx$MES <- str_to_title(tolower(mesadrn.gmx$MES))
mesadrn.gmx$ADRN <- str_to_title(tolower(mesadrn.gmx$ADRN))
write.table(mesadrn.gmx, file = "GSEA/alldata/mesadrn.gmx", sep = "\t")
