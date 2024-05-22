library(Seurat)
library(monocle3)

NB_DBHiCre <- readRDS(file = "/Volumes/clusterhome/rchapple/scRNAseq/mouse/analyses/NMF/cNMF/input/alldata/unfiltered_data/umap/NB_DBHiCre.RDS")


expr <- NB_DBHiCre@assays$RNA@counts
meta <- NB_DBHiCre@meta.data
gene_ann <- as.data.frame(rownames(NB_DBHiCre@assays$RNA@counts))
rownames(gene_ann) <- gene_ann[, 1]
colnames(gene_ann) <- "gene_short_name"

NB_DBHiCre_monocle <- new_cell_data_set(expression_data = as(expr, "sparseMatrix"), cell_metadata = meta, gene_metadata = gene_ann)

#Transfer Seurat embeddings to monocle cell_data_set
reducedDim(NB_DBHiCre_monocle, type = "PCA") <- NB_DBHiCre@reductions$umap@cell.embeddings
#NB_DBHiCre_monocle@preprocess_aux$prop_var_expl <- NB_DBHiCre@reductions$pca@stdev
NB_DBHiCre_monocle@int_colData@listData$reducedDims$UMAP <- NB_DBHiCre@reductions$umap@cell.embeddings
#NB_DBHiCre_monocle@clusters$UMAP_seurat$clusters <- NB_DBHiCre@meta.data$seurat_clusters
NB_DBHiCre_monocle <- cluster_cells(NB_DBHiCre_monocle, reduction_method = "UMAP", resolution = 1e-20)

saveRDS(NB_DBHiCre_monocle, file = "NB_DBHiCre_monocle.RDS")

NB_DBHiCre_monocle <- learn_graph(NB_DBHiCre_monocle, verbose = T)
NB_DBHiCre_monocle <- order_cells(NB_DBHiCre_monocle)

plot_cells(NB_DBHiCre, color_cells_by = "pseudotime")

#temp <- ProjectDim(NB_DBHiCre, reduction = "pca")

################################
#Subset out Neuroblastoma cells
##############################
#NB_DBHiCre_NB <- subset(x = NB_DBHiCre, subset = SingleR.main.labels == "Neurons" | SingleR.main.labels == "Neuroepithelial_cell")
NB_DBHiCre_NB <- subset(x = NB_DBHiCre_NB, subset = SingleR.main.labels == "Neurons")

saveRDS(NB_DBHiCre_NB, file = '/Volumes/clusterhome/rchapple/scRNAseq/mouse/analyses/Monocle3/NB_DBHiCre_NB.RDS')

expr <- NB_DBHiCre_NB@assays$RNA@counts
meta <- NB_DBHiCre_NB@meta.data
gene_ann <- as.data.frame(rownames(NB_DBHiCre_NB@assays$RNA@counts))
rownames(gene_ann) <- gene_ann[, 1]
colnames(gene_ann) <- "gene_short_name"

NB_DBHiCre_NB_monocle <- new_cell_data_set(expression_data = as(expr, "sparseMatrix"), cell_metadata = meta, gene_metadata = gene_ann)

#Transfer Seurat embeddings to monocle cell_data_set
reducedDim(NB_DBHiCre_NB_monocle, type = "PCA") <- NB_DBHiCre_NB@reductions$umap@cell.embeddings
NB_DBHiCre_NB_monocle@int_colData@listData$reducedDims$UMAP <- NB_DBHiCre_NB@reductions$umap@cell.embeddings
NB_DBHiCre_NB_monocle <- cluster_cells(NB_DBHiCre_NB_monocle, reduction_method = "UMAP", resolution = 1e-20)



NB_DBHiCre_NB_monocle <- learn_graph(NB_DBHiCre_NB_monocle, verbose = T)
NB_DBHiCre_NB_monocle <- order_cells(NB_DBHiCre_NB_monocle)

saveRDS(NB_DBHiCre_NB_monocle, file = "/Volumes/clusterhome/rchapple/scRNAseq/mouse/analyses/Monocle3/NB_DBHiCre_NB_monocle.RDS")

plot_cells(NB_DBHiCre_NB_monocle)
NB_DBHiCre_NB_monocle <- order_cells(NB_DBHiCre_NB_monocle)
plot_cells(NB_DBHiCre_NB_monocle, color_cells_by = "pseudotime")

pseudotime_genes <- graph_test(NB_DBHiCre_NB_monocle, neighbor_graph = "principal_graph", cores = 4)
sig_pseudotime_genes <- row.names(subset(pseudotime_genes, q_value == 0))


#MYCN genes
mycgenes <- c("MYCN", "CRABP1", "THUMPD3")
myc_cds <- NB_DBHiCre_NB_monocle[rowData(NB_DBHiCre_NB_monocle)$gene_short_name %in% mycgenes, ]
plot_genes_in_pseudotime(myc_cds, color_cells_by = "drug", min_expr = 0.5)

#Mid genes
midgenes <- c("PHOX2A", "CHGB", "DBH", "DDC", "NEFL", "NNAT")
mid_cds <- NB_DBHiCre_NB_monocle[rowData(NB_DBHiCre_NB_monocle)$gene_short_name %in% midgenes, ]
plot_genes_in_pseudotime(mid_cds, color_cells_by = "drug", min_expr = 0.5)


#Drug response genes
drgenes <- c("TWIST1", "SNCA", "PTN")
dr_cds <- NB_DBHiCre_NB_monocle[rowData(NB_DBHiCre_NB_monocle)$gene_short_name %in% drgenes, ]
plot_genes_in_pseudotime(dr_cds, color_cells_by = "drug", min_expr = 0.5)
