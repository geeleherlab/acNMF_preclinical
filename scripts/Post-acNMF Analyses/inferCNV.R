library(Seurat)
library(SeuratDisk)
library(infercnv)
library(dplyr)

input.dir <- "/home/rchapple/scRNAseq/human/analyses/cNMF/input/GOSH/"

# load GOSH dataset
GOSH_readcounts <- readRDS(paste0(input.dir, "GOSH.readcounts.RDS"))
GOSH_meta <- read.delim2(paste0(input.dir, "GOSH_metadata.csv"), as.is = T, sep = ",")

# load .gtf file with gene ordering (downloaded from Ensembl)
gtf <- readRDS("../../annotation/annotation.gtf")

# remove duplicates
duplicate <- which(duplicated(gtf$gene_name))
new_gtf <- gtf[-duplicate, ]

# convert gene ordering into .txt file
gene_order_df <- data.frame(new_gtf$gene_name, new_gtf$chr, new_gtf$start, new_gtf$end)
write.table(gene_order_df, "gene_order.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
gene_order <- paste0(getwd(), "/gene_order.txt")

sample_names <- c(unique(GOSH_meta$orig.ident))

infercnv_objs <- list()
infercnv_results <- list()

for(name in sample_names){
  sample_ann <- GOSH_meta %>% filter(orig.ident == name & !(SingleR.main.labels %in% c("B_cell", "T_cells", "Macrophage", "Endothelial_cells")))
  anntable <- table(sample_ann$SingleR.main.labels)
  
  #check for instances of annotations with only one cell
  rem_ann <- NULL
  for(i in 1:length(anntable)){
    if(anntable[i] == 1){
      rem_ann <- c(rem_ann, names(anntable[i]))
    }
  }
  
  if(!(is.null(rem_ann))){
    GOSH_meta_filtered <- GOSH_meta %>% filter(!(SingleR.main.labels %in% rem_ann))
    GOSH_sub <- subset(GOSH_meta_filtered, orig.ident == name | SingleR.main.labels %in% c("Endothelial_cells", "Macrophage", "T_cells", "B_cell"))
  }
  else {
    GOSH_sub <- subset(GOSH_meta, SampleID == name | SingleR.main.labels %in% c("Endothelial_cells", "Macrophage", "T_cells", "B_cell"))
  }
  
  GOSH_ann_df <- data.frame(row.names(GOSH_sub), GOSH_sub$SingleR.main.labels)
  
  write.table(GOSH_ann_df, paste0("GOSH_", name, "_annotations.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
  GOSH_annotations <- paste0(getwd(), "/GOSH_", name, "_annotations.txt")
  
  infercnvobj = CreateInfercnvObject(raw_counts_matrix = GOSH_readcounts,
                                     annotations_file = GOSH_annotations,
                                     delim = "\t",
                                     gene_order_file = gene_order,
                                     ref_group_names = c("Endothelial_cells", "Macrophage", "T_cells", "B_cell"))
  infercnv_objs[[name]] <- infercnvobj
  
  infercnvres = infercnv::run(infercnvobj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = paste0("infercnv_GOSH_", name, "_0.3"),  # dir is auto-created for storing outputs
                              cluster_by_groups = TRUE,   # cluster
                              denoise = TRUE,
                              HMM = TRUE,
                              num_threads = 16,
                              BayesMaxPNormal = 0.3,
                              analysis_mode = "subclusters")
  
  infercnv_results[[name]] <- readRDS(paste0("infercnv_GOSH_", name, "_0.3/run.final.infercnv_obj"))
}

saveRDS(infercnv_objs, file = "infercnv_objs.RDS")
saveRDS(infercnv_results, file = "infercnv_results.RDS")
