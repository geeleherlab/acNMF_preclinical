library(Seurat)
library(SeuratDisk)
library(dplyr)

#Load Seurat file
NB_cell_line <- readRDS(file = "/home/rchapple/scRNAseq/NB_cell_line/data/NB_cell_line.merged.RDS")
#Obtain read counts
counts <- as.data.frame(NB_cell_line@assays$RNA@counts)
dim(counts)

#Prepare cNMF input as TSV files
#_______________________________________
countmatrix <- as.matrix(NB_cell_line@assays$RNA@counts)
dim(countmatrix)

#Processing for cNMF compatability
cnmf_input <- t(countmatrix)
dim(cnmf_input)
head(rownames(cnmf_input))
head(colnames(cnmf_input))

saveRDS(cnmf_input, file = "cnmf_input_NB_cell_line.RDS")
write.table(cnmf_input, file = "cnmf_input_NB_cell_line.tsv", quote = FALSE, sep = "\t", col.names = NA)

#Benchmark

#Add split groups
splits <- sample(c("split0", "split1"), size = length(colnames(NB_cell_line)), replace = T)
names(splits) <- colnames(NB_cell_line)
NB_cell_line <- AddMetaData(object = NB_cell_line, metadata = splits, col.name = "split.group")
split.list <- SplitObject(NB_cell_line, split.by = "split.group")
split.list

NB_cell_line.split0 <- split.list$split0
NB_cell_line.split1 <- split.list$split1
dim(NB_cell_line.split0)
dim(NB_cell_line.split1)

metadata <- NB_cell_line@meta.data
write.table(metadata, file = "NB_metadata.csv", quote = FALSE, sep = ",")

write.table(NB_cell_line.split0@meta.data, file = "NB_metadata_split0.csv", quote = FALSE, sep = ",")
write.table(NB_cell_line.split1@meta.data, file = "NB_metadata_split1.csv", quote = FALSE, sep = ",")

#Prepare format for H5AD input to cNMF including splits

#Convert to H5AD
SaveH5Seurat(NB_cell_line, filename = "NB_cell_line.h5Seurat")
Convert("NB_cell_line.h5Seurat", dest = "h5ad", assay = "RNA")

SaveH5Seurat(NB_cell_line.split0, filename = "NB_cell_line.split0.h5Seurat")
Convert("NB_cell_line.split0.h5Seurat", dest = "h5ad", assay = "RNA")

SaveH5Seurat(NB_cell_line.split1, filename = "NB_cell_line.split1.h5Seurat")
Convert("NB_cell_line.split1.h5Seurat", dest = "h5ad", assay = "RNA")

#The H5AD files have to be processed in python before they are compatible with the cNMF package
#The following commands are performed on all of the H5AD files produced from this script

# python
# > import scanpy as sc
# > adata = sc.read_h5ad("filename.h5ad")
# > del adata.raw
# > adata.write("filename.edited.h5ad")
