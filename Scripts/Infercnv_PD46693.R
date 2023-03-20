library(infercnv)
library(readr)
library(tidyverse)
# load data
split0 <- readRDS(paste0(getwd(),"/GOSH_data/cnmf_input_benchmark_split0.RDS"))
split1 <- readRDS(paste0(getwd(), "/GOSH_data/cnmf_input_benchmark_split1.RDS"))
allData <- rbind(split0, split1)
annotation <- read_tsv("GOSH_data/GOSH_metadata.tsv")

annotation <- annotation[annotation$SampleName == "PD46693",]

annot <- annotation[,c(1,3)]
annot <- annot %>% column_to_rownames(.,var = "...1")

gene_order <- read.table("GOSH_data/gencode_v19_gene_pos.txt", 
                        header = F, row.names = 1, sep = '\t')

if (file.exists("GOSH_data/Infercnv_GOSH_PD46693")){
  
} else {
  dir.create(file.path("GOSH_data/Infercnv_GOSH_PD46693"))
}

cell_ls <- rownames(annot)# sample specific cells
train_normal <- annotation[annotation$Annotation == "Leukocytes" & annotation$SampleName == "PD46693",]
train_normal_ls <- train_normal$...1 # all immune cells of PD46693
valid_normal <- annotation[annotation$Annotation == "Endothelium" & annotation$SampleName == "PD46693",]
valid_normal_ls <- valid_normal$...1 # all endothelial cells of PD46693

cancers <- annotation[annotation$Annotation %in% c("Tumour cluster 1", "Tumour cluster 2", 
                                                  "Tumour cluster 3") & annotation$SampleName == "PD46693",]
cancer_ls <- cancers$...1 # all cancer cells of PD46693 

set.seed(1334)
proportion <- 0.75
split <- sample(c(rep(0, (1-proportion)*length(cancer_ls)), rep(1,length(cancer_ls)*proportion)))
train_cancer_ls <- cancer_ls[split==1]
valid_cancer_ls <- cancer_ls[split==0]

ambiguous <- annotation[annotation$Annotation == "Mesenchyme" & annotation$SampleName == "PD46693",]
ambiguous_ls <- ambiguous$...1

sub_Data <- allData[rownames(allData) %in% rownames(annot),]
sub_Data_t <- t(sub_Data)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = sub_Data_t,
  annotations_file = annot,
  delim = "\t",
  gene_order_file = gene_order,
  ref_group_names = c("Leukocytes")
)

infercnv_obj1 <- infercnv::run(infercnv_obj,
                             cutoff = .1,
                             out_dir = paste0("GOSH_data/Infercnv_GOSH_PD46693/prelim"),
                             cluster_by_groups = T,
                             denoise = T,
                             HMM = F,
                             num_threads = 32,
                             analysis_mode = "samples",
                             # diagnostics = T,
                             resume_mode = T,
                             noise_logistic = T
)


infercnv_obj2 <- infercnv::run(infercnv_obj,
                             cutoff = .1,
                             out_dir = paste0("GOSH_data/Infercnv_GOSH_PD46693/HMM"),
                             cluster_by_groups = T,
                             denoise = T,
                             HMM = T,
                             num_threads = 32,
                             BayesMaxPNormal = .3,
                             diagnostics = T,
                             analysis_mode = "samples",
                             resume_mode = T,
                             noise_logistic = T
)
saveRDS(list(train_normal_ls, train_cancer_ls, valid_normal_ls, valid_cancer_ls, ambiguous_ls), 
        paste0("GOSH_data/cellnames.RDS"))
