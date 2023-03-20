# Load data and libraries
library(infercnv)
library(dplyr)

## Dong data
Dong = readRDS("Dong_data/GSE137804.readcounts.RDS")
annotation_Dong = read.csv("Dong_data/GSE137804_metadata.csv",
                            header = TRUE, sep = ",")

GEP_ind = readRDS("Dong_data/GSE137804_gepindex.RDS")
cNMF = read.table("Dong_data/cNMF_split0.usages.k_59.dt_0_10.consensus.txt",
                  header = T, row.names = 1, sep = '\t')

# load group ids
normal_GEPs = c(1,4,5,18,20,26,32,35) # GEPs 1 and 18 are endothelial
cancer_GEPs = c(2,10,11,17,19,24,31,36,37 )
unknown_GEPs = c(6,7,25,27,30,34)

# find highly expressed cells
cell_GEP = data.frame(matrix(data = 0, nrow = nrow(cNMF), ncol = length(c(normal_GEPs,cancer_GEPs,unknown_GEPs))))
colnames(cell_GEP) = paste0("GEP_",c(normal_GEPs,cancer_GEPs,unknown_GEPs))
for (i in c(normal_GEPs,cancer_GEPs,unknown_GEPs)){
  if (length(GEP_ind[[i]]$Split0)>1){
    cell_GEP[,paste0("GEP_",i)] = 
      rowMeans(cNMF[,GEP_ind[[i]]$Split0])
  }
  else {
    cell_GEP[,paste0("GEP_",i)] = cNMF[,GEP_ind[[i]]$Split0]
  }
}



cell_GEP = cell_GEP %>% rowwise() %>% 
  mutate(m = max(c_across(GEP_1:GEP_34)), id = names(.)[which.max(c_across(GEP_1:GEP_34))]) %>%
  select(m,id)

cell_GEP = data.frame(cell_GEP)
rownames(cell_GEP) = rownames(cNMF)
cell_GEP = cell_GEP[cell_GEP$m >100,]   

# subsetting and preparing for inferCNV

colnames(cell_GEP) = c("h_mean","Annotation")
gene_order = "Dong_data/gencode_v19_gene_pos.txt"

# Self-processed data ------------------------------
annot_Dong = annotation_Dong[rownames(annotation_Dong) %in% rownames(cell_GEP),]
annot = data.frame(cell_GEP$Annotation)
rownames(annot) = rownames(cell_GEP)

for (name in c("T69", "T200", "T214", "T230")){
  print(name)
  set.seed(1334)	 
  cell_ls = rownames(annot_Dong)[annot_Dong$SampleID == name] # sample specific cells
  vehicle_ls = rownames(annot_Dong)[(annot_Dong$orig.ident %in% vehicle)] # vehicle treated cells
  normal_ls = rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", normal_GEPs)] # normal cells
  model_ls = unique(c(cell_ls,intersect(vehicle_ls, normal_ls)))
  model_ls = cell_ls
   
  train_normal_ls = model_ls[model_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% 
                                                                c(paste0("GEP_", c(4,5,20,26,32,35)))]]# all immune cells from all same-type mice
  valid_normal_ls = model_ls[model_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% 
                                                                c(paste0("GEP_", c(1,18)))]]# all endothelial cells from all same-type mice
  
  cancer_ls = cell_ls[cell_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", cancer_GEPs)]]
  proportion = 0.75
  split = sample(c(rep(0, (1-proportion)*length(cancer_ls)), rep(1,length(cancer_ls)*proportion)))
  train_cancer_ls = cancer_ls[split==1]
  valid_cancer_ls = cancer_ls[split==0]
  
  ambiguous_ls = cell_ls[cell_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", unknown_GEPs)]]
  train = Dong[,colnames(Dong) %in% c(train_normal_ls, train_cancer_ls)]
  
  annot_Dong_GEP = data.frame(cell_GEP[rownames(cell_GEP) %in% colnames(train),"Annotation"])
  rownames(annot_Dong_GEP) = colnames(train)
  annotable = table(annot_Dong_GEP$cell_GEP.rownames.cell_GEP...in..colnames.train....Annotation..) 
  for (i in 1:length(annotable)){
    if (annotable[i] == 1){
      annot_Dong_GEP = subset(annot_Dong_GEP, cell_GEP.rownames.cell_GEP...in..colnames.train....Annotation..
                               != names(annotable[i]))
    }
  }
  
  annotable = table(annot_Dong_GEP$cell_GEP.rownames.cell_GEP...in..colnames.train....Annotation..) 
  refs = intersect(paste0("GEP_", c(4,5,20,26,32,35)), names(annotable))
  
  if (file.exists(paste0("Dong_data/Infercnv_Dong_by_sample_",name))){
  } else {
    
    # create a new sub directory inside
    # the main path
    dir.create(file.path(paste0("Dong_data/Infercnv_Dong_by_sample_",name)))
    
  }
  
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix = train,
    annotations_file = annot_Dong_GEP,
    delim = "\t",
    gene_order_file = gene_order,
    ref_group_names = refs 
  )
  
  infercnv_obj1 = infercnv::run(infercnv_obj,
                                cutoff = .1,
                                out_dir = paste0("Dong_data/Infercnv_Dong_by_sample_", name,"/prelim"),
                                cluster_by_groups = T,
                                denoise = T,
                                HMM = F,
                                num_threads = 32,
                                analysis_mode = "samples",
                                BayesMaxPNormal = .3,
                                diagnostics = T,
                                resume_mode = F,
                                #sd_amplifier = 1,
                                noise_logistic = T
  )
  
  infercnv_obj2 = infercnv::run(infercnv_obj,
                                cutoff = .1,
                                out_dir = paste0("Dong_data/Infercnv_Dong_by_sample_", name,"/HMM"),
                                cluster_by_groups = T,
                                denoise = T,
                                HMM = T,
                                num_threads = 32,
                                analysis_mode = "samples",
                                BayesMaxPNormal = .3,
                                diagnostics = T,
                                resume_mode = F,
                                #sd_amplifier = 1,
                                noise_logistic = T
  )
  saveRDS(list(train_normal_ls, train_cancer_ls, valid_normal_ls, valid_cancer_ls, ambiguous_ls), paste0("Dong_data/Infercnv_Dong_by_sample_", name, "/cellnames.RDS"))
}

