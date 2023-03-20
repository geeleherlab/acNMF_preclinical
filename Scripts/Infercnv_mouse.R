# Load data and libraries
library(infercnv)
library(dplyr)


## GEMM data

mouse <- readRDS("mouse_data/NB_DBHiCre.readcounts.RDS")
annotation_mouse <- read.csv("mouse_data/NB_DBHiCre_metadata.csv",
                            header = TRUE, sep = ",")

GEP_ind <- readRDS("mouse_data/mouse_gepindex.RDS")
cNMF <- read.table("mouse_data/cNMF_split0.k_59.dt_0_10.consensus.txt",
                  header = T, row.names = 1, sep = '\t')

# load group ids
normal_GEPs <- c(6,7,8,11,18,19,20,23,28,29,35,38,41,45) # GEPs 7, 38 and 41 are endothelial
cancer_GEPs <- c(5,15,17,26,27,30,34,42,44)
unknown_GEPs <- c(2,3,4,9,10,12,13,21,24,31,32,33,36,37,39)

# find highly expressed cells
cell_GEP <- data.frame(matrix(data = 0, nrow = nrow(cNMF), ncol = length(c(normal_GEPs,cancer_GEPs,unknown_GEPs))))
colnames(cell_GEP) <- paste0("GEP_",c(normal_GEPs,cancer_GEPs,unknown_GEPs))
for (i in c(normal_GEPs,cancer_GEPs,unknown_GEPs)){
  if (length(GEP_ind[[i]]$Split0)>1){
    cell_GEP[,paste0("GEP_",i)] <- 
      rowMeans(cNMF[,GEP_ind[[i]]$Split0])
  }
  else {
    cell_GEP[,paste0("GEP_",i)] <- cNMF[,GEP_ind[[i]]$Split0]
  }
}



cell_GEP <- cell_GEP %>% rowwise() %>% 
  mutate(m = max(c_across(GEP_6:GEP_39)), id = names(.)[which.max(c_across(GEP_6:GEP_39))]) %>%
  select(m,id)

cell_GEP <- data.frame(cell_GEP)
rownames(cell_GEP) <- rownames(cNMF)
cell_GEP <- cell_GEP[cell_GEP$m >100,]   

# subsetting and preparing for inferCNV

colnames(cell_GEP) <- c("h_mean","Annotation")
gene_order <- read.table("mouse_data/mouse_gene_order.txt", 
                        header = F, row.names = 1, sep = '\t')


vehicle <- c("NB831", "NB837",  "NB856", "NB864")
drug_treated <- c("NB839", "NB847", "NB849", "NB853", "NB883", "NB887")

# Self-processed data ------------------------------
annot_mouse <- annotation_mouse[rownames(annotation_mouse) %in% rownames(cell_GEP),]
annot <- data.frame(cell_GEP$Annotation)
rownames(annot) <- rownames(cell_GEP)

for (name in vehicle){
  print(name)
  set.seed(1334)	 
  if (name %in% vehicle){
    cell_ls <- rownames(annot_mouse)[annot_mouse$orig.ident == name] # sample specific cells
    vehicle_ls <- rownames(annot_mouse)[(annot_mouse$orig.ident %in% vehicle)] # vehicle treated cells
    normal_ls <- rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", normal_GEPs)] # normal cells
    model_ls <- unique(c(cell_ls,intersect(vehicle_ls, normal_ls)))
    
   
    train_normal_ls <- model_ls[model_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% 
                                                                                c(paste0("GEP_", c(6,8,11,18,19,20,23,28,29,35,45)))]]# all immune cells from all same-type mice
    valid_normal_ls <- model_ls[model_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% 
                                                                  c(paste0("GEP_", c(7,38,41)))]]# all endothelial cells from all same-type mice

    cancer_ls <- cell_ls[cell_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", cancer_GEPs)]]
    proportion <- 0.75
    split <- sample(c(rep(0, (1-proportion)*length(cancer_ls)), rep(1,length(cancer_ls)*proportion)))
    train_cancer_ls <- cancer_ls[split==1]
    valid_cancer_ls <- cancer_ls[split==0]
    
    ambiguous_ls <- cell_ls[cell_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", unknown_GEPs)]]
    train <- mouse[,colnames(mouse) %in% model_ls]
  } else{ # drug treated
    cell_ls <- rownames(annot_mouse)[annot_mouse$orig.ident == name] # sample specific cells
    vehicle_ls <- rownames(annot_mouse)[(annot_mouse$orig.ident %in% drug_treated)] # drug treated cells
    normal_ls <- rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", normal_GEPs)] # normal cells
    model_ls <- unique(c(cell_ls,intersect(vehicle_ls, normal_ls)))
    
   
    train_normal_ls <- model_ls[model_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% 
                                                                  c(paste0("GEP_", c(6,8,11,18,19,20,23,28,29,35,45)))]]# all immune cells from all same-type mice
    valid_normal_ls <- model_ls[model_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% 
                                                                  c(paste0("GEP_", c(7,38,41)))]]# all endothelial cells from all same-type mice
    
    cancer_ls <- cell_ls[cell_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", cancer_GEPs)]]
    proportion <- 0.75
    split <- sample(c(rep(0, (1-proportion)*length(cancer_ls)), rep(1,length(cancer_ls)*proportion)))
    train_cancer_ls <- cancer_ls[split==1]
    valid_cancer_ls <- cancer_ls[split==0]
    
    ambiguous_ls <- cell_ls[cell_ls %in% rownames(cell_GEP)[cell_GEP$Annotation %in% paste0("GEP_", unknown_GEPs)]]
    train <- mouse[,colnames(mouse) %in% model_ls]
  }
    
    annot_mouse_GEP <- data.frame(cell_GEP[rownames(cell_GEP) %in% colnames(train),"Annotation"])
    rownames(annot_mouse_GEP) <- colnames(train)
    annotable <- table(annot_mouse_GEP$cell_GEP.rownames.cell_GEP...in..colnames.train....Annotation..) 
    for (i in 1:length(annotable)){
      if (annotable[i] == 1){
        annot_mouse_GEP <- subset(annot_mouse_GEP, cell_GEP.rownames.cell_GEP...in..colnames.train....Annotation..
                                 != names(annotable[i]))
      }
    }
    
    cell_GEP_filtered <- cell_GEP[rownames(cell_GEP) %in% rownames(annot_mouse_GEP),]
    annotable <- table(annot_mouse_GEP$cell_GEP.rownames.cell_GEP...in..colnames.train....Annotation..) 
    refs <- intersect(paste0("GEP_", c(6,8,11,18,19,20,23,28,29,35,45)), names(annotable))
    
if (file.exists(paste0("mouse_data/Infercnv_mouse_by_sample_",name))){
 
} else {
  # create a new sub directory inside
  # the main path
  dir.create(file.path(paste0("mouse_data/Infercnv_mouse_by_sample_",name)))
}

    infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = train,
    annotations_file = annot_mouse_GEP,
    delim = "\t",
    gene_order_file = gene_order,
    ref_group_names = refs 
  )
  
    infercnv_obj1 <- infercnv::run(infercnv_obj,
                                 cutoff = .1,
                                 out_dir = paste0("mouse_data/Infercnv_mouse_by_sample_", name,"/prelim"),
                                 cluster_by_groups = T,
                                 denoise = T,
                                 HMM = F,
                                 num_threads = 32,
                                 analysis_mode = "samples",
                                 BayesMaxPNormal = .3,
                                 diagnostics = T,
                                 resume_mode = T,
                                 noise_logistic = T
                                 )
                                 
    infercnv_obj2 <- infercnv::run(infercnv_obj,
                               cutoff = .1,
                               out_dir = paste0("mouse_data/Infercnv_mouse_by_sample_", name,"/HMM"),
                               cluster_by_groups = T,
                               denoise = T,
                               HMM = T,
                               num_threads = 32,
                               analysis_mode = "samples",
                               BayesMaxPNormal = .3,
                               diagnostics = T,
                               resume_mode = T,
                               noise_logistic = T
    )
saveRDS(cell_GEP_filtered, paste0("mouse_data/Infercnv_mouse_by_sample_", name, "/cell_GEP_annot_filtered.RDS"))
saveRDS(list(train_normal_ls, train_cancer_ls, valid_normal_ls, valid_cancer_ls, ambiguous_ls), paste0("mouse_data/Infercnv_mouse_by_sample_", name, "/cellnames.RDS"))
}
  
