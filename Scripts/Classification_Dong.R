# load packages
library(data.table)
library(readr)
library(tidyverse)
library(dplyr)
library(caret)
library(randomForest)
library(e1071)
library(car)
library(glmnet)
library(ggplot2)
library(svglite)
library(RColorBrewer)
# set working directory
setwd("")
sample_ls <- c("T230", "T214", "T69", "T230")
# TPM
geneLength <- colMeans(data_sample)
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
  
}

for (sample in sample_ls){
  # load cell list
  cells <- readRDS(paste0("Dong_data/Infercnv_Dong_by_sample_", sample, "/cellnames.RDS"))
  train_normal_ls <- cells[[1]]
  train_cancer_ls <- cells[[2]]
  valid_normal_ls <- cells[[3]]
  valid_cancer_ls <- cells[[4]] 
  ambiguous_ls <- cells[[5]]

  # load InferCNV processed data and predicted regions----
  CNV_regions <- read.table(paste0("Dong_data/Infercnv_Dong_by_sample_", sample,
                                  "/HMM/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.3.pred_cnv_regions.dat"),
                           header = T)
  infercnv_obj <- readRDS(paste0("Dong_data/Infercnv_Dong_by_sample_",sample,"/prelim/run.final.infercnv_obj"))
  CNV_pred <- infercnv_obj@expr.data
  gene_order <- read.table("Dong_data/gencode_v19_gene_pos.txt", 
                          header = F, row.names = 1, sep = '\t')

  # filter out cells that are removed from inferCNV
  train_normal_ls <- intersect(train_normal_ls, colnames(CNV_pred))
  train_cancer_ls <- intersect(train_cancer_ls, colnames(CNV_pred))
  valid_normal_ls <- intersect(valid_normal_ls, colnames(CNV_pred))
  valid_cancer_ls <- intersect(valid_cancer_ls, colnames(CNV_pred))
  ambiguous_ls <- intersect(ambiguous_ls, colnames(CNV_pred))

  CNV_pred <- as.data.frame(CNV_pred)
  CNV_pred <- cbind(gene_order$V2[rownames(gene_order) %in% rownames(CNV_pred)], CNV_pred)
  colnames(CNV_pred)[1] <- c("chr")


  inferCNV_grouped <- data.frame(matrix(nrow = ncol(CNV_pred)-1, ncol = 0))
  rownames(inferCNV_grouped) <- colnames(CNV_pred)[-1]
  for(i in 1:22) {
    ss_region <- CNV_regions[CNV_regions$chr==paste0("chr",i),c("cnv_name","start","end")]
    ss_prob <- CNV_pred[CNV_pred$chr==paste0("chr",i),]
    if (nrow(ss_region)){
      for (j in 1:nrow(ss_region)){
        ls <- rownames(gene_order)[gene_order$V2 == paste0("chr",i) & 
                                    gene_order$V3>=ss_region$start[j] & 
                                    gene_order$V4 <=ss_region$end[j]]
        tmp <- data.frame(colMeans(ss_prob[rownames(ss_prob)%in%ls,names(ss_prob) != "chr"]))
        rownames(tmp) <- colnames(ss_prob)[-1]
        colnames(tmp) <- ss_region$cnv_name[j]
        inferCNV_grouped[ss_region$cnv_name[j]] <- tmp
      }
    }
  }
  # remove chr6
  inferCNV_grouped <- inferCNV_grouped %>% select(!starts_with("chr_6"))
  # reformat col names
  colnames(inferCNV_grouped) <- sapply(colnames(inferCNV_grouped), str_replace, '-', '_')

  # load raw data and create binning features of 50 genes----
  Dong <- readRDS("Dong_data/GSE137804.readcounts.RDS")

  data_sample <- Dong[,colnames(Dong) %in% c(train_normal_ls,train_cancer_ls, valid_normal_ls,
                                                 valid_cancer_ls, ambiguous_ls)]
  data_sample <- data_sample[rownames(data_sample) %in% rownames(gene_order),]

  tpm <- data.frame(matrix(nrow = nrow(data_sample), ncol = ncol(data_sample)))
  rownames(tpm) <- rownames(data_sample)
  colnames(tpm) <- colnames(data_sample)

  for (i in 1:ncol(data_sample)){
    tpm[,i] <- countToTpm(data_sample[,i], geneLength[i])
  }

  # Z-score normalization
  tpm_std <- as.data.frame(tpm) %>% mutate_all(~(scale(.) %>% as.vector))

  # binning
  tpm_std <- cbind(gene_order[rownames(gene_order) %in% rownames(tpm_std),], tpm_std)
  colnames(tpm_std)[1] <- c("chr")
  tpm_std <- tpm_std[,!(colnames(tpm_std) %in% c("V3","V4"))]

  interval <- 50
  raw_grouped <- data.frame(matrix(nrow = ncol(data_sample), ncol = 0))
  rownames(raw_grouped) <- colnames(tpm_std)[-1]
  for(i in 1:22) {
    ss <- tpm_std[tpm_std$chr==paste0("chr",i),]
    if (nrow(ss)){
      my_cuts <- cut(1:nrow(ss),
                    breaks = unique(c(seq(0,nrow(ss),interval), nrow(ss))))
      for(j in 1:length(levels(my_cuts))) {
        tmp <- data.frame(colMeans(ss[my_cuts == levels(my_cuts)[j], -1]))
        rownames(tmp) <- colnames(ss)[-1]
        colnames(tmp) <- paste0("chr_",i,"_region_",j)
        raw_grouped[paste0("chr_",i,"_region_",j)] <- tmp
      }
    }
  }
  raw_grouped <- raw_grouped %>% select(!starts_with("chr_6"))


  # classification----
  df <- data.frame()
  # binning features----
  train_cancers <- raw_grouped[rownames(raw_grouped) %in% train_cancer_ls,]
  train_cancers$response <- 1
  valid_cancers <- raw_grouped[rownames(raw_grouped) %in% valid_cancer_ls,]
  valid_cancers$response <- 1

  train_normals <- raw_grouped[rownames(raw_grouped) %in% train_normal_ls,]
  train_normals$response <- 0
  valid_normals <- raw_grouped[rownames(raw_grouped) %in% valid_normal_ls,]
  valid_normals$response <- 0
  # test set
  unknown <- raw_grouped[rownames(raw_grouped) %in% ambiguous_ls,]

  train <- rbind(train_cancers, train_normals)
  valid <- rbind(valid_cancers, valid_normals)
  # downsampling (library Caret)
  train <- downSample(train, as.factor(train$response))
  train <- train[,names(train) != "Class"]

  ## fitting glm - manual ----
  mylogit <- glmnet(train[, !(names(train)%in% c("response"))], train$response,
                    family = "binomial")
  glm_pred_manual <- predict(mylogit,as.matrix(valid[,!(names(valid) %in% c("response"))]),
                            s = 0.01, type = "response")
  glm_classes <- ifelse(glm_pred_manual > 0.5, 1, 0)
  glm_confusion <- confusionMatrix(factor(glm_classes), factor(valid$response), positive='1') 
  # saveRDS(mylogit, file = paste0("Models/glm_",sample,"_manual.rda"))
  
  ## fitting random forest
  forest <- randomForest(factor(response)~., data=train)
  forest_pred <- predict(forest, valid)
  rf_confusion <- confusionMatrix(factor(forest_pred), factor(valid$response))
  # saveRDS(forest, file = paste0("Models/rf_",sample,"_manual.rda"))
  
  
  ## fitting support vector machine
  svmfit <- svm(factor(response) ~ ., data = train,
               kernel = "linear",  cost = 5, scale = TRUE)
  ypred <- predict(svmfit, valid)
  svm_confusion <- confusionMatrix(factor(ypred), factor(valid$response), positive='1')
  # saveRDS(svmfit, file = paste0("Models/svm_",sample,"_manual.rda"))
  
  d <- data.frame(rbind(c("GLM","Manual",sample,as.numeric(glm_confusion$byClass["Balanced Accuracy"])),
             c("RF","Manual",sample,as.numeric(rf_confusion$byClass["Balanced Accuracy"])), 
             c("SVM","Manual", sample,as.numeric(svm_confusion$byClass["Balanced Accuracy"]))))
  df <- rbind(df,d)

  # prediction
  glm_test_manual <- predict(mylogit, as.matrix(unknown),
                     s = 0.01, type = "response")
  glm_class <- ifelse(glm_test_manual > 0.5, 1, 0)
  
  
  # InferCNV features----
  train_cancers <- inferCNV_grouped[rownames(inferCNV_grouped) %in% train_cancer_ls,]
  train_cancers$response <- 1
  valid_cancers <- inferCNV_grouped[rownames(inferCNV_grouped) %in% valid_cancer_ls,]
  valid_cancers$response <- 1
  
  train_normals <- inferCNV_grouped[rownames(inferCNV_grouped) %in% train_normal_ls,]
  train_normals$response <- 0
  valid_normals <- inferCNV_grouped[rownames(inferCNV_grouped) %in% valid_normal_ls,]
  valid_normals$response <- 0
  # test set
  unknown <- inferCNV_grouped[rownames(inferCNV_grouped) %in% ambiguous_ls,]
  
  train <- rbind(train_cancers, train_normals)
  valid <- rbind(valid_cancers, valid_normals)
  # downsampling
  train <- downSample(train, as.factor(train$response))
  train <- train[,names(train) != "Class"]
  

  ## fitting glm - infercnv----
  mylogit <- glmnet(train[, !(names(train)%in% c("response"))], train$response,
                    family = "binomial")
  glm_pred_infercnv <- predict(mylogit,as.matrix(valid[,!(names(valid) %in% c("response"))]),
                              s = 0.01, type = "response")
  glm_classes <- ifelse(glm_pred_infercnv > 0.5, 1, 0)
  glm_confusion <- confusionMatrix(factor(glm_classes), factor(valid$response), positive='1')
  # saveRDS(mylogit, file = paste0("Models/glm_",sample,"_infercnv.rda"))
  
  
  ## fitting random forest
  forest <- randomForest(factor(response)~., data=train)
  forest_pred <- predict(forest, valid)
  rf_confusion <- confusionMatrix(factor(forest_pred), factor(valid$response))
  # saveRDS(forest, file = paste0("Models/rf_",sample,"_infercnv.rda"))
  
  
  ## fitting support vector machine
  svmfit <- svm(factor(response) ~ ., data = train,
               kernel = "linear",  cost = 5, scale = TRUE)
  ypred <- predict(svmfit, valid)
  svm_confusion <- confusionMatrix(factor(ypred), factor(valid$response), positive='1')
  # saveRDS(svmfit, file = paste0("Models/svm_",sample,"_infercnv.rda"))
  
  d <- data.frame(rbind(c("GLM","InferCNV",sample,as.numeric(glm_confusion$byClass["Balanced Accuracy"])),
                       c("RF","InferCNV",sample,as.numeric(rf_confusion$byClass["Balanced Accuracy"])), 
                       c("SVM","InferCNV", sample,as.numeric(svm_confusion$byClass["Balanced Accuracy"]))))
  df <- rbind(df,d)


  # prediction
  glm_test_infercnv <- predict(mylogit, as.matrix(unknown),
                     s = 0.01, type = "response")
  glm_class <- ifelse(glm_test_infercnv > 0.5, 1, 0)

  summ <- df %>%
    group_by(X1,X2) %>%
    summarize(Accuracy = mean(as.numeric(X4)))
  summ <- data.frame(summ, c(rep(sample,6)))
  # write_csv(summ, paste0("Figure_human/model_acc.csv"), append= TRUE)


  # plot
  df <- data.frame(
    Sample = rep(sample, length(glm_pred_manual)*2+ length(glm_test_manual)*2),
    Feature = c(rep("Manual", length(glm_pred_manual)), rep("InferCNV", length(glm_pred_infercnv)),
                rep("Manual", length(glm_test_manual)), rep("InferCNV", length(glm_test_infercnv))),
    Cell = c(valid$response, valid$response, rep("MES-like", length(glm_test_manual)*2)),
    Prediction = rbind(glm_pred_manual, glm_pred_infercnv, glm_test_manual, glm_test_infercnv)
  )
  colnames(df) <- c("Sample", "Feature", "Celltype", "Prediction")
  # write_csv(df, paste0("Figure/GLM_score_human_" sample, ".csv"), append = TRUE)

  
  df$Celltype <- factor(df$Celltype, levels = c("MES-like", "Cancer", "Normal"))
  ggplot(df, aes(x = Feature, y = as.numeric(Prediction), fill = Celltype)) +
    geom_boxplot(outlier.size = 0.2, lwd=0.2) + 
    xlab("Feature selection approach") + ylab("GLM score") + ggtitle("") +
    scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#1f78b4")) +
    # annotation_custom(g, xmin = 0.8, xmax = 1.5, ymin = 0., ymax = 0.5) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.key = element_rect(colour = NA, fill = NA),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)
    )

  ggsave(paste0("Figure_mouse/glm_score_by_feature_", sample, ".svg"), width = 3, height = 3, units = "in")

}
