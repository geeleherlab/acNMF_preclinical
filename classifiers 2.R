setwd("/Users/xliu38/Documents/")
# load packages -----------------------------------------------------------
library(data.table)
library(infercnv)
library(readr)
library(tidyverse)
library(dplyr)
library(caret)
library(rpart)
library(rpart.plot)
library(tree)
library(randomForest)
library(xgboost)
library(e1071)
library(car)
library(glmnet)
split0 <- readRDS(paste0(getwd(),"/GOSH_data/cnmf_input_benchmark_split0.RDS"))
split1 <- readRDS(paste0(getwd(), "/GOSH_data/cnmf_input_benchmark_split1.RDS"))
allData = rbind(split0, split1)

allData_t = t(allData)
rownames(allData_t) = colnames(allData)
colnames(allData_t) = rownames(allData)

annotation = read_tsv(paste0(getwd(),"/GOSH_data/GOSH_metadata.tsv"))
# remove PD42184 relevant cells from gene expression data
# allData_filtered = allData_t[,colnames(allData_t) %in% annotation$...1[annotation$SampleName!="PD42184"]]
# annotation = annotation[annotation$SampleName != "PD42184",]
annot = annotation[,c(1,3)]
annot = annot %>% column_to_rownames(.,var = "...1")

gene_order = read.table(paste0(getwd(),"/GOSH_data/gencode_v19_gene_pos.txt"), 
                        header = F, row.names = 1, sep = '\t')
allData_t = allData_t[rowSums(allData_t) != 0,]
allData_t = allData_t[rownames(allData_t) %in% rownames(gene_order),]

# load GOSH data ----------------------------------------------------------
split0 <- readRDS(paste0(getwd(),"/GOSH_data/cnmf_input_benchmark_split0.RDS"))
split1 <- readRDS(paste0(getwd(), "/GOSH_data/cnmf_input_benchmark_split1.RDS"))
allData = rbind(split0, split1)

allData_t = t(allData)
rownames(allData_t) = colnames(allData)
colnames(allData_t) = rownames(allData)

annotation = read_tsv(paste0(getwd(),"/GOSH_data/GOSH_metadata.tsv"))
# remove PD42184 relevant cells from gene expression data
# allData_filtered = allData_t[,colnames(allData_t) %in% annotation$...1[annotation$SampleName!="PD42184"]]
# annotation = annotation[annotation$SampleName != "PD42184",]
annot = annotation[,c(1,3)]
annot = annot %>% column_to_rownames(.,var = "...1")

gene_order = read.table(paste0(getwd(),"/GOSH_data/gencode_v19_gene_pos.txt"), 
                        header = F, row.names = 1, sep = '\t')


# Self-processed data ------------------------------
allData_t = allData_t[rowSums(allData_t) != 0,]
allData_t = allData_t[rownames(allData_t) %in% rownames(gene_order),]
# TPM
geneLength = colMeans(allData_t)
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
  
}
tpm = data.frame(matrix(nrow = nrow(allData_t), ncol = ncol(allData_t)))
colnames(tpm) = colnames(allData_t)
for (i in 1:ncol(allData_t)){
  tpm[,i] = countToTpm(allData_t[,i], geneLength[i])
}

rownames(tpm) = rownames(allData_t)
colnames(tpm) = colnames(allData_t)
# Z-score normalization
tpm_std <- as.data.frame(tpm) %>% mutate_all(~(scale(.) %>% as.vector))

# binning
# ss_gene = gene_order[rownames(gene_order) %in% rownames(tpm_std),]
# tpm_std = cbind(ss_gene$V2,tpm_std)
tpm_std = cbind(gene_order[rownames(gene_order) %in% rownames(tpm_std),], tpm_std)
colnames(tpm_std)[1] = c("chr")
tpm_std = tpm_std[,!(colnames(tpm_std) %in% c("V3","V4"))]

interval = 50
raw_grouped = data.frame(matrix(nrow = 6442, ncol = 0))
rownames(raw_grouped) = colnames(tpm_std)[-1]
for(i in 1:22) {
  ss = tpm_std[tpm_std$chr==paste0("chr",i),]
  my_cuts = cut(1:nrow(ss),
                breaks = unique(c(seq(0,nrow(ss),interval), nrow(ss))))
  for(j in 1:length(levels(my_cuts))) {
    tmp = data.frame(colMeans(ss[my_cuts == levels(my_cuts)[j], -1]))
    rownames(tmp) = colnames(ss)[-1]
    colnames(tmp) = paste0("chr_",i,"_region_",j)
    raw_grouped[paste0("chr_",i,"_region_",j)] = tmp
  }
}
raw_grouped = raw_grouped %>% select(!starts_with("chr_6"))


# train-test split --------------------------------------------------------
# training set - positive
cancers = raw_grouped[(annot[rownames(raw_grouped),] %in%
                              c("Tumour cluster 1", "Tumour cluster 2", "Tumour cluster 3")),]
cancers$sample = annotation$SampleName[annotation$...1 %in% rownames(cancers)]
cancers$response = 1
# test set
unknown = raw_grouped[annot[rownames(raw_grouped),] == "Mesenchyme",]
unknown$sample = annotation$SampleName[annotation$...1 %in% rownames(unknown)]

# training set - negative
normals = raw_grouped[(annot[rownames(raw_grouped),] %in%
                              c("Endothelium", "Leukocytes")),]
normals$sample = annotation$SampleName[annotation$...1 %in% rownames(normals)]
normals$response = 0

# set validation
# set.seed(1334)
result = data.frame(matrix(nrow = 100, ncol = 3))
for (rep in 1:100){
  proportion = 0.75 # proportion of training set
  split1 = sample(c(rep(0, (1-proportion)*nrow(normals)), rep(1,nrow(normals)*proportion)))
  normals_train = normals[split1==1,]
  normals_valid = normals[split1==0,]
  split2 = sample(c(rep(0, (1-proportion)*nrow(cancers)), rep(1,nrow(cancers)*proportion)))
  cancers_train = cancers[split2==1,]
  cancers_valid = cancers[split2==0,]
  train = rbind(cancers_train, normals_train)
  valid = rbind(cancers_valid, normals_valid)
  # downsampling
  train = downSample(train, as.factor(train$response))
  train = train[,names(train) != "Class"]
  
  
  ## fitting glm
  mylogit <- glmnet(train[, !(names(train)%in% c("sample","response"))], train$response,
                     family = "binomial")
  glm_pred = predict(mylogit,as.matrix(valid[,!(names(valid)%in%c("response","sample"))]),
                      s = 0.01, type = "response")
  glm_classes <- ifelse(glm_pred > 0.5, 1, 0)
  glm_confusion = confusionMatrix(factor(glm_classes), factor(valid$response), positive='1')
  
 
  
  ## fitting random forest
  forest = randomForest(factor(response)~., data=train[,names(train)!= "sample"])
 
  forest_pred = predict(forest, valid[,names(valid)!="sample"])
  rf_confusion = confusionMatrix(factor(forest_pred), factor(valid$response))
  
  
 
  
  
  ## fitting support vector machine
  svmfit = svm(factor(response) ~ ., data = train[, names(train)!= "sample"],
                kernel = "linear",  cost = 5, scale = TRUE)
  ypred = predict(svmfit, valid[,names(valid)!="sample"])
  svm_confusion = confusionMatrix(factor(ypred), factor(valid$response), positive='1')
  
  
  result[rep,] =  c(glm_confusion$overall["Accuracy"],
                    rf_confusion$overall["Accuracy"],svm_confusion$overall["Accuracy"])
  
}
# inferCNV processed data ------------------------
CNV_regions = read.table(paste0(getwd(), "/GOSH_r/infercnv_GOSH_1018_.3/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.3.pred_cnv_regions.dat"),
                         header = T)
CNV_pred = read.table(paste0(getwd(), "/GOSH_r/infercnv_GOSH_1024_prelim/expr.infercnv.dat"),
                      header = T)
CNV_pred = cbind(gene_order$V2[rownames(gene_order) %in% rownames(CNV_pred)], CNV_pred)
colnames(CNV_pred)[1] = c("chr")

inferCNV_grouped = data.frame(matrix(nrow = 6442, ncol = 0))
rownames(inferCNV_grouped) = colnames(tpm_std)[-1]
for(i in 1:22) {
  ss_region = CNV_regions[CNV_regions$chr==paste0("chr",i),c("cnv_name","start","end")]
  ss_prob = CNV_pred[CNV_pred$chr==paste0("chr",i),]
  if (nrow(ss_region)){
    for (j in 1:nrow(ss_region)){
      ls = rownames(gene_order)[gene_order$V2 == paste0("chr",i) & 
                                  gene_order$V3>=ss_region$start[j] & 
                                  gene_order$V4 <=ss_region$end[j]]
      tmp = data.frame(colMeans(ss_prob[rownames(ss_prob)%in%ls,names(ss_prob) != "chr"]))
      rownames(tmp) = colnames(ss_prob)[-1]
      colnames(tmp) = ss_region$cnv_name[j]
      inferCNV_grouped[ss_region$cnv_name[j]] = tmp
    }
  }
}
# remove chr6
inferCNV_grouped = inferCNV_grouped %>% select(!starts_with("chr_6"))
# reformat col names
colnames(inferCNV_grouped) = sapply(colnames(inferCNV_grouped), str_replace, '-', '_')
# train-test split --------------------------------------------------------
# training set - positive
cancers = inferCNV_grouped[(annot[rownames(inferCNV_grouped),] %in%
                         c("Tumour cluster 1", "Tumour cluster 2", "Tumour cluster 3")),]
cancers$sample = annotation$SampleName[annotation$...1 %in% rownames(cancers)]
cancers$response = 1
# test set
unknown = inferCNV_grouped[annot[rownames(inferCNV_grouped),] == "Mesenchyme",]
unknown$sample = annotation$SampleName[annotation$...1 %in% rownames(unknown)]

# training set - negative
normals = inferCNV_grouped[(annot[rownames(inferCNV_grouped),] %in%
                         c("Endothelium", "Leukocytes")),]
normals$sample = annotation$SampleName[annotation$...1 %in% rownames(normals)]
normals$response = 0

# set validation
# set.seed(1334)
result2 = data.frame(matrix(nrow = 100, ncol = 3))
for (rep in 1:100){
  proportion = 0.75 # proportion of training set
  split1 = sample(c(rep(0, (1-proportion)*nrow(normals)), rep(1,nrow(normals)*proportion)))
  normals_train = normals[split1==1,]
  normals_valid = normals[split1==0,]
  split2 = sample(c(rep(0, (1-proportion)*nrow(cancers)), rep(1,nrow(cancers)*proportion)))
  cancers_train = cancers[split2==1,]
  cancers_valid = cancers[split2==0,]
  train = rbind(cancers_train, normals_train)
  valid = rbind(cancers_valid, normals_valid)
  # downsampling
  train = downSample(train, as.factor(train$response))
  train = train[,names(train) != "Class"]
  
  
  ## fitting glm
  mylogit <- glmnet(train[, !(names(train)%in% c("sample","response"))], train$response,
                    family = "binomial")
  glm_pred = predict(mylogit,as.matrix(valid[,!(names(valid)%in%c("response","sample"))]),
                     s = 0.01, type = "response")
  glm_classes <- ifelse(glm_pred > 0.5, 1, 0)
  glm_confusion = confusionMatrix(factor(glm_classes), factor(valid$response), positive='1')
  
  
  
  ## fitting random forest
  forest = randomForest(factor(response)~., data=train[,names(train)!= "sample"])
  
  forest_pred = predict(forest, valid[,names(valid)!="sample"])
  rf_confusion = confusionMatrix(factor(forest_pred), factor(valid$response))
  
  
  
  
  
  ## fitting support vector machine
  svmfit = svm(factor(response) ~ ., data = train[, names(train)!= "sample"],
               kernel = "linear",  cost = 5, scale = TRUE)
  ypred = predict(svmfit, valid[,names(valid)!="sample"])
  svm_confusion = confusionMatrix(factor(ypred), factor(valid$response), positive='1')
  
  
  result2[rep,] =  c(glm_confusion$overall["Accuracy"],
                     rf_confusion$overall["Accuracy"],svm_confusion$overall["Accuracy"])
  
}



# load Jansky data-----------------------------------------------
Data = readRDS(paste0(getwd(), "/Jansky/Jansky.readcounts.RDS"))
annotation = read_csv(paste0(getwd(),"/Jansky/Jansky_metadata.csv"))
# load Dong data-------------------------------------------------
Dong = readRDS(paste0(getwd(), "/Dong/GSE137804.readcounts.RDS"))
annotation_Dong = read_csv(paste0(getwd(), "/Dong/GSE137804_metadata.csv"))


# Statistical test --------------------------------------------------------
# list of regions not significantly different
ls_insignificant = c()

for (i in 1:(ncol(cancers)-2)){
  tt = t.test(cancers[,i],normals[,i])
  if (tt$p.value > 0.05){
    ls_insignificant = c(ls_insignificant, colnames(cancers)[i])
  }
}
# Entire dataset ----------------------------------------------------------
set.seed(1334)
proportion = 0.75 # proportion of training set
split1 = sample(c(rep(0, (1-proportion)*nrow(cancers)), rep(1,nrow(cancers)*proportion)))
cancers_train = cancers[split1==1,]
cancers_valid = cancers[split1==0,]
split2 = sample(c(rep(0, (1-proportion)*nrow(normals_filtered)), rep(1,nrow(normals_filtered)*proportion)))
normals_train = normals_filtered[split2==1,]
normals_valid = normals_filtered[split2==0,]
train = rbind(cancers_train, normals_train)
valid = rbind(cancers_valid, normals_valid)

# downsampling
train = downSample(train, as.factor(train$response))
train = train[,names(train) != "Class"]

# logistic regression
mylogit <- glmnet(train[, !(names(train)%in% c("sample","response"))], train$response,
                  family = "binomial")
summary(mylogit)
glm_pred = predict(mylogit,as.matrix(valid[,!(names(valid)%in% c("sample","response"))]),
                   s = 0.001, type = "response")
glm_classes <- as.numeric(glm_pred > 0.5)
confusionMatrix(factor(glm_classes), factor(valid$response), positive='1')

# trees
train$response <- factor(train$response)
tree <- rpart(response ~ ., data = train[, names(train)!= "sample"], 
              method = "class", control = list(maxdepth = 6))
rpart.plot(tree)
pred <- predict(tree, valid[,names(valid)!="sample"], type = "class")
confusionMatrix(factor(pred), factor(valid$response), positive='1')

# forest
forest = randomForest(response~., data=train[,names(train)!= "sample"],
                      control = list(maxdepth = 4),
                      proximity = TRUE)
forest_pred = predict(forest, valid[,names(valid)!="sample"])
confusionMatrix(factor(forest_pred), factor(valid$response))

ytest = predict(forest, unknown[,names(unknown)!="sample"], positive = "1")


# svm - wins
svmfit = svm(factor(response) ~ ., data = train[, names(train)!= "sample"],
             kernel = "linear", cost = 5, scale = TRUE)
ypred = predict(svmfit, valid[,names(valid)!="sample"])
confusionMatrix(factor(ypred), factor(valid$response), positive='1')

ytest = predict(svmfit, unknown[,names(unknown)!="sample"], positive = "1")


# On PD46693 --------------------------------------------------------------
cancers_ss = cancers[cancers$sample == "PD46693",]
normals_ss = normals[normals$sample == "PD46693",]
set.seed(1334)
proportion = 0.75 # proportion of training set
split1 = sample(c(rep(0, (1-proportion)*nrow(cancers_ss)), rep(1,nrow(cancers_ss)*proportion)))
cancers_ss_train = cancers_ss[split1==1,]
cancers_ss_valid = cancers_ss[split1==0,]
split2 = sample(c(rep(0, (1-proportion)*nrow(normals_ss)), rep(1,nrow(normals_ss)*proportion)))
normals_ss_train = normals_ss[split2==1,]
normals_ss_valid = normals_ss[split2==0,]
train = rbind(cancers_ss_train, normals_ss_train)
valid = rbind(cancers_ss_valid, normals_ss_valid)

# downsampling
train = downSample(train, as.factor(train$response))
train = train[,names(train) != "Class"]

# logistic regression
mylogit <- glmnet(response ~ ., data = train[, names(train)!= "sample"], family = "quasibinomial")
summary(mylogit)
glm_pred = predict(mylogit,as.matrix(valid[,!(names(valid)%in%c("response","sample"))]),
                   s = 0.01, type = "response")
glm_classes <- ifelse(glm_pred > 0.5, 1, 0)
confusionMatrix(factor(glm_classes), factor(valid$response), positive='1')

# trees
train$response <- factor(train$response)
tree <- rpart(response ~ ., data = train[, names(train)!= "sample"], 
              method = "class")
rpart.plot(tree)
pred <- predict(tree, valid[,names(valid)!="sample"], type = "class")
confusionMatrix(factor(pred), factor(valid$response), positive='1')

# forest
forest = randomForest(response~., data=train[,names(train)!= "sample"],
                      proximity = TRUE)
forest_pred = predict(forest, valid[,names(valid)!="sample"])
confusionMatrix(factor(forest_pred), factor(valid$response))

# svm
svmfit = svm(factor(response) ~ ., data = train[, names(train)!= "sample"],
             kernel = "radial", cost = 5, scale = TRUE)
ypred = predict(svmfit, valid[,names(valid)!="sample"])
confusionMatrix(factor(ypred), factor(valid$response), positive='1')

ytest = predict(svmfit, unknown[unknown$sample=="PD46693",names(unknown)!="sample"], positive = "1")



# On PD43255 --------------------------------------------------------------
cancers_ss = cancers[cancers$sample == "PD43255",]
normals_ss = normals[normals$sample == "PD43255",]
set.seed(1334)
proportion = 0.75 # proportion of training set
split1 = sample(c(rep(0, (1-proportion)*nrow(cancers_ss)), rep(1,nrow(cancers_ss)*proportion)))
cancers_ss_train = cancers_ss[split1==1,]
cancers_ss_valid = cancers_ss[split1==0,]
split2 = sample(c(rep(0, (1-proportion)*nrow(normals_ss)), rep(1,nrow(normals_ss)*proportion)))
normals_ss_train = normals_ss[split2==1,]
normals_ss_valid = normals_ss[split2==0,]
train = rbind(cancers_ss_train, normals_ss_train)
valid = rbind(cancers_ss_valid, normals_ss_valid)

# logistic regression
mylogit <- glm(response ~ ., data = train[, names(train)!= "sample"], family = "binomial")
summary(mylogit)
glm_pred = predict(mylogit,as.matrix(valid[,!(names(valid)%in%c("response","sample"))]),
                   s = 0.01, type = "response")
glm_classes <- ifelse(glm_pred > 0.5, 1, 0)
confusionMatrix(factor(glm_classes), factor(valid$response), positive='1')

# trees
train$response <- factor(train$response)
tree <- rpart(response ~ ., data = train[, names(train)!= "sample"], 
              method = "class")
rpart.plot(tree)
pred <- predict(tree, valid[,names(valid)!="sample"], type = "class")
confusionMatrix(factor(pred), factor(valid$response), positive='1')

# forest - wins
forest = randomForest(response~., data=train[,names(train)!= "sample"],
                      proximity = TRUE)
forest_pred = predict(forest, valid[,names(valid)!="sample"])
confusionMatrix(factor(forest_pred), factor(valid$response))

# svm
svmfit = svm(factor(response) ~ ., data = train[, names(train)!= "sample"],
             kernel = "radial", cost = 1, scale = TRUE)
ypred = predict(svmfit, valid[,names(valid)!="sample"])
confusionMatrix(factor(ypred), factor(valid$response), positive='1')

ytest = predict(svmfit, unknown[unknown$sample=="PD46693",names(unknown)!="sample"], positive = "1")

