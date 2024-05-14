library(jaccard)
library(dplyr)
library(igraph)
library(vroom)
library(elbow)
library(ggplot2)

set.seed(1234)

#############################################################################################3
#Interpret the Jaccard Test results for individual dataset benchmarking study
#Obtain the Rank and Jaccard length that maximize the number of communitites in network graph
#Once these paramaters are identified, extract GEPs from cNMF output (in next code block)
###############################################################################################
dataset.names = c("Kildisiute_GOSH", "GSE137804", "Kildisiute_PMC", "Kameneva_human", "Kameneva_mouse", "Jansky", "NB_cell_line", "mouse", "PDX")

#JaccardTest result directories

#Human
gosh.jaccard.dir <- "/home/rchapple/scRNAseq/human/analyses/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/Kildisiute_GOSH/h5ad/results/"
dong.jaccard.dir <- "/home/rchapple/scRNAseq/human/analyses/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/GSE137804/results/"
pmc.jaccard.dir <- "/home/rchapple/scRNAseq/human/analyses/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/Kildisiute_PMC/h5ad/results/"
jansky.jaccard.dir <- "/home/rchapple/scRNAseq/human/analyses/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/Jansky/results/"
#NB_cell_line
cellline.jaccard.dir <- "/home/rchapple/scRNAseq/NB_cell_line/analyses/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/results/"
#Mouse
mouse.jaccard.dir <- "/home/rchapple/scRNAseq/mouse/analyses/NMF/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/results/"
#PDX
pdx.jaccard.dir <- "/home/rchapple/scRNAseq/PDX/analyses/cNMF/meta_analysis/benchmark/jaccard/unfiltered_data/results/"
#Fetal
kameneva.human.jaccard.dir <- "/home/rchapple/scRNAseq/fetal/analyses/cNMF/meta_analysis/unfiltered_data/Kameneva_human/results/"
kameneva.mouse.jaccard.dir <- "/home/rchapple/scRNAseq/fetal/analyses/cNMF/meta_analysis/unfiltered_data/Kameneva_mouse/results/"

jaccard.dirs <- list(gosh.jaccard.dir, dong.jaccard.dir, pmc.jaccard.dir, jansky.jaccard.dir, cellline.jaccard.dir, mouse.jaccard.dir, pdx.jaccard.dir, kameneva.human.jaccard.dir, kameneva.mouse.jaccard.dir)
names(jaccard.dirs) <- dataset.names

ranks = c(seq(2,59), seq(70,150,10))
jaccard.val <- seq(10,100,10)

#This will store the parameters (jaccard length and rank) that maximizes community number 
optimized.community.params <- list()
#This will store the community membership label information
com.labels <- list()

for(d in 1:length(jaccard.dirs)){
    
    jaccardtest_results <- list()
    
    for(k in ranks){
        #These are the results of the calcSplitStats function, which calculates the FWER-corrected jaccard significance test 
        splitstats <- readRDS(file = paste0(jaccard.dirs[d], "jaccardtest_results_K", k, ".RDS"))
        jaccardtest_results[[k]] <- splitstats
    }

    jaccard.val <- as.character(jaccard.val)

    com.df <- NULL

    for(x in ranks){
        u <- NULL
        #print(paste0("Rank: ", x))
        for(y in jaccard.val){
            u <- c(u, jaccardtest_results[[x]][[y]][["jaccard.com"]]) 
        }
        
        u <- cbind(u, x)
        u <- data.frame(u, jaccard.val)
        com.df <- rbind(com.df, u)
    }   

    colnames(com.df) <- c("Community_Number", "Rank", "Jaccard_Length")
    saveRDS(com.df, file = paste0(names(jaccard.dirs[d]), "_jaccardtest_communities.RDS"))

    #Insert plotting function here
    
    #This block of code evaluates which curve has the best inflection point
    jl.eval <- data.frame()
  
    for(x in jaccard.val){
        x <- as.numeric(x)  
        jl.df <- com.df %>% filter(Jaccard_Length == x)
        #Run a loess regression on JL data to smooth the curve
        sl <- jl.df %>% ggplot(aes(Rank, Community_Number)) + geom_point() + geom_smooth(method = "loess", span = 0.3, method.args = list(degree = 1))
        gbuild <- ggplot_build(sl)
        elbow.df <- data.frame(X = gbuild$data[[2]]$x, Y = gbuild$data[[2]]$y)
        #Find inflection point
        ipoint <- elbow(data = elbow.df, plot = F)
        #Since the data was run on smoothed data, retreive the closest rank to the x intercept calculated by elbow package
        optimum.rank <- ranks[which.min(abs(ranks - ipoint$X_selected))]
        #Run linear regression on cost values (X values after inflection point) and determine slope
        ipoint.filtered <- ipoint$data %>% filter(X >= ipoint$X_selected)
        reg <- lm(ipoint.filtered$benefits ~ ipoint.filtered$X)
        jl.eval <- rbind(jl.eval, c(x, optimum.rank, reg$coefficients[2]))
    }
  
    colnames(jl.eval) <- c("JL", "Rank", "Slope")
    #This will determine which JL curve gives the steepest slope for the cost values
    max.slope.params <- jl.eval %>% filter(JL == jl.eval$JL[jl.eval$Slope == min(jl.eval$Slope)])
    comnum <- com.df %>% filter(Jaccard_Length == max.slope.params$JL & Rank == max.slope.params$Rank)
    params <- as.data.frame(cbind(comnum$Community_Number, max.slope.params$Rank, max.slope.params$JL))
    colnames(params) <- c("Community_Number", "Rank", "Jaccard_Length")
    
    optimized.community.params[[names(jaccard.dirs[d])]] <- params   

    #Get the significant matches(indeces) from the Jaccard test
    
    sig.labels <- as.data.frame(jaccardtest_results[[params$Rank]][[as.character(params$Jaccard_Length)]]$sigFwerJaccard)
    sig.labels$row <- paste0("Split0_GEP", sig.labels$row)
    sig.labels$col <- paste0("Split1_GEP", sig.labels$col)
    com.labels[[names(jaccard.dirs[d])]] <- sig.labels
}

saveRDS(optimized.community.params, file = "optimized_community_paramters.RDS")
saveRDS(com.labels, file = "community_labels.RDS")

###############################################################################
#Use network (igraph) results to obtain mean spectra score for each community 
##############################################################################3

dataset.com.graphs <- list()
community_geps <- list()

for(dat in dataset.names){
    print(paste0("Calculating ", dat, " averaged community-wide GEPs"))
    #Build network graph from splits
    dat.graph <- graph_from_data_frame(com.labels[[dat]], directed = F)
    ceb <- cluster_edge_betweenness(dat.graph)
    membership <- sort(membership(ceb))

    dataset.com.graphs[[dat]] <- list(Graph = dat.graph, CEB = ceb)

    #Obtain the labels for the community membership

    gep.index <- list()
    x <- numeric()
    y <- numeric()

    for(i in 1:max(membership)){
        z <- names(membership[membership == i])
        for(j in z){
            if(grepl("Split0", j)){
                temp <- as.numeric(substr(j,11,nchar(j)))
                x <- c(x, temp)
                }else{
                    temp <- as.numeric(substr(j,11,nchar(j)))
                    y <- c(temp, y)
                }
            }
            gep.index[[i]] <- list("Split0" = x, "Split1" = y)
            x <- NULL
            y <- NULL
    }
    
    saveRDS(gep.index, file = paste0(dat, "_gepindex.RDS"))

    #cNMF output
    
    K = optimized.community.params[[dat]]$Rank
    max.com.num = optimized.community.params[[dat]]$Community_Number

    if(dat == "NB_cell_line"){
        data.dir <- "/home/rchapple/scRNAseq/NB_cell_line/analyses/cNMF/output/benchmark/unfiltered_data/"
    }else if(dat == "mouse"){
        data.dir <- "/home/rchapple/scRNAseq/mouse/analyses/NMF/cNMF/output/benchmark/unfiltered_data/"
    }else if(dat == "PDX"){
        data.dir <- "/home/rchapple/scRNAseq/PDX/analyses/cNMF/output/benchmark/unfiltered_data/"
    }else{
        data.dir = paste0("/home/rchapple/scRNAseq/human/analyses/cNMF/output/benchmark/unfiltered_data/", dat, "/")
    }

    split0 <- t(vroom(file = paste0(data.dir, "cNMF_split0/cNMF_split0.gene_spectra_score.k_", K, ".dt_0_10.txt"), col_select = c(-...1), show_col_types = FALSE))
    split1 <- t(vroom(file = paste0(data.dir, "cNMF_split1/cNMF_split1.gene_spectra_score.k_", K, ".dt_0_10.txt"), col_select = c(-...1), show_col_types = FALSE))

    community_vectors <- list()
    
    #Average Z-score ranked GEP for each community
    for(i in 1:length(gep.index)){
        
        #These are the node names in each community
        #These names correspond to columns in the cNMF output
        s0 <- unlist(gep.index[i][[1]][1])
        s1 <- unlist(gep.index[i][[1]][2])

        if(length(s0) == 0 | length(s1) == 0){
            if(length(s0) > length(s1)){
                community_vectors[[i]] <- split0[, s0]
            }else{
                community_vectors[[i]] <- split1[, s1]
            }
        }else{
        
        #Compute the average score for the community
        mean.df.temp <- NULL
        mean.df.temp2 <- NULL
        if(length(s0) == 1){
            #print("Length Split0 = 1")
            if(length(s1) == 1){
                #print("Length Split1 = 1")
                community_vectors[[i]] <- rowMeans(cbind(split0[, s0], split1[, s1]), na.rm = T)
            }else{
                #print("Length Split1 > 1")
                for(j in s1){
                    mean.df.temp <- cbind(mean.df.temp, split1[, j])
                }
                mean.df.temp <- rowMeans(mean.df.temp, na.rm = T)
                community_vectors[[i]] <- rowMeans(cbind(split0[, s0], mean.df.temp), na.rm = T)
            }
        }else if(length(s1) == 1){
            #print("Length Split0 > 1")
            #print("Length Split1 = 1")
            for(j in s0){
                mean.df.temp <- cbind(mean.df.temp, split0[, j])
            }
            mean.df.temp <- rowMeans(mean.df.temp, na.rm = T)
            community_vectors[[i]] <- rowMeans(cbind(mean.df.temp, split1[, s1]), na.rm = T)
        }else{
            #print("Length Split1 > 1")
            for(j in s0){
                mean.df.temp <- cbind(mean.df.temp, split0[, j])
            }
            for(k in s1){
                mean.df.temp2 <- cbind(mean.df.temp2, split1[, k])
            }
    
        mean.df.temp <- rowMeans(mean.df.temp, na.rm = T)
        mean.df.temp2 <- rowMeans(mean.df.temp2, na.rm = T)
        community_vectors[[i]] <- rowMeans(cbind(mean.df.temp, mean.df.temp2), na.rm = T)
        }
        }
    }

    community_vector_mat <- matrix(as.numeric(unlist(community_vectors)), ncol = max.com.num, nrow = nrow(split0))
    rownames(community_vector_mat) <- rownames(split0)
    community_vector_mat_filtered <- community_vector_mat[complete.cases(community_vector_mat), ]

    community_geps[[dat]] <- community_vector_mat_filtered
        
}

saveRDS(dataset.com.graphs, file = "indv_dataset_community_graphs.RDS")

##################################
#Parse PDX GEPs into human or mouse
###################################

pdx.human <- NULL
pdx.mouse <- NULL

for(i in 1:ncol(community_vector_mat_filtered)){
    pdx.sorted <- head(names(sort(community_vector_mat_filtered[, i], decreasing = T)), 25L)
    topgenes.mouse.pctg <- sum(startsWith(pdx.sorted, "mm10") == TRUE) / 25
    if(topgenes.mouse.pctg > 0.9){
        pdx.mouse <- cbind(pdx.mouse, community_vector_mat_filtered[, i])
    }else{
        pdx.human <- cbind(pdx.human, community_vector_mat_filtered[, i])
    }
}

mouse.indx <- which(startsWith(rownames(pdx.human), "mm10") == TRUE)
human.indx <- which(startsWith(rownames(pdx.human), "GRCh38") == TRUE)

pdx.mouse <- pdx.mouse[mouse.indx, ]
pdx.human <- pdx.human[human.indx, ] 

#trim prefix off of gene name
mse.names <- rownames(pdx.mouse)
rownames(pdx.mouse) <- toupper(substr(mse.names,8,nchar(mse.names)))
hum.names <- rownames(pdx.human)
rownames(pdx.human) <- toupper(substr(hum.names,8,nchar(hum.names)))

community_geps[["PDX_Human"]] <- pdx.human
community_geps[["PDX_Mouse"]] <- pdx.mouse
community_geps[["PDX"]] <- NULL

#This is the object that contains all of the GEPs representative of each community for each dataset
saveRDS(community_geps, file = "community_geps.RDS")

#####################################################################################################
#Obtain the shared genes across all datasets, and filter community-wide GEPs based on the common genes
#Needed because the Jaccard test has to be performed on symmetrical vectors
######################################################################################################

allgenes <- NULL
for(i in 1:length(community_geps)){
    allgenes <- c(allgenes, rownames(community_geps[[i]]))
}

duplicate.genes <- function(x){
    dup.genes <- data.frame(table(x))
    dup.genes <- filter(dup.genes, Freq == max(dup.genes$Freq))
    dup.genes <- as.vector(dup.genes[, 1])
}

shared.genes <- duplicate.genes(allgenes)
saveRDS(shared.genes, file = "shared.genes.RDS")

community_geps_filtered <- NULL
for(i in 1:length(community_geps)){
    filt <- community_geps[[i]]
    filt <- filt[rownames(filt) %in% shared.genes, ]
    community_geps_filtered[[i]] <- filt[order(rownames(filt)), ] 
}

names(community_geps_filtered) <- names(community_geps)

saveRDS(community_geps_filtered, file = "community_geps_filtered.RDS")


###############################################################
#Meta-analysis for all pairwise comparisons across all datasets
##############################################################3

#JaccardTest function
calcJaccardStats <- function(dataset_a, dataset_b, jval_a, jval_b){
    jaccardInd <- numeric()
    jaccardIndPs <- numeric()
    
    #perform Jaccard test across all GEPs for the two datasets
    for(i in 1:ncol(community_geps_filtered[[dataset_a]])){
        for(j in 1:ncol(community_geps_filtered[[dataset_b]])){
            GEP_A <- as.numeric(rank(-community_geps_filtered[[dataset_a]][, i]) < jval_a)
            GEP_B <- as.numeric(rank(-community_geps_filtered[[dataset_b]][, j]) < jval_b)

            jaccardTest <- jaccard.test.mca(GEP_A, GEP_B)
            jaccardInd <- c(jaccardInd, jaccardTest$statistics)
            jaccardIndPs <- c(jaccardIndPs, jaccardTest$pvalue)
        }
    }
    
    #multiple testing correction
    fwerJaccard <- p.adjust(jaccardIndPs, method="bonferroni")
    #fwerJaccard <- t(matrix(fwerJaccard, nrow = ncol(community_geps_filtered[[dataset_a]]), ncol = ncol(community_geps_filtered[[dataset_b]])))
    fwerJaccard <- t(matrix(fwerJaccard, nrow = ncol(community_geps_filtered[[dataset_b]]), ncol = ncol(community_geps_filtered[[dataset_a]])))

    #pull out significant hits
    #more stringent FWER threshold (used 0.05 / (n comparisons) ~= 0.001)
    sigFwerJaccard <- as.data.frame(which(fwerJaccard < 0.05, arr.ind = T))
    
    if(nrow(sigFwerJaccard) == 0){
        sigFwerJaccard <- NULL
        return(sigFwerJaccard)
    }else{
        #formatting and adding adjusted pvalues for edge weights
        for(i in 1:nrow(sigFwerJaccard)){
            sigFwerJaccard[i,"pval.adj"] <- fwerJaccard[sigFwerJaccard[i, "row"], sigFwerJaccard[i, "col"]]
        }
        sigFwerJaccard$row <- paste0(dataset_a, "_", sigFwerJaccard$row)
        sigFwerJaccard$col <- paste0(dataset_b, "_", sigFwerJaccard$col)

        return(sigFwerJaccard)
    }
}


#Set up all pairwise comparisons
pairwise.comb <- combn(x = names(community_geps_filtered), 2)

connectivities <- NULL

for(i in 1:ncol(pairwise.comb)){
  
    dataset_a <- pairwise.comb[1, i]
    dataset_b <- pairwise.comb[2, i]
    jval_a <- optimized.community.params[[dataset_a]]$Jaccard_Length
    jval_b <- optimized.community.params[[dataset_b]]$Jaccard_Length
  
    if(startsWith(dataset_a, "PDX")){
        jval_a <- optimized.community.params[["PDX"]]$Jaccard_Length
    }
    if(startsWith(dataset_b, "PDX")){
        jval_b <- optimized.community.params[["PDX"]]$Jaccard_Length
    } 
  
    print(paste0("Comaparison: ", i))
    print(paste0("Dataset_A: ", dataset_a, " JL: ", jval_a))
    print(paste0("Dataset_B: ", dataset_b, " JL: ", jval_b))
    
    sigFwerJaccard <- calcJaccardStats(dataset_a, dataset_b, jval_a, jval_b)
    connectivities <- rbind(connectivities, sigFwerJaccard)

}

saveRDS(connectivities, 'connectivities.RDS')
