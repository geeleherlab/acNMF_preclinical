library(jaccard)
library(igraph)
library(vroom)

set.seed(1234)

#This function computes the Jaccard index and significance test for all programs in each benchmark split
#Once this is completed, a network graph is built using igraph package
#Community detection of the network graph is then calculated using betweeness scores

calcSplitStats <- function(split0, split1, k, jval){
    print("calcSplitStats")
    # for each row of split0, compare it to every row of split1
    jaccardInd <- numeric()
    jaccardIndPs <- numeric()
    
    for(i in 1:ncol(split0)){
        for(j in 1:ncol(split1)){
            GEP_split0 = as.numeric(rank(-split0[, i]) < jval)
            GEP_split1 = as.numeric(rank(-split1[, j]) < jval)
            jaccardTest <- jaccard.test.mca(GEP_split0, GEP_split1)
            jaccardInd <- c(jaccardInd, jaccardTest$statistics)
            jaccardIndPs <- c(jaccardIndPs, jaccardTest$pvalue)
        }
    }
  
    # Calculate FWERs
    k <- as.numeric(k) 
    fwerJaccard <- p.adjust(jaccardIndPs, method="bonferroni")
    fwerJaccard <- t(matrix(fwerJaccard, nrow = k, ncol = k))
  
    # Calculate number of significant findings.
    sigFwerJaccard <- which(fwerJaccard < 0.05, arr.ind = T)
  
    #Calculate community number
    if(nrow(sigFwerJaccard) > 0){
        community.df <- as.data.frame(sigFwerJaccard)
        community.df$row <- paste0("Split0_GEP", community.df$row)
        community.df$col <- paste0("Split1_GEP", community.df$col)
        community.graph <- graph_from_data_frame(community.df, directed = F)
        ceb <- cluster_edge_betweenness(community.graph)
        community_number <- length(communities(ceb))
    }else{
        community_number <- 0
    }
    return(list(jaccardIndPs=jaccardIndPs, jaccardInd=jaccardInd, sigFwerJaccard=sigFwerJaccard, jaccard.com=community_number))
}

#Calculate Jaccard similarity index over range of values
ranks = c(seq(2,60), seq(70,200,10))
data.dir <- "/<path>/<to>/<cNMF_output>"
jaccard.val <- seq(10,100,10)

for(k in ranks){
  #Results list
  #Gets written to RDS file at end of R script
  jaccardtest_results <- list()

  #Read in NMF gene spectra scores  
  print(paste0("Rank: ", k))
  split0 <- t(vroom(file = paste0(data.dir, "cNMF_split0/cNMF_split0.gene_spectra_score.k_", k, ".dt_0_10.txt"), col_select = c(-...1), show_col_types = FALSE))
  split1 <- t(vroom(file = paste0(data.dir, "cNMF_split1/cNMF_split1.gene_spectra_score.k_", k, ".dt_0_10.txt"), col_select = c(-...1), show_col_types = FALSE))
  
#Compute jaccard index for different jaccard lengths
  for(jval in jaccard.val){
      print(paste0("Jaccard Length: ", jval))
      splitstatlist <- list()
      splitstatlist <- calcSplitStats(split0, split1, k, jval) 
      jval.name <- as.character(jval)
      jaccardtest_results[[jval.name]] <- splitstatlist
  }      
  
  #This list contains the computations across all jaccard lengths for the rank provided on the command line
  saveRDS(jaccardtest_results, file = paste0("jaccardtest_results_K", k, ".RDS"))
}

#Plot results
library(ggplot2)
library(dplyr)

data.dir <- "/<path>/<to>/<jaccardtest_results>/"
data.df <- NULL

jaccardtest_results <- list()

for(k in ranks){
  splitstats <- readRDS(file = paste0(data.dir, seed, "/results/jaccardtest_results_K", k, ".RDS"))
  jaccardtest_results[[k]] <- splitstats
}

jaccard.val <- as.character(jaccard.val)

#for plotting
com.df <- NULL

for(x in ranks){
  u <- NULL
  for(y in jaccard.val){
    u <- c(u, jaccardtest_results[[x]][[y]][["jaccard.com"]]) 
  }
  u <- cbind(u, x)
  u <- data.frame(u, jaccard.val)
  com.df <- rbind(com.df, u)
}
    
colnames(com.df) <- c("Community_Number", "Rank", "Jaccard_Length")

#rank on x axis
com.df$Rank <- as.integer(com.df$Rank)
com.df$Jaccard_Length <- as.character(com.df$Jaccard_Length)
com.df$Jaccard_Length <- factor(com.df$Jaccard_Length, levels = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))


plot <- ggplot(com.df, aes(x=Rank,y=Community_Number, group = Jaccard_Length, color = Jaccard_Length)) +
                stat_smooth(method="loess", span=0.2, se=TRUE, aes(fill=Jaccard_Length), alpha=0.3) +
                theme_linedraw()
plot
