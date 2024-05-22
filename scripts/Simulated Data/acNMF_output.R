#This script passes two arguments to the jaccard.R script
#The first argument is the rank for which the jaccard test should be calculated 
ranks = c(seq(2,60), seq(70,200,10))

#The second is the data directory in which the gene spectra scores from cNMF are held
data.dir = "/home/rchapple/scRNAseq/human/analyses/cNMF/output/benchmark/cNMF_simulated_data/results_replicated/deloc_1.00/Seed_13657/"

#This loop parses out LSF jobs such that jaccard significance tests are calculated on single ranks across a range of jaccard lengths
#The total number of jobs submitted to LSF will be equivalent to the number of ranks (established above)
#All output from the jaccard.R script will be directed to outfile_<K>.txt

for(k in ranks){
    outfile = paste0("outfile_K", k)
    cmd = paste0('bsub -P RC -J cNMF_sim_jaccard -R "rusage[mem=1GB]" -oo ', outfile, ".out -eo ", outfile, ".err 'Rscript /home/rchapple/software/jaccardtest/jaccard.R ", k, " ", data.dir, " > ", outfile, ".txt'")
    system(cmd)
}

#Set up the i/o log directory
wait_cmd = "bwait -w 'ended(cNMF_sim_jaccard)'"
system(wait_cmd)
Sys.sleep(20)
system('mkdir logdir/')
system('mv outfile* logdir/')

#Set up the results directory
system('mkdir results/')
system('mv jaccardtest_results* results/')

################################################################################
#After the jaccard results are generated, the following code will plot the results
##############################################################################

library(jaccard)
library(igraph)
set.seed(1234)

#same directory as above
simdata.dir <- "/Volumes/clusterhome/rchapple/scRNAseq/simulation/analyses/cNMF/meta_analyis/benchmark/jaccard/cNMF_simulated_data/replicated_results/deloc_1.00/Seed_"

ranks = c(seq(2,59), seq(70,200,10))
jaccard.val <- seq(10,100,10)
seeds <- c(13657, 18839, 19052) 
           22856, 9485)

sim.df <- NULL

for(seed in seeds){
  benchmark <- list()

  for(k in ranks){
    splitstats <- readRDS(file = paste0(simdata.dir, seed, "/results/jaccardtest_results_K", k, ".RDS"))
    benchmark[[k]] <- splitstats
  }

  jaccard.val <- as.character(jaccard.val)

  #for plotting
  sim.data.df <- NULL

  for(x in ranks){
    u <- NULL
    print(paste0("Rank: ", x))
    for(y in jaccard.val){
      print(y)
      u <- c(u, benchmark[[x]][[y]][["jaccard.com"]]) 
    }
    u <- cbind(u, x)
    u <- data.frame(u, jaccard.val)
    sim.data.df <- rbind(sim.data.df, u)
  }
  
  sim.df <- rbind(sim.df, sim.data.df)
}
    
colnames(sim.df) <- c("Community_Number", "Rank", "Jaccard_Length")

#rank on x axis
sim.df$Rank <- as.integer(sim.df$Rank)
sim.df$Jaccard_Length <- as.character(sim.df$Jaccard_Length)
sim.df$Jaccard_Length <- factor(sim.df$Jaccard_Length, levels = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100", "110", "120", "130", "140", "150"))


simplot <- ggplot(sim.df, aes(x=Rank,y=Community_Number, group = Jaccard_Length, color = Jaccard_Length)) +
  stat_smooth(method="loess", span=0.2, se=TRUE, aes(fill=Jaccard_Length), alpha=0.3) +
  theme_linedraw()

ggsave(file = "/Users/rchapple/Desktop/simplot.svg", plot = simplot)
