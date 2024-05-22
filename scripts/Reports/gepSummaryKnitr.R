#+ setTheProgNum, echo=FALSE, include=FALSE

# Load the genereric data / code for the GOSH dataset (this code is also used in creating Dataset and Community level reports).
source("commonCodeForSummaries.R")

#+ printPageTitle, echo=FALSE, results='asis'
# Print the page title.
thisProgName <- programNames[thisProgNumber]
pageTitle <- paste("## Program: ", thisProgName, ".", sep="")
cat(pageTitle)

cat("\n<br><hr>")

cat(paste0("<br><em>Submit a comment on this gene expression program's interpretation:</em> <a href='mailto:paul.geeleher@stjude.org?subject=Comments on program ", thisProgName, "'>CLICK</a><br>"))

# Load the gene sets for GSA, GMT format, provide as a list for convenience below (can then programatically add as many databases as we want....).
#+ thenamehere, echo=FALSE, message = FALSE, warnings = FALSE, include = FALSE

HallmarkSets <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/h.all.v7.4.symbols.gmt")
CelltypeSets <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/c8.all.v7.4.symbols.gmt")
KeggPathways <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/c2.cp.kegg.v7.4.symbols.gmt")
PositionalGeneSets <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/c1.all.v7.4.symbols.gmt")
TranscriptionFactorTragets <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/c3.tft.v7.4.symbols.gmt")
GoBPs <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/c5.go.bp.v7.4.symbols.gmt")
ImmunologicGeneSets <- GSA.read.gmt("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/c7.immunesigdb.v7.4.symbols.gmt")

geneSetDatabases <- list("mSigDB Cell Types Gene Set"=CelltypeSets, 
                         "mSigDB Hallmark Gene Sets"=HallmarkSets,
                         "KEGG Pathways"=KeggPathways,
                         "CHR Positional Gene Sets"=PositionalGeneSets,
                         "Transcription Factor Targets"=TranscriptionFactorTragets,
                         "GO Biological Processes"= GoBPs,
                         "Immunological Gene Sets"=ImmunologicGeneSets
) # Note, the names in this list will be used in the heading in the resulting HTML.

# Human transcription factors from this Cell paper: https://www.sciencedirect.com/science/article/pii/S0092867418301065

humanTfs <- read.csv("/Users/rchapple/Documents/scRNAseq/cNMF_meta_analysis/data/jaccard/ontology/DatabaseExtract_v_1.01.csv", as.is=TRUE, row.names = 3)
humanTfs <- humanTfs[, -1]

##################################################
## Program specific code starts below.         ##
##################################################

#+ makeplot, echo=FALSE, message = FALSE, warnings = FALSE
i <- 1:length(commsSorted)
fi <- (i - 0.5) / length(commsSorted)
xNorm <- qnorm(fi)

qqlineCords <- function (y, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75), qtype = 7) 
{
  stopifnot(length(probs) == 2, is.function(distribution))
  y <- quantile(y, probs, names = FALSE, type = qtype, na.rm = TRUE)
  x <- distribution(probs)
  if (datax) {
    slope <- diff(x)/diff(y)
    int <- x[1L] - slope * y[1L]
  }
  else {
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
  }
  return(c(int, slope))
}
lineCoOrds <- qqlineCords(commsSorted)


qqDf <- data.frame(commsSorted, xNorm, theNames=commsSorted, theText=paste(length(commsSorted):1, ": ", names(commsSorted), sep=""))
p <- ggplot(qqDf, aes(xNorm, commsSorted, text=theText)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() +
  xlab("Theoretical Normal Distribution (Z-score)") + 
  ylab("Observed Loadings") +
  geom_abline(intercept = lineCoOrds[1], slope=lineCoOrds[2])

#+ r, echo=FALSE
widgetQQ <- suppressWarnings(plotly::ggplotly(p, tooltip="text"))
qqPlotFileName <- paste("qqPlot_", thisProgNumber, ".html", sep="")
htmlwidgets::saveWidget(widgetQQ, qqPlotFileName)

#' **Interactive QQ-plot of gene loadings, averaged over both independent splits of the data<br>**
#' This plot highlights the relative contribution of each gene to the GEP
#+ qqPlotiFrame, results="asis", echo=FALSE
cat(paste("\n<iframe scrolling='no' seamless='seamless' style='border:none' src='", qqPlotFileName, "' width='1200' height='750'></iframe>\n", sep=""))

#+ r2, echo=FALSE, message = FALSE
dfSorted <- qqDf[order(qqDf[,1], decreasing=T), c(3,1)]
gtexString <- "https://gtexportal.org/home/gene/"
dfSorted[,3] <- paste("[", dfSorted[,1], "]", "(", gtexString, dfSorted[,1], ")", sep="")
colnames(dfSorted)[1:2] <- c("Gene", "Loading")

theSymbs <- as.character(rownames(dfSorted))
cols <- c("SYMBOL", "GENENAME", "ENSEMBL") ## /addedCode, ENSEMBL bit (used to query decartes)

symTab <- AnnotationDbi::select(org.Hs.eg.db, keys=theSymbs, columns=cols, keytype="SYMBOL")

symTabNoDups <- symTab[which(!duplicated(symTab[,1])), ]
dftoDisplay <- cbind(dfSorted[, c(3,2)], symTabNoDups[,2], rownames(dfSorted))

## \addedCode Add code for gene counts in genes scoring > 50 on this program.
progNumInSplit0 <- communitiesOther[[thisProgNumber]]$Split0
activityList <- list()
for(j in 1:length(progNumInSplit0))
{
  activityList[[j]] <- split0_H[,progNumInSplit0[j]] # need to provide sampleNames (i.e. patient IDs and Annotation, i.e. Cancer/normal/lymphocte etc.)
}
activityMatS0 <- do.call(cbind, activityList)
rownames(activityMatS0) <- rownames(split0_H)
activityMaxesS0 <- apply(activityMatS0, 1, max)
cellsHGreaterThan50 <- names(activityMaxesS0)[activityMaxesS0 > 50]
meanCountsHOver50 <- round(apply(allData[cellsHGreaterThan50, ], 2, mean), digits=2)
meanTpmHOver50 <- round(apply(allDataTpm[cellsHGreaterThan50, ], 2, mean), digits=2)

colnames(dftoDisplay) <- c("Gene Symbol", "Loading", "Gene Name", "Gene")
rownames(dftoDisplay) <- as.character(nrow(dftoDisplay):1)
gtexLinks <- paste("<a href='https://gtexportal.org/home/gene/", dftoDisplay[,4], "'>GTEx</a>", sep="")
depmapLinks <- paste("<a href='https://depmap.org/portal/gene/", dftoDisplay[,4], "?tab=dependency'>DepMap</a>", sep="")
geneCardsLink <- paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", dftoDisplay[,4], "'>", dftoDisplay[,4], "</a>", sep="")

decartesLink <- paste("<a href='https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/gene/", symTabNoDups[,3], "/in/adrenal'>Descartes</a>", sep="") ## /addedCode

countHighScoringCells <- meanCountsHOver50[as.character(dftoDisplay$"Gene")]
tpmHighScoringCells <- meanTpmHOver50[as.character(dftoDisplay$"Gene")]

#sets rownames to numbers
rn <- seq(1:nrow(dftoDisplay))

dftForKable <- data.frame("Gene"=geneCardsLink, "Loading"=dftoDisplay[,2], "Gene Name"=dftoDisplay[,3], "GTEx"=gtexLinks, "DepMap"=depmapLinks, "Descartes"=decartesLink, "Mean Counts"=countHighScoringCells, "Mean Tpm"=tpmHighScoringCells) ## /addedCode
rownames(dftForKable) <- rn

#' **Top genes driving this program.<br>** 
#' Note: Decartes website is buggy, try refreshing. Also, Decartes fetal adrenal data have been collected at specific time points (89-122 days), all possible cell types of interest may not be represented, do not overinterpret.<br> 
#' The Mean Count column shows the mean read count in cells scoring highly (H > 50) on this gene expression program. <br>
#+ echo=FALSE

kable(dftForKable[1:50, ], row.names=TRUE)
fileNameForGepTable <- paste("GEP_table", thisProgNumber, ".csv", sep="")
write.csv(cbind(dftoDisplay[,c(4, 2, 3)], Counts=countHighScoringCells, TPM=tpmHighScoringCells), row.names=FALSE, file=fileNameForGepTable) ## /addedCode (to include tpms and counts in this file)

#+ topProgramGenesTable, results="asis", echo=FALSE
cat("\n<br><a href='", fileNameForGepTable, "'>Dowload full table</a><br>")


cat("\n<br><hr>")

#+ echo=FALSE


#+ umapTitle, echo=FALSE, results='asis'
umapTitle <- paste0("<strong>UMAP plots showing activity of gene expression program identified in GEP ", thisProgName, ":", "</strong>")
cat(umapTitle)

#+ umapBlock_split0, echo=FALSE, results='asis'
for(j in 1:(length(progNumInSplit0)))
{
  # Subset the UMAP to the cells we want to show
  umapData <- cbind(umapOut[rownames(split0_H), ], split0_H[, progNumInSplit0[j]])
  uMatDataOrd <- umapData[order(umapData[,3]), ] # Need to order these lowest to highest, so the high scoring cells don't get buried and obscured.
  
  # Add the Activity Scores to the already prepared generic Hover Text
  theseCellsHoverText <- paste(genericCellHoverText[rownames(uMatDataOrd)], "\nActivity Score: ", round(uMatDataOrd[, 3], digits=2), sep="")
 
  uMatDataOrd <- cbind.data.frame(uMatDataOrd, theseCellsHoverText)
  colnames(uMatDataOrd) <- c("x", "y", "Activity", "Text")
  
  # Make the plot
  plot <- ggplot(uMatDataOrd, aes(x, y, colour=Activity, text=Text)) + 
    geom_point() +
    theme_bw() +
    labs(title=paste("Split0 ", thisProgName, "_", j, sep=""), x = "UMAP 1", y = "UMAP 2") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")))
  
  # Write the plot out to to a HTML file
  widgetUmap <- plotly::ggplotly(plot, tooltip = "Text")
  htmlFileLoc <- paste("umap_", thisProgNumber, "_split0_", j, ".html", sep="")
  htmlwidgets::saveWidget(widgetUmap, htmlFileLoc)
  
  # Create an iframe in the HTML that will display the UMAP.
  thestr <- paste("<iframe scrolling='no' seamless='seamless' style='border:none' src='", htmlFileLoc , "' width='1200' height='750'></iframe>", sep="")
  cat(thestr)
}

#+ boxplotBlock_split0, echo=FALSE, results='asis'
progNumInSplit0 <- communitiesOther[[thisProgNumber]]$Split0 # get the progam number...
for(j in 1:length(progNumInSplit0))
{
  theseCellsHoverText <- paste(genericCellHoverText[rownames(split0_H)], "\nActivity Score: ", round(split0_H[,progNumInSplit0[j]], digits=2), sep="")
  
  boxplotDf <- data.frame(SampleNames=boxplotCategories[rownames(split0_H), "SampleName"], Categories=boxplotCategories[rownames(split0_H), "Annotation"], Activity=split0_H[,progNumInSplit0[j]], Text=theseCellsHoverText) # need to provide sampleNames (i.e. patient IDs and Annotation, i.e. Cancer/normal/lymphocte etc.)
  
  rownames(boxplotDf) <- rownames(split0_H)
  plot <- ggplot(boxplotDf, aes(x=Categories, y=Activity, color=SampleNames, text=Text)) +
    labs(title=paste("Split0 ", thisProgName, "_", j, sep=""), x ="Original Author Cell Annotation", y = "Activity Score") + 
    geom_boxplot(color="black", outlier.shape = NA, outlier.color = NA) +
    geom_jitter(position=position_jitter(0.2), alpha = 1/2)
  
  # Write the plot out to to a HTML file
  widgetUmap <- plotly::ggplotly(plot, tooltip = "Text")
  htmlFileLoc <- paste("boxplot_", thisProgNumber, "_split0_", j, ".html", sep="")
  htmlwidgets::saveWidget(widgetUmap, htmlFileLoc)
  
  # Create an iframe in the HTML that will display the UMAP.
  thestr <- paste("<iframe scrolling='no' seamless='seamless' style='border:none' src='", htmlFileLoc , "' width='1200' height='750'></iframe>", sep="")
  cat(thestr)
}


cat("\n<br><hr>")


#' **CNV Data procured from inferCNV.<br>**
#' Outer tracks are putative CNV regions (gains = red, losses = blue) for each patient<br>
#' Inner track is expression data representing:<br>
#' &emsp;&emsp;&emsp; The top cells expressing this GEP (purple)<br>
#' &emsp;&emsp;&emsp; Random cells (n =50) from the reference set used in inferCNV (orange)
#+ echo=FALSE

#Read in inferCNV results
infercnv_datadir <- "/Volumes/groups/geelegrp/home/common/Papers/Rich_NBpreclinical/CNVs_Maggie/inferCNV/by_sample_with_ref/GOSH/gosh_0.3/"
#sample.names <- c("PD42752-1", "PD42752-2", "PD43255", "PD46693")
sample.names <- unique(GOSH@meta.data$SampleName)

regionsdf <- NULL

for(name in sample.names){
  
  regions <- read.delim2(file = paste0(infercnv_datadir, "infercnv_gosh_", name, "_0.3/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.3.pred_cnv_regions.dat")) 
  
  #Filter out cells that were previously labeled as leukocytes
  #Only needed for GOSH?????
  regions <- regions %>% filter(!startsWith(cell_group_name, "Leukocytes"))
  regions$sample_id <- name
  regionsdf <- rbind(regionsdf, regions)
}

# Remove losses on Chr6 in test cells
# This is a result of using immune cells as the reference due to HLA copy number differences in these cells
chr6_ind <- NULL
for(i in 1:nrow(regionsdf)){
  if(regionsdf$chr[i] == "chr6" && regionsdf$state[i] < 3){
    chr6_ind <- c(chr6_ind, i)
  }
}
regionsdf <- regionsdf[-chr6_ind, ]
#Add a numeric chromosome column (remove the 'chr' prefix) - this is needed for BioCircos plotting
regionsdf$numeric_chr <- as.numeric(substr(regionsdf$chr, 4, nchar(regionsdf$chr)))

#################################
#######ARC TRACKS for CNV REGIONS
##################################

#Define colors for arcs
#Losses (i.e. states 1 and 2) are blue
#Gains (i.e. states > 3) are red

cols <- NULL
for(i in 1:nrow(regionsdf)){
  if(regionsdf$state[i] == 1){
    cols <- c(cols, "#2166AC")
  }else if(regionsdf$state[i] == 2){
    cols <- c(cols, "#92C5DE")
  }else if(regionsdf$state[i] == 3){
    cols <- c(cols, "#F7F7F7")
  }else if(regionsdf$state[i] == 4){
    cols <- c(cols, "#F4A582")
  }else if(regionsdf$state[i] == 5){
    cols <- c(cols, "#D6604D")
  }else if(regionsdf$state[i] == 6){
    cols <- c(cols, "#B2182B")
  }
}
regionsdf$cols <- cols

#Opacities for arcs
#Greater CN gains/loss are less opaque
ops <- NULL
for(i in 1:nrow(regionsdf)){
  if(regionsdf$state[i] == 1){
    ops <- c(ops, 1)
  }else if(regionsdf$state[i] == 2){
    ops <- c(ops, 0.5)
  }else if(regionsdf$state[i] == 3){
    ops <- c(ops, 0.3)
  }else if(regionsdf$state[i] == 4){
    ops <- c(ops, 0.5)
  }else if(regionsdf$state[i] == 5){
    ops <- c(ops, 0.7)
  }else if(regionsdf$state[i] == 6){
    ops <- c(ops, 1)
  }
}
regionsdf$opacity <- ops

#Define ARC Tracklists

tracklist <- NULL

#Sets the position for the first track
minRadius = 1.16

sample.names.new <- sample.names[2:5]
for(name in sample.names.new){
  
  #Make one track per sample_id
  sample_regions <- regionsdf %>% filter(sample_id == name)
  
  arcs_chromosomes <- as.character(sample_regions$numeric_chr)
  arcs_begin <- sample_regions$start
  arcs_end <- sample_regions$end
  
  maxRadius <- minRadius + 0.04 #width of arc, must be determined empirically for each dataset so all of the tracks fit on the plot
  
  track <- BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end,
                             minRadius = minRadius, maxRadius = maxRadius, 
                             colors =  sample_regions$cols, opacities = sample_regions$opacity,
                             labels = name)
  tracklist <- tracklist + track
  minRadius <- maxRadius + 0.01 #adds gap between tracks, stays the same for all datasets
}


###################################################3
#######SNP TRACKS - Repurposed for Expression Values
####################################################

#Obtain list of genes (filtered by those used in inferCNV analysis), along with their chromosomal positions

#Full list of genes (~55k)
gene_order <- read.delim(paste0(infercnv_datadir, "gene_order.txt"), header = F)
names(gene_order) <- c("Gene", "Chr", "Start", "End")

allgene_orders <- NULL

for(name in sample.names){
  infercnv_obj <- readRDS(paste0(infercnv_datadir, "infercnv_gosh_", name, "_0.3/run.final.infercnv_obj"))  
  allgene_orders <- c(allgene_orders, rownames(infercnv_obj@gene_order))
}

allgenes <- as.data.frame(table(allgene_orders))
allgenes_filt <- allgenes %>% filter(Freq == max(allgenes$Freq))

#This contains the common genes across all samples in chromosomal order (includes coordinates)
gene_order <- infercnv_obj@gene_order[rownames(infercnv_obj@gene_order) %in% allgenes_filt$allgene_orders, ]
gene_order$chr <- as.character(gene_order$chr)
gene_order$chr_num_char <- substr(gene_order$chr, 4, nchar(gene_order$chr))
gene_order$chr_num <- as.numeric(substr(gene_order$chr, 4, nchar(gene_order$chr)))

#Remove genes from HLA region of Chromosome 6
chr6_gene_ind <- gene_order %>% filter(chr == "chr6" & start > 26055967 & stop < 56951643)
gene_order <- gene_order[!(rownames(gene_order) %in% rownames(chr6_gene_ind)), ]

#Expression data (from inferCNV) 

refcell_expr <- NULL
obscell_expr <- NULL

#sample.names <- c("PD42184", sample.names)
for(name in sample.names){
  infercnv_obj <- readRDS(paste0(infercnv_datadir, "infercnv_gosh_", name, "_0.3/run.final.infercnv_obj"))  
  
  ref_ind <- unlist(infercnv_obj@reference_grouped_cell_indices)
  obs_ind <- unlist(infercnv_obj@observation_grouped_cell_indices)
  
  exprmat <- infercnv_obj@expr.data 
  
  refcell <- exprmat[rownames(exprmat) %in% rownames(gene_order), ref_ind]
  obscell <- exprmat[rownames(exprmat) %in% rownames(gene_order), obs_ind]
  
  #refcell_expr <- cbind(refcell_expr, refcell)
  obscell_expr <- cbind(obscell_expr, obscell)
}

exprmat_filt <- cbind(refcell, obscell_expr)

#Obtain cell information for those cells expressing the program and random cells (sampled from the reference)
allAct <- boxplotDf[, "Activity"]
names(allAct) <- rownames(boxplotDf)

set.seed(1)
randomCells <- sample(colnames(refcell), 50)
sorted_act <- sort(allAct, decreasing = T)[1:50]

#Use an activity threshold of 50 to filter out cells that don't express the program 
nonexp <- NULL
for(i in 1:length(sorted_act)){
  if(sorted_act[i] <= 50){
    nonexp <- c(nonexp, i)
  }
}

if(is.null(nonexp)){
  topProgCells <- names(sorted_act)
}else{
  topProgCells <- sorted_act[-nonexp]
  topProgCells <- names(topProgCells)
}


topProgCell_expr <- t(exprmat_filt[, topProgCells])
randomCell_expr <- t(exprmat_filt[, randomCells])

#Obtain average expression for each gene across test and reference cells
topProgCell_expr_vals <- apply(topProgCell_expr, 2, mean)
randomCell_expr_vals <- apply(randomCell_expr, 2, mean)

#Dataframe for plotting 
expr_vals <- as.data.frame(cbind(topProgCell_expr_vals, randomCell_expr_vals))

expr_vals$topProgcols <- "darkmagenta"
expr_vals$toprandomcols <- "goldenrod1"

chrs <- c(gene_order$chr_num_char, gene_order$chr_num_char)
starts <- c(gene_order$start, (gene_order$start + 1))
expr_vals_filt <- c(expr_vals$topProgCell_expr_vals, expr_vals$randomCell_expr_vals)
snpcols <- c(expr_vals$topProgcols, expr_vals$toprandomcols)
snpopacities <- c(rep(1, nrow(expr_vals)), rep(0.15, nrow(expr_vals)))

#this will add two additional points that will not be shown on the plot
#sets the range - makes consistent across all plots

chrs <- c("1", "1", chrs)
starts <- c(1, 2, starts)
expr_vals_filt <- c(1.2, 0.8, expr_vals_filt)
spcols <- c("#F7F7F7", "#F7F7F7", snpcols)
snpopacities <- c(0, 0, snpopacities)


tracklist <- tracklist + BioCircosSNPTrack("expr", chrs, starts, expr_vals_filt, 
                                           maxRadius = 0.95, minRadius = 0.65,
                                           color = snpcols, opacities = snpopacities)
tracklist <- tracklist + BioCircosBackgroundTrack("myBackgroundTrack_testcell", 
                                                  maxRadius = 0.95, minRadius = 0.65, 
                                                  borderColors = "#E8FDE8", borderSize = 0.3, fillColors = "#F7F7F7")  


widgetBioCircos <- BioCircos(tracklist, yChr = FALSE,
                             genomeFillColor = magma(23, direction = -1),
                             chrPad = 0.02,
                             displayGenomeBorder = F, genomeTicksDisplay = F, genomeLabelDisplay = F, zoom = T, 
                             ARCMouseOverDisplay = T, ARCMouseOverTooltipsHtml04 = "<br/>Patient ID: ")

#+ echo=FALSE

fnameBioCircos <- paste0("BioCircos_plot_", thisProgNumber, ".html") 
htmlwidgets::saveWidget(widgetBioCircos, fnameBioCircos)

# #' Interactive BioCircos plot:
#+ BioCircosPlotiFrame, results="asis", echo=FALSE
cat(paste("\n<iframe scrolling='no' seamless='seamless' style='border:none' src='", fnameBioCircos, "' width='800' height='700'></iframe>\n", sep=""))

cat("\n<br><hr>")


# mSigDB enrichments for cell type, biological processes, scored for the top 50 genes in the program.
#+ echo=FALSE
getGeneSetRes <- function(geneSetDb)
{
  totalNumGenes <- length(commsSorted)
  theSums <- numeric()
  pVals <- numeric()
  theProp <- numeric() # the proportion of genes in the top 50 from the gene expression program featuring in a gene set.
  theOddsRatio <- numeric()
  oddsRatioCILowerBound <- numeric()
  genesFoundAsText <- character() # store the character vector of the genes found in each gene set.
  
  for(i in 1:length(geneSetDb$genesets))
  {
    genesFound <- intersect(names(sort(commsSorted, decreasing=T)[1:50]), geneSetDb$genesets[[i]])
    genesFoundAsText[i] <- paste(genesFound, collapse=", ")
    theSums[i] <- sum(names(sort(commsSorted, decreasing=T)[1:50]) %in% geneSetDb$genesets[[i]])
    
    # calculate Odds ratio
    inGset_inTop50 <- theSums[i]
    inGset_notInTop50 <- length(geneSetDb$genesets[[i]]) - inGset_inTop50
    notInGset_inTop50 <- 50 - inGset_inTop50
    notInGset_notInTop50 <- totalNumGenes - notInGset_inTop50
    
    fisherOut <- fisher.test(matrix(c(inGset_inTop50, inGset_notInTop50, notInGset_inTop50, notInGset_notInTop50), nrow=2))
    theOddsRatio[i] <- fisherOut$estimate
    pVals[i] <- fisherOut$p.value
    oddsRatioCILowerBound[i] <- fisherOut$conf.int[1]
  }
  
  # add names and format as a data frame.
  fwerList <- p.adjust(pVals, method="bonf")
  fdrList <- p.adjust(pVals, method="BH")
  names(theSums) <- geneSetDb$geneset.names
  names(pVals) <- geneSetDb$geneset.names
  names(theOddsRatio) <- geneSetDb$geneset.names
  names(oddsRatioCILowerBound) <- geneSetDb$geneset.names
  names(fwerList) <- geneSetDb$geneset.names
  names(fdrList) <- geneSetDb$geneset.names
  totalGeneSetSizes <- sapply(geneSetDb$genesets, length)
  names(genesFoundAsText) <- geneSetDb$geneset.names
  
  # Format the data frame that contains the results and info.
  resOut <- data.frame(pVals, theOddsRatio, oddsRatioCILowerBound, fdrList, fwerList, theSums, totalGeneSetSizes, genesFoundAsText)
  resOrd <- resOut[order(resOut[,3], decreasing=T), ] # rank by the lower bound of 95% confidence interval of OR
  
  colnames(resOrd) <- c("P-value", "OR", "Lower 95% CI", "FDR", "FWER", "Genes Found", "Gene Set Size", "Specific genes identified")
  
  return(resOrd)
}

#+ cssText, results="asis", echo=FALSE
cssText <- "<style>
.tooltip {
  position: relative;
  display: inline-block;
  border-bottom: 1px dotted black; 
}

.tooltip .tooltiptext {
  visibility: hidden;
  width: 120px;
  background-color: black;
  color: #fff;
  text-align: center;
  padding: 5px 0;
  border-radius: 6px;
 
  position: absolute;
  z-index: 1;
}

.tooltip:hover .tooltiptext {
  visibility: visible;
}
</style>"
cat(cssText)





#' **Gene set Enrichments for this program, caculated from top 50 genes**
#+ tablesBlock, results="asis", echo=FALSE
for(i in 1:length(geneSetDatabases))
{
  cat(paste("<br><strong>", names(geneSetDatabases)[i], ":</strong>", sep=""))
  thisRes <- getGeneSetRes(geneSetDatabases[[i]])
  
  # write out the CSV file
  fileName <- paste( "Progam", thisProgNumber, "_", thisProgName, "_", gsub(' ', '', names(geneSetDatabases)[i]), ".csv" , sep="") # should be unique to this program, and this gene set.
  write.csv(thisRes, file=paste0(getwd(), "/", fileName))
  
  # Code below to format this table for display in HTML, including links, hover etc.
  thisRes[,c(1, 4, 5)] <- formatC(data.matrix(thisRes[,c(1, 4, 5)]), format = "e", digits = 2)
  thisRes[,c(2)] <- round(as.numeric(thisRes[,c(2)]), digits = 2)
  thisRes[,c(3)] <- round(as.numeric(thisRes[,c(3)]), digits = 2)
  
  # Make the rownames be linkable.
  namesToDescription <- geneSetDatabases[[i]]$geneset.descriptions
  names(namesToDescription) <- geneSetDatabases[[i]]$geneset.names
  linkableNames <- paste("<a href='", namesToDescription[rownames(thisRes)], "'>", rownames(thisRes), "</a>", sep="")
  rownames(thisRes) <- linkableNames
  
  thisRes[, 6] <- paste("<div class='tooltip'>", thisRes[,6], "<span class='tooltiptext'>", thisRes[,8], "</span></div>", sep="") # uses the CSS code defined above to show the genes upon hover.
  
  print(kable(thisRes[1:20, 1:7], row.names=TRUE, table.attr = "style='width:150%;'"))
  cat("\n")
  
  # put the link to the CSV file at the bottom
  cat("<a href='", fileName, "'>Dowload full table</a><br>")
}


#' **Top Ranked Transcription Factors for this Gene Expression Program:**
#+ results="asis", echo=FALSE
ranksVector <- rank(-sort(commsSorted, decreasing=T))
measuredTfs <- names(ranksVector)[names(ranksVector) %in% rownames(humanTfs)]
theTopTfs <- ranksVector[measuredTfs][1:20]
topEnsTfs <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theTopTfs), 3] ## /addedCode
names(topEnsTfs) <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theTopTfs), 1] ## /addedCode
decartesTfLink <- paste("<a href='https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/gene/", topEnsTfs[names(theTopTfs)], "/in/adrenal'>Descartes</a>", sep="") ## /addedCode
gtexLinks <- paste("<a href='https://gtexportal.org/home/gene/", names(theTopTfs), "'>GTEx</a>", sep="")
depmapLinks <- paste("<a href='https://depmap.org/portal/gene/", names(theTopTfs), "?tab=dependency'>DepMap</a>", sep="")
theTfTableOut <- data.frame("TFRank"=theTopTfs, humanTfs[names(theTopTfs), ], GtexLink=gtexLinks, depmapLink=depmapLinks, decartesTfLink=decartesTfLink) ## /addedCode
theTfTableOut[,2] <- paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", names(theTopTfs), "'>", names(theTopTfs), "</a>", sep="")
colnames(theTfTableOut) <- c("TF Rank", "Gene Symbol", "ID", "DNA Binding Domain", "Motif Status", "IUPAC PWM", "GTEx", "DepMap", "Decartes") ## /addedCode
print(kable(theTfTableOut[, c(2, 1, 4, 5, 6, 7, 8, 9)], row.names=F, table.attr = "style='width:150%;'")) ## /addedCode

cat("\n<br><hr>")

#' **QQ Plot showing correlations with other GEPs in this dataset, calculated by Spearman correlation:**
#+ echo=FALSE, warnings=FALSE, message=FALSE, include=FALSE
progNumInSplit0 <- communitiesOther[[thisProgNumber]]$Split0 # get the progam number...
thisProgH <- split0_H[,progNumInSplit0[1]] # just go with the first split of the data for now (this should really be averaged over the splits).
theHCors <- numeric()
theHPsCors <- numeric()
for(i in 1:length(communitiesOther))
{
  progNumInSplit0 <- communitiesOther[[i]]$Split0 # get the progam number...
  otherProgH <- split0_H[,progNumInSplit0[1]]
  corOut <- cor.test(thisProgH, otherProgH, method="spearman")
  theHCors[i] <- corOut$estimate
  theHPsCors[i] <- corOut$p.value
  print(i)
}
names(theHCors) <- programNames
names(theHPsCors) <- programNames
theHFdrCors <- p.adjust(theHPsCors, method="BH")

# Make interactive QQplot...
qqData <- sort(theHCors[-thisProgNumber])
qqPs <- theHPsCors[names(qqData)]
qqFdrs <- theHFdrCors[names(qqData)]
i <- 1:length(qqData)
fi <- (i - 0.5) / length(qqData)
xNorm <- qnorm(fi)

lineCoOrds <- qqlineCords(qqData)

qqDf <- data.frame(qqData, xNorm, theNames=qqData, theText=paste("Program ", names(qqData), "\nSpearman Correlation: ", round(qqData, digits=2), "\nP-value: ", formatC(qqPs, format = "e", digits = 2), "\nFDR: ", formatC(qqFdrs, format = "e", digits = 2), sep=""))
p <- ggplot(qqDf, aes(xNorm, qqData, text=theText)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() +
  xlab("Theoretical Normal Distribution (Z)") + 
  ylab("Observed Spearman Correlation") +
  geom_abline(intercept = lineCoOrds[1], slope=lineCoOrds[2])

#+ echo=FALSE
widgetQQ <- plotly::ggplotly(p, tooltip="text")
fnameQQcov <- paste("covariationQQPlot", thisProgNumber, ".html", sep="") # This needs to be UPDATED so this is filename is program specific.
htmlwidgets::saveWidget(widgetQQ, fnameQQcov)

#' Interactive QQ-plot of gene loadings:
#+ correlationPlotiFrame, results="asis", echo=FALSE
cat(paste("\n<iframe scrolling='no' seamless='seamless' style='border:none' src='", fnameQQcov, "' width='1200' height='750'></iframe>\n", sep=""))

cat("\n<br><hr>")

#' **A similar QQ-plot as above, but only for instances where the H value is e.g. > 25, i.e. we are confident that the expression program is active above noise. Agreemenet between these binary vectors is tested using the Jaccard Index, with the P-values calculated by an exact test:**
#+ echo=FALSE, warnings=FALSE, message=FALSE, include=FALSE 
progNumInSplit0 <- communitiesOther[[thisProgNumber]]$Split0 # get the progam number...
thisProgH <- split0_H[,progNumInSplit0[1]] # just go with the first split of the data for now (this should really be averaged over the splits).
thisProgH_over25 <- as.numeric(thisProgH > 25)
jaccardInd25 <- numeric()
jaccardInd25Ps <- numeric()
for(i in 1:length(communitiesOther))
{
  progNumInSplit0 <- communitiesOther[[i]]$Split0 # get the progam number...
  otherProgH_over25 <- as.numeric(split0_H[,progNumInSplit0[1]] > 25)
  
  jaccardTest25 <- jaccard.test.mca(thisProgH_over25, otherProgH_over25)
  jaccardInd25[i] <- jaccardTest25$statistics
  jaccardInd25Ps[i] <- jaccardTest25$pvalue
  
  print(i)    
  
}
names(jaccardInd25) <- programNames
names(jaccardInd25Ps) <- programNames
jaccardInd25Fdrs <- p.adjust(jaccardInd25Ps, method="BH")


# Make interactive QQplot...
qqData <- sort(jaccardInd25[-thisProgNumber])
qqPs <- jaccardInd25Ps[names(qqData)]
qqFdrs <- jaccardInd25Fdrs[names(qqData)]
i <- 1:length(qqData)
fi <- (i - 0.5) / length(qqData)
xNorm <- qnorm(fi)

lineCoOrds <- qqlineCords(qqData)

qqDf <- data.frame(qqData, xNorm, theNames=qqData, theText=paste("Program ", names(qqData), "\nJaccard Index: ", round(qqData, digits=2), "\nP-value: ", formatC(qqPs, format = "e", digits = 2), "\nFDR: ", formatC(qqFdrs, format = "e", digits = 2), sep=""))
p <- ggplot(qqDf, aes(xNorm, qqData, text=theText)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() +
  xlab("Theoretical Normal Distribution (Z)") + 
  ylab("Observed Jaccard Index") +
  geom_abline(intercept = lineCoOrds[1], slope=lineCoOrds[2])

#+ echo=FALSE
widgetQQ <- plotly::ggplotly(p, tooltip="text")
fnameQQcov <- paste("jaccardAgreementActivityOver25QQPlot", thisProgNumber, ".html", sep="") # This needs to be UPDATED so this is filename is program specific.
htmlwidgets::saveWidget(widgetQQ, fnameQQcov)

#' Interactive QQ-plot:
#+ jaccardAgreementPlotiFrame, results="asis", echo=FALSE
cat(paste("\n<iframe scrolling='no' seamless='seamless' style='border:none' src='", fnameQQcov, "' width='1200' height='750'></iframe>\n", sep=""))


cat("\n<br><hr>")


#' **Singler cell type annotations for the top 50 cells on this program.**
#+ SinglerTableBlock, results="asis", echo=FALSE
# First make a string to display the top cell types scoring on this program, we'll show this on hover.

topCellsAnnotMat <- singlerCellLabels[topProgCells,] # pull out the top 50 cells
theMaxesString <- character()

singlerval <- character()
for(i in 1:nrow(topCellsAnnotMat))
{
  scores <- topCellsAnnotMat[i, ]
  names(scores) <- colnames(topCellsAnnotMat)
  topscores <- sort(scores, decreasing = T)[1: 10]
  topnames <- names(sort(scores, decreasing = T)[1:10])
  singlerval <- c(singlerval, topnames[1])
  
  thisstring <- paste(topnames[1], ": ", round(topscores[1], digits=2), sep="")
  for(j in 2:length(topscores))
  {
    thisstring <- paste(thisstring, ", ", topnames[j], ": ", round(topscores[j], digits=2), sep="")
  }
  theMaxesString[i] <- thisstring
}


# Create the table to show the singler results
singlerTooltip <- paste("<div class='tooltip'>Raw Scores<span class='tooltiptext'>", theMaxesString, "</span></div>", sep="") # uses the CSS code defined above to show the genes upon hover.
colnames(singlerDeltaValues) <- c("Singler label", "Singler Delta")
singlerDeltaValues[,2] <- round(singlerDeltaValues[,2], digits=2)
singlerTable <- data.frame("Cell ID"=topProgCells, singlerDeltaValues[topProgCells, ], "Activity Score"=round(allAct[topProgCells], digits=2), "Top Singler Scores"=singlerTooltip)
colnames(singlerTable) <- c("Cell ID", "Singler label", "Singler Delta", "Activity Score", "Top Singler Raw Scores")
print(kable(singlerTable, row.names=FALSE, table.attr = "style='width:150%;'"))

cat("\n<br><hr>")

source("literatureCuratedGenes.R") # The file where these genes are in a structured format, in 2 objects literatureMarkers_smallGeneSets and literatureMarkers_largeGeneSets

## /addedBlock
# Figure out the ordering of these literature curated gene sets (order by mean rank.).
meanRanks <- numeric()
wilcoxPVals_smallSets <- numeric()
for(i in 1:length(literatureMarkers_smallGeneSets))
{
  theseGenes <- commsSorted[names(commsSorted) %in% literatureMarkers_smallGeneSets[[i]][-1]]
  
  if(length(theseGenes > 0)) # are any of these genes in the dataset.
  {
    geneRanks <- 1:length(commsSorted)
    names(geneRanks) <- names(rev(commsSorted))
    theseGenesRanks <- geneRanks[names(theseGenes)]
    
    meanRanks[i] <- mean(theseGenesRanks) # this will be used to rank these gene sets.
    notTheseGenes <- names(commsSorted)[!names(commsSorted) %in% names(theseGenes)]
    wilcoxPVals_smallSets[i] <- wilcox.test(commsSorted[notTheseGenes], commsSorted[names(theseGenes)], alternative="less")$p.value        
  }
}
## end /addedBlock

#' **Below shows the significant enrichments of this GEP for literature curated gene lists<br>**
#' &emsp;&emsp;This data was procured from existing single cell RNA-seq maps of neuroblastoma or related relevant data.<br> 
#' &emsp;&emsp;High ranks indicate this gene is a driver of this GEP.<br> 
#' &emsp;&emsp;These curated gene list are ranked by the P-value (on this GEP) of their constituent genes.<br>
#' &emsp;&emsp;The Mean Count column shows the mean read count in cells scoring highly (H > 50) on this gene expression program.<br>
#+ LiteratureGenesTablesBlock, results="asis", echo=FALSE

literatureMarkers_smallGeneSets_ord <- literatureMarkers_smallGeneSets[order(wilcoxPVals_smallSets)] ## /addedCode, I've also changed literatureMarkers_smallGeneSets to literatureMarkers_smallGeneSets_ord in the rest of this block
wilcoxPs_ordSmall <- wilcoxPVals_smallSets[order(wilcoxPVals_smallSets)]
meanRanks_ordSmall <- meanRanks[order(wilcoxPVals_smallSets)]


#+ jsText, results="asis", echo=FALSE
jsText <- "<script>
  function ShowAndHide() {
    var x = document.getElementById(\'smallset\');
    if (x.style.display == \'none\') {
      x.style.display = \'block\';
    } else {
      x.style.display = \'none\';
    }
  }
</script>"
cat(jsText)


for(i in 1:length(literatureMarkers_smallGeneSets_ord))  
{
  if(i == 4)
  {
    cat("<br><button onclick=\"ShowAndHide()\">Click to show/hide remaining gene sets</button><br><div id=\"smallset\" style=\"display:none\">")
  }
    
  cat(paste("<br><strong>", names(literatureMarkers_smallGeneSets_ord)[i], "</strong>", sep="")) # Print the name of this small gene set
  
  cat(paste("<br><em>", literatureMarkers_smallGeneSets_ord[[i]][1], ":</em>", sep="")) # Print the descruotuib of this small gene set
  
  theseGenes <- commsSorted[names(commsSorted) %in% literatureMarkers_smallGeneSets_ord[[i]][-1]]
  
  if(length(theseGenes > 0)) # are any of these genes in the dataset.
  {
    
    cat(paste("<br>Wilcoxon ranksum test P-value for gene set overrepresentation: ", formatC(wilcoxPs_ordSmall[i], format = "e", digits = 2), "<br>Mean rank of genes in gene set: ", round(meanRanks_ordSmall[i], digits=2), "<br>Rank on gene expression program of genes in gene set:", sep=""))
    
    
    geneRanks <- 1:length(commsSorted)
    names(geneRanks) <- names(rev(commsSorted))
    theseGenesRanks <- geneRanks[names(theseGenes)]
    
    # Include gene cards, GTEx and DepMap links...
    geneCardsLinkSmallGsets <- paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", names(theseGenes), "'>", names(theseGenes), "</a>", sep="")
    gtexLinksSmallGsets <- paste("<a href='https://gtexportal.org/home/gene/", names(theseGenes), "'>GTEx</a>", sep="")
    depmapLinksSmallGsets <- paste("<a href='https://depmap.org/portal/gene/", names(theseGenes), "?tab=dependency'>DepMap</a>", sep="")
    
    topEnsGenes <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theseGenes), 3] ## /addedCode
    names(topEnsGenes) <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theseGenes), 1] ## /addedCode
    topEnsGenes_ord <- topEnsGenes[names(theseGenes)] ## /addedCode
    decartesGenesLink <- paste("<a href='https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/gene/", topEnsGenes_ord[names(topEnsGenes_ord)], "/in/adrenal'>Descartes</a>", sep="") ## /addedCode
    
    outTab <- data.frame(Genes=geneCardsLinkSmallGsets, Weight=theseGenes, Rank=theseGenesRanks, GTEx=gtexLinksSmallGsets, DepMap=depmapLinksSmallGsets, Decartes=decartesGenesLink, "Mean Counts"=countHighScoringCells[names(theseGenes)], "Mean TPM"=meanTpmHOver50[names(theseGenes)])
    
    print(kable(outTab[nrow(outTab):1, ], row.names=FALSE, table.attr = "style='width:150%;'"))
    cat("\n")
    
    # print genes with no detectable expression
    noDetGenes <- noDetectableExpression[noDetectableExpression %in% literatureMarkers_smallGeneSets_ord[[i]][-1]]
    if(length(noDetGenes) > 0)
    {
      cat(paste("<br><em>No detectable expression in this dataset: </em>", paste(noDetGenes, collapse=" "), "<br>", sep="")) # Print the descruotuib of this small gene set
    }
    cat("\n")
  }
}

cat("</div>")


# A matrix of genes x gene set name would be interesting... Allow to be downloaded???? I.e. you could see gene -> gene sets too, rather than just gene sets -> gene...
# Would need to drop all element 1s before making this table....
allGenes <- unique(unlist(literatureMarkers_smallGeneSets))
markerMatFull <- numeric(length(literatureMarkers_smallGeneSets)*length(allGenes))
dim(markerMatFull) <- c(length(literatureMarkers_smallGeneSets), length(allGenes))

theseGenes <- commsSorted[names(commsSorted) %in% literatureMarkers_largeGeneSets[[1]][-1]]
geneRanks <- 1:length(commsSorted)
names(geneRanks) <- names(rev(commsSorted))
theseGenesRanks <- geneRanks[names(theseGenes)]

cat("\n<br><hr>")

#' **Below shows ranks on this GEP for literature curated gene lists for large gene sets <br>**
#' &emsp;&emsp;These include those reported as mesenchymal/adrenergic by Van Groningen et al.<br> 
#' &emsp;&emsp;High ranks indicate this gene is a driver of this GEP (note these results are not ordered).<br> 
#' &emsp;&emsp;The Mean Count column shows the mean read count in cells scoring highly (H > 50) on this gene expression program.<br>
#+ LiteratureGenesTablesBlockLARGE, results="asis", echo=FALSE
for(i in 1:length(literatureMarkers_largeGeneSets))
{
  cat(paste("<br><strong>", names(literatureMarkers_largeGeneSets)[i], "</strong>", sep="")) # Print the name of this small gene set
  cat(paste("<br><em>", literatureMarkers_largeGeneSets[[i]][1], "</em>", sep="")) # Print the name of this small gene set
  
  theseGenes <- commsSorted[names(commsSorted) %in% literatureMarkers_largeGeneSets[[i]][-1]] # get the gene set genes
  geneRanks <- 1:length(commsSorted)
  names(geneRanks) <- names(rev(commsSorted))
  theseGenesRanks <- geneRanks[names(theseGenes)]
  
  # put a little summary here, Wilcox P-value, average and median ranks...
  notTheseGenes <- names(commsSorted)[!names(commsSorted) %in% names(theseGenes)]
  wilcoxPVal <- wilcox.test(commsSorted[notTheseGenes], commsSorted[names(theseGenes)], alternative="less")$p.value
  
  cat(paste("<br>Wilcoxon ranksum test P-value for gene set overrepresentation: ", formatC(wilcoxPVal, format = "e", digits = 2), "<br>Mean rank of genes in gene set: ", round(mean(theseGenesRanks), digits=2), "<br>Median rank of genes in gene set: ", median(theseGenesRanks), "<br>Rank on gene expression program of top 30 genes in gene set:", sep=""))
  
  # Include gene cards, GTEx and DepMap, decartes links...
  geneCardsLinkSmallGsets <- paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", names(theseGenes), "'>", names(theseGenes), "</a>", sep="")
  gtexLinksSmallGsets <- paste("<a href='https://gtexportal.org/home/gene/", names(theseGenes), "'>GTEx</a>", sep="")
  depmapLinksSmallGsets <- paste("<a href='https://depmap.org/portal/gene/", names(theseGenes), "?tab=dependency'>DepMap</a>", sep="")
  
  topEnsGenes <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theseGenes), 3] ## /addedCode
  names(topEnsGenes) <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theseGenes), 1] ## /addedCode
  topEnsGenes_ord <- topEnsGenes[names(theseGenes)] ## /addedCode
  decartesGenesLink <- paste("<a href='https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/gene/", topEnsGenes_ord[names(topEnsGenes_ord)], "/in/adrenal'>Descartes</a>", sep="") ## /addedCode
  
  outTab <- data.frame(Genes=geneCardsLinkSmallGsets, Weight=theseGenes, Rank=theseGenesRanks, GTEx=gtexLinksSmallGsets, DepMap=depmapLinksSmallGsets, Descartes=decartesGenesLink, "Mean Counts"=countHighScoringCells[names(theseGenes)], "Mean TPM"=meanTpmHOver50[names(theseGenes)])
  
  print(kable(outTab[nrow(outTab):1, ][1:30, ], row.names=FALSE, table.attr = "style='width:150%;'"))
  cat("\n")
}


## /addedBlock (ALL new IMMUNE marker gene data from CellTypist is BELOW)
## IMMUNE GENE MAKERES HERE
cat("\n<br><hr>")

# Figure out the ordering of these literature curated gene sets (order by mean rank.).
meanRanks <- numeric()
wilcoxPVals_smallSets <- numeric()
for(i in 1:length(literatureMarkers_celltypistImmuneGeneSets))
{
  theseGenes <- commsSorted[names(commsSorted) %in% literatureMarkers_celltypistImmuneGeneSets[[i]][-1]]
  
  if(length(theseGenes > 0)) # are any of these genes in the dataset.
  {
    geneRanks <- 1:length(commsSorted)
    names(geneRanks) <- names(rev(commsSorted))
    theseGenesRanks <- geneRanks[names(theseGenes)]
    
    meanRanks[i] <- mean(theseGenesRanks) # this will be used to rank these gene sets.
    notTheseGenes <- names(commsSorted)[!names(commsSorted) %in% names(theseGenes)]
    wilcoxPVals_smallSets[i] <- wilcox.test(commsSorted[notTheseGenes], commsSorted[names(theseGenes)], alternative="less")$p.value        
  }
}

#' **Below shows the ranks on this GEP for immune marker gene lists curated by CellTypist (https://www.celltypist.org/encyclopedia).<br>**
#' &emsp;&emsp;High ranks indicate this gene is a driver of this GEP.<br> 
#' &emsp;&emsp;These curated gene list are ranked by P-value (on this GEP) of their constituent genes.<br>
#' &emsp;&emsp;The Mean Count column shows the mean read count in cells scoring highly (H > 50) on this gene expression program.<br>
#+ CellTypistGenesTablesBlock, results="asis", echo=FALSE
literatureMarkers_celltypistImmuneGeneSets_ord <- literatureMarkers_celltypistImmuneGeneSets[order(wilcoxPVals_smallSets)]
wilcoxPs_ordSmall <- wilcoxPVals_smallSets[order(wilcoxPVals_smallSets)]
meanRanks_ordSmall <- meanRanks[order(wilcoxPVals_smallSets)]


#+ jsText2, results="asis", echo=FALSE
jsText2 <- "<script>
  function ShowAndHide_celltypist() {
    var x = document.getElementById(\'celltypist\');
    if (x.style.display == \'none\') {
      x.style.display = \'block\';
    } else {
      x.style.display = \'none\';
    }
  }
</script>"
cat(jsText2)

for(i in 1:length(literatureMarkers_celltypistImmuneGeneSets_ord))
{
  if(i == 4)
  {
    cat("<br><button onclick=\"ShowAndHide_celltypist()\">Click to show/hide remaining gene sets</button><br><div id=\"celltypist\" style=\"display:none\">")
  }
  
  cat(paste("<br><strong>", names(literatureMarkers_celltypistImmuneGeneSets_ord)[i], "</strong>", sep="")) # Print the name of this small gene set
  
  cat(paste("<br><em>", literatureMarkers_celltypistImmuneGeneSets_ord[[i]][1], ":</em>", sep="")) # Print the descruotuib of this small gene set
  
  theseGenes <- commsSorted[names(commsSorted) %in% literatureMarkers_celltypistImmuneGeneSets_ord[[i]][-1]]
  
  if(length(theseGenes > 0)) # are any of these genes in the dataset.
  {
    
    cat(paste("<br>Wilcoxon ranksum test P-value for gene set overrepresentation: ", formatC(wilcoxPs_ordSmall[i], format = "e", digits = 2), "<br>Mean rank of genes in gene set: ", round(meanRanks_ordSmall[i], digits=2), "<br>Rank on gene expression program of genes in gene set:", sep=""))
    
    
    geneRanks <- 1:length(commsSorted)
    names(geneRanks) <- names(rev(commsSorted))
    theseGenesRanks <- geneRanks[names(theseGenes)]
    
    # Include gene cards, GTEx and DepMap links...
    geneCardsLinkSmallGsets <- paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", names(theseGenes), "'>", names(theseGenes), "</a>", sep="")
    gtexLinksSmallGsets <- paste("<a href='https://gtexportal.org/home/gene/", names(theseGenes), "'>GTEx</a>", sep="")
    depmapLinksSmallGsets <- paste("<a href='https://depmap.org/portal/gene/", names(theseGenes), "?tab=dependency'>DepMap</a>", sep="")
    
    topEnsGenes <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theseGenes), 3] 
    names(topEnsGenes) <- symTabNoDups[symTabNoDups$SYMBOL %in% names(theseGenes), 1] 
    topEnsGenes_ord <- topEnsGenes[names(theseGenes)] 
    decartesGenesLink <- paste("<a href='https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/gene/", topEnsGenes_ord[names(topEnsGenes_ord)], "/in/adrenal'>Descartes</a>", sep="") 
    
    outTab <- data.frame(Genes=geneCardsLinkSmallGsets, Weight=theseGenes, Rank=theseGenesRanks, GTEx=gtexLinksSmallGsets, DepMap=depmapLinksSmallGsets, Decartes=decartesGenesLink, "Mean Counts"=countHighScoringCells[names(theseGenes)], "Mean TPM"=meanTpmHOver50[names(theseGenes)])
    
    print(kable(outTab[nrow(outTab):1, ], row.names=FALSE, table.attr = "style='width:150%;'"))
    cat("\n")
    
    # print genes with no detectable expression
    noDetGenes <- noDetectableExpression[noDetectableExpression %in% literatureMarkers_celltypistImmuneGeneSets_ord[[i]][-1]]
    if(length(noDetGenes) > 0)
    {
      cat(paste("<br><em>No detectable expression in this dataset: </em>", paste(noDetGenes, collapse=" "), "<br>", sep="")) # Print the descruotuib of this small gene set
    }
    cat("\n")
  }
}

cat("</div>")


