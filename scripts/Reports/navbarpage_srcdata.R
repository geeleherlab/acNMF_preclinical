#This file is needed to create the tables and figures that are included in the NB_meta_shiny app

library(igraph)
library(networkD3)
library(RColorBrewer)
library(viridisLite)
library(stringr)

#####################
#METADATA TABLES

#Human
#####################

gosh.table <- readRDS(file = "www/metadata_tables/gosh.table.RDS")
pmc.table <- readRDS(file = "www/metadata_tables/pmc.table.RDS")
jansky.table <- readRDS(file = "www/metadata_tables/jansky.table.RDS")
dong.table <- readRDS(file = "www/metadata_tables/dong.table.RDS")
kameneva.table <- readRDS(file = "www/metadata_tables/kameneva.table.RDS")

#Preclinical
#####################
pdx.table <- readRDS(file = "www/metadata_tables/pdx.table.RDS")
mouse.table <- readRDS(file = "www/metadata_tables/mouse.table.RDS")
cellline.table <- readRDS(file = "www/metadata_tables/cellline.table.RDS")

set.seed(1234)

#############
#IGRAPH
#probably could be added from previous code
###########3

#meta.graph <- graph_from_data_frame(connectivities, directed = F)
meta.graph <- readRDS(file = 'www/igraph/meta.graph.RDS')
edges <- E(meta.graph)
vertices <- as.character(names(V(meta.graph)))
meta.degree <- degree(meta.graph, v = vertices, mode = "total")

layout <- readRDS(file = 'www/igraph/layout.RDS')

#establishes node colors for both dataset and community levels
#also generates the weblinks directory for hyperlinks
h.col <- brewer.pal(n = 6, name = "Blues")

db.vertex.name <- NULL
db.vertex.col <- NULL
weblinks <- NULL

for(i in 1:length(vertices)){
  splt <- str_split(vertices[i], fixed("_"))
  splt_ind <- as.numeric(splt[[1]][length(splt[[1]])])
  
  if(startsWith(vertices[i], "Kildisiute_GOSH")){
    db.vertex.col[i] <- h.col[2]
    db.vertex.name[i] <- "GOSH"
    wblink.dir <- paste0("GEP_reports/human/GOSH/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "Kildisiute_PMC")){
    db.vertex.col[i] <- h.col[3]
    db.vertex.name[i] <- "PMC"
    wblink.dir <- paste0("GEP_reports/human/PMC/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "GSE137804")){
    db.vertex.col[i] <- h.col[4]
    db.vertex.name[i] <- "Dong"
    wblink.dir <- paste0("GEP_reports/human/GSE137804/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "Jansky")){
    db.vertex.col[i] <- h.col[5]
    db.vertex.name[i] <- "Jansky"
    wblink.dir <- paste0("GEP_reports/human/Jansky/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "Kameneva")){
    db.vertex.col[i] <- h.col[6]
    db.vertex.name[i] <- "Verhoeven/Kameneva"
    wblink.dir <- paste0("GEP_reports/human/Kameneva/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "NB_cell_line")){
    db.vertex.col[i] <- "#CC0000"
    db.vertex.name[i] <- "NB_cellline"
    wblink.dir <- paste0("GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "mouse")){
    db.vertex.col[i] <- "#FDCC0D"
    db.vertex.name[i] <- "Mouse"
    wblink.dir <- paste0("GEP_reports/preclinical/Mouse/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "PDX_Human")){
    db.vertex.col[i] <- "#336600"
    db.vertex.name[i] <- "PDX_Human"
    wblink.dir <- paste0("GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_", splt_ind, ".html")
  }else if(startsWith(vertices[i], "PDX_Mouse")){
    db.vertex.col[i] <- "#00CC99"
    db.vertex.name[i] <- "PDX_Mouse"
    wblink.dir <- paste0("GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_", splt_ind, ".html")
  }
  
  weblinks <- c(weblinks, wblink.dir)
}
hyperlink <- weblinks

#Community detection via cluster walktrap
wt <- cluster_walktrap(meta.graph, weights = E(meta.graph)$weight, steps = 50)
wt.communities <- communities(wt)
community.group <- membership(wt)


#renaming community groups
community.names <- c("Cancer Associated Fibroblast (Other)",
                     "Plasmacytoid Dendritic Cells",
                     "T Cells",
                     "B Cells",
                     "Neuroblastoma: Unknown #1", 
                     "Neuroblastoma: MYCN", 
                     "Neuroblastoma: Adrenergic",
                     "Myeloid",
                     "Endothelial",
                     "Cancer Associated Fibroblast (Myofibroblastic)",
                     "Schwannian Stromal Cells",
                     "Cell Cycle",
                     "Neuroblastoma: ALK",
                     "Erythroblasts",
                     "Plasma Cells",
                     "Adrenal Cortex Cells",
                     "Mast Cells",
                     "Neuroblastoma: Unknown #2",
                     "Hypoxia")

for(i in 1:length(community.group)){
  community.group[i] <- community.names[as.numeric(community.group[i])]
}
V(meta.graph)$community <- community.group

#Community colors
vcol <- adjustcolor(rainbow(length(wt.communities)), alpha.f = 0.7)
community.colors <- vcol[V(meta.graph)$community]
