library(tidyverse)
library(phyloseq)
library(NetCoMi)

rm(list = ls())

iniDir <- "C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/Red1"
setwd(iniDir)

# Data import
physeq.Celiac <- readRDS("Celiac-net-diego.rds")
physeq.a4f <- readRDS("AI-net-diego.rds")

# Helper function to export networks 
exportNet <- function(net, physeq, props, filename){
  # EDGES
  # Create edge object from the edge list exported by netConstruct()
  edges <- net$edgelist1 %>%
    select(v1, v2, asso) %>%
    rename(Source = v1,
           Target = v2,
           Weight = asso) %>%
    mutate(Type = "Undirected",
           Sign = if_else(Weight < 0, "Negative", "Positive"))
  
  # NODES
  # get taxonomic info:
  taxa <- data.frame(physeq@tax_table) %>%
    rownames_to_column(var = "Label")
  
  # get mean abundances:
  mean_ab <- data.frame(
    "Abundance" = apply(data.frame(physeq@otu_table), 1, mean)) %>%
    rownames_to_column(var = "Label")
  
  # get clusters:
  clusters <- data.frame("Cluster" =  props$clustering$clust1) %>%
    rownames_to_column(var = "Label")
  
  # get hubs:
  hubs <- data.frame(isHub = rep(1, length(props$hubs$hubs1)),
                     "Label" = props$hubs$hubs1)
  
  # join all of them together:
  metadata <- taxa %>%
    left_join(mean_ab, by = "Label") %>%
    left_join(clusters, by = "Label") %>%
    left_join(hubs, by = "Label") %>%
    mutate(isHub = if_else(is.na(isHub), 0, 1))
  
  
  # WRITE CSV FILES:
  write_csv(edges,
            file = paste0(filename, "_edges.csv"),
            quote = "needed")
  write_csv(metadata,
            paste0(filename, "_metadata.csv"),
            quote = "needed")
}

# Celiac
net.Celiac <- netConstruct(data = physeq.Celiac,
                        dataType = "counts",
                        measure = "spieceasi",
                        sparsMethod = "none",
                        nboot = 1000,
                        cores = 4L)
saveRDS(net.Celiac, "net_celiac.rds")

props.Celiac <- netAnalyze(net.Celiac,
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = c("degree", "betweenness"),
                        hubQuant = 0.9,
                        weightDeg = FALSE, normDeg = FALSE)
saveRDS(props.Celiac, "props_celiac.rds")

# exportar red a Cytoscape:
exportNet(net.Celiac, physeq.Celiac, props.Celiac, "celiac")

# AI4Food
net.a4f <- netConstruct(data = physeq.a4f,
                           dataType = "counts",
                           measure = "spieceasi",
                           sparsMethod = "none",
                           nboot = 1000,
                           cores = 4L)
saveRDS(net.a4f, "net_a4f.rds")

props.a4f <- netAnalyze(net.a4f,
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "betweenness"),
                           hubQuant = 0.9,
                           weightDeg = FALSE, normDeg = FALSE)
saveRDS(props.a4f, "props_a4f.rds")

# Cytoscape:
exportNet(net.a4f, physeq.a4f, props.a4f, "a4f")
