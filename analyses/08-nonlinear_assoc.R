## .................................................................................
## Purpose: Network and Nonlinear module assessment
##
## Author: Nick Brazeau
##
## Date: 04 December, 2025
##
## Notes:
## .................................................................................
library(tidyverse)
library(plotly)
library(igraph)
library(tidygraph)
library(WGCNA)
library(ape)
library(ggtree)
#++++++++++++++++++++++++++++++++++++++++++
### Section 0: Import Data        ####
#++++++++++++++++++++++++++++++++++++++++++
#......................
# metadata
#......................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")
#......................
# genomic data
#......................
# panaroo
panPA <- readRDS(file ="results/panaroo_out.rds")
panPAmat <- panPA$panaroo_pa_mat
panPAmatshell <- panPA$panaroo_pa_mat[, colnames(panPA$panaroo_pa_mat) %in% panPA$genomes$shellgenome]


#++++++++++++++++++++++++++++++++++++++++++
### Section 1: Make Distances        ####
#++++++++++++++++++++++++++++++++++++++++++
# make a jaccard distnace

#
  tidygraph::as_tbl_graph(., directed = F)
#++++++++++++++++++++++++++++++++++++++++++
### Section 2: Gene Matrix/Network     ####
#++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++
### Section 999: Mash Distances and Tree      ####
#++++++++++++++++++++++++++++++++++++++++++
mashdist <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/mashtree-20251118-225618/mashtree/mashtree.tsv")
mashdistNJtree <- ape::read.tree("data/derived/prealm_ret/bactopia-runs/mashtree-20251118-225618/mashtree/mashtree.dnd")
plot(mashdistNJtree)
ggtree


#++++++++++++++++++++++++++++++++++++++++++
### Section 3: ***    ####
#++++++++++++++++++++++++++++++++++++++++++
# install.packages("WGCNA")
# BiocManager::install("WGCNA")
# library(WGCNA)
# https://link.springer.com/article/10.1186/1471-2105-9-559
# https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

#++++++++++++++++++++++++++++++++++++++++++
### Section 4: Bipartite Network       ####
#++++++++++++++++++++++++++++++++++++++++++
igraph bipartite


