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
library(jaccard)
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
panPAmatshell <- readRDS("data/derived/panPAmatshells.RDS")$panPAmatshell
masked_panPAmatshell <- readRDS("data/derived/panPAmatshells.RDS")$masked_panPAmatshell


#......................
# QC Shorten Search for reasonable NxM (samples x features)
#......................
# we are going to prune the network here. We want at least 20-80% variation for a signal
prn <- apply(panPAmatshell, 2, mean) # pop allele freq
prn <- ( prn <= 0.9 & prn >= 0.1 )
panPAmatshell <- panPAmatshell[, prn]

#++++++++++++++++++++++++++++++++++++++++++
### Section 1: Make Distances        ####
#++++++++++++++++++++++++++++++++++++++++++
# make a jaccard distance
# assume rows are samples, N x M format
jacc_dist_wrap <- function(binarymat, type){
  if (type == "N") { # sample level
    # hold distmat
    distmat <- matrix(NA, nrow = nrow(binarymat), ncol = nrow(binarymat))
    colnames(distmat) <- rownames(distmat) <- rownames(binarymat)
    # get upper tri
    uppertri <- combn(nrow(binarymat), 2)
    # iter through
    for (i in 1:ncol(uppertri)) {
      distmat[ uppertri[1,i], uppertri[2,i] ] <- jaccard::jaccard(
        x = binarymat[uppertri[1,i], ],
        y = binarymat[uppertri[2,i], ],
      )
    }

  } else if (type == "M") { # gene level.
    # hold distmat
    distmat <- matrix(NA, nrow = ncol(binarymat), ncol = ncol(binarymat))
    colnames(distmat) <- rownames(distmat) <- colnames(binarymat)
    # get upper tri
    uppertri <- combn(ncol(binarymat), 2)
    # iter through
    for (i in 1:ncol(uppertri)) {
      distmat[ uppertri[1,i], uppertri[2,i] ] <- jaccard::jaccard(
        x = binarymat[, uppertri[1,i]],
        y = binarymat[, uppertri[2,i]],
      )
    } # end for loop
  } else {
    stop("error")
  }

  # tidy up
  diag(distmat) <- 0
  distmat[lower.tri(distmat)] <- t(distmat)[lower.tri(distmat)]
  return(as.dist(distmat))
}

# make matrix for shell
shell_sample_jcc <- jacc_dist_wrap(panPAmatshell, type = "N")
shell_gene_jcc <- jacc_dist_wrap(panPAmatshell, type = "M")

# make matrix for masked
shell_sample_jcc_masked <- jacc_dist_wrap(masked_panPAmatshell, type = "N")
shell_gene_jcc_masked <- jacc_dist_wrap(masked_panPAmatshell, type = "M")

#++++++++++++++++++++++++++++++++++++++++++
### Section 2: Shell Sample Matrix/Network     ####
#++++++++++++++++++++++++++++++++++++++++++
# shell long
shell_sample_jcc_long <- broom::tidy(shell_sample_jcc) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "jaccard")) %>%
  dplyr::mutate(
    smpl1 = as.character(smpl1),
    smpl2 = as.character(smpl2)
  )

# make an network and pull in mtdt
shell_sample_jcc_net <- shell_sample_jcc_long %>%
  tidygraph::as_tbl_graph(., directed = F) %>%
  tidygraph::activate("nodes") %>%
  dplyr::mutate(community = as.factor(tidygraph::group_louvain(weights = jaccard))) %>%
  dplyr::rename(S_No = name) %>%
  dplyr::left_join(., mtdt, by = "S_No")

p2 <- shell_sample_jcc_net %>%
  tidygraph::activate("edges") %>%
  dplyr::filter(jaccard > 0.5) %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(color = jaccard)) +
  ggraph::geom_node_point(aes(color = Phenotype),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  ggraph::scale_edge_color_viridis("Jaccard", values = c(0,1), option = "plasma") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PanGenome-SHELL \n Sample-Level Network") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        legend.position = "bottom")



#++++++++++++++++++++++++++++++++++++++++++
### Section 3: Shell Gene Matrix/Network     ####
#++++++++++++++++++++++++++++++++++++++++++
shell_gene_jcc_long <- broom::tidy(shell_gene_jcc) %>%
  magrittr::set_colnames(c("gene1", "gene2", "jaccard")) %>%
  dplyr::mutate(
    gene1 = as.character(gene1),
    gene2 = as.character(gene2)
  )

# make an network and pull in mtdt
shell_gene_jcc_long_net <- shell_gene_jcc_long %>%
  tidygraph::as_tbl_graph(., directed = F) %>%
  tidygraph::activate("nodes") %>%
  dplyr::mutate(Community = as.factor(tidygraph::group_louvain(weights = jaccard))) %>%
  dplyr::rename(gene = name)

p3 <- shell_gene_jcc_long_net %>%
  tidygraph::activate("edges") %>%
  dplyr::filter(jaccard > 0.8) %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(color = jaccard)) +
  ggraph::geom_node_point(aes(color = Community),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  ggraph::scale_edge_color_viridis("Jaccard", values = c(0,1), option = "plasma") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PanGenome-SHELL \n Gene-Level Network") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        legend.position = "bottom")

shell_net_plot <- cowplot::plot_grid(
  p2, p3, nrow = 1
)

#++++++++++++++++++++++++++++++++++++++++++
### Section 4: Masked Sample Matrix/Network     ####
#++++++++++++++++++++++++++++++++++++++++++
# shell long
shell_sample_jcc_maskedlong <- broom::tidy(shell_sample_jcc_masked) %>%
  magrittr::set_colnames(c("smpl1", "smpl2", "jaccard")) %>%
  dplyr::mutate(
    smpl1 = as.character(smpl1),
    smpl2 = as.character(smpl2)
  )

# make an network and pull in mtdt
shell_sample_jcc_maskednet <- shell_sample_jcc_maskedlong %>%
  tidygraph::as_tbl_graph(., directed = F) %>%
  tidygraph::activate("nodes") %>%
  dplyr::mutate(community = as.factor(tidygraph::group_louvain(weights = jaccard))) %>%
  dplyr::rename(S_No = name) %>%
  dplyr::left_join(., mtdt, by = "S_No")

p4 <- shell_sample_jcc_maskednet %>%
  tidygraph::activate("edges") %>%
  dplyr::filter(jaccard > 0.5) %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(color = jaccard)) +
  ggraph::geom_node_point(aes(color = Phenotype),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  ggraph::scale_edge_color_viridis("Jaccard", values = c(0,1), option = "plasma") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PanGenome-SHELL-MASKED \n Sample-Level Network") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        legend.position = "bottom")



#++++++++++++++++++++++++++++++++++++++++++
### Section 5: Masked Gene Matrix/Network     ####
#++++++++++++++++++++++++++++++++++++++++++
shell_gene_jcc_maskedlong <- broom::tidy(shell_gene_jcc_masked) %>%
  magrittr::set_colnames(c("gene1", "gene2", "jaccard")) %>%
  dplyr::mutate(
    gene1 = as.character(gene1),
    gene2 = as.character(gene2)
  )

# make an network and pull in mtdt
shell_gene_jcc_maskedlong_net <- shell_gene_jcc_maskedlong %>%
  tidygraph::as_tbl_graph(., directed = F) %>%
  tidygraph::activate("nodes") %>%
  dplyr::mutate(Community = as.factor(tidygraph::group_louvain(weights = jaccard))) %>%
  dplyr::rename(gene = name)

p5 <- shell_gene_jcc_maskedlong_net %>%
  tidygraph::activate("edges") %>%
  dplyr::filter(jaccard > 0.8) %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(color = jaccard)) +
  ggraph::geom_node_point(aes(color = Community),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  ggraph::scale_edge_color_viridis("Jaccard", values = c(0,1), option = "plasma") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PanGenome-SHELL-MASKED \n Gene-Level Network") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        legend.position = "bottom")

# bring together plot
shell_masked_net_plot <- cowplot::plot_grid(
  p4, p5, nrow = 1
)


#++++++++++++++++++++++++++++++++++++++++++
### Section 999: Mash Distances and Tree      ####
#++++++++++++++++++++++++++++++++++++++++++
mashdist <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/mashtree-20251118-225618/mashtree/mashtree.tsv")
mashdistNJtree <- ape::read.tree("data/derived/prealm_ret/bactopia-runs/mashtree-20251118-225618/mashtree/mashtree.dnd")
plot(mashdistNJtree)
ggtree



#++++++++++++++++++++++++++++++++++++++++++
### out ####
#++++++++++++++++++++++++++++++++++++++++++
shell_sample_jcc_net_comm <- shell_sample_jcc_net %>%
  tidygraph::activate("nodes") %>%
  tibble::as_tibble()
out <- list(
  comm_shell_sample_jcc_net = shell_sample_jcc_net_comm,
  shell_net_plot = shell_net_plot,
  shell_masked_net_plot = shell_masked_net_plot
)

saveRDS(out, file = "results/network_out.RDS")

#++++++++++++++++++++++++++++++++++++++++++
### parking lot  ####
#++++++++++++++++++++++++++++++++++++++++++
# install.packages("WGCNA")
# BiocManager::install("WGCNA")
# library(WGCNA)
# https://link.springer.com/article/10.1186/1471-2105-9-559
# https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0

# Bipartite Network  --> igraph
