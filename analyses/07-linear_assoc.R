## .................................................................................
## Purpose: Identify signatures of structure in data from presence-absence matrices
##
## Author: Nick Brazeau
##
## Date: 26 November, 2025
##
## Notes:
## .................................................................................
library(tidyverse)
library(plotly)
library(vcfR)
source("R/pca.R")

#++++++++++++++++++++++++++++++++++++++++++
### Section 0: Setup      ####
#++++++++++++++++++++++++++++++++++++++++++
#............................................................
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")


#............................................................
# read in genomic data
#...........................................................
# panaroo
panPA <- readRDS(file ="results/panaroo_out.rds")
panPAmat <- panPA$panaroo_pa_mat
panPAmatshell <- panPA$panaroo_pa_mat[, colnames(panPA$panaroo_pa_mat) %in% panPA$genomes$shellgenome]
vcf <- readRDS(file = "results/vcf_out.RDS")
vcfmat <- t(vcf[,3:ncol(vcf)])
rownames(vcfmat) <- colnames(vcf)[3:ncol(vcf)]
colnames(vcfmat) <- paste0(vcf$CHROM, vcf$POS)

#++++++++++++++++++++++++++++++++++++++++++
### Section 1: PCA       ####
#++++++++++++++++++++++++++++++++++++++++++
#............................................................
# Structure
# NB, PCA expects samples in rows and columns as loci given SVD
#...........................................................
panPA_pca <- pca_pa(panPAmat)
panPAshell_pca <- pca_pa(panPAmatshell)
vcf_pca <- pca_pa(vcfmat)

#++++++++++++++++++++++++++++++++++++++++++
### Section 2: PCA Visualization        ####
#++++++++++++++++++++++++++++++++++++++++++
panPA_pcaPlotObj <- plot_PCA_3d(panPA_pca, mtdt)
panPAshell_pcaPlotObj <- plot_PCA_3d(panPAshell_pca, mtdt)
vcf_pcaPlotObj <- plot_PCA_3d(vcf_pca, mtdt)
#......................
# loci contributions
#......................
# pangemome
panPA_pca_loci_cont_PlotObj <- plot_loci_contributions(panPA_pca, num_components = 1, contmin = 0) +
  labs(title = "PanGenome Loci/Gene Contributions",
       x = "PCA Contribution",
       y = "Loci") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        axis.text.x = element_text(family = "Helvetica", hjust = 0.5, size = 6, angle = 45),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))
# shell
shellPA_pca_loci_cont_PlotObj <- plot_loci_contributions(panPAshell_pca, num_components = 1, contmin = 0) +
  labs(title = "PanGenome-SHELL Loci/Gene Contributions",
       x = "PCA Contribution",
       y = "Loci") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        axis.text.x = element_text(family = "Helvetica", hjust = 0.5, size = 6, angle = 45),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


vcf_pca_loci_cont_PlotObj <- plot_loci_contributions(vcf_pca, num_components = 1, contmin = 0) +
  labs(title = "SNP Contributions",
       x = "PCA Contribution",
       y = "SNP") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        axis.text.x = element_text(family = "Helvetica", hjust = 0.5, size = 6, angle = 45),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


#............................................................
# Out
#...........................................................

out <- list(
  panPA_pcaPlotObj = panPA_pcaPlotObj,
  panPAshell_pcaPlotObj = panPAshell_pcaPlotObj,
  vcf_pcaPlotObj = vcf_pcaPlotObj,
  panlocicont = panPA_pca_loci_cont_PlotObj,
  shelllocicont = shellPA_pca_loci_cont_PlotObj,
  vcflocicont = vcf_pca_loci_cont_PlotObj
)
saveRDS(out, "results/linear_assoc.RDS")
