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
library(pheatmap)
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
# read in gene rosetta
#...........................................................
rosetta <- readr::read_csv("~/Documents/Github/P-REALM/data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_presence_absence.csv")
rosetta <- rosetta[,1:3]

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
### Section 2: PCA Visualization All      ####
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
        axis.text.x = element_text(family = "Helvetica", hjust = 0.8, size = 6, angle = 45),
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
        axis.text.x = element_text(family = "Helvetica", hjust = 0.8, size = 6, angle = 45),
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
        axis.text.x = element_text(family = "Helvetica", hjust = 0.8, size = 6, angle = 45),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))



#++++++++++++++++++++++++++++++++++++++++++
### Section 3: Masekd PCA Visualization      ####
#++++++++++++++++++++++++++++++++++++++++++
# there are PCA loadings with perfect collinearity, this means they are
# "LD"/all being conserved on the same genetic unit
# PC1 is explaining ~53% of variation, this is likely structural not phenotypic
# let's mask and reapply
pc1con <- table( round(panPAshell_pca$pca$contribution[,1], digits = 4) ) # NB losing some float here, going to push it down even more
pc1con_val <- as.numeric(as.character(names(pc1con)))
pc1con_val <- pc1con_val[pc1con > 5]

# extract names of loci that we want to exclude
collinear_loci <- names(
  panPAshell_pca$pca$contribution[,1][
    round(panPAshell_pca$pca$contribution[,1], 4) %in% pc1con_val
  ]
)

#...........
# heatmap blocks
#...........
collinear_loci_mat <- panPAmatshell[,colnames(panPAmatshell) %in% collinear_loci]
hclustfig <- pheatmap(mat = t(collinear_loci_mat),
                      fontsize = 6,
                      legend_breaks = c(0,1),
                      legend_labels = c("Absent", "Present"))

#......................
# Let's Look at what is in our Extracted Genes
#......................
gncl <- left_join(tibble::tibble(Gene = collinear_loci),
                  rosetta, by = "Gene")


#...........
# wordcloud figure
#...........
words <- sapply(gncl$Annotation, stringr::str_split, " |;")
wordstab <- table(unlist(words))
wordsdf <- tibble::tibble(words = names(wordstab),
                          freq = as.numeric(wordstab))

wordcloud <- wordsdf %>%
  dplyr::filter(freq > 2) %>%
  dplyr::filter(!(words %in% c("hypothetical", "protein", "family", "domain-containing", "subunit"))) %>%
  ggplot() +
  ggwordcloud::geom_text_wordcloud(aes(label = words,
                                       size = freq,
                                       color = freq)) +
  theme_minimal() +
  scale_color_viridis_c()



#............................................................
# masked PCA
#...........................................................
# make df
panPAmatshelldf <- tibble::as_tibble(panPAmatshell)
panPAmatshelldf <- dplyr::bind_cols(
  tibble::tibble(S_No = rownames(panPAmatshell)),
  panPAmatshelldf) %>%
  dplyr::left_join(., mtdt, by = "S_No")
# split into resolve and persist
panPA_res_pa <- panPAmatshelldf %>%
  dplyr::filter(Phenotype == "Resolving") %>%
  dplyr::select(-c("Phenotype", "Sample_ID", "nm")) %>%
  tidyr::pivot_longer(., cols = c(-"S_No"), names_to = "Gene", values_to = "pa") %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    mn = mean(pa),
  )

panPA_per_pa <- panPAmatshelldf %>%
  dplyr::filter(Phenotype == "Persistent") %>%
  dplyr::select(-c("Phenotype", "Sample_ID", "nm")) %>%
  tidyr::pivot_longer(., cols = c(-"S_No"), names_to = "Gene", values_to = "pa") %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    mn = mean(pa),
  )

# PAF genes
pafgenes <- dplyr::bind_rows(panPA_per_pa, panPA_res_pa) %>%
  dplyr::filter(mn >= 0.4 & mn <= 0.6) %>%
  dplyr::pull(Gene)

#......................
# NOW MASK
#......................
masked_panPAmatshell <- panPAmatshell[,!colnames(panPAmatshell) %in% pafgenes]

# PCA for masked shell
masked_panPAshell_pca <- pca_pa(masked_panPAmatshell)
# PCA plot for masked shell
masked_panPAshell_pcaPlotObj <- plot_PCA_3d(masked_panPAshell_pca, mtdt)

# loci contributions for masked shell
masked_shellPA_pca_loci_cont_PlotObj <- plot_loci_contributions(panPAshell_pca, num_components = 1, contmin = 0) +
  labs(title = "PanGenome-SHELL-MASKED Loci/Gene Contributions",
       x = "PCA Contribution",
       y = "Loci") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        axis.text.x = element_text(family = "Helvetica", hjust = 0.8, size = 6, angle = 45),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


hclustfig2 <- pheatmap(mat = t(masked_panPAmatshell),
                       fontsize = 6,
                       legend_breaks = c(0,1),
                       legend_labels = c("Absent", "Present"))




#............................................................
# Out
#...........................................................
# long list
out <- list(
  panPA_pca = panPA_pca,
  panPAshell_pca = panPAshell_pca,
  vcf_pca = vcf_pca,
  panPA_pcaPlotObj = panPA_pcaPlotObj,
  panPAshell_pcaPlotObj = panPAshell_pcaPlotObj,
  vcf_pcaPlotObj = vcf_pcaPlotObj,
  panlocicont = panPA_pca_loci_cont_PlotObj,
  shelllocicont = shellPA_pca_loci_cont_PlotObj,
  vcflocicont = vcf_pca_loci_cont_PlotObj,
  hclustfig2 = hclustfig2,
  gncl = gncl,
  wordcloud = wordcloud,
  masked_panPAshell_pca = masked_panPAshell_pca,
  masked_panPAshell_pcaPlotObj = masked_panPAshell_pcaPlotObj,
  masked_shellPA_pca_loci_cont_PlotObj
)
saveRDS(out, "results/linear_assoc.RDS")

new_mats <- list(
  masked_panPAmatshell = masked_panPAmatshell,
  panPAmatshell = panPAmatshell
)
saveRDS(new_mats, "data/derived/panPAmatshells.RDS")




#++++++++++++++++++++++++++++++++++++++++++
### Section ZZZ: BATCH EFFECT       ####
#++++++++++++++++++++++++++++++++++++++++++
batch <- readRDS("data/derived/batch_sample_identifiers.RDS")
plot_PCA_3d_batch <- function(pca_run, mtdt){
  col_palette <- viridis::plasma(n = 3)
  dat <- as.data.frame(pca_run$pca$x[, 1:3])
  dat$S_No <- rownames(dat)
  dat <- dplyr::left_join(x = dat, y = mtdt, by = "S_No")
  plotly::plot_ly(dat, x = ~PC1, y = ~PC2, z = ~PC3,
                  color = ~ instrument,
                  Phenotype = ~ instrument,
                  type = "scatter3d", colors = col_palette,
                  mode = "markers", marker = list(size = 5),
                  alpha = 0.8)
}

# plots out
panPA_pcaPlotObj_batch <- plot_PCA_3d_batch(panPA_pca,  mtdt = batch)
panPAshell_pcaPlotObj_batch <- plot_PCA_3d_batch(panPAshell_pca,  mtdt = batch)
vcf_pcaPlotObj_batch <- plot_PCA_3d_batch(vcf_pca,  mtdt = batch)

batchout <- list(
  panPA_pcaPlotObj_batch = panPA_pcaPlotObj_batch,
  panPAshell_pcaPlotObj_batch = panPAshell_pcaPlotObj_batch,
  vcf_pcaPlotObj_batch = vcf_pcaPlotObj_batch
)

saveRDS(batchout, "results/bacth_qc.RDS")

