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
source("R/pca.R")

#............................................................
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")


#............................................................
# read in genomic data
#...........................................................
panPAmat <- readRDS(file ="results/panaroo_out.rds")$panaroo_pa_mat
panPAmatdiff <- panPAmat[,colSums(panPAmat) != nrow(panPAmat)]
MOBPAmat <- readRDS(file ="results/MOB_out.rds")$MOBPAmatrix
MOBPAmatdiff <- MOBPAmat[,colSums(MOBPAmat) != nrow(MOBPAmat)]
AMRPAmat <- readRDS(file ="results/AMR_out.rds")$AMRPAmat
AMRPAmatdiff <- AMRPAmat[,colSums(AMRPAmat) != nrow(AMRPAmat)]

# dims
dim(panPAmat)
dim(panPAmatdiff)

dim(MOBPAmat)
dim(MOBPAmatdiff)

dim(AMRPAmat)
dim(AMRPAmatdiff)

# image
image(panPAmatdiff)
image(MOBPAmatdiff)
image(AMRPAmatdiff)

#............................................................
# Structure
#...........................................................
panPA_pca <- pca_pa(panPAmatdiff)
MOBPA_pca <- pca_pa(MOBPAmatdiff)
AMRPA_pca <- pca_pa(AMRPAmatdiff)


plot_PCA_3d(panPA_pca, mtdt)
plot_PCA_3d(MOBPA_pca, mtdt)
plot_PCA_3d(AMRPA_pca, mtdt)

p1 <- plot_loci_contributions(panPA_pca, num_components = 1, contmin = 0)
p1
p1 +
  labs(title = "My Plot Title",
       subtitle = "A descriptive subtitle",
       x = "X-axis Label",
       y = "Y-axis Label",
       caption = "Source: Data source") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        legend.position = "right",
        legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))

plot_loci_contributions(MOBPA_pca)
plot_loci_contributions(AMRPA_pca)





#............................................................
# ENTROPY
#...........................................................
panaroo_entropy_data <- readr::read_csv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/alignment_entropy.csv")

Entropy is a measure of how informative this pattern is across genomes.
Mathematically: 𝐻 𝑖 = − 𝑝 𝑖 log ⁡ 2 ( 𝑝 𝑖 ) − ( 1 − 𝑝 𝑖 ) log ⁡ 2 ( 1 − 𝑝 𝑖 ) H i ​ =−p i ​ log 2 ​ (p i ​ )−(1−p i ​ )log 2 ​ (1−p i ​ )
ent <- -p*log2(p) - (1-p)*log2(1-p)

plot(ent, abs(loadings))
