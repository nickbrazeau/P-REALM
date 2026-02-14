## .................................................................................
## Purpose: Panaroo organization and analyses
##
## Author: Nick Brazeau
##
## Date: 25 November, 2025
##
## Notes:
## .................................................................................
library(tidyverse)
source("R/themes.R")
#++++++++++++++++++++++++++++++++++++++++++
### Section 0: Import Data     ####
#++++++++++++++++++++++++++++++++++++++++++
#............................................................
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")

#............................................................
# read in panaroo
#...........................................................
#......................
# presence absence data
#......................
panaroo_pa_data <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_presence_absence.Rtab")
# reshape so it is more like a VCF
panaroo_pa_data <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = c(-"Gene"), names_to = "S_No", values_to = "pa") %>%
  tidyr::pivot_wider(., names_from = "Gene", values_from = "pa")
# drop samples not wanted
panaroo_pa_data <- panaroo_pa_data %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))
# Presence Absence matrix
panaroo_pa_mat <- as.matrix(panaroo_pa_data[,2:ncol(panaroo_pa_data)])
rownames(panaroo_pa_mat) <- panaroo_pa_data$S_No

#......................
# csv data
#......................
panaroo_pa_csv_data <- readr::read_csv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_presence_absence.csv")


#++++++++++++++++++++++++++++++++++++++++++
### Section 1: Basic Presence Absence Table    ####
#++++++++++++++++++++++++++++++++++++++++++
paFig <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  ggplot() +
  geom_tile(aes(x = Gene, y = S_No, fill = pa)) +
  map_theme +
  labs(
    title = "Presence-Absence Gene-Sample Matrix",
    x = "Genes",
    y = "Samples") +
  scale_fill_viridis_c() +
  theme(
    plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"
  )

#++++++++++++++++++++++++++++++++++++++++++
### Section 2: Pangenome Statistics        ####
#++++++++++++++++++++++++++++++++++++++++++
# core genome: >=95%; cloud genome <15%; shell genome everything else
coregenome <- colnames(panaroo_pa_mat)[colSums(panaroo_pa_mat)/nrow(panaroo_pa_data) >= 0.95]
cloudgenome <- colnames(panaroo_pa_mat)[colSums(panaroo_pa_mat)/nrow(panaroo_pa_data) < 0.15]
shellgenome <- colnames(panaroo_pa_mat)[!(colnames(panaroo_pa_mat) %in% c(coregenome, cloudgenome))]

#......................
# tidy up coregenome for any phylogenetics
#......................
entropy <- readr::read_csv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/alignment_entropy.csv",
                           col_names = F) %>%
  magrittr::set_colnames(c("loci", "shannon_entropy")) %>%
  dplyr::mutate(
    loci = purrr::map_chr(loci, function(x){ stringr::str_replace(x, ".aln", "") })
  )
hist(entropy$shannon_entropy)
summary(entropy$shannon_entropy)
# make a cut to be conservative
clngns <- entropy$loci[entropy$shannon_entropy < 0.025]
coregenome_clean <- coregenome[ coregenome %in% clngns ]

#++++++++++++++++++++++++++++++++++++++++++
### Section 3: Pangenome Visualization        ####
#++++++++++++++++++++++++++++++++++++++++++
#......................
# Shell Figure
#......................
paFigShell <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::filter(Gene %in% shellgenome) %>%
  ggplot() +
  geom_tile(aes(x = Gene, y = S_No, fill = pa)) +
  map_theme +
  labs(
    title = "Presence-Absence Gene-Sample Matrix \n for Shell Genome",
    x = "Genes",
    y = "Samples") +
  scale_fill_viridis_c() +
  theme(
    plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"
  )

# Shell Figure by phenotype
paFigShellresolve <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::filter(Gene %in% shellgenome) %>%
  dplyr::inner_join(mtdt, ., by = "S_No") %>%
  dplyr::filter(Phenotype == "Resolving") %>%
  ggplot() +
  geom_tile(aes(x = Gene, y = S_No, fill = pa)) +
  map_theme +
  labs(
    title = "Presence-Absence Gene-Sample Matrix \n for Cloud Genome amongst Resolving",
    x = "Genes",
    y = "Samples") +
  scale_fill_viridis_c() +
  theme(
    plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"
  )

paFigShellpersist <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::filter(Gene %in% shellgenome) %>%
  dplyr::inner_join(mtdt, ., by = "S_No") %>%
  dplyr::filter(Phenotype == "Persistent") %>%
  ggplot() +
  geom_tile(aes(x = Gene, y = S_No, fill = pa)) +
  map_theme +
  labs(
    title = "Presence-Absence Gene-Sample Matrix \n for Cloud Genome amongst Persistent",
    x = "Genes",
    y = "Samples") +
  scale_fill_viridis_c() +
  theme(
    plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"
  )


#......................
# hierarchical cluster
#......................
panselldf <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::filter(Gene %in% shellgenome) %>%
  tidyr::pivot_wider(., names_from = "Gene", values_from = "pa")
pansellmat <- as.matrix( panselldf[,2:ncol(panselldf)] )
rownames(pansellmat) <- panselldf$S_No

hclustfig <- pheatmap(mat = t(pansellmat),
                      fontsize = 6,
                      legend_breaks = c(0,1),
                      legend_labels = c("Absent", "Present"))




# combine
paFigShell_facet <- cowplot::plot_grid(paFigShellresolve, paFigShellpersist)



#......................
# Invert-Diff Plot
#......................
pers_pa <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::inner_join(mtdt, ., by = "S_No") %>%
  dplyr::filter(Phenotype == "Persistent") %>%
  dplyr::group_by(Gene, Phenotype) %>%
  dplyr::summarise(
    GeneSum = sum(pa)
  )
res_pa <- panaroo_pa_data %>%
  tidyr::pivot_longer(., cols = -c("S_No"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::inner_join(mtdt, ., by = "S_No") %>%
  dplyr::filter(Phenotype == "Resolving") %>%
  dplyr::group_by(Gene, Phenotype) %>%
  dplyr::summarise(
    GeneSum = sum(pa)
  )

# find genes of interest for smaller figure
gnint <- dplyr::full_join(pers_pa, res_pa, by = "Gene") %>%
  dplyr::mutate(L1diff = abs(GeneSum.x - GeneSum.y)) %>%
  dplyr::filter(L1diff > 10) %>%
  dplyr::pull("Gene")


GeneCountPlot_bypheno <- dplyr::bind_rows(pers_pa, res_pa) %>%
  dplyr::filter(Gene %in% gnint) %>%
  ggplot() +
  geom_bar(aes(x = Gene, y = GeneSum,
               fill = Phenotype, color = NA), stat="identity", color="black", position = position_dodge()) +
  scale_fill_viridis_d("Phenotype") +
  labs(
    title = "Sample Genes Counts by Phenotype in Shell Genome with reasonable Diff",
    x = "Genes",
    y = "Number of Samples with Gene") +
  plot_theme +
  theme(
    plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
    axis.title = element_blank(),
    axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 6, angle = 45),
    legend.position = "none"
  )

#......................
# Entropy considerations
#......................
entropyhigh <- entropy[entropy$shannon_entropy > 0.3, ]
panaroo_pa_data[,colnames(panaroo_pa_data) %in% c("S_No", entropyhigh$loci)] %>%
  dplyr::left_join(., mtdt, by = "S_No") %>%
  tidyr::pivot_longer(., cols = -c("S_No", "Sample_ID", "Phenotype", "nm"),
                      names_to = "Gene", values_to = "pa") %>%
  dplyr::group_by(Gene, Phenotype) %>%
  dplyr::summarise(
    GenePA = sum(pa)
  )


#............................................................
# out
#...........................................................
out <-list(
  panaroo_pa_data = panaroo_pa_data,
  panaroo_pa_mat = panaroo_pa_mat,
  panaroo_gen_annot = panaroo_pa_csv_data[,1:3],
  genomes = list(coregenome = coregenome,
                 cloudgenome = cloudgenome,
                 shellgenome = shellgenome,
                 coregenomelowent = coregenome_clean),
  paFigShell = paFigShell,
  paFigShell_facet = paFigShell_facet,
  GeneCountPlot_bypheno = GeneCountPlot_bypheno,
  hclustfig = hclustfig
)
saveRDS(out, file = "results/panaroo_out.rds")



