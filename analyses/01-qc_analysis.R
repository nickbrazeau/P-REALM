## .................................................................................
## Purpose: QC Analysis
##
## Author: Nick Brazeau
##
## Date: 21 November, 2025
##
## Notes:
## .................................................................................
library(tidyverse)
source("R/themes.R")

#............................................................
# storage for out
#...........................................................
out <- list()

#............................................................
# read in data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv")
qc <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/bactopia-20251118-102712/merged-results/assembly-scan.tsv") %>%
  dplyr::rename(S_No = sample)
# quality control
qc <- dplyr::left_join(mtdt, qc, by = "S_No")

#............................................................
# Overview
#...........................................................
# table
qctbl <- qc %>%
  dplyr::select(-c("S_No", "Sample_ID")) %>%
  DT::datatable(.,
                rownames = F,
                extensions='Buttons',
                options = list(
                  searching = T,
                  pageLength = 10,
                  dom = 'Bfrtip',
                  autoWidth = TRUE,
                  buttons = c('csv')))

# append out
out <- append(out, list(qctbl = qctbl))

#............................................................
# Hardcoded excludes samples with >X contigs
#   this was based on reasonable cutoffs
#...........................................................
excludedsmpl <- qc %>%
  dplyr::select(c("S_No", "Sample_ID", "nm", "Phenotype", "total_contig")) %>%
  dplyr::filter(total_contig >= 200)
# append out
out <- append(out, list(excludedsmpl = excludedsmpl))
# write out
readr::write_csv(excludedsmpl, file = "data/derived/excluded_smpls.csv")

# datatable
excludedsmpl %>%
  DT::datatable(.,
                rownames = F,
                extensions='Buttons',
                options = list(
                  searching = T,
                  pageLength = 10,
                  dom = 'Bfrtip',
                  autoWidth = TRUE,
                  buttons = c('csv')))

# plots
# p1
(p1 <- qc %>%
    dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
    ggplot() +
    geom_point(aes(x = n50_contig_length, y = total_contig,
                   shape = Phenotype, color = Sample_ID),
               size = 2, alpha = 0.5) +
    plot_theme +
    theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, angle = 45, size = 11)) +
    theme(legend.position = "none"))

# p2
(p2 <- qc %>%
    dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
    ggplot() +
    geom_point(aes(x = min_contig_length, y = max_contig_length,
                   shape = Phenotype, color = Sample_ID),
               size = 2, alpha = 0.5) +
    plot_theme +
    theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, angle = 45, size = 11)) +
    theme(legend.position = "none"))

# p3
(p3 <- qc %>%
    dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
    ggplot() +
    geom_point(aes(x = mean_contig_length, y = median_contig_length,
                   shape = Phenotype, color = Sample_ID),
               size = 2, alpha = 0.5) +
    plot_theme +
    theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, angle = 45, size = 11)) +
    theme(legend.position = "none"))


# p4
(p4 <- qc %>%
    dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
    ggplot() +
    geom_point(aes(x = total_contig, y = total_contig_length,
                   shape = Phenotype, color = Sample_ID),
               size = 2, alpha = 0.5) +
    plot_theme +
    theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, angle = 45, size = 11)) +
    theme(legend.position = "none"))


# append out
qcplots <- cowplot::plot_grid(p1,p2,p3,p4, nrow = 2)



# p5
(gccont_plot <- qc %>%
    dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
    dplyr::mutate(gc_content = (contig_percent_g + contig_percent_c)) %>%
    ggplot() +
    geom_point(aes(x = n50_contig_length, y = gc_content,
                   shape = Phenotype, color = Sample_ID),
               size = 2, alpha = 0.5) +
    plot_theme +
    theme(legend.position = "none"))

# append out
out <- append(out, list(qcplots = qcplots))
out <- append(out, list(gccont_plot = gccont_plot))

#......................
# out
#......................
saveRDS(out, file = "results/qc_ret_out.rds")



