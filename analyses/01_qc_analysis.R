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

#............................................................
# Hardcoded excludes samples with >150 contigs
#   this was based on reasonable cutoffs
#...........................................................
excludedsmpl <- qc %>%
  dplyr::select(c("S_No", "Sample_ID", "nm", "Phenotype", "total_contig")) %>%
  dplyr::filter(total_contig >= 1500)
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
qc %>%
  dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
  ggplot() +
  geom_point(aes(x = n50_contig_length, y = total_contig,
                 shape = Phenotype, color = Sample_ID),
             size = 2, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")

# p2
qc %>%
  dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
  ggplot() +
  geom_point(aes(x = min_contig_length, y = max_contig_length,
                 shape = Phenotype, color = Sample_ID),
             size = 2, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")

# p3
qc %>%
  dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
  ggplot() +
  geom_point(aes(x = mean_contig_length, y = median_contig_length,
                 shape = Phenotype, color = Sample_ID),
             size = 2, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")


# p3
qc %>%
  dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
  ggplot() +
  geom_point(aes(x = total_contig, y = total_contig_length,
                 shape = Phenotype, color = Sample_ID),
             size = 2, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")



# p4
qc %>%
  dplyr::filter(!c(S_No %in% excludedsmpl$S_No)) %>%
  dplyr::mutate(gc_content = (contig_percent_g + contig_percent_c)/2) %>%
  ggplot() +
  geom_point(aes(x = n50_contig_length, y = gc_content,
                 shape = Phenotype, color = Sample_ID),
             size = 2, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")





