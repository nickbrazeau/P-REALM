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
panaroo_pa_data <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_presence_absence.Rtab")
keepcols <- colnames(panaroo_pa_data)[ !(colnames(panaroo_pa_data) %in% excludedsmpl$S_No) ]
panaroo_pa_data <- panaroo_pa_data %>%
  dplyr::select(c("Gene", keepcols))
panaroo_pa_csv_data <- readr::read_csv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_presence_absence.csv")

#%>%
  dplyr::select(c("Gene", "Non-unique Gene name", "Annotation", keepcols))
panaroo_gene_data <- readr::read_csv("data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_data.csv.gz") %>%
  dplyr::rename(S_No = gff_file) %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))

#++++++++++++++++++++++++++++++++++++++++++
### Section 1: Basic Presence Absence Table    ####
#++++++++++++++++++++++++++++++++++++++++++

#......................
# Presence Abscence natrix
#.....................
panaroo_pa_mat <- panaroo_pa_data %>%
  dplyr::select(-c("Gene")) %>%
  as.matrix()
rownames(panaroo_pa_mat) <- panaroo_pa_data$Gene
panaroo_pa_mat <- t(panaroo_pa_mat)


#++++++++++++++++++++++++++++++++++++++++++
### Section 2: Pangenome Statistics        ####
#++++++++++++++++++++++++++++++++++++++++++

# https://github.com/gtonkinhill/panstripe/


#............................................................
# out
#...........................................................
out <-list(
  panaroo_pa_data = panaroo_pa_data,
  panaroo_pa_mat = panaroo_pa_mat,
  panaroo_gene_data = panaroo_gene_data
)
saveRDS(out, file = "results/panaroo_out.rds")



