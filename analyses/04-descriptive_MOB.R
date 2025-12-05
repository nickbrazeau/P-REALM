## .................................................................................
## Purpose: Basic descriptives from Bactopia
##
## Author: Nick Brazeau
##
## Date: 21 November, 2025
##
## Notes:
## .................................................................................
library(tidyverse)
source("R/clean_accessions.R")
#............................................................
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")


#............................................................
# Mobile Units
#...........................................................
#..................
# plasmidfinder
#..................
# contig provided (need to clean)
# start end sense (based on order of start/stop) is encoded in
# Position_in_contig
plasmidfinder <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/plasmidfinder-20251118-163037/merged-results/plasmidfinder.tsv")
colnames(plasmidfinder) <- paste0("pf_", gsub(colnames(plasmidfinder), pattern = " ", replacement = "_"))
# tidy
plasmidfinder <- plasmidfinder %>%
  dplyr::rename(S_No = pf_Sample) %>%
  dplyr::filter(pf_Identity > 99) %>%
  dplyr::rename(acc = pf_Accession_number) %>%
  dplyr::mutate(acc = clean_accessions(acc, source = "plasmidfinder")) %>%
  dplyr::filter(!is.na(acc)) %>%
  dplyr::select(c("S_No", "pf_Database", "pf_Plasmid", "acc")) %>%
  dplyr::left_join(., mtdt, by = "S_No") %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))

#.................
# mobsuite
#..................
mobsuite <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/mobsuite-20251118-160413/merged-results/mobsuite.tsv")
colnames(mobsuite) <- paste0("mb_", gsub(colnames(mobsuite), pattern = " ", replacement = "_"))
mobsuite <- mobsuite %>%
  dplyr::mutate(S_No = stringr::str_split_fixed(mb_sample_id, ":", n = 2)[,1]) %>%
  dplyr::rename(acc = `mb_rep_type_accession(s)`) %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))

# clean accessions
mobsuite <- process_mobsuite_accessions(mobsuite)


# drop to relevant
mobsuite <- mobsuite %>%
  dplyr::select(c("S_No", "acc", "mb_primary_cluster_id", "mb_predicted_mobility",
                  "mb_mpf_type", "mb_relaxase_type(s)", "mb_orit_type(s)"))

#............................................................
# come together
#...........................................................
pfmobcomb <- dplyr::full_join(plasmidfinder, mobsuite,
                              by = c("S_No", "acc"))
# confusion mat
mobconfusionmatrix <- tibble::tibble(
  concordant = sum(!is.na(pfmobcomb$pf_Plasmid) & !is.na(pfmobcomb$mb_primary_cluster_id)),
  MOBmissedPFdetect = sum(is.na(pfmobcomb$mb_primary_cluster_id)),
  PFmissedMOBdetect = sum(is.na(pfmobcomb$pf_Plasmid))
)


#......................
# Presence Abscence Table
#.....................
MOBPA <- pfmobcomb %>%
  dplyr::mutate(mobunit = case_when(
    !is.na(pf_Plasmid) & !is.na(mb_primary_cluster_id) ~ mb_primary_cluster_id,
    is.na(mb_primary_cluster_id) ~ pf_Plasmid,
    is.na(pf_Plasmid) ~ mb_primary_cluster_id,
  )
  ) %>%
  dplyr::select(c("S_No", "mobunit")) %>%
  dplyr::mutate(p = 1) %>%
  tidyr::pivot_wider(., names_from = "mobunit",
                     values_from = p)

MOBPAmatrix <- matrix(data = NA, nrow = nrow(MOBPA), ncol = ncol(MOBPA)-1)
rownames(MOBPAmatrix) <- MOBPA$S_No
colnames(MOBPAmatrix) <- colnames(MOBPA)[2:ncol(MOBPA)]
for(i in 2:ncol(MOBPA)) {
  MOBPAmatrix[,i-1] <- as.numeric(sapply(MOBPA[i][[1]], function(x){is.null(x)}))
}


#......................
# out
#......................
out <-list(
  plasmidfinder = plasmidfinder,
  mobsuite = mobsuite,
  pfmobcomb = pfmobcomb,
  mobconfusionmatrix = mobconfusionmatrix,
  MOBPA = MOBPA,
  MOBPAmatrix = MOBPAmatrix
)
saveRDS(out, file = "results/MOB_out.rds")
