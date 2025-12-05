## .................................................................................
## Purpose: Basic descriptives from Bactopia
##
## Author: Nick Brazeau
##
## Date: 21 November, 2025
##
## Notes: Intersecting AMRFinder and RGI by gene names. Cannot easily do genomic
## coordinates bc RGI use protein fasta from contig faa and AMR finder is using bakta gff3
## still imperfect off of gene names because of typos
## .................................................................................
library(tidyverse)
source("R/amrfinder_fuzzy_join.R")
#............................................................
# storage for out
#...........................................................
out <- list()

#............................................................
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")

#............................................................
# AMRFinder plus
#............................................................
amrfinder <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/bactopia-20251118-102712/merged-results/amrfinderplus.tsv") %>%
  dplyr::rename(S_No = Name)
colnames(amrfinder) <- gsub(colnames(amrfinder), pattern = " ", replacement = "_")
amrfinder <- amrfinder %>%
  dplyr::select(c("S_No", "Element_symbol", "Element_name", "%_Identity_to_reference", "%_Coverage_of_reference", "Start", "Stop")) %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))
# remove duplicates - assume was from subclass
amrfinder <- amrfinder[!duplicated(amrfinder), ]
colnames(amrfinder)[2:ncol(amrfinder)] <- paste0("AMRF_", colnames(amrfinder)[2:ncol(amrfinder)])

#............................................................
# rgi
#...........................................................
rgi <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/rgi-20251118-161628/merged-results/rgi.tsv")
colnames(rgi) <- gsub(colnames(rgi), pattern = " ", replacement = "_")
rgi <- rgi %>%
  dplyr::mutate(S_No = purrr::map_chr(Contig, function(x){
    stringr::str_split_fixed(x, "_0", n= 2)[,1]})) %>%
  dplyr::select("S_No", "Best_Hit_ARO", "AMR_Gene_Family", "Drug_Class", "Pass_Bitscore", "Hit_Start", "Hit_End") %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))
colnames(rgi)[2:ncol(rgi)] <- paste0("RGI_", colnames(rgi)[2:ncol(rgi)])


#............................................................
# Fuzzy Join
#...........................................................
# unique hit from samples for query seq
uniRGIARO <- unique(rgi$RGI_Best_Hit_ARO)
# unique hit from samples for query seq
uniAMR <- unique(amrfinder$AMRF_Element_symbol)


amrfinder_key <- amrfinder %>%
  dplyr::mutate(RGI_Best_Hit_ARO_joiner = purrr::map_chr(AMRF_Element_symbol,
                                                         amrfuzzjoin,
                                                         RGI_Best_Hit_ARO = uniRGIARO))


#............................................................
# now can bring together
#...........................................................
rgi <- rgi %>%
  dplyr::mutate(RGI_Best_Hit_ARO_joiner = stringr::str_to_lower(RGI_Best_Hit_ARO))

# full database
AMRRGIcomb <- dplyr::left_join(amrfinder_key, amrfinder) %>%
  dplyr::full_join(., rgi, by = c("S_No", "RGI_Best_Hit_ARO_joiner"),
                   relationship = "many-to-many") %>%
  dplyr::select(-c("RGI_Best_Hit_ARO_joiner")) %>%
  dplyr::mutate(
    length_comp = abs(AMRF_Stop - AMRF_Start) / abs(RGI_Hit_End - RGI_Hit_Start)
  )

# now tidy out
AMRRGIcombtidy <- AMRRGIcomb %>%
  dplyr::select(c("S_No", "AMRF_Element_symbol", "RGI_Best_Hit_ARO")) %>%
  dplyr::mutate(RGI_Best_Hit_ARO = ifelse(RGI_Best_Hit_ARO == "", NA, RGI_Best_Hit_ARO))

#......................
# confusion matrix
# this considers sample level agreement
#......................
smplconfusionmatrix <- tibble::tibble(
  concordant = sum(!is.na(AMRRGIcombtidy$AMRF_Element_symbol) & !is.na(AMRRGIcombtidy$RGI_Best_Hit_ARO)),
  AMRmissedRGIdetect = sum(is.na(AMRRGIcombtidy$AMRF_Element_symbol)),
  RGImissedAMRdetect = sum(is.na(AMRRGIcombtidy$RGI_Best_Hit_ARO))
)

gntmp <- AMRRGIcombtidy %>%
  dplyr::select(-c("S_No")) %>%
  dplyr::filter(!duplicated(.))
# see the csv gntmp_annoted for apparent mismatch

# confusion mat
geneconfusionmatrix <- tibble::tibble(
  concordant = sum(!is.na(gntmp$AMRF_Element_symbol) & !is.na(gntmp$RGI_Best_Hit_ARO)),
  AMRmissedRGIdetect = sum(is.na(gntmp$AMRF_Element_symbol)),
  RGImissedAMRdetect = sum(is.na(gntmp$RGI_Best_Hit_ARO))
)


#......................
# Presence Abscence Table
#.....................
# TODO why are there duplicates
AMRPA <- AMRRGIcombtidy %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(gene = case_when(
    !is.na(RGI_Best_Hit_ARO) & !is.na(AMRF_Element_symbol) ~ AMRF_Element_symbol,
    is.na(RGI_Best_Hit_ARO) ~ AMRF_Element_symbol,
    is.na(AMRF_Element_symbol) ~ RGI_Best_Hit_ARO,
  )
  ) %>%
  dplyr::select(c("S_No", "gene")) %>%
  dplyr::mutate(p = 1) %>%
  tidyr::pivot_wider(., names_from = "gene",
                     values_from = p)

# make mat
AMRPAmat <- AMRPA %>%
  dplyr::select(-c(S_No)) %>%
  as.matrix()
rownames(AMRPAmat) <- AMRPA$S_No
AMRPAmat[is.na(AMRPAmat)] <- 0

#......................
# out
#......................
out <-list(
  amrfinder = amrfinder,
  rgi = rgi,
  AMRRGIcomb = AMRRGIcomb,
  AMRRGIcombtidy = AMRRGIcombtidy,
  smplconfusionmatrix = smplconfusionmatrix,
  gntmp = gntmp,
  geneconfusionmatrix = geneconfusionmatrix,
  AMRPA = AMRPA,
  AMRPAmat = AMRPAmat
)
saveRDS(out, file = "results/AMR_out.rds")
