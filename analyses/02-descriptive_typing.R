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
# Genomic Typing
#...........................................................
#......................
# MLST
#......................
mlst <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/bactopia-20251118-102712/merged-results/mlst.tsv",
                        col_names = F)
colnames(mlst) <- c("sample", "species", "ST", "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL")
mlst <- mlst %>%
  dplyr::mutate(S_No = stringr::str_split_fixed(sample, pattern = "\\.", n = 3)[,1]) %>%
  dplyr::left_join(., mtdt, by = "S_No")


# utility in MLST is mostly in the ST
mlstsummout <- mlst %>%
  dplyr::group_by(Phenotype, ST) %>%
  dplyr::summarize( n = n() ) %>%
  dplyr::filter(ST != "-") %>%
  dplyr::filter(Phenotype != "unknown") %>%
  tidyr::pivot_wider(., names_from = "Phenotype",
                     values_from = n) %>%
  dplyr::mutate(diff = Persistent - Resolving)

# append out
out <- append(out, list(mlstsummout = mlstsummout))

#................
# staphtyper
#.................
# agr liftover
agrkey <- tibble::tibble(agr_group = c(paste0("gp", 1:4)),
          agr = c("I", "II", "III", "IV"))

agrstaphtyper <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphtyper-20251118-154929/merged-results/agrvate.tsv") %>%
  dplyr::rename(S_No = `#filename`) %>%
  dplyr::left_join(., agrkey, by = "agr_group")
mecstaphtyper <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphtyper-20251118-154929/merged-results/sccmec.tsv")  %>%
  dplyr::rename(S_No = sample,
                mec_type = type,
                mec_subtype = subtype) %>%
  dplyr::left_join(., mtdt, by = "S_No")
spastaphtyper <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphtyper-20251118-154929/merged-results/spatyper.tsv") %>%
  dplyr::mutate(S_No = purrr::map_chr(`Sequence name`, function(x){
    stringr::str_split_fixed(x, "_0", n= 2)[,1]})) %>%
  dplyr::left_join(., mtdt, by = "S_No") %>%
  dplyr::rename(spa = Type)


#............................................................
# come together
#...........................................................
comb <- dplyr::full_join(mlst, agrstaphtyper, by = "S_No") %>%
  dplyr::full_join(., mecstaphtyper, by = "S_No") %>%
  dplyr::full_join(., spastaphtyper, by = "S_No") %>%
  dplyr::select(c("nm", "S_No", "Phenotype",
                  "agr",
                  "mec_type", "mec_subtype",
                  "spa")) %>%
  dplyr::filter(Phenotype != "unknown") %>%
  dplyr::filter(!(S_No %in% excludedsmpl$S_No))


#......................
# out
#......................
out <- append(out, list(typing_summ = comb))
saveRDS(out, file = "results/typing_out.rds")


