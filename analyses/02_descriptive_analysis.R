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
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))

#............................................................
# HOSPTIAL EPI
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

mlstout <- mlst %>%
  dplyr::group_by(Phenotype, ST) %>%
  dplyr::summarize(
    n = n()
  )

# utility in MLST is mostly in the ST
mlstout <- mlstout %>%
  dplyr::filter(ST != "-") %>%
  dplyr::filter(Phenotype != "unknown") %>%
  tidyr::pivot_wider(., names_from = "Phenotype",
                     values_from = n) %>%
  dplyr::mutate(diff = Persistent - Resolving)
saveRDS(object = mlstout, file = "results/mlst_out.RDS")

#..................
# staphopiasccmec
#..................
staphmec <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphopiasccmec-20251118-160029/merged-results/staphopiasccmec.tsv")

#................
# staphtyper
#.................
agrstaphtyper <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphtyper-20251118-154929/merged-results/agrvate.tsv")
mecstaphtyper <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphtyper-20251118-154929/merged-results/sccmec.tsv")
spastaphtyper <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/staphtyper-20251118-154929/merged-results/spatyper.tsv")


#
# | Sample   | MLST  | agr    | SCCmec     | spa                    |
#   | -------- | ----- | ------ | ---------- | ---------------------- |
#   | S1       | ST105 | agr II | SCCmec IIa | t062                   |
#   | S10      | ST8   | agr I  | SCCmec IVa | t008                   |
#   | S11_S137 | ST105 | agr II | SCCmec IIa | t002                   |
#   | S13      | …     | agr I  | SCCmec IV… | t9121 / 34-24-… (etc.) |
#







#----


#............................................................
# AMR
#...........................................................
amrfinder <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/bactopia-20251118-102712/merged-results/amrfinderplus.tsv") %>%
  dplyr::rename(S_No = Name) %>%
  dplyr::left_join(., mtdt, by = "S_No")

# pivot wider
amrfinder_presence_absence <- amrfinder %>%
  dplyr::filter(`% Coverage of reference` > 0.95) %>%
  dplyr::filter(`% Identity to reference` > 0.975) %>%
  dplyr::filter(!is.na(`Protein id`)) %>%
  dplyr::mutate(geneid = paste0("p", `Protein id`, "_es", `Element symbol`, "_b", Start, "_e", Stop, "_s", Strand)) %>%
  dplyr::select(c("S_No", "geneid")) %>%
  dplyr::mutate(presence = 1) %>%
  tidyr::pivot_wider(., names_from = "geneid",
                     values_from = presence)


amrfinder_out <- amrfinder %>%
  dplyr::filter(`% Coverage of reference` > 0.95) %>%
  dplyr::filter(`% Identity to reference` > 0.975) %>%
  dplyr::filter(!is.na(`Protein id`)) %>%
  dplyr::mutate(geneid = paste0("p", `Protein id`, "_es", `Element symbol`, "_b", Start, "_e", Stop, "_s", Strand)) %>%
  dplyr::select(c("Phenotype", "S_No", "geneid"))

# out
amrfinder_presence_absence[is.na(amrfinder_presence_absence)] <- 0
saveRDS(object = amrfinder_presence_absence,
        file = "results/amr_presence_absence.RDS")
saveRDS(object = amrfinder_out,
        file = "results/amrfinder_out.RDS")

#............................................................
# plasmidfinder
#...........................................................
plasmidfinder <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/plasmidfinder-20251118-163037/merged-results/plasmidfinder.tsv") %>%
  dplyr::rename(S_No = Sample) %>%
  dplyr::left_join(., mtdt, by = "S_No") %>%
  dplyr::filter(Identity > 99)


#............................................................
# mobsuite
#...........................................................
mobsuite <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/mobsuite-20251118-160413/merged-results/mobsuite.tsv")

#............................................................
# rgi
#...........................................................
rgi <- readr::read_tsv("data/derived/prealm_ret/bactopia-runs/rgi-20251118-161628/merged-results/rgi.tsv")

