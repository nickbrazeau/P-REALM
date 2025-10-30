## .................................................................................
## Purpose: Metadata sorter and symlink writer
##
## Author: Nick Brazeau
##
## Date: 30 October, 2025
## .................................................................................
library(tidyverse)

#............................................................
# Read in Raw Data
#...........................................................
fq_raw <- readr::read_csv("data/meta/fastq_filenames.csv")
mrsa_phen <- readr::read_csv("data/meta/MRSA_phenotypes.csv")

#............................................................
# Parse Raw Data
#...........................................................

fq_raw <- fq_raw %>%
  dplyr::mutate(S_No = purrr::map_chr(fastqfiles, function(x){
    stringr::str_extract(x, "(?<=BS-|4988-).*?(?=_R[12]|_L007)")

  }))

fq_raw_R1 <- fq_raw[stringr::str_detect(fq_raw$fastqfiles, "R1"), ]  %>%
  dplyr::rename(R1 = fastqfiles)
fq_raw_R2 <- fq_raw[stringr::str_detect(fq_raw$fastqfiles, "R2"), ]  %>%
  dplyr::rename(R2 = fastqfiles)

fq_join <- dplyr::left_join(fq_raw_R1, fq_raw_R2, by = "S_No") %>%
  dplyr::select(S_No, R1, R2)

# now mrsa
mrsa_phen <- readr::read_csv("~/Desktop/MRSA_phenotypes.csv") %>%
  dplyr::select(c("S_No", "Sample_ID_Phenotype", "Clin"))

# bring together
combined <- dplyr::full_join(mrsa_phen, fq_join, by = "S_No")

# write out
readr::write_csv(combined, file = "data/meta/combined_mtdt.csv")

#............................................................
# Make Symlink
#...........................................................
arch <- combined %>%
  dplyr::select(c("S_No", "R1", "R2")) %>%
  tidyr::pivot_longer(., cols = c("R1", "R2"),
                      names_to = "drop",
                      values_to = "reads") %>%
  dplyr::select(-c("drop")) %>%
  dplyr::mutate(reads = paste0("/hpc/group/taylorlab/raw_sequencing_reads/nfbwork/staphaureus/",
                               reads),
                out = paste0("/hpc/group/taylorlab/users/nfb/projects/P-REALM/fastq/",
                             S_No, ifelse(stringr::str_detect(reads, "R1"),
                                          "_R1.fastq.gz", "_R2.fastq.gz"))
  )

# out
readr::write_csv(arch, file = "analyses/00_organization/symlink_architecture.csv")


