## .................................................................................
## Purpose: Pull in and tidy up SNPs
##
## Author: Nick Brazeau
##
## Date: 26 November, 2025
##
## Notes:
## .................................................................................
library(tidyverse)
library(plotly)
library(vcfR)
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

# snippy
vcf <- vcfR::read.vcfR("data/derived/prealm_ret/bactopia-runs/snippy-20251118-163505/snippy-core/core-snp.vcf.gz")
POSCHROM <- vcf@fix[,1:2]
vcfgt <- vcf %>%
  vcfR::extract.gt(., element = "GT", as.numeric = T)

#++++++++++++++++++++++++++++++++++++++++++
### Section 1: Filter       ####
#++++++++++++++++++++++++++++++++++++++++++
# drop excluded samples
vcfgt <- vcfgt[, !(colnames(vcfgt) %in% excludedsmpl$S_No)]

# filter loci based on heterozygosity
locihet <- rowMeans(vcfgt)
# assume a 5% erorr rate
locikeep <- locihet > 0.05 &  locihet < 0.95
# bring togeth
vcfclean <- dplyr::bind_cols(POSCHROM, vcfgt)[locikeep, ]

#++++++++++++++++++++++++++++++++++++++++++
### out    ####
#++++++++++++++++++++++++++++++++++++++++++
saveRDS(vcfclean, file = "results/vcf_out.RDS")
