## .................................................................................
## Purpose: GWAS controlling for lineage effects
##
## Author: Nick Brazeau
##
## Date: 10 Feb 2026
##
## Notes:
##    file:///Users/nbrazeau/Documents/Github/P-REALM/ignore/brms-dashboard/index.html
##    https://search.r-project.org/CRAN/refmans/spaMM/html/MaternCorr.html
##    https://www.image.ucar.edu/GSP/Software/Fields/Help/matern.cov.html
## .................................................................................
library(tidyverse)
library(spaMM)
library(pheatmap)
library(brms)

#++++++++++++++++++++++++++++++++++++++++++
### Section 0: Setup      ####
#++++++++++++++++++++++++++++++++++++++++++
#............................................................
# read in meta-data
#...........................................................
mtdt <- readr::read_csv("data/meta/combined_mtdt.csv") %>%
  dplyr::select(c("S_No", "Sample_ID", "Phenotype", "nm"))
excludedsmpl <- readr::read_csv("data/derived/excluded_smpls.csv")

#............................................................
# read in Panaroo Data
#...........................................................
# panaroo
panPA <- readRDS(file ="results/panaroo_out.rds")
panPAmat <- panPA$panaroo_pa_mat
panPAmatshell <- panPA$panaroo_pa_mat[, colnames(panPA$panaroo_pa_mat) %in% panPA$genomes$shellgenome]
panPAmatshelldf <- dplyr::bind_cols(S_No = rownames(panPAmatshell),
                                    panPAmatshell) %>%
  dplyr::left_join(., mtdt, by = "S_No") %>%
  dplyr::select(-c("Sample_ID", "nm")) %>%
  dplyr::select(c("S_No", "Phenotype", dplyr::everything())) %>%
  dplyr::mutate(S_No = factor(S_No)) %>%
  dplyr::filter(Phenotype %in% c("Persistent", "Resolving")) %>%
  dplyr::mutate(Phenotype = factor(Phenotype))

#............................................................
# Read in Structure Data to Control for
#...........................................................
# pca
pca <- readRDS("results/linear_assoc.RDS")$panPAshell_pca
# network
comm <- readRDS("results/network_out.RDS")$comm_shell_sample_jcc_net %>%
  dplyr::select(c("S_No", "community"))
# bring together
panPAmatshelldf_net <- dplyr::left_join(panPAmatshelldf, comm, by = "S_No") %>%
  dplyr::select(c("S_No", "Phenotype", "community", dplyr::everything()))

#++++++++++++++++++++++++++++++++++++++++++
### Section 1: Find our Covariances   ####
#++++++++++++++++++++++++++++++++++++++++++
# \Sigma = W \Lambda W^T
lambda <- pca$pca$sdev^2
lambda <- diag( lambda[1:10] )
W <- sweep(pca$pca$x, 2, pca$pca$sdev, "/")
W <- W[,1:10]
K <- W %*% lambda %*% t(W)
K_scaled <- K / mean(diag(K))
image(K_scaled)
pheatmap::pheatmap(K_scaled)
diag(K_scaled)
# work arounds
K_scaled_nearPD <- Matrix::nearPD(K_scaled)$mat
all(eigen(K_scaled_nearPD)$values > 0)
image(as.matrix(K_scaled_nearPD))
pheatmap::pheatmap(K_scaled_nearPD)
summary(as.vector(K_scaled - K_scaled_nearPD))
hist(as.vector(K_scaled - K_scaled_nearPD))

#++++++++++++++++++++++++++++++++++++++++++
### Section 2: Regression with Fixed PC Effects      ####
#++++++++++++++++++++++++++++++++++++++++++
# storage
modelret_pca <- tibble::tibble(loci = colnames(panPAmatshelldf[,3:ncol(panPAmatshelldf)]),
                               model = NA,
                               pvalue = NA)

# regression
for (i in 3:ncol(panPAmatshelldf)) {
  dat <- panPAmatshelldf[,c(1,2,i)]
  colnames(dat) <- c("S_No", "Phenotype", "gene_i")

  fit_pca_gwas <- spaMM::fitme(
      Phenotype ~ gene_i + corrMatrix(1 | S_No),
      data = dat,
      corrMatrix = K_scaled_nearPD,
      family = binomial(link = "logit")
    )
  # extract out infro
  modelret_pca$model[i-2] <- fit_pca_gwas
  modelret_pca$pvalue[i-2] <- 2*pnorm(abs(summary(fit_pca_gwas)$beta_table[2,3]), lower.tail = F)

}


#++++++++++++++++++++++++++++++++++++++++++
### Section 3: Regression with Fixed Network Effects      ####
#++++++++++++++++++++++++++++++++++++++++++
# storage
modelret_net <- tibble::tibble(loci = colnames(panPAmatshelldf_net[,4:ncol(panPAmatshelldf_net)]),
                               model = NA,
                               pvalue = NA)

# regression
for (i in 4:ncol(panPAmatshelldf_net)) {
  dat <- panPAmatshelldf[,c(1,2,3,i)]
  colnames(dat) <- c("S_No", "Phenotype", "community", "gene_i")

  fit_pca_gwas <- spaMM::fitme(
    Phenotype ~ gene_i + corrMatrix(1 | S_No),
    data = dat,
    corrMatrix = K_scaled_nearPD,
    family = binomial(link = "logit")
  )
  # extract out infro
  modelret_net$model[i-3] <- fit_pca_gwas
  modelret_net$pvalue[i-3] <- 2*pnorm(abs(summary(fit_pca_gwas)$beta_table[2,3]), lower.tail = F)

}


#++++++++++++++++++++++++++++++++++++++++++
### Section 4: Investigate further     ####
#++++++++++++++++++++++++++++++++++++++++++
rosetta <- readr::read_csv("~/Documents/Github/P-REALM/data/derived/prealm_ret/bactopia-runs/pangenome-20251118-225922/panaroo/gene_presence_absence.csv")
rosetta <- rosetta[,1:3] %>%
  magrittr::set_colnames(c("loci", "Non-unique Gene name", "Annotation"))
modelret_pca_rosetta <- dplyr::left_join(modelret_pca, rosetta, by = "loci")
hist(modelret_pca$pvalue)
filtmodelret_pca_rosetta <- modelret_pca_rosetta %>%
  dplyr::filter(pvalue < 0.05)
# filt
keeploci <- colnames(panPAmatshelldf) %in% filtmodelret_pca_rosetta$loci
keeploci[1:2] <- T
filtpanPAmatshelldf <- panPAmatshelldf[, keeploci]
pheatmap(filtpanPAmatshelldf[3:ncol(filtpanPAmatshelldf)])


#............................................................
# Bayesian Dig into those Loci
#...........................................................
bayesmodelret_pca <- tibble::tibble(loci = colnames(filtpanPAmatshelldf[,3:ncol(filtpanPAmatshelldf)]),
                               model = NA,
                               LCI = NA,
                               UCI = NA)

for (i in 1:(ncol(filtpanPAmatshelldf) - 2)) {
  dat <- filtpanPAmatshelldf[,c(1,2,(i+2))]
  colnames(dat) <- c("S_No", "Phenotype", "gene_i")

  fit_pca_bgwas <- brms::brm(
    Phenotype ~ gene_i + (1 | gr(S_No, cov = K_scaled_nearPD)),
    data   = dat,
    data2  = list(K_scaled_nearPD = K_scaled_nearPD),
    family = bernoulli(link = "logit"),
    prior  = c(
      prior(normal(0, 2.5), class = "b"),
      prior(normal(0, 5), class = "Intercept"),
      prior(cauchy(0, 5), class = "sd")
    ),
    chains = 3,
    iter   = 1e3,
    cores  = 2
  )
  # extract out infro
  bayesmodelret_pca$model[i] <- fit_pca_bgwas
  bayesmodelret_pca$LCI[i] <- summary(fit_pca_bgwas)$fixed$`l-95% CI`[2]
  bayesmodelret_pca$UCI[i] <- summary(fit_pca_bgwas)$fixed$`u-95% CI`[2]


}





#............................................................
# out
#...........................................................
out <- list(
  gwas_modelret_pca = modelret_pca,
  gwas_modelret_net = modelret_net,
  bayesmodelret_pca = bayesmodelret_pca
)

saveRDS(out, file = "results/mod_gwas.RDS")
