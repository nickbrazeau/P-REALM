#'
#' @description Conduct principal components analysis (PCA) on a presence/absence matrix.
#'    Output includes the raw components, the variance in the
#'   data explained by each component, and the loadings of each component also
#'   returned.
#'
#' @details Contributions of each variable are computed from the loading values
#'   (stored as "rotation" within the \code{prcomp} object). The percent
#'   contribution of a variable is defined as the absolute loading value for
#'   this variable, divided by the sum of loadings over all variables and
#'   multiplied by 100.
#'
#' @return Invisibly returns a list of class `prcomp` with the following
#'   components
#'   \itemize{
#'       \item{"sdev"}{ the standard deviations of the principal components
#'       (i.e., the square roots of the eigenvalues of the
#'       covariance/correlation matrix, though the calculation is actually done
#'       with the singular values of the data matrix).} \item{"rotation"}{ the
#'       matrix of variable loadings (i.e., a matrix whose columns contain the
#'       eigenvectors). The function \code{princomp} returns this in the element
#'       \code{loadings}.} \item{"center, scale"}{ the centering and scaling
#'       used.} \item{"x"}{ the value of the rotated data (the centred data
#'       multiplied by the rotation matrix). Hence, \code{cov(x)} is the
#'       diagonal matrix \code{diag(sdev^2)}.}
#'       \item{"var"}{ the variance in the data explained by each component.}
#'       \item{"contribution"}{ the percent contribution of a variable (i.e. a
#'       locus) to the overall variation.}
#'       }
#' @importFrom stats prcomp
#' @export

pca_pa <- function(x, num_components = 10) {
  #............................................................
  # computations
  #...........................................................
  pca <- prcomp(x)

  # compute variance explained
  pca$var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2) * 100

  # compute percent contribution of each locus
  pca$contribution <- abs(pca$rotation)
  pca$contribution <- sweep(pca$contribution, 2, colSums(pca$contribution), "/") * 100

  #............................................................
  # variance of components plot
  #............................................................
  # potential components
  nc <- min(length(colnames(pca$x)), num_components)
  # x has to be factored in order to preserve PC order
  compfct <- factor(colnames(pca$x)[seq_len(nc)],
                    levels = colnames(pca$x)[seq_len(nc)])
  outvarplot <- tibble::tibble(compfct = compfct, var = pca$var[seq_len(nc)]) %>%
    ggplot(., aes(x = compfct, y = var)) +
    geom_col() +
    theme_minimal() +
    labs(
      x = "Principal Component",
      y = "Percentage of variance explained",
      title = "Percentage of total variance explaned by PCA"
    )

  #............................................................
  # out
  #...........................................................
  out <- list(
    pca = pca,
    outvarplot = outvarplot
  )

}


#' @title 3d plotly pca
plot_PCA_3d <- function(pca_run, mtdt){
  col_palette <- viridis::plasma(n = 3)
  dat <- as.data.frame(pca_run$pca$x[, 1:3])
  dat$S_No <- rownames(dat)
  dat <- dplyr::left_join(x = dat, y = mtdt, by = "S_No")
  plotly::plot_ly(dat, x = ~PC1, y = ~PC2, z = ~PC3,
                  color = ~ Phenotype,
                  smpl = ~nm,
                  Phenotype = ~ Phenotype,
                  type = "scatter3d", colors = col_palette,
                  mode = "markers", marker = list(size = 5),
                  alpha = 0.8)
}

#' @title Look at PCA contributions

plot_loci_contributions <- function(pca_run, num_components = 3, contmin = 0.01) {
  # tidy
  dat <- tibble::as_tibble(pca_run$pca$contribution[,1:num_components])
  dat$locus <- names(pca_run$pca$contribution[,1])
  dat <- dat %>%
    tidyr::pivot_longer(., cols = -locus, names_to = "components", values_to = "contribution")

  # filter
  dat <- dat %>%
    dplyr::filter(contribution > contmin)

  # plot
  dat %>%
    ggplot() +
    geom_point(aes(x = locus, y = contribution),
             stat = "identity", position = position_dodge(),
             color = "#9ecae1") +
    theme_minimal() +
    theme(axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 10),
          axis.text.x = element_text(family = "Helvetica", angle = 90, hjust = 0.5, size = 8))
}


