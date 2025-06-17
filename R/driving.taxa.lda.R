#' Identify Driving Taxa Using Mantel Correlation and LEfSe-Style LDA Scores
#'
#' This function quantifies the importance of taxa in driving sample-level ordination patterns
#' by integrating Mantel correlations between taxa profiles and ordination embeddings,
#' LDA-based separation scores, and ordination loadings.
#'
#' @param OTU A numeric matrix of OTU/taxa abundances with taxa as rows and samples as columns.
#' @param CORDE_res A list output from the \code{CORDE} method containing elements:
#'   \code{Loadings} (taxa loadings), \code{Embeddings} (sample coordinates),
#'   \code{Coordinates} (low-dimensional representations), and \code{Y} (class labels).
#' @param dist_method Character string specifying the dissimilarity metric used,
#'   consistent with the one used in the \code{CORDE} method (e.g., "bray").
#' @param n_cores Number of CPU cores to use for parallel computation (default = 4).
#'
#' @return A named numeric vector of driving scores for each taxon.
#'   Higher scores indicate stronger contributions to the observed separation in the ordination space.
#'
#' @details
#' The driving score integrates three components:
#' \enumerate{
#'   \item \strong{Mantel correlation} between individual taxa and ordination axes.
#'   \item \strong{LDA-style group separation scores} from sample ordination coordinates.
#'   \item \strong{Normalized taxa loadings} indicating taxa contributions to axes.
#' }
#' Each component is normalized, and scores are combined via weighted matrix multiplication,
#' scaled by 100 for interpretability.
#'
#'
#' @seealso \code{\link{get_mantel_C_matrix}}, \code{\link{compute_lefse_style_score}}
#'
#' @export
driving.taxa.lda <- function(
    OTU, # OTU matrix, taxa as rows
    CORDE_res, # CORDE output
    dist_method, # dissimilarity measure (the same as CORDE)
    n_cores=4
){

  W <- CORDE_res$Loadings

  E <- CORDE_res$Embeddings
  C_mat <- get_mantel_C_matrix(t(OTU), E, dist_method, n_perm = 0, n_cores=n_cores)

  X <- CORDE_res$Coordinates
  Y <- CORDE_res$Y
  lda_score <- compute_lefse_style_score(X,Y)

  W_norm <- apply(W, 2, function(x) x / sqrt(sum(x^2)))
  C_norm <- t(apply(C_mat, 1, function(x) x / sqrt(sum(x^2))))

  Score <- 100*as.vector(C_norm %*% W_norm %*% lda_score)
  names(Score) <- rownames(OTU)

  return(Score)
}
