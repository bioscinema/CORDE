#’ Identify Driving Taxa Using Mantel Correlation and LEfSe‐Style LDA Scores
#’
#’ Quantify each taxon’s contribution to sample‐level ordination patterns by combining:
#’ 1) Mantel correlations between taxa abundance profiles and ordination embeddings,
#’ 2) LDA‐style group separation scores on the sample coordinates, and
#’ 3) Taxa loadings on the ordination axes.
#’
#’ @param OTU Numeric matrix of taxa abundances (rows = taxa, columns = samples).
#’ @param CORDE_res List output from the \code{CORDE} workflow, with components:
#’   \describe{
#’     \item{\code{Loadings}}{Numeric matrix of taxa loadings (taxa × axes).}
#’     \item{\code{Embeddings}}{Numeric matrix of sample coordinates (samples × axes).}
#’     \item{\code{Coordinates}}{Numeric matrix used for LDA (samples × features).}
#’     \item{\code{Y}}{Factor or character vector of sample‐level class labels.}
#’   }
#’
#’ @return A named list with two components:
#’   \describe{
#’     \item{\code{Coord_score}}{Numeric vector of per‐axis LDA separation scores, normalized as \eqn{|lda|/\sqrt{1+lda^2}}.}
#’     \item{\code{Taxon_score}}{Numeric vector of overall driving scores for each taxon (length = number of rows in \code{OTU}), normalized as \eqn{Score/\sqrt{1+Score^2}}.}
#’   }
#’   Names of both vectors correspond to \code{rownames(OTU)}.
#’
#’ @details
#’ The driving score is computed as:
#’ \enumerate{
#’   \item Compute the cosine‐similarity matrix \eqn{C} between each taxa profile (columns of \code{t(OTU)}) and each ordination axis (\code{Embeddings}).
#’   \item Obtain per‐axis LDA‐style separation scores (\code{lda_score}) on \code{Coordinates} and \code{Y}.
#’   \item Multiply \eqn{C} by the loadings matrix (\code{Loadings}) to get taxa contributions to each axis.
#’   \item Project these onto \code{lda_score} and normalize to yield \code{Taxon_score}.
#’ }
#’
#’ @seealso
#’ \code{\link{get_mantel_C_matrix}}, \code{\link{compute_lefse_style_score}}
#’
#’
#’ @export
driving.taxa.lda <- function(
    OTU, # OTU matrix, taxa as rows
    CORDE_res
){

  W <- CORDE_res$Loadings

  E <- CORDE_res$Embeddings
  C_mat <- get_mantel_C_matrix(t(OTU), E)

  X <- CORDE_res$Coordinates
  Y <- CORDE_res$Y
  lda_score <- compute_lefse_style_score(X,Y)

  M <- C_mat %*% W

  Score <- as.vector(M %*% lda_score)

  Score <- Score/sqrt(1+Score^2)
  names(Score) <- rownames(OTU)

  return(list(Coord_score = abs(lda_score)/sqrt(1+lda_score^2), Taxon_score = Score))
}
