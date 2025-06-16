#' Evaluate Ordination Coordinates Across Multiple Metrics
#'
#' This function evaluates the quality of a low-dimensional ordination (e.g., CORDE or PCoA)
#' based on clustering, classification, and neighborhood preservation metrics.
#'
#' @param D Original high-dimensional dissimilarity matrix (samples × samples).
#' @param Coordinates Ordination coordinates (samples × ncomp).
#' @param expl_var Vector of explained variance per component (usually from PCoA or sPLS).
#' @param Y Outcome vector (factor or numeric).
#' @param X Confounder matrix (samples × covariates), used to compute confounder-associated R².
#' @param ncomp Number of ordination components (default: 2).
#'
#' @return A named list of evaluation metrics, including:
#' \describe{
#'   \item{ICR}{Information Compression Rate: total variance explained by first two components.}
#'   \item{TW}{Trustworthiness: neighbor preservation from low- to high-dimensional space.}
#'   \item{Cont}{Continuity: neighbor preservation from high- to low-dimensional space.}
#'   \item{Rand}{Rand index of clustering agreement.}
#'   \item{ARand}{Adjusted Rand index.}
#'   \item{ASil}{Average silhouette width from PAM clustering.}
#'   \item{CARI}{Confounder-adjusted R² improvement (McFadden).}
#'   \item{rCMI}{Relative conditional mutual information between clustering and outcome.}
#' }
#' @export
Ordination_Eval <- function(
    D, Coordinates, expl_var, Y, X, ncomp = 2
) {
  # Load packages safely
  stopifnot(requireNamespace("fossil", quietly = TRUE))
  stopifnot(requireNamespace("factoextra", quietly = TRUE))
  stopifnot(requireNamespace("mclust", quietly = TRUE))
  stopifnot(requireNamespace("cluster", quietly = TRUE))
  stopifnot(requireNamespace("DescTools", quietly = TRUE))
  stopifnot(requireNamespace("infotheo", quietly = TRUE))

  measures <- list(
    ICR = NA, TW = NA, Cont = NA,
    Rand = NA, ARand = NA, ASil = NA,
    CARI = NA, rCMI = NA
  )

  # Information Compression Rate
  measures$ICR <- sum(expl_var[1:min(2, length(expl_var))])

  # Neighborhood preservation: Trustworthiness & Continuity
  D_low <- as.matrix(dist(Coordinates[, 1:2]))
  low_rank <- apply(D_low, 1, order)[-1, ]
  high_rank <- apply(D, 1, order)[-1, ]
  n <- nrow(Coordinates)
  k <- if (n <= 50) 5 else if (n <= 500) ceiling(n / 10) else 50

  trust <- 0
  cont <- 0
  for (i in 1:n) {
    u <- low_rank[1:k, i]
    ranks_high <- match(u, high_rank[, i])
    trust <- trust + sum((ranks_high[ranks_high > k] - k), na.rm = TRUE)

    v <- high_rank[1:k, i]
    ranks_low <- match(v, low_rank[, i])
    cont <- cont + sum((ranks_low[ranks_low > k] - k), na.rm = TRUE)
  }
  norm <- n * k * (2 * n - 3 * k - 1) / 2
  measures$TW <- 1 - (2 / norm) * trust
  measures$Cont <- 1 - (2 / norm) * cont

  # Clustering via PAM
  pam_res <- cluster::pam(D_low, k = nlevels(as.factor(Y)))
  cluster_labels <- pam_res$clustering
  measures$Rand <- fossil::rand.index(cluster_labels, as.integer(as.factor(Y)))
  measures$ARand <- mclust::adjustedRandIndex(cluster_labels, Y)
  measures$ASil <- pam_res$silinfo$avg.width

  # Confounder-Associated R² Improvement (CARI)
  model1 <- glm(Y ~ cluster_labels, family = binomial)
  model2 <- glm(Y ~ cluster_labels + X, family = binomial)
  r2_yc <- DescTools::PseudoR2(model1, which = "McFadden")
  r2_ycz <- DescTools::PseudoR2(model2, which = "McFadden")
  measures$CARI <- (r2_ycz - r2_yc) / r2_yc

  # Relative Conditional Mutual Information (rCMI)
  mi_cy <- infotheo::mutinformation(as.factor(cluster_labels), as.factor(Y))
  mi_cy_x <- infotheo::condinformation(as.factor(cluster_labels), as.factor(Y), as.data.frame(X))
  measures$rCMI <- mi_cy_x / mi_cy

  return(measures)
}
