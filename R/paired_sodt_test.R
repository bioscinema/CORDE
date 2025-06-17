#' Paired SoDT Test for Group-Level Changes in Dissimilarity Matrices
#'
#' Performs a paired Sum-of-Distances Test (SoDT) to compare group-level
#' dispersion and separation across two dissimilarity matrices (e.g., before and after treatment),
#' while accounting for within-group variances and providing permutation-based p-values.
#'
#' Internally, this function uses a high-performance C++ implementation
#' (via Rcpp and Armadillo) to accelerate permutation-based inference.
#'
#' @param D0 A symmetric dissimilarity matrix for condition 0 (e.g., pre-treatment).
#' @param D1 A symmetric dissimilarity matrix for condition 1 (e.g., post-treatment).
#' @param Y A vector of group labels (e.g., treatment groups or clusters).
#' @param nperm Number of permutations to compute empirical p-values (default = 999).
#' @param seed Integer seed for reproducibility (default = 2025).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Source}{The source of variation (e.g., Between-group, Within-group).}
#'   \item{Delta}{The observed change in each statistic from \code{D0} to \code{D1}.}
#'   \item{P_value}{The permutation-based two-sided p-value for each statistic.}
#' }
#'
#' @examples
#' \dontrun{
#' D0 <- as.matrix(vegan::vegdist(iris[1:50, 1:4]))
#' D1 <- as.matrix(vegan::vegdist(iris[51:100, 1:4]))
#' Y <- rep(c("A", "B"), each = 25)
#' paired_sodt_test(D0, D1, Y)
#' }
#'
#' @useDynLib CORDE, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
paired_sodt_test <- function(D0, D1, Y, nperm = 999, seed = 2025) {
  set.seed(seed)

  if (!is.matrix(D0) || !all(D0 == t(D0))) stop("D0 must be symmetric.")
  if (!is.matrix(D1) || !all(D1 == t(D1))) stop("D1 must be symmetric.")
  if (any(diag(D0) != 0) || any(diag(D1) != 0)) stop("D must have 0 diagonal.")

  n <- nrow(D0)
  groups <- unique(Y)
  d <- length(groups)

  compute_gram <- function(D) {
    D2 <- D^2
    J <- diag(n) - matrix(1, n, n) / n
    -0.5 * J %*% D2 %*% J
  }

  hat_matrix <- function(Yvec) {
    Y_mat <- model.matrix(~ 0 + factor(Yvec))
    Y_mat %*% solve(t(Y_mat) %*% Y_mat) %*% t(Y_mat)
  }

  compute_sodt_simple <- function(G, H_Y) {
    trace_G <- sum(diag(G))
    SSB <- sum(H_Y * G)
    SSW <- trace_G - SSB
    F_stat <- (SSB / (d - 1)) / (SSW / (n - d))
    R <- diag(n) - H_Y
    GR <- R %*% G
    SSW_g <- sapply(groups, function(g) {
      idx <- which(Y == g)
      sum(diag(GR[idx, ] %*% t(R[idx, , drop = FALSE])))
    })
    list(SSB = SSB, SSW = SSW, SSW_g = SSW_g, F = F_stat)
  }

  # Compute Gram matrices and observed stats
  G0 <- compute_gram(D0)
  G1 <- compute_gram(D1)
  H_Y <- hat_matrix(Y)
  sodt0 <- compute_sodt_simple(G0, H_Y)
  sodt1 <- compute_sodt_simple(G1, H_Y)

  delta_SSB <- sodt1$SSB - sodt0$SSB
  delta_SSW <- sodt1$SSW - sodt0$SSW
  delta_SSW_g <- sodt1$SSW_g - sodt0$SSW_g
  delta_F <- sodt1$F - sodt0$F

  # Prepare permutation inputs
  G0_list <- replicate(nperm, G0, simplify = FALSE)
  G1_list <- replicate(nperm, G1, simplify = FALSE)
  HY_list <- lapply(seq_len(nperm), function(i) hat_matrix(sample(Y)))

  # Call Rcpp backend
  perm_res <- sodt_permutation_loop(G0_list, G1_list, HY_list)

  # Compute empirical p-values
  pval <- function(delta, dist) {
    (1 + sum(abs(dist) >= abs(delta))) / (length(dist) + 1)
  }

  result_table <- data.frame(
    Source = c(
      "Between-group", "Within-group", "PERMANOVA F",
      paste0("Within-group (", groups, ")")
    ),
    Delta = c(delta_SSB, delta_SSW, delta_F, delta_SSW_g),
    P_value = c(
      pval(delta_SSB, perm_res$delta_SSB),
      pval(delta_SSW, perm_res$delta_SSW),
      pval(delta_F, perm_res$delta_F),
      sapply(1:d, function(j) pval(delta_SSW_g[j], perm_res$delta_SSW_g[, j]))
    ),
    stringsAsFactors = FALSE
  )

  return(result_table)
}


