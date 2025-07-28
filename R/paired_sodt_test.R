#’ Paired Sum‐of‐Distances Test (SoDT) for Paired Dissimilarity Matrices
#’
#’ Performs a paired Sum‐of‐Distances Test to compare group‐level dispersion and separation
#’ between two symmetric dissimilarity matrices (e.g. pre‐ vs. post‐treatment), while
#’ accounting for within‐group variability.  Empirical p‐values are obtained via two
#’ permutation schemes—group‐label shuffling (between‐group) and within‐sample swapping.
#’
#’ @param D0 Numeric matrix.  A symmetric dissimilarity matrix for condition 0
#’   (e.g. “pre” measurements).  Must have zeros on the diagonal.
#’ @param D1 Numeric matrix.  A symmetric dissimilarity matrix for condition 1
#’   (e.g. “post” measurements).  Must have zeros on the diagonal.
#’ @param Y  Factor or character vector of length n, giving the group label for each
#’   of the n samples.  Groups define the “within‐group” strata for testing.
#’ @param nperm Integer.  Number of permutations to use when computing empirical
#’   p‐values (default = 999).
#’ @param seed Integer.  Random seed for reproducibility of the permutation schemes
#’   (default = 2025).
#’
#’ @return A `data.frame` with one row per test statistic:
#’   \itemize{
#’     \item **Between‐group**: change in between‐group sum‐of‐squares (ΔSSB).
#’     \item **Within‐group**: change in within‐group sum‐of‐squares (ΔSSW).
#’     \item **PERMANOVA F**: change in the classic F‐statistic.
#’     \item **Within‐group (G)**: per‐group ΔSSW for each level of `Y`.
#’   }
#’   Columns are:
#’   \describe{
#’     \item{Source}{Name of the statistic.}
#’     \item{Delta}{Observed difference \(\Delta\) = statistic\(_{D1}\) – statistic\(_{D0}\).}
#’     \item{P_value}{Two‐sided empirical p‐value from the corresponding permutation null.}
#’   }
#’
#’ @details
#’ Internally, this function:
#’ \enumerate{
#’   \item Computes centering‐based Gram matrices \(G = -\tfrac12\,J D^2 J\) for each input.
#’   \item Calculates observed SSB, SSW, and F by projecting onto the group hat‐matrix.
#’   \item Uses a high‐performance C++ routine (via Rcpp/Armadillo) to accelerate
#’     permutation inference in two schemes:
#’     \itemize{
#’       \item *Between‐group*: shuffle the labels in `Y` and recompute the hat‐matrix.
#’       \item *Paired‐swap*: randomly swap entire rows/columns of \(`D0`,`D1`\) within samples.
#’     }
#’   \item Returns empirical p‐values with the standard \((1 + \#\{\lvert\Delta_{\rm perm}\rvert \ge \lvert\Delta_{\rm obs}\rvert\})/(n_{\rm perm}+1)\) adjustment.
#’ }
#’

#’
#’ @useDynLib CORDE, .registration = TRUE
#’ @importFrom Rcpp evalCpp
#’ @export
paired_sodt_test <- function(D0, D1, Y, nperm = 999, seed = 2025) {
  # ... your function body ...
}

paired_sodt_test <- function(D0, D1, Y, nperm = 999, seed = 2025) {
  set.seed(seed)

  # 1) sanity checks -------------------------------------------------------
  if (!is.matrix(D0) || !all(D0 == t(D0)))
    stop("D0 must be symmetric.")
  if (!is.matrix(D1) || !all(D1 == t(D1)))
    stop("D1 must be symmetric.")
  if (any(diag(D0) != 0) || any(diag(D1) != 0))
    stop("D0/D1 must have zero diagonals.")

  n      <- nrow(D0)
  groups <- unique(Y)
  d      <- length(groups)

  # 2) helper closures -----------------------------------------------------
  compute_gram <- function(D) {
    D2 <- D^2
    J  <- diag(n) - matrix(1, n, n) / n
    -0.5 * J %*% D2 %*% J
  }

  hat_matrix <- function(Yvec) {
    Y_mat <- model.matrix(~ 0 + factor(Yvec))
    H     <- Y_mat %*% solve(crossprod(Y_mat)) %*% t(Y_mat)
    H
  }

  compute_sodt_simple <- function(G, H) {
    trace_G <- sum(diag(G))
    SSB     <- sum(H * G)
    SSW     <- trace_G - SSB
    F_stat  <- (SSB / (d - 1)) / (SSW / (n - d))

    R       <- diag(n) - H
    GR      <- R %*% G
    SSW_g   <- sapply(groups, function(g) {
      idx <- which(Y == g)
      sum(diag(GR[idx, ] %*% t(R[idx, , drop = FALSE])))
    })

    list(SSB = SSB, SSW = SSW, SSW_g = SSW_g, F = F_stat)
  }

  # 3) observed statistics -------------------------------------------------
  G0   <- compute_gram(D0)
  G1   <- compute_gram(D1)
  H_Y  <- hat_matrix(Y)

  sodt0 <- compute_sodt_simple(G0, H_Y)
  sodt1 <- compute_sodt_simple(G1, H_Y)

  delta_SSB   <- sodt1$SSB   - sodt0$SSB
  delta_SSW   <- sodt1$SSW   - sodt0$SSW
  delta_SSW_g <- sodt1$SSW_g - sodt0$SSW_g
  delta_F     <- sodt1$F     - sodt0$F

  # 4) between‐group permutations (shuffle Y only) ------------------------
  H_list  <- replicate(nperm, {
    H_perm <- hat_matrix(sample(Y))
    H_perm
  }, simplify = FALSE)
  G0_list <- replicate(nperm, G0, simplify = FALSE)
  G1_list <- replicate(nperm, G1, simplify = FALSE)

  perm1 <- sodt_permutation_loop(G0_list, G1_list, H_list)
  # unpack
  delta_SSB_perm <- perm1$delta_SSB
  delta_SSW_perm <- perm1$delta_SSW
  delta_F_perm   <- perm1$delta_F

  # 5) paired‐swap permutations (swap rows/cols of D0/D1) -----------------
  H_const  <- H_Y
  HY_list2 <- replicate(nperm, H_const, simplify = FALSE)

  G0b_list <- vector("list", nperm)
  G1b_list <- vector("list", nperm)
  for (b in seq_len(nperm)) {
    swap_idx <- which(rbinom(n, 1, 0.5) == 1)

    D0_b <- D0; D1_b <- D1
    D0_b[swap_idx, ]    <- D1[swap_idx, ]
    D0_b[, swap_idx]    <- D1[, swap_idx]
    D1_b[swap_idx, ]    <- D0[swap_idx, ]
    D1_b[, swap_idx]    <- D0[, swap_idx]

    G0b_list[[b]] <- compute_gram(D0_b)
    G1b_list[[b]] <- compute_gram(D1_b)
  }

  perm2 <- sodt_permutation_loop(G0b_list, G1b_list, HY_list2)
  delta_SSW_g_perm <- perm2$delta_SSW_g

  # 6) empirical p‐value helper --------------------------------------------
  pval_emp <- function(obs, dist) {
    (1 + sum(abs(dist) >= abs(obs))) / (length(dist) + 1)
  }

  # 7) assemble results ----------------------------------------------------
  result_table <- data.frame(
    Source  = c(
      "Between-group",
      "Within-group",
      "PERMANOVA F",
      paste0("Within-group (", groups, ")")
    ),
    Delta   = c(delta_SSB, delta_SSW, delta_F,      delta_SSW_g),
    P_value = c(
      pval_emp(delta_SSB,   delta_SSB_perm),
      pval_emp(delta_SSW,   delta_SSW_perm),
      pval_emp(delta_F,     delta_F_perm),
      sapply(seq_along(groups), function(j)
        pval_emp(delta_SSW_g[j], delta_SSW_g_perm[, j])
      )
    ),
    stringsAsFactors = FALSE
  )

  return(result_table)
}

