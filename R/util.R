#' Compute Sum-of-Distances Test (SoDT)
#'
#' Performs the SoDT (Sum-of-Distances Test), a dissimilarity-based ANOVA-like test
#' for group differences, with optional permutation testing for significance.
#'
#' @param D A symmetric dissimilarity matrix (e.g., Bray-Curtis, UniFrac).
#' @param Y A factor or vector indicating group membership.
#' @param nperm Number of permutations for significance testing (default = 999).
#' @param seed Random seed for reproducibility (default = 2025).
#'
#' @return A list containing:
#' \describe{
#'   \item{table}{A data frame with source of variation, sum of squares, degrees of freedom, and notes.}
#'   \item{F_stat}{The SoDT pseudo-F statistic.}
#'   \item{p_value}{The permutation p-value.}
#' }
#'
#' @importFrom stats glm binomial model.matrix var cmdscale as.dist wilcox.test
#' @importFrom utils combn
#' @importFrom grDevices chull
#' @examples
#' \dontrun{
#' dist_mat <- as.matrix(vegan::vegdist(iris[,1:4]))
#' compute_sodt(dist_mat, iris$Species)
#' }
#' @export
compute_sodt <- function(D, Y, nperm = 999, seed = 2025) {
  set.seed(seed)

  # Input checks
  if (!is.matrix(D) || !all(D == t(D))) stop("D must be a symmetric matrix.")
  if (any(diag(D) != 0)) stop("D must have 0 diagonal.")
  if (length(Y) != nrow(D)) stop("Y length must match number of rows in D.")

  n <- nrow(D)
  groups <- unique(Y)
  d <- length(groups)
  group_sizes <- table(Y)

  # Step 1: Compute centered Gram matrix G
  D2 <- D^2
  J <- diag(n) - matrix(1, n, n) / n
  G <- -0.5 * J %*% D2 %*% J

  # Step 2: Hat matrix from design matrix
  Y_mat <- model.matrix(~ 0 + factor(Y))
  H_Y <- Y_mat %*% solve(t(Y_mat) %*% Y_mat) %*% t(Y_mat)
  I_n <- diag(n)

  # Step 3: Compute SoDT statistics
  SST <- sum(diag(G))
  SSB <- sum(diag(H_Y %*% G %*% H_Y))
  SSW <- sum(diag((I_n - H_Y) %*% G %*% (I_n - H_Y)))
  F_stat <- (SSB / (d - 1)) / (SSW / (n - d))

  # Step 4: Compute group-specific SSW
  SSW_g <- sapply(groups, function(g) {
    idx <- which(Y == g)
    P_g <- matrix(0, nrow = length(idx), ncol = n)
    P_g[cbind(1:length(idx), idx)] <- 1
    sum(diag(P_g %*% (I_n - H_Y) %*% G %*% (I_n - H_Y) %*% t(P_g)))
  })

  # Step 5: Permutation test
  F_perm <- replicate(nperm, {
    Y_perm <- sample(Y)
    Y_perm_mat <- model.matrix(~ 0 + factor(Y_perm))
    H_Y_perm <- Y_perm_mat %*% solve(t(Y_perm_mat) %*% Y_perm_mat) %*% t(Y_perm_mat)
    SSB_perm <- sum(diag(H_Y_perm %*% G %*% H_Y_perm))
    SSW_perm <- sum(diag((I_n - H_Y_perm) %*% G %*% (I_n - H_Y_perm)))
    (SSB_perm / (d - 1)) / (SSW_perm / (n - d))
  })
  p_val <- (1 + sum(F_perm >= F_stat)) / (nperm + 1)

  # Step 6: Format output table
  df_table <- data.frame(
    Source = c("Total", "Between-group", "Within-group", paste0("Within-group (", groups, ")")),
    SumOfSquares = c(SST, SSB, SSW, SSW_g),
    df = c(n - 1, d - 1, n - d, as.numeric(group_sizes) - 1),
    Notes = c(
      "Total dissimilarity",
      "Group separation",
      "Residual dispersion",
      paste0("Group-specific dispersion: ", groups)
    ),
    stringsAsFactors = FALSE
  )

  return(list(
    table = df_table,
    F_stat = F_stat,
    p_value = p_val
  ))
}


#' Fast Sum-of-Distances Test (SoDT) for Group Differences
#'
#' Performs an efficient version of the SoDT to assess whether groups defined by `Y`
#' differ in their distribution in a dissimilarity matrix `D`.
#' This implementation is optimized for speed and includes permutation-based significance testing.
#'
#' @param D A symmetric dissimilarity matrix with zeros on the diagonal.
#' @param Y A grouping factor or vector of class labels.
#' @param nperm Number of permutations for empirical p-value computation (default = 999).
#' @param seed Integer seed for reproducibility (default = 2025).
#'
#' @return A list containing:
#' \describe{
#'   \item{table}{A data frame summarizing total, between-group, and within-group sums of squares.}
#'   \item{F_stat}{The observed SoDT F-statistic.}
#'   \item{p_value}{Permutation-based p-value for the F-statistic.}
#' }
#'
#' @examples
#' \dontrun{
#' dist_mat <- as.matrix(vegan::vegdist(iris[1:50, 1:4]))
#' groups <- rep(c("A", "B"), each = 25)
#' result <- compute_sodt_fast(dist_mat, groups)
#' print(result$table)
#' }
#'
#' @export
compute_sodt_fast <- function(D, Y, nperm = 999, seed = 2025) {
  set.seed(seed)
  if (!is.matrix(D) || !all(D == t(D))) stop("D must be symmetric.")
  if (any(diag(D) != 0)) stop("D must have 0 diagonal.")

  n <- nrow(D)
  groups <- unique(Y)
  d <- length(groups)
  group_sizes <- table(Y)

  # Compute Gram matrix
  D2 <- D^2
  J <- diag(n) - matrix(1, n, n) / n
  G <- -0.5 * J %*% D2 %*% J
  trace_G <- sum(diag(G))

  # Hat matrix
  Y_mat <- model.matrix(~ 0 + factor(Y))
  H_Y <- Y_mat %*% solve(t(Y_mat) %*% Y_mat) %*% t(Y_mat)

  # Compute F-statistic
  SSB <- sum(H_Y * G)
  SSW <- trace_G - SSB
  SST <- trace_G
  F_stat <- (SSB / (d - 1)) / (SSW / (n - d))

  # Compute group-specific residual dispersion
  R <- diag(n) - H_Y
  GR <- R %*% G
  SSW_g <- sapply(groups, function(g) {
    idx <- which(Y == g)
    sum(diag(GR[idx, ] %*% t(R[idx, , drop = FALSE])))
  })

  # Permutation test
  F_perm <- replicate(nperm, {
    Y_perm <- sample(Y)
    Y_mat_perm <- model.matrix(~ 0 + factor(Y_perm))
    H_Y_perm <- Y_mat_perm %*% solve(t(Y_mat_perm) %*% Y_mat_perm) %*% t(Y_mat_perm)
    SSB_perm <- sum(H_Y_perm * G)
    SSW_perm <- trace_G - SSB_perm
    (SSB_perm / (d - 1)) / (SSW_perm / (n - d))
  })

  p_val <- (1 + sum(F_perm >= F_stat)) / (nperm + 1)

  df_table <- data.frame(
    Source = c("Total", "Between-group", "Within-group", paste0("Within-group (", groups, ")")),
    SumOfSquares = c(SST, SSB, SSW, SSW_g),
    df = c(n - 1, d - 1, n - d, as.numeric(group_sizes) - 1),
    Notes = c(
      "Total dissimilarity",
      "Group separation",
      "Residual dispersion",
      paste0("Group-specific dispersion: ", groups)
    ),
    stringsAsFactors = FALSE
  )

  return(list(
    table = df_table,
    F_stat = F_stat,
    p_value = p_val
  ))
}


#' Compute Mantel-based Correlation Matrix Between OTUs and Environmental Axes
#'
#' This function computes a Mantel statistic matrix (`C matrix`) measuring the correlation
#' between dissimilarity profiles of each OTU and each axis of a given environmental or
#' ordination matrix.
#'
#' @param otu_mat A matrix of OTU/ASV abundances (samples × features).
#' @param E A numeric matrix of environmental or ordination axes (samples × axes).
#' @param dist_method Character string indicating the distance metric for OTU dissimilarity
#'   (default is `"bray"` for Bray-Curtis; passed to \code{\link[vegan]{vegdist}}).
#' @param n_perm Integer, number of permutations for the Mantel test (default = 0, i.e., no permutation test).
#'
#' @return A numeric matrix with dimensions (number of OTUs × number of axes) containing
#'   Mantel correlation coefficients (Mantel r) for each OTU-axis pair.
#'
#' @details A small pseudocount (0.5) is added to each OTU profile to avoid zero-distance artifacts.
#'   The function uses \code{\link[vegan]{mantel}} to compute the Mantel statistic.
#'
#' @importFrom vegan vegdist mantel
#' @importFrom stats dist
#'
#' @examples
#' \dontrun{
#' # Example OTU table (samples × OTUs) and ordination axes
#' otu_mat <- matrix(rpois(100, lambda = 5), nrow = 10)
#' colnames(otu_mat) <- paste0("OTU", 1:10)
#' E <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' colnames(E) <- c("Axis1", "Axis2")
#'
#' C_mat <- get_mantel_C_matrix(otu_mat, E, dist_method = "bray", n_perm = 99)
#' print(C_mat)
#' }
#'
#' @export
get_mantel_C_matrix <- function(otu_mat, E, dist_method = "bray", n_perm = 0) {
  n_otus <- ncol(otu_mat)
  n_axes <- ncol(E)
  C_mat <- matrix(NA, nrow = n_otus, ncol = n_axes)
  rownames(C_mat) <- colnames(otu_mat)
  colnames(C_mat) <- colnames(E)

  for (j in 1:n_otus) {
    otu_j <- otu_mat[, j, drop = FALSE] + 0.5 # add pseudo counts
    dist_j <- vegdist(otu_j, method = dist_method)

    for (l in 1:n_axes) {
      e_l <- E[, l, drop = FALSE]
      dist_e <- dist(e_l)
      mantel_res <- mantel(dist_j, dist_e, permutations = n_perm)
      C_mat[j, l] <- mantel_res$statistic  # Mantel r
    }
  }
  return(C_mat)
}


#' Compute LEfSe-Style Effect Size Scores for Features
#'
#' This function computes a LEfSe-style score for each feature, based on standardized group
#' differences. For two-class problems, it returns signed log-transformed effect sizes. For
#' multi-class problems, it returns the maximum absolute pairwise standardized difference.
#'
#' @param X A numeric matrix or data frame of features (samples × features), such as OTU abundances or component scores.
#' @param Y A factor vector of class labels corresponding to rows of \code{X}. Must contain at least two levels.
#' @param eps A small number added to denominators to avoid division by zero or log(0). Default is \code{1e-6}.
#'
#' @return A named numeric vector of LEfSe-style scores (length = number of features).
#'   For binary classes, the scores are signed \code{log10}-transformed effect sizes.
#'   For multi-class, the scores are unsigned maximum standardized pairwise differences.
#'
#' @details
#' The function mimics the scoring strategy in the LEfSe tool by computing the difference in means
#' between groups standardized by the pooled standard deviation. For two groups, this is the
#' standardized effect size, log-transformed for interpretability. For multiple groups, the
#' maximum pairwise absolute standardized difference is returned.
#'
#'
#' @export
compute_lefse_style_score <- function(X, Y, eps = 1e-6) {
  # X: feature matrix (samples × features), e.g. sPLS components or OTU scores
  # Y: class labels (factor), with 2 or more levels
  # eps: small number to avoid division by zero and log(0)

  stopifnot(is.factor(Y))
  groups <- levels(Y)
  p <- ncol(X)

  lda_score <- numeric(p)
  names(lda_score) <- colnames(X)

  for (j in 1:p) {
    xj <- X[, j]
    mean_g <- tapply(xj, Y, mean)
    var_g  <- tapply(xj, Y, var)
    n_g    <- table(Y)

    pooled_var <- sum((n_g - 1) * var_g) / sum(n_g - 1)
    pooled_sd <- sqrt(pooled_var)

    if (length(groups) == 2) {
      effect <- (mean_g[2] - mean_g[1]) / (pooled_sd + eps)
      lda_score[j] <- sign(effect) * log10(abs(effect) + eps)
    } else {
      # Optional: still return unsigned max pairwise effect size
      combs <- combn(groups, 2)
      effect <- max(abs(mean_g[combs[1,]] - mean_g[combs[2,]]) / (pooled_sd + eps))
      lda_score[j] <- effect + eps
    }
  }

  return(lda_score)
}


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
    dist_method # dissimilarity measure (the same as CORDE)
){

  W <- CORDE_res$Loadings

  E <- CORDE_res$Embeddings
  C_mat <- get_mantel_C_matrix(t(OTU), E, dist_method, n_perm = 0)

  X <- CORDE_res$Coordinates
  Y <- CORDE_res$Y
  lda_score <- compute_lefse_style_score(X,Y)

  W_norm <- apply(W, 2, function(x) x / sqrt(sum(x^2)))
  C_norm <- t(apply(C_mat, 1, function(x) x / sqrt(sum(x^2))))

  Score <- 100*as.vector(C_norm %*% W_norm %*% lda_score)
  names(Score) <- rownames(OTU)

  return(Score)
}

#' Hadza gut microbiome dataset
#'
#' A sample phyloseq object from the Hadza gut microbiome, used as a demo.
#'
#' @format A `phyloseq` object named \code{physeq}.
#' @source \url{https://www.nature.com/articles/ncomms3654}
#' @name hadza
#' @docType data
"hadza"


