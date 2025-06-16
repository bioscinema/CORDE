#' Semi-Supervised Ordination (CORDE or ORDE)
#'
#' Performs outcome-guided, confounder-adjusted ordination using sPLS on PCoA embeddings.
#'
#' @param dist A dissimilarity matrix (e.g., from Bray-Curtis or UniFrac).
#' @param Y A vector of outcome labels (numeric or factor).
#' @param X A matrix of confounders (optional, default = NULL).
#' @param ncomp Number of components to extract (default = 2).
#' @param eta Proportion of PCoA variance to retain in embedding (default = 0.99).
#' @param method Either "CORDE" (adjust for confounders) or "ORDE" (no adjustment).
#' @param seed Integer seed for reproducibility.
#'
#' @return A list with:
#' \item{Embeddings}{Embedding matrix after confounder adjustment (or raw PCoA if ORDE).}
#' \item{Coordinates}{sPLS component scores.}
#' \item{Loadings}{sPLS loadings.}
#' \item{expl_var}{Variance explained by each component.}
#' \item{Y}{Outcome vector used in model.}
#'
#' @examples
#' \dontrun{
#' SSordination(dist = bray_mat, Y = factor(group), X = meta[, c("age", "sex")], method = "CORDE", seed = 42)
#' }
#' @export
SSordination <- function(dist, Y, X = NULL, ncomp = 2, eta = 0.99, method = c("CORDE", "ORDE"), seed = 1) {
  requireNamespace("ape", quietly = TRUE)
  requireNamespace("mixOmics", quietly = TRUE)

  method <- match.arg(method)

  # Step 1: PCoA embedding
  pcoa_res <- ape::pcoa(dist)
  eigenvalues <- pcoa_res$values$Eigenvalues
  var_cumsum <- cumsum(eigenvalues) / sum(eigenvalues)
  k_optimal <- which(var_cumsum >= eta)[1]
  Z <- pcoa_res$vectors[, 1:k_optimal]

  # Step 2: Confounder adjustment if CORDE
  if (method == "CORDE") {
    if (is.null(X)) stop("Confounder matrix X must be provided for CORDE method.")
    Hx <- X %*% solve(t(X) %*% X) %*% t(X)  # projection matrix
    E <- (diag(nrow(Hx)) - Hx) %*% Z       # residualized embedding
  } else {
    E <- Z  # ORDE uses raw embedding
  }

  # Step 3: Run sPLS-DA
  set.seed(seed)
  tune_res <- mixOmics::tune.splsda(E, Y, ncomp = ncomp, validation = "Mfold", folds = 5)
  spls_model <- mixOmics::splsda(E, Y, ncomp = ncomp, keepX = tune_res$choice.keepX)

  # Step 4: Extract and sort components by explained variance
  expl_var <- mixOmics::explained_variance(E, spls_model$variates$X, ncomp)
  ranked_idx <- order(expl_var, decreasing = TRUE)

  return(list(
    Embeddings = E,
    Coordinates = spls_model$variates$X[, ranked_idx, drop = FALSE],
    Loadings = spls_model$loadings$X[, ranked_idx, drop = FALSE],
    expl_var = expl_var[ranked_idx],
    Y = Y
  ))
}
