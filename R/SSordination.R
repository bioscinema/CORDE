#' Semi Supervised Ordination with Confounder Adjustment
#'
#' Performs PCoA followed by sPLS-DA, with optional confounder adjustment (CORDE).
#'
#' @param dist A distance matrix (class 'dist') for PCoA.
#' @param Y A factor vector of outcome labels.
#' @param X Optional confounder matrix. Required for method = "CORDE".
#' @param ncomp Number of components to extract via sPLS-DA.
#' @param eta Proportion of variance to retain from PCoA (default: 0.99).
#' @param method Either "CORDE" (confounder-adjusted) or "ORDE" (raw).
#' @param seed Random seed for reproducibility (default: 1).
#'
#' @return A list containing:
#'   \item{Embeddings}{PCoA embeddings (residualized if CORDE)}
#'   \item{Coordinates}{Coordinates from sPLS-DA (PC1, PC2, ...)}
#'   \item{Loadings}{Taxa/component loadings from sPLS-DA}
#'   \item{expl_var}{Explained variance for selected components}
#'   \item{Y}{Input outcome vector}
#'
#' @importFrom ape pcoa
#' @importFrom mixOmics tune.splsda splsda explained_variance
#' @export
SSordination <- function(dist, Y, X = NULL, ncomp = 2, eta = 0.99, method = c("CORDE", "ORDE"), seed = 1) {
  method <- match.arg(method)

  # Step 1: PCoA
  pcoa_res <- ape::pcoa(dist)
  eigenvalues <- pcoa_res$values$Eigenvalues
  var_cumsum <- cumsum(eigenvalues) / sum(eigenvalues)
  k_optimal <- which(var_cumsum >= eta)[1]
  Z <- pcoa_res$vectors[, 1:k_optimal]

  # Step 2: Residualize if CORDE
  if (method == "CORDE") {
    if (is.null(X)) stop("Confounder matrix X must be provided for CORDE method.")
    Hx <- X %*% solve(t(X) %*% X) %*% t(X)
    E <- (diag(nrow(Hx)) - Hx) %*% Z
  } else {
    E <- Z
  }

  # Step 3: sPLS-DA
  set.seed(seed)
  tune_res <- mixOmics::tune.splsda(E, Y, ncomp = ncomp, validation = "Mfold", folds = 5)
  spls_model <- mixOmics::splsda(E, Y, ncomp = ncomp, keepX = tune_res$choice.keepX)

  # Step 4: Variance and return
  expl_var <- mixOmics::explained_variance(E, spls_model$variates$X, ncomp)
  ranked_idx <- order(expl_var, decreasing = TRUE)

  return(list(
    Embeddings  = E,
    Coordinates = spls_model$variates$X[, ranked_idx, drop = FALSE],
    Loadings    = spls_model$loadings$X[, ranked_idx, drop = FALSE],
    expl_var    = expl_var[ranked_idx],
    Y           = Y
  ))
}
