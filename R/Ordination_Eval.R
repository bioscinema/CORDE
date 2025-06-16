#' Evaluate Ordination Quality Using Multiple Metrics
#'
#' This function evaluates the quality of an ordination (e.g., from PCA, t-SNE, CORDE)
#' based on trustworthiness, continuity, clustering accuracy, silhouette width, and
#' confounder association.
#'
#' @param D A dissimilarity matrix (e.g., Bray-Curtis, UniFrac).
#' @param Coordinates A numeric matrix or data frame containing the ordination coordinates.
#' @param expl_var A numeric vector of explained variance for the axes (e.g., from PCA).
#' @param Y A factor or vector indicating the true group labels for each sample (e.g., phenotype).
#' @param X A matrix or data frame of confounder variables to adjust for (e.g., batch, age).
#' @param ncomp Integer specifying the number of components (default = 2).
#'
#' @return A named list containing the following evaluation metrics:
#' \describe{
#'   \item{ICR}{Information compression rate (sum of explained variance from first `ncomp` components).}
#'   \item{TW}{Trustworthiness (neighborhood preservation from low-dim to high-dim).}
#'   \item{Cont}{Continuity (neighborhood preservation from high-dim to low-dim).}
#'   \item{Rand}{Rand index between clustering from PAM and true labels.}
#'   \item{ARand}{Adjusted Rand index (ARI).}
#'   \item{ASil}{Average silhouette width of PAM clustering in 2D space.}
#'   \item{CARI}{Confounder-associated R² improvement (from McFadden's pseudo-R²).}
#'   \item{rCMI}{Relative conditional mutual information between clustering and outcome, conditioned on confounders.}
#' }
#'
#' @importFrom cluster pam
#' @importFrom mclust adjustedRandIndex
#' @importFrom DescTools PseudoR2
#' @importFrom infotheo mutinformation condinformation
#' @importFrom fossil rand.index
#' @examples
#' \dontrun{
#' # Example usage:
#' ord_res <- Ordination_Eval(D = dist_mat, Coordinates = pcoa_coords,
#'                             expl_var = eigenvals[1:2], Y = group_labels, X = batch_data)
#' print(ord_res)
#' }
#' @export


### function to evaluate ordination methods
Ordination_Eval <- function(
    D, # input dissimilarity matrix
    Coordinates, # input CORDE coordinates
    expl_var, # information explained by coordinates
    Y, # input outcome matrix (numeric)
    X, # input confounder matrix (numeric)
    ncomp = 2) # input number of features to be extracted
{
  require(fossil)
  require(factoextra)
  require(mclust)
  require(cluster)
  require(DescTools)
  require(ggforce)
  require(infotheo)

  ### measures to calculate
  measures <- list(
    ICR = NA, # information compression rate
    TW = NA, # Trustworthiness
    Cont = NA, # Continuity
    Rand = NA, # Rand index
    ARand = NA, # Adjusted Rand index
    ASil = NA, # Average silhouette score
    CARI = NA, # Confounder-Associated R-square Improvement
    rCMI = NA # relative conditional mutual information
  )

  # ICR
  measures$ICR = expl_var[1] + expl_var[2]

  # Trustworthiness and Continuity
  D_low <- as.matrix(dist(Coordinates[,1:2]))
  low_rank <- apply(D_low, 1, order)[-1, ]
  high_rank <- apply(D, 1, order)[-1, ]

  n <- nrow(Coordinates)
  if(n<=50){k = 5}
  if(n>50 & n<= 500){k=ceiling(n/10)}
  if(n>500){k=50}
  trust <- 0
  cont <- 0
  for (i in 1:n) {
    # Trustworthiness: neighbors in low-dim also in high-dim, 1 - type 1 error
    u <- low_rank[1:k, i]
    ranks_high <- match(u, high_rank[, i])
    partial_sum <- sum((ranks_high[ranks_high > k] - k), na.rm = TRUE)
    trust <- trust + partial_sum

    # Continuity: neighbors in high-dim also in low-dim, power
    v <- high_rank[1:k, i]
    ranks_low <- match(v, low_rank[, i])
    cont <- cont + sum((ranks_low[ranks_low > k] - k), na.rm = TRUE)
  }

  norm <- n * k * (2 * n - 3 * k - 1) / 2
  trustworthiness <- 1 - (2 / norm) * trust
  continuity <- 1 - (2 / norm) * cont

  measures$TW <- trustworthiness
  measures$Cont <- continuity

  # 2-medoids and evaluation
  kmed <- pam(D_low, k = nlevels(Y))
  ri <- rand.index(kmed$clustering, as.integer(Y))
  ari <- adjustedRandIndex(kmed$clustering, Y)

  measures$Rand <- ri
  measures$ARand <- ari
  measures$ASil <- kmed$silinfo$avg.width

  # adjusted clustering evaluation
  C <- kmed$clustering
  model1 <- glm(Y ~ C, family = binomial)
  model2 <- glm(Y ~ C + X, family = binomial)

  r2_yc <- PseudoR2(model1, which = "McFadden")
  r2_ycz <- PseudoR2(model2, which = "McFadden")

  CARI <- (r2_ycz - r2_yc)/r2_yc
  measures$CARI <- as.numeric(CARI)

  # rCMI
  C <- as.factor(C)
  Y <- as.factor(Y)
  X <- as.factor(X)

  mi_cy <- mutinformation(C, Y)
  mi_cy_x <- condinformation(C, Y, X)
  rCMI <- mi_cy_x / mi_cy

  measures$rCMI <- rCMI

  return(measures)
}
