#' Plot Vector Field After Removing a Driving OTU
#'
#' This function visualizes the effect of removing a selected OTU on the low-dimensional ordination
#' space derived from sPLSDA. It overlays vector arrows indicating the change in sample coordinates
#' after OTU removal, helping to interpret the contribution of a taxon to group separation.
#'
#' @param otu_df A matrix or data frame of OTU abundances (samples × taxa).
#' @param otu_to_remove Character vector of one or more OTU (column) names to be removed.
#' @param diss Character string specifying the dissimilarity metric (e.g., \code{"bray"}).
#' @param corde_df A data frame returned from a CORDE-like analysis, containing columns named \code{"SampleID"}, \code{"PC1"}, \code{"PC2"}, \code{"Cluster"}, and \code{"Group"}.
#' @param X A design matrix of confounders (samples × covariates).
#' @param Y A factor vector of group labels.
#' @param alpha Numeric scalar controlling arrow length in the plot (e.g., 0.5).
#'
#' @return A \code{ggplot} object showing the vector field of sample movement after OTU removal.
#' Polygon hulls are plotted around re-labeled clusters, and arrow segments illustrate the shift direction.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Removes the selected OTU(s) from the input.
#'   \item Recomputes a PCoA embedding of the new OTU table.
#'   \item Removes confounding effects using a projection matrix.
#'   \item Applies \code{splsda} to extract new components.
#'   \item Plots arrows from the original component scores to the new ones.
#'   \item Reassigns clusters using the Kuhn-Munkres algorithm for interpretability.
#' }
#'
#' @import vegan
#' @import ape
#' @import mixOmics
#' @import clue
#' @import dplyr
#' @import ggplot2
#' @importFrom grid unit
#'
#' @examples
#' \dontrun{
#' plot_driving_field(
#'   otu_df = otu_matrix,
#'   otu_to_remove = "OTU_42",
#'   diss = "bray",
#'   corde_df = corde_output,
#'   X = confounder_matrix,
#'   Y = sample_labels,
#'   alpha = 0.7
#' )
#' }
#'
#' @export
plot_driving_field <- function(otu_df, otu_to_remove, diss, corde_df, X, Y, alpha) {
  cat("\033[32mRemoving OTU:", otu_to_remove, "\033[0m\n")

  ## Step 1: Remove selected OTU
  otu_df_removed <- otu_df[, !(colnames(otu_df) %in% otu_to_remove)]

  ## Step 2: Compute distance & PCoA
  dist_removed <- vegdist(otu_df_removed, method=diss)
  pcoa_removed <- ape::pcoa(dist_removed)
  Z_removed <- pcoa_removed$vectors[, 1:ncol(corde_df[, c("PC1", "PC2")])]

  ## Step 3: Remove confounder
  HatX <- X %*% solve(t(X) %*% X) %*% t(X)
  E_removed <- (diag(nrow(HatX)) - HatX) %*% Z_removed

  ## Step 4: Project using sPLS with same parameters
  splsda_removed <- splsda(E_removed, Y, ncomp = 2)
  P_removed <- splsda_removed$variates$X

  ## Step 5: Compute delta vectors
  delta_mat <- P_removed - as.matrix(corde_df[, c("PC1", "PC2")])
  delta_norm <- sqrt(rowSums(delta_mat^2))
  unit_vec <- delta_mat / delta_norm

  ## Step 6: Assemble dataframe for plotting
  vec_df <- as.data.frame(corde_df) %>%
    dplyr::select(SampleID, PC1, PC2, Cluster, Group) %>%
    mutate(dx = delta_mat[,1],
           dy = delta_mat[,2],
           dx_unit = unit_vec[,1],
           dy_unit = unit_vec[,2],
           vector_length = delta_norm)

  # Step 6.1: Relabel the clusters using Kuhn-Munkres algorithm
  true_labels <- vec_df$Group
  cluster_labels <- vec_df$Cluster
  perm <- clue::solve_LSAP(table(cluster_labels, true_labels), maximum = TRUE)
  vec_df$Cluster <- factor(perm[cluster_labels],
                                      levels = 1:length(levels(true_labels)),
                                      labels = levels(true_labels))

  ## Step 7: Plot
  hull_points <- vec_df %>%
    group_by(Cluster) %>%
    slice(chull(PC1, PC2))

  p <- ggplot(vec_df, aes(x = PC1, y = PC2)) +
    geom_polygon(data = hull_points, aes(fill = Cluster, group = Cluster), alpha = 0.1, color = "gray70") +
    geom_segment(
      aes(xend = PC1 - alpha*dx_unit, yend = PC2 - alpha*dy_unit),
      arrow = arrow(length = unit(0.15, "inches")),
      color = "black", linewidth = 0.3, alpha = 0.7
    ) +
    geom_point(aes(color = Group), size = 3, alpha = 0.9) +
    theme_minimal() +
    coord_fixed() +
    labs(
      x = "PCo 1",
      y = "PCo 2"
    )
  return(p)
}
