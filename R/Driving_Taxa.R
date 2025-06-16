#' Evaluate the Effect of Removing a Taxon on Dissimilarity Structure
#'
#' This function quantifies the impact of removing one or more taxa from an OTU table
#' on community structure, using dissimilarity matrices, silhouette scores, Mantel tests,
#' and Procrustes analysis.
#'
#' @param otu_df A numeric matrix or data frame of OTU abundances (rows = taxa, columns = samples).
#' @param remove_otu A character vector of taxon (row) names to be removed from \code{otu_df}.
#' @param group A factor or character vector indicating group membership for each sample (length must match number of columns in \code{otu_df}).
#' @param diss Character string specifying the dissimilarity metric used by \code{vegdist} (e.g., \code{"bray"}, \code{"euclidean"}).
#'
#' @return A data frame summarizing the effect of OTU removal, with the following columns:
#' \describe{
#'   \item{comparison}{Comparison level: overall ("between group") and per group.}
#'   \item{mantel_p}{Mantel test p-value comparing dissimilarities before and after OTU removal.}
#'   \item{proc_p}{Procrustes test p-value comparing ordinations before and after OTU removal.}
#'   \item{silhouette_difference}{Change in silhouette width (overall or medoid distance).}
#'   \item{silhouette_difference_pvalue}{Wilcoxon test p-value comparing silhouette widths or medoid distances.}
#' }
#'
#' @details
#' This function performs the following:
#' \itemize{
#'   \item Computes baseline dissimilarity and silhouette scores.
#'   \item Identifies group-specific medoids using PAM clustering.
#'   \item Recomputes dissimilarity matrix after removing specified taxa.
#'   \item Assesses structural shifts via Mantel and Procrustes tests at group and global levels.
#'   \item Compares silhouette and medoid distances using paired Wilcoxon tests.
#' }
#'
#' @import vegan cluster
#' @seealso \code{\link[vegan]{mantel}}, \code{\link[vegan]{protest}}, \code{\link[cluster]{silhouette}}, \code{\link[cluster]{pam}}
#'
#' @examples
#' \dontrun{
#' otu_mat <- matrix(rpois(100, 5), nrow = 10)
#' rownames(otu_mat) <- paste0("OTU", 1:10)
#' colnames(otu_mat) <- paste0("Sample", 1:10)
#' group <- rep(c("A", "B"), each = 5)
#' evaluate_taxon_removal_effect(otu_mat, remove_otu = "OTU5", group = group, diss = "bray")
#' }
#'
#' @export
evaluate_taxon_removal_effect <- function(otu_df, remove_otu, group, diss) {

  # Step 1: Compute baseline dissimilarity matrix
  n <- ncol(otu_df)
  dist1 <- vegdist(t(otu_df), method = diss)  # transpose so samples are rows
  sil1 <- silhouette(as.numeric(as.factor(group)), dist1)
  sil_score1 <- mean(sil1[, "sil_width"])  # average silhouette width

  # Step 1.5: Compute within-group dissimilarity to medoids
  group_levels <- unique(group)
  group1 <- group_levels[1]
  group2 <- group_levels[2]

  dmat1 <- as.matrix(dist1)

  group1_samples <- which(group == group1)
  pam1 <- pam(dmat1[group1_samples, group1_samples], k = 1)
  medoid1_index <- group1_samples[pam1$id.med]
  avg_d1 <- mean(dmat1[group1_samples, medoid1_index])

  group2_samples <- which(group == group2)
  pam2 <- pam(dmat1[group2_samples, group2_samples], k = 1)
  medoid2_index <- group2_samples[pam2$id.med]
  avg_d2 <- mean(dmat1[group2_samples, medoid2_index])

  # Step 2: Remove specified OTU(s) and recompute dissimilarity
  otu_df_new <- otu_df[!(rownames(otu_df) %in% remove_otu), ]
  dist2 <- vegdist(t(otu_df_new), method = diss)

  # Step 3: Matrix-level comparison using Mantel and PERMANOVA
  mantel_p <- rep(0,3)
  proc_p <- rep(0,3)

  # between group
  mantel_p[1] <- mantel(dist1, dist2)$signif
  pcoa1 <- cmdscale(dist1, k = floor(0.5*n))
  pcoa2 <- cmdscale(dist2, k = floor(0.5*n))
  proc <- protest(pcoa1, pcoa2, permutations = 999)
  proc_p[1] <- proc$signif

  # group 1
  mantel_p[2] <- mantel(as.dist(as.matrix(dist1)[group1_samples,group1_samples]),
                        as.dist(as.matrix(dist2)[group1_samples,group1_samples]))$signif
  pcoa1 <- cmdscale(as.dist(as.matrix(dist1)[group1_samples,group1_samples]), k = floor(0.5*length(group1_samples)))
  pcoa2 <- cmdscale(as.dist(as.matrix(dist2)[group1_samples,group1_samples]), k = floor(0.5*length(group1_samples)))
  proc <- protest(pcoa1, pcoa2, permutations = 999)
  proc_p[2] <- proc$signif

  # group 2
  mantel_p[3] <- mantel(as.dist(as.matrix(dist1)[group2_samples,group2_samples]),
                        as.dist(as.matrix(dist2)[group2_samples,group2_samples]))$signif
  pcoa1 <- cmdscale(as.dist(as.matrix(dist1)[group2_samples,group2_samples]), k = floor(0.5*length(group2_samples)))
  pcoa2 <- cmdscale(as.dist(as.matrix(dist2)[group2_samples,group2_samples]), k = floor(0.5*length(group2_samples)))
  proc <- protest(pcoa1, pcoa2, permutations = 999)
  proc_p[3] <- proc$signif

  # Step 4: Recompute silhouette and medoid distances
  sil2 <- silhouette(as.numeric(as.factor(group)), dist2)
  sil_score2 <- mean(sil2[, "sil_width"])

  dmat2 <- as.matrix(dist2)

  pam1 <- pam(dmat2[group1_samples, group1_samples], k = 1)
  medoid1_new_index <- group1_samples[pam1$id.med]
  avg_d1_new <- mean(dmat2[group1_samples, medoid1_new_index])

  pam2 <- pam(dmat2[group2_samples, group2_samples], k = 1)
  medoid2_new_index <- group2_samples[pam2$id.med]
  avg_d2_new <- mean(dmat2[group2_samples, medoid2_new_index])

  # Step 5: Compute differences before and after OTU removal
  sil_diff <- sil_score1 - sil_score2
  d1_diff <- avg_d1 - avg_d1_new
  d2_diff <- avg_d2 - avg_d2_new

  sil_diff_wil <- wilcox.test(sil1[, "sil_width"],sil2[, "sil_width"])$p.value
  d1_diff_wil <- wilcox.test(dmat1[group1_samples, medoid1_index],dmat2[group1_samples, medoid1_new_index])$p.value
  d2_diff_wil <- wilcox.test(dmat1[group2_samples, medoid2_index],dmat2[group2_samples, medoid2_new_index])$p.value

  # Step 6: Return results as a named list
  return(data.frame(
    comparison = c("between group",as.character(group1),as.character(group2)),
    mantel_p = mantel_p, # correlation before and after removing
    proc_p = proc_p, # ordination structure change
    silhouette_difference = c(sil_diff,d1_diff,d2_diff),
    silhouette_difference_pvalue = c(sil_diff_wil,d1_diff_wil,d2_diff_wil)
  ))
}
