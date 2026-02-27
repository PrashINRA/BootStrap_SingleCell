#' Plot Co-assignment Probability Heatmap
#'
#' Visualizes the co-assignment probability matrix as a heatmap using
#' \code{pheatmap}. Diagonal values represent cluster self-stability.
#' Off-diagonal values indicate the degree of overlap between cluster pairs.
#'
#' @param coassign_matrix A numeric matrix of co-assignment probabilities,
#'   as returned by \code{\link{bootstrap_clusters}}.
#' @param cluster_rows Logical. Whether to hierarchically cluster the rows.
#'   Default is \code{FALSE} to preserve original cluster ordering.
#' @param cluster_cols Logical. Whether to hierarchically cluster the columns.
#'   Default is \code{FALSE}.
#' @param title Character string for the heatmap title. Default is
#'   \code{"Co-assignment Probabilities"}.
#' @param color Color palette vector. Default uses \code{viridis::magma(100)}
#'   reversed.
#' @param display_numbers Logical. If \code{TRUE}, displays probability values
#'   in each cell. Default is \code{FALSE}.
#' @param number_format Character string for \code{sprintf} formatting of
#'   displayed numbers. Default is \code{"\%.2f"}.
#' @param ... Additional arguments passed to \code{pheatmap::pheatmap}.
#'
#' @return Invisibly returns the pheatmap object.
#'
#' @examples
#' \dontrun{
#' coassign <- bootstrap_clusters(seurat_obj, reduction = "pca")
#' plot_coassignment(coassign)
#'
#' # With numbers displayed and custom title
#' plot_coassignment(coassign, display_numbers = TRUE,
#'                   title = "My Dataset - Cluster Stability")
#'
#' # With hierarchical clustering of rows/columns
#' plot_coassignment(coassign, cluster_rows = TRUE, cluster_cols = TRUE)
#' }
#'
#' @export
plot_coassignment <- function(coassign_matrix,
                               cluster_rows = FALSE,
                               cluster_cols = FALSE,
                               title = "Co-assignment Probabilities",
                               color = rev(viridis::magma(100)),
                               display_numbers = FALSE,
                               number_format = "%.2f",
                               ...) {

  if (!is.matrix(coassign_matrix) || !is.numeric(coassign_matrix)) {
    stop("'coassign_matrix' must be a numeric matrix.")
  }

  p <- pheatmap::pheatmap(coassign_matrix,
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           main = title,
                           angle_col = 45,
                           color = color,
                           display_numbers = display_numbers,
                           number_format = number_format,
                           na_col = "white",
                           ...)
  invisible(p)
}
