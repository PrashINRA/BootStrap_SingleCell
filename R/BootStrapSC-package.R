#' @keywords internal
"_PACKAGE"

#' BootStrapSC: Bootstrap Cluster Stability for Single-Cell Data
#'
#' Evaluates single-cell cluster stability using bootstrap resampling and
#' co-assignment probability analysis. Cells are repeatedly resampled with
#' replacement and reclustered to measure how consistently they reassemble
#' into the same groups. Works directly with Seurat objects and supports
#' multiple dimensionality reduction methods.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{bootstrap_clusters}}: Run bootstrap stability analysis
#'   \item \code{\link{plot_coassignment}}: Visualize co-assignment heatmap
#' }
#'
#' @docType _PACKAGE
#' @name BootStrapSC-package
NULL
