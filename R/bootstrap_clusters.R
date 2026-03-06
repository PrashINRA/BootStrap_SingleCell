#' Bootstrap Cluster Stability Analysis
#'
#' Evaluates cluster stability in a Seurat object by bootstrap resampling
#' and computing co-assignment probabilities. In each iteration, cells are
#' resampled with replacement and reclustered independently. Co-assignment
#' is computed per iteration using a label-agnostic pairwise approach
#' (Monti et al. 2003): for each pair of original clusters, we measure the
#' fraction of cell pairs that land in the same recluster, then average
#' across iterations.
#'
#' @param seurat_obj A processed Seurat object with dimensionality reduction
#'   already computed (e.g., PCA, NMF).
#' @param clusters A character or factor vector of cluster labels for each cell.
#'   Defaults to \code{Idents(seurat_obj)} if not provided.
#' @param reduction Character string specifying the dimensionality reduction to
#'   use for neighbor finding (e.g., \code{"pca"}, \code{"nmf"}, \code{"mnn"}).
#'   Default is \code{"pca"}.
#' @param dims Integer vector of dimensions to use from the reduction. Default
#'   is \code{1:30}.
#' @param resolution Numeric value for the clustering resolution parameter
#'   passed to \code{Seurat::FindClusters}. Default is \code{0.3}.
#' @param algorithm Integer specifying the clustering algorithm.
#'   1 = Louvain (default), 2 = Louvain with multilevel refinement,
#'   3 = SLM, 4 = Leiden (requires leidenalg Python module).
#' @param n.cells Integer specifying the number of cells to sample per
#'   iteration. Default is the total number of cells in the object.
#' @param iterations Integer specifying the number of bootstrap iterations.
#'   Default is \code{50}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'   Default is \code{TRUE}.
#'
#' @return A numeric matrix of co-assignment probabilities. Rows and columns
#'   correspond to original cluster identities. Values range from 0 to 1.
#'   Diagonal values represent within-cluster stability (fraction of same-cluster
#'   cell pairs that co-cluster). Off-diagonal values represent cross-cluster
#'   mixing.
#'
#' @details
#' **Bootstrap resampling:** Each iteration samples \code{n.cells} cells with
#' replacement from the Seurat object. With full-size resampling
#' (\code{n.cells = ncol(seurat_obj)}), approximately 63.2\% of unique cells
#' appear per resample.
#'
#' **Pairwise co-assignment (Monti et al. 2003):** For each bootstrap iteration,
#' a contingency table is built between original cluster labels and new recluster
#' assignments. For a given pair of original clusters (A, B):
#' \itemize{
#'   \item **Self-stability (A == B):** The fraction of within-cluster cell pairs
#'     that land in the same recluster, computed as
#'     \code{sum(choose(n_per_recluster, 2)) / choose(n_total, 2)}.
#'   \item **Cross-assignment (A != B):** The fraction of cross-cluster cell pairs
#'     that land in the same recluster, computed as
#'     \code{sum(n_A_per_recluster * n_B_per_recluster) / (n_A * n_B)}.
#' }
#' This approach is label-agnostic: it does not depend on recluster label
#' identity across iterations, avoiding the label permutation problem that
#' affects spread-vector methods.
#'
#' **Missing clusters:** Small clusters may be absent from some bootstrap
#' samples. The function tracks valid iterations per cluster pair and averages
#' only over iterations where both clusters were represented.
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J. & Golub, T. (2003). Consensus Clustering:
#' A Resampling-Based Method for Class Discovery and Visualization of Gene
#' Expression Microarray Data. \emph{Machine Learning}, 52, 91-118.
#'
#' @examples
#' \dontrun{
#' # Basic usage with PCA
#' coassign <- bootstrap_clusters(seurat_obj, reduction = "pca",
#'                                dims = 1:20, resolution = 0.3)
#'
#' # Using NMF reduction with custom clusters
#' coassign <- bootstrap_clusters(seurat_obj,
#'                                clusters = seurat_obj$cell_type,
#'                                reduction = "nmf", dims = 1:15,
#'                                resolution = 0.5, iterations = 100)
#'
#' # Quick test with fewer iterations
#' coassign <- bootstrap_clusters(seurat_obj, iterations = 10)
#' }
#'
#' @export
bootstrap_clusters <- function(seurat_obj,
                               clusters = NULL,
                               reduction = "pca",
                               dims = 1:30,
                               resolution = 0.3,
                               algorithm = 1,
                               n.cells = NULL,
                               iterations = 50,
                               verbose = TRUE) {
  
  # --- Input validation ---
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.")
  }
  
  # Check that the specified reduction exists
  available_reductions <- Seurat::Reductions(seurat_obj)
  if (!reduction %in% available_reductions) {
    stop(paste0("Reduction '", reduction, "' not found in the Seurat object.\n",
                "Available reductions: ", paste(available_reductions, collapse = ", ")))
  }
  
  # Check dims do not exceed available dimensions
  max_dim <- ncol(Seurat::Embeddings(seurat_obj, reduction = reduction))
  if (max(dims) > max_dim) {
    stop(paste0("Requested dims go up to ", max(dims),
                " but reduction '", reduction, "' only has ", max_dim, " dimensions."))
  }
  
  # Handle clusters
  if (is.null(clusters)) {
    clusters <- Seurat::Idents(seurat_obj)
    if (verbose) message("Using active Idents as cluster labels.")
  }
  clusters <- as.character(clusters)
  
  if (length(clusters) != ncol(seurat_obj)) {
    stop("Length of 'clusters' must match the number of cells in the Seurat object.")
  }
  
  # Name clusters by cell barcodes
  names(clusters) <- colnames(seurat_obj)
  cluster.ids <- sort(unique(clusters))
  n.clust <- length(cluster.ids)
  
  if (n.clust < 2) {
    stop("Need at least 2 clusters to compute co-assignment probabilities.")
  }
  
  if (is.null(n.cells)) {
    n.cells <- ncol(seurat_obj)
  }
  
  if (iterations < 1) {
    stop("'iterations' must be at least 1.")
  }
  
  if (verbose) {
    message(paste0("BootStrapSC: ", n.clust, " clusters, ", ncol(seurat_obj),
                   " cells, ", iterations, " iterations"))
    message(paste0("Reduction: ", reduction, ", Dims: 1:", max(dims),
                   ", Resolution: ", resolution))
  }
  
  # --- Initialize output matrices ---
  coassign_sum <- matrix(0, n.clust, n.clust,
                         dimnames = list(cluster.ids, cluster.ids))
  valid_count  <- matrix(0, n.clust, n.clust,
                         dimnames = list(cluster.ids, cluster.ids))
  
  # --- Bootstrap loop ---
  for (i in seq_len(iterations)) {
    if (verbose) message(paste0("  Iteration ", i, "/", iterations))
    
    # Sample with replacement
    chosen <- sample(ncol(seurat_obj), n.cells, replace = TRUE)
    unique_idx <- unique(chosen)
    
    # Subset and recluster only unique cells
    sub <- seurat_obj[, unique_idx]
    sub <- Seurat::FindNeighbors(sub, reduction = reduction,
                                 dims = dims, verbose = FALSE)
    sub <- Seurat::FindClusters(sub, resolution = resolution,
                                algorithm = algorithm, verbose = FALSE)
    
    recluster <- as.character(sub$seurat_clusters)
    orig <- clusters[colnames(sub)]
    
    # Contingency table: original_cluster x recluster
    cont <- table(orig, recluster)
    present <- intersect(cluster.ids, rownames(cont))
    
    # Compute pairwise co-assignment for this iteration
    for (a in present) {
      row_a <- cont[a, ]
      n_a <- sum(row_a)
      
      for (b in present) {
        row_b <- cont[b, ]
        n_b <- sum(row_b)
        
        if (a == b) {
          # Self: fraction of within-cluster pairs in same recluster
          if (n_a < 2) next
          val <- sum(choose(row_a, 2)) / choose(n_a, 2)
        } else {
          # Cross: fraction of cross-pairs in same recluster
          val <- sum(row_a * row_b) / (n_a * n_b)
        }
        
        coassign_sum[a, b] <- coassign_sum[a, b] + val
        valid_count[a, b]  <- valid_count[a, b] + 1
      }
    }
  }
  
  # --- Compute final probabilities ---
  valid_count[valid_count == 0] <- NA
  output <- coassign_sum / valid_count
  
  if (verbose) message("Done.")
  return(output)
}