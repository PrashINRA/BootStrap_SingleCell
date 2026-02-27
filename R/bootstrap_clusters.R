#' Bootstrap Cluster Stability Analysis
#'
#' Evaluates cluster stability in a Seurat object by bootstrap resampling
#' and computing co-assignment probabilities. In each iteration, cells are
#' resampled with replacement and reclustered independently. The co-assignment
#' probability between two clusters is the dot product of their normalized
#' spread vectors across bootstrap reclusters, averaged over all iterations.
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
#'   The lower triangle is set to \code{NA} (matrix is symmetric).
#'
#' @details
#' **Bootstrap resampling:** Each iteration samples \code{n.cells} cells with
#' replacement from the Seurat object. With full-size resampling
#' (\code{n.cells = ncol(seurat_obj)}), approximately 63.2\% of unique cells
#' appear per resample.
#'
#' **Handling Seurat duplicates:** Since Seurat objects cannot hold cells with
#' identical barcodes, only unique cells are subset for reclustering. Results
#' are mapped back to the full bootstrap sample via index matching, preserving
#' correct co-assignment frequencies.
#'
#' **Missing clusters:** Small clusters may be absent from some bootstrap
#' samples. The function tracks valid iterations per cluster pair and averages
#' only over iterations where both clusters were represented.
#'
#' **Probability bounds:** Co-assignment scores are dot products of vectors
#' that each sum to 1, guaranteeing values in \[0, 1\] by the Cauchy-Schwarz
#' inequality.
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
  output <- matrix(0, n.clust, n.clust)
  dimnames(output) <- list(cluster.ids, cluster.ids)
  count_matrix <- matrix(0, n.clust, n.clust)

  # --- Bootstrap loop ---
  for (i in seq_len(iterations)) {
    if (verbose) message(paste0("  Iteration ", i, "/", iterations))

    # Sample with replacement
    chosen <- sample(ncol(seurat_obj), n.cells, replace = TRUE)
    unique_chosen <- unique(chosen)

    # Subset and recluster only unique cells
    resampled <- seurat_obj[, unique_chosen]
    reclusters_unique <- .recluster(resampled, reduction, dims, resolution, algorithm)

    # Map back to full bootstrap sample
    idx_map <- match(chosen, unique_chosen)
    reclusters <- reclusters_unique[idx_map]

    # Contingency table: original labels vs new labels
    tab <- table(clusters[chosen], reclusters)

    # Build normalized spread vectors
    spreads <- list()
    valid <- character(0)

    for (cid in cluster.ids) {
      if (cid %in% rownames(tab)) {
        row_counts <- tab[cid, ]
        row_sum <- sum(row_counts)
        if (row_sum > 0) {
          spreads[[cid]] <- row_counts / row_sum
          valid <- c(valid, cid)
        }
      }
    }

    # Accumulate co-assignment scores for valid pairs
    for (j1 in seq_along(cluster.ids)) {
      cid1 <- cluster.ids[j1]
      if (!(cid1 %in% valid)) next
      for (j2 in seq_len(j1)) {
        cid2 <- cluster.ids[j2]
        if (!(cid2 %in% valid)) next

        all_cols <- union(names(spreads[[cid1]]), names(spreads[[cid2]]))
        s1 <- spreads[[cid1]][all_cols]
        s2 <- spreads[[cid2]][all_cols]
        s1[is.na(s1)] <- 0
        s2[is.na(s2)] <- 0

        output[j2, j1] <- output[j2, j1] + sum(s1 * s2)
        count_matrix[j2, j1] <- count_matrix[j2, j1] + 1
      }
    }
  }

  # --- Compute final probabilities ---
  count_matrix[count_matrix == 0] <- NA
  output <- output / count_matrix

  # Fill symmetric lower triangle for complete matrix storage
  for (j1 in seq_len(n.clust)) {
    for (j2 in seq_len(j1 - 1)) {
      output[j1, j2] <- output[j2, j1]
    }
  }

  if (verbose) message("Done.")
  return(output)
}


#' Internal: Recluster a Seurat Subset
#'
#' @param seurat_sub A subset Seurat object.
#' @param reduction Dimensionality reduction to use.
#' @param dims Dimensions to use.
#' @param resolution Clustering resolution.
#' @param algorithm Clustering algorithm.
#' @return Character vector of new cluster assignments.
#' @keywords internal
#' @noRd
.recluster <- function(seurat_sub, reduction, dims, resolution, algorithm) {
  seurat_sub <- Seurat::FindNeighbors(seurat_sub, reduction = reduction,
                                       dims = dims, verbose = FALSE)
  seurat_sub <- Seurat::FindClusters(seurat_sub, resolution = resolution,
                                      algorithm = algorithm, verbose = FALSE)
  as.character(seurat_sub$seurat_clusters)
}
