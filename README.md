# BootStrap_SingleCell to Rub BootStrappin on Seurat (single cell object) to evaluate cluster stability.
bootstrap_myclusters <- function(x, FUN, clusters=NULL, transposed=FALSE, n.cells=5000, 
                                 iterations=30, ...) {
  if (is.null(clusters)) {
    clusters <- FUN(x, ...)
  }
  
  cluster.ids <- as.character(sort(unique(clusters)))
  output <- matrix(0, length(cluster.ids), length(cluster.ids))
  output[lower.tri(output)] <- NA_real_
  dimnames(output) <- list(cluster.ids, cluster.ids)
  
  for (i in seq_len(iterations)) {
    if (transposed) {
      chosen <- sample(nrow(x), n.cells)
      resampled <- x[chosen,]
    } else {
      chosen <- sample(ncol(x), n.cells)
      resampled <- x[,chosen]
    }
    
    reclusters <- FUN(resampled, ...)
    tab <- table(clusters[chosen], reclusters)
    
    for (j1 in seq_along(cluster.ids)) {
      spread1 <- tab[cluster.ids[j1],]
      spread1 <- spread1/sum(spread1)
      
      for (j2 in seq_len(j1)) {
        spread2 <- tab[cluster.ids[j2],]
        spread2 <- spread2/sum(spread2)
        output[j2,j1] <- output[j2,j1] + sum(spread1 * spread2)/iterations
      }
    }
  }
  
  output
}

#
