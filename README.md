# BootStrapping on Seurat (single cell object) to evaluate cluster stability.

**Step1:Load function to sample iteratively from previously loaded Seurat object**

```{r}
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
```


**Step 2 function to Run the clustering iteratively**

```{r}
myknn_FUN <- function(x) {
  g <- FindNeighbors(x, verbose = F, reduction='pca', dims=1:30 ) 
  g <- FindClusters(g, verbose = F, resolution = 0.2) #Use the resloution of your choice (I prefer optimized via clustree function)
  as.numeric(g$seurat_clusters)}
```
  
**Step3:Run BootStrap**
```{r}
originals<- seurat$seurat_clusters #This is the cluster or cluster information, you already have stored in Seurat object
coassign <- bootstrap_myclusters(seurat, clusters = originals, FUN = myknn_FUN, 
                                n.cells = ncol(seurat)-1, iterations = 30) #You can choose n.cells and iterations of your choice

#Plot heatmap of coassignmnet probabilities
pheatmap(coassign, cluster_row=F, cluster_col=F, main= "Coassignment probabilities", angle_col = 45,
         color=rev(viridis::magma(100)))
```

