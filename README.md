# BootStrapping on Seurat (single cell object) to evaluate cluster stability.
Single-cell clusters are mathematical constructs while cell types are biological truth. There must be a consensus between biology and mathematics to interpret single-cell clusters. **BootStrapping** could be a nice way to evaluate the quality/stability of your SC clusters. Overlapping clusters may also indicate cellular hierarchies (e.g- Hemetopoitic Stem cells could show an overlap with progenitor cells) and **BootStrapping** could be a nice tool to unravel these connections as well.

**Step1:Load your seurat object, perform QC and clustering then load the function bellow to sample iteratively from previously loaded Seurat object**

```{r}
bootstrap_myclusters <- function(x, FUN, clusters=NULL, transposed=FALSE, n.cells=5000, 
                                 iterations=50, ...) {
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


**Step 2 Load function to Run the clustering iteratively**

```{r}
myknn_FUN <- function(x) {
  g <- FindNeighbors(x, verbose = F, reduction='pca', dims=1:30 ) 
  #Choose the appropriate dimension reduction method, e.g- 'mnn' if dataset was integrated with MNN batch correction method
  g <- FindClusters(g, verbose = F, resolution = 0.2) #Use the resloution of your choice (I prefer optimized via clustree function)
  as.numeric(g$seurat_clusters)}
```
  
**Step3:Run the BootStrap**
```{r}
originals<- seurat$seurat_clusters #This is the cluster or CellType information, you already have stored in Seurat object
coassign <- bootstrap_myclusters(seurat, clusters = originals, FUN = myknn_FUN, 
                                n.cells = ncol(seurat)-1, iterations = 50) #You can choose n.cells and iterations of your choice

#Plot heatmap of coassignmnet probabilities
pheatmap(coassign, cluster_row=F, cluster_col=F, main= "Coassignment probabilities", angle_col = 45,
         color=rev(viridis::magma(100)))
```
This will generate a heatmap of **Coassignment Probabilities (CP)** ranging from 0 to ~1. Higher the CP, Higher is the stability of your cluster. If the CP overlaps with other, check the cell type info of the cluster if they are biologically related then its probably a true biology otherwise re-evaluate your clustering. 
