library(Seurat)
library(dplyr)
library(ggplot2)

# If your expression matrix contains a mix of different pre-determined cell clusters or labels, it may be wise to run OscoNet separately for each cell type.
# To do so specify the cell label you would like to isolate and run below. Remember to modify your output file names (counts, clusters) to reflect this.

# cellType <- 'her4_les' 

recluster <- function(ds) {
  
  # If you are running OscoNet on only a subpopulation then uncomment and use the code below, otherwise proceed as normal.
  # cols <- colnames(ds)
  # cellsList <- c()
  # for (I in 1:length(cols)) {if (startsWith(cols[I], 'MCF7_RM_CD44L') == T) {cellsList <- c(cellsList, cols[I])}}
  # df = subset(ds, select = cellsList)
  
  df = ds
  
  pbmc_MGH <- CreateSeuratObject(counts = df)
  
  # CHANGED FOR SEURAT V5: Reordered workflow. NormalizeData and FindVariableFeatures 
  # are now executed BEFORE ScaleData, which aligns with standard v5 best practices.
  pbmc_MGH <- NormalizeData(object = pbmc_MGH)
  pbmc_MGH <- FindVariableFeatures(object = pbmc_MGH)
  pbmc_MGH <- ScaleData(object = pbmc_MGH)
  
  pbmc_MGH <- RunPCA(object = pbmc_MGH)
  
  # CHANGED FOR SEURAT V5: Explicitly added dims and reduction parameters to FindNeighbors, 
  # and dims to RunTSNE to ensure proper data usage across v5 assay structures.
  pbmc_MGH <- FindNeighbors(object = pbmc_MGH, reduction = "pca", dims = 1:10)
  pbmc_MGH <- FindClusters(object = pbmc_MGH, resolution = 0.5)
  pbmc_MGH <- RunTSNE(object = pbmc_MGH, dims = 1:10)
  pbmc_MGH <- RunUMAP(object = pbmc_MGH, dims = 1:10)
  
  DimPlot(pbmc_MGH, reduction = "umap")
  DimPlot(pbmc_MGH, reduction = "tsne")
  
  # CHANGED FOR SEURAT V5: Replaced legacy slot access (pbmc_MGH@assays$RNA@counts) 
  # with the recommended Seurat v5 layer extraction syntax.
  counts <- GetAssayData(pbmc_MGH, assay = "RNA", layer = "counts")
  clusters <- pbmc_MGH$seurat_clusters
  names(clusters) <- colnames(pbmc_MGH)
  
  out <- list("counts" = counts, "clusters" = clusters, "shape" = "round")
  return(out) 
}

# Select appropriate line to match file extension of your matrix 
# ds = read.csv("/OscoNet/YourDataset/Data/expressionMatrix.csv", sep=',')
# ds = read.csv("/OscoNet/YourDataset/Data/expressionMatrix.txt", sep='\t')
ds = read.csv("C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/her4_les.tsv", sep='\t')

row.names(ds) <- ds$X
ds$X <- NULL  # Remove the 'X' column

out = recluster(ds)

# Modify these file names below if only running on a subset of this data.
write.csv(out$counts, file="C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/her4_les_counts.csv")
write.csv(out$clusters, file="C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/her4_les_clusters.csv")
rm(ds)