library(Seurat)
library(dplyr)

# pipeline and parameters based on report done by Lucy


####################################################################
#                           EDIT SECTION                           #

# working directory
workDir <- "~/osconet"

# 'case' name, carried through to script 2. Ends up in the OsconetInput
# filenames as filter<case>D<cluster>_<alpha>.csv, so case='WT' gives
# filterWTD0_0.2.csv
case <- "WT"

# input: folder, and the count file inside it
# .h5   -> Read10X_h5      (Cigliola)
# .rds  -> readRDS         (Matson, split by 0-SplitMatson.R)
# other -> ReadMtx prefix  (Yao, Saraswathy) e.g. "GSM6122390_ctrl"
inDir    <- "processed_raw_files/Cigliola_2023"
dataFile <- "GSM6586114_sample3_filtered_matrix.h5"

# output directory, where the seurat and counts will be stored for this sample 
outDir <- "processed_raw_files/Cigliola_2023/seurat_object_and_counts"

####################################################################


setwd(workDir)

# read count files
if (grepl("\\.h5$", dataFile)) {
  data <- Read10X_h5(file.path(inDir, dataFile))
} else if (grepl("\\.rds$", dataFile)) {
  data <- readRDS(file.path(inDir, dataFile))
} else {
  data <- ReadMtx(mtx      = file.path(inDir, paste0(dataFile, "_matrix.mtx.gz")),
                  features = file.path(inDir, paste0(dataFile, "_features.tsv.gz")),
                  cells    = file.path(inDir, paste0(dataFile, "_barcodes.tsv.gz")))
}

# check data
dim(data)
data[1:5, 1:5]

# create suerat object
obj_t <- CreateSeuratObject(counts = data, project = case)
head(obj_t)

# adding percentage of genes that are mitochondrial to metadata
obj_t[["percent.mt"]] <- PercentageFeatureSet(obj_t, pattern = "^mt-")
head(obj_t@meta.data, 5) # viewing QC metrics

# visualise QC with violin plot
VlnPlot(obj_t, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filtering cells
# NOTE: comment out for Matson - those counts are the authors' already
# filtered cell set, so this would be filtering twice
obj_t <- subset(obj_t, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                  percent.mt < 20)

# Normalise Data
# default param normalization.method = "LogNormalize", scale.factor = 10000
obj_t <- NormalizeData(obj_t)

# Identification of highly variable features
obj_t <- FindVariableFeatures(obj_t, nfeatures = 2000)

# Scaling the data
obj_t <- ScaleData(obj_t)

# Principal component analysis (PCA) on highly variable genes
obj_t <- RunPCA(obj_t, features = VariableFeatures(object = obj_t))

ElbowPlot(obj_t) # ranking of principle components based on the percentage of variance

# Clustering of cells
# K-nearest neighbor graph of first 10 principal components
obj_t <- FindNeighbors(obj_t, dims = 1:10)
# cluster with resolution of 0.8 to improve homogeneity
obj_t <- FindClusters(obj_t, resolution = 0.8)

# UMAP visualization
obj_t <- RunUMAP(obj_t, dims = 1:10)
DimPlot(obj_t, reduction = "umap")

# take raw counts and cluster assignment
counts <- obj_t[["RNA"]]$counts
clusters <- obj_t$seurat_clusters

#Attaches cell barcodes as names
names(clusters) <- colnames(obj_t)

# write outputs
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

# input to script 2: point its 'filename' and 'metafile' at these,
# and set its 'case' to the same value used above
write.csv(as.matrix(counts), file = file.path(outDir, paste0(case, "_counts.csv")))
write.csv(clusters,          file = file.path(outDir, paste0(case, "_clusters.csv")))

# saving Seurat object
saveRDS(obj_t, file = file.path(outDir, paste0(case, ".rds")))
