setwd("scSpatialReconstructCompare-Paper")

#########
## Load and process tabular muris data, downloaded from figshare.
# https://figshare.com/articles/Robject_files_for_tissues_processed_by_Seurat/5821263
# Download droplet_Liver_seurat_tiss.Robj
#########

#########
# Downlaoded the data from Figshare:
#########
library(Seurat)
load("DATA/droplet_Liver_seurat_tiss.Robj")

#########
## Subset just the labeled hepatocytes
#########
droplet_liver_data <- as.matrix(tiss@raw.data[,which(tiss@meta.data$subtissue == "hepatocytes")])
head(droplet_liver_data)
dim(droplet_liver_data)


#########
## QC on the cells:
#########
hepatocytes <- CreateSeuratObject(counts = droplet_liver_data, project = "hepatocytes", 
                                  min.cells = 3, min.features = 200)
FeatureScatter(hepatocytes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
hepatocytes <- subset(hepatocytes, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
hepatocytes <- NormalizeData(hepatocytes, normalization.method = "LogNormalize", scale.factor = 10000)
hepatocytes <- FindVariableFeatures(hepatocytes, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(hepatocytes)
# Save this...will use this to get the variable genes later on for ordering!
hvg <- HVFInfo(hepatocytes)

# Check if data all cluster together:
hepatocytes <- ScaleData(hepatocytes, features = rownames(hepatocytes))
hepatocytes <- RunPCA(hepatocytes, features = VariableFeatures(object = hepatocytes))
ElbowPlot(hepatocytes)
hepatocytes <- FindNeighbors(hepatocytes, dims = 1:20)
hepatocytes <- FindClusters(hepatocytes, resolution = 0.5)
hepatocytes <- RunUMAP(hepatocytes, dims = 1:20)
DimPlot(hepatocytes, reduction = "umap")

# Use the cells that are in the larger group area only:
use.cells <- rownames(subset(hepatocytes@meta.data, seurat_clusters %in% c(0,2,3,5)))
length(use.cells) # 606 cells after QC

############## Subset the data with only these cells:
droplet_liver_data <- droplet_liver_data[,use.cells]


#########
## Data needs to be normalized: Use scran for UMI data.
#########

library(scran)
sce <- SingleCellExperiment(list(counts=droplet_liver_data))
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

sce <- normalize(sce, return_log=FALSE)

library(SingleCellExperiment)
droplet.norm <- normcounts(sce)


####################################

save(droplet_liver_data, droplet.norm, hvg,
     file="RDATA/dataReady_DropletPaper.RData")

