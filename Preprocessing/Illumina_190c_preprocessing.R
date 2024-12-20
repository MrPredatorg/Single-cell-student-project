###################################################190 illumina

library(Seurat)

#reading matrix and SeuratObject

counts <- Read10X(data.dir = "sc_illumina_190/")

seurat <- CreateSeuratObject(counts, project="DS1", min.features = 0, min.cells = 0)

#checking mitochondrial genes and RNA count
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt[-\\.]")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, , pt.size = 2)

plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)

#Normalizing 

seurat <- NormalizeData(seurat)

seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)

plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)

seurat <- ScaleData(seurat)

#PCA

seurat <- RunPCA(seurat, npcs = 50)

ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

#tSNE

seurat <- RunTSNE(seurat, dims = 1:8)


plot1 <- TSNEPlot(seurat)

seurat <- FindNeighbors(seurat, dims = 1:8)

seurat <- FindClusters(seurat, resolution = 1)

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)

plot1

#saving_data

saveRDS(seurat, file="sc_illumina_190/seurat_i_190")

rm(list=ls())