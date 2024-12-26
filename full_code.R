install.packages("R.utils")
install.packages("Matrix")
install.packages("dplyr")
install.packages("Seurat")


library("R.utils")

untar("GSE130708_RAW.tar")

dir.create("sc_illumina_190")
dir.create("sc_illumina_951")
dir.create("sc_nanopore_190")
dir.create("sc_nanopore_951")


file.copy(from = "GSM3748086_190c_barcodes.tsv.gz",
          to   = "sc_illumina_190/barcodes.tsv.gz")

file.copy(from = "GSM3748086_190c_genes.tsv.gz",
          to   = "sc_illumina_190/features.tsv.gz")

file.copy(from = "GSM3748086_190c_matrix.mtx.gz",
          to   = "sc_illumina_190/matrix.mtx.gz")

file.copy(from = "GSM3748088_951c_barcodes.tsv.gz",
          to   = "sc_illumina_951/barcodes.tsv.gz")

file.copy(from = "GSM3748088_951c_genes.tsv.gz",
          to   = "sc_illumina_951/features.tsv.gz")

file.copy(from = "GSM3748088_951c_matrix.mtx.gz",
          to   = "sc_illumina_951/matrix.mtx.gz")

file.copy(from = "GSM3748087_190c.isoforms.matrix.txt.gz",
          to   = "sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt.gz")

file.copy(from = "GSM3748089_951c.isoforms.matrix.txt.gz",
          to   = "sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt.gz")

files_list <- list.files(pattern = ".gz")
files_list
file.remove(files_list)

rm(list=ls())



##############Making Nanopore 190c matrix for analisys

library("R.utils")
library('Matrix')
library("dplyr")

gunzip("sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt.gz", remove=TRUE)

df <-read.table("sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt",header=TRUE,sep="")

df <- df %>%
  group_by(geneId) %>%
  summarize(across(CAACTAGAGCTGTTCA:TTCTTAGTCTGTTGAG, sum))

df2 <- df[, -1]

rownames(df2) <- df$geneId

# save sparse matrix

gene_matrix <- Matrix(as.matrix(df2) , sparse = T )

head(gene_matrix)

writeMM(obj = gene_matrix, file="sc_nanopore_190/matrix.mtx")

# save genes and cells names

write.table(x = rownames(gene_matrix), file = "sc_nanopore_190/genes.tsv", col.names = FALSE, sep = ",")

write.table(x = colnames(gene_matrix), file = "sc_nanopore_190/barcodes.tsv", col.names = FALSE, sep = ",")

file.remove("sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt")


###############Making Nanopore  951c matrix for analisys

gunzip("sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt.gz", remove=TRUE)

df <-read.table("sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt",header=TRUE,sep="")

df <- df %>%
  group_by(geneId) %>%
  summarize(across(GGCGACTAGGCTCATT:GGCTCGAGTACCGAGA, sum))

df2 <- df[, -1]

rownames(df2) <- df$geneId

# save sparse matrix
gene_matrix <- Matrix(as.matrix(df2) , sparse = T )

head(gene_matrix)

writeMM(obj = gene_matrix, file="sc_nanopore_951/matrix.mtx")

# save genes and cells names

write.table(x = rownames(gene_matrix), file = "sc_nanopore_951/genes.tsv", col.names = FALSE, sep = ",")

write.table(x = colnames(gene_matrix), file = "sc_nanopore_951/barcodes.tsv", col.names = FALSE, sep = ",")

rm(list=ls())

file.remove("sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt")




#####################################################
###################################################
###################################################951 illumina

library(Seurat)

#reading matrix and SeuratObject

counts <- Read10X(data.dir = "sc_illumina_951/")

seurat <- CreateSeuratObject(counts, project="DS1", min.features = 0, min.cells = 0)

#checking mitochondrial genes and RNA count
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt[-\\.]")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, , pt.size = 2)

plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

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

saveRDS(seurat, file="sc_illumina_951/seurat_i_951")

rm(list=ls())


###################################################190 nanopore

library(Seurat)

#reading matrix and SeuratObject

counts <- Read10X(data.dir = "sc_nanopore_190/")

seurat <- CreateSeuratObject(counts, project="DS1", min.features = 0, min.cells = 0)

#checking mitochondrial genes and RNA count
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt[-\\.]")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, , pt.size = 2)

plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

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

saveRDS(seurat, file="sc_nanopore_190/seurat_N_190")

rm(list=ls())


###################################################951 nanopore

library(Seurat)

#reading matrix and SeuratObject

counts <- Read10X(data.dir = "sc_nanopore_951/")

seurat <- CreateSeuratObject(counts, project="DS1", min.features = 0, min.cells = 0)

#checking mitochondrial genes and RNA count
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt[-\\.]")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, , pt.size = 2)

plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

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

saveRDS(seurat, file="sc_nanopore_951/seurat_N_951")

rm(list=ls())


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

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

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


################seurat's merging

library('Seurat')
library('dplyr')

seurat_DS1 <- readRDS("sc_nanopore_190/seurat_N_190")

seurat_DS2 <- readRDS("sc_nanopore_951/seurat_N_951")

seurat_DS1$orig.ident <- '190'

seurat_DS2$orig.ident <- '951'

nanopore_ser <- merge(seurat_DS1, seurat_DS2) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:8)

plot1 <- DimPlot(nanopore_ser, group.by = "orig.ident", cols = c('190' ='green', '951'='purple'))

seurat_DS3 <- readRDS("sc_illumina_190/seurat_i_190")

seurat_DS4 <- readRDS("sc_illumina_951/seurat_i_951")

seurat_DS3$orig.ident <- '190'

seurat_DS4$orig.ident <- '951'

illumina_ser <- merge(seurat_DS3, seurat_DS4) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:8)

plot2 <- DimPlot(illumina_ser, group.by = "orig.ident", cols = c('190' ='green', '951'='purple'))

plot1 + plot2

dir.create("Final_fold")

illumina_ser <- merge(seurat_DS1, seurat_DS2, merge.dr = TRUE, merge.data = TRUE)

saveRDS(nanopore_ser, file="Final_fold/nanopore_ser")

illumina_ser <- merge(seurat_DS3, seurat_DS4, merge.dr = TRUE, merge.data = TRUE)

saveRDS(illumina_ser, file="Final_fold/illumina_ser")

unlink("sc_illumina_190/", recursive = TRUE)
unlink("sc_illumina_951/", recursive = TRUE)
unlink("sc_nanopore_190/", recursive = TRUE)
unlink("sc_nanopore_951/", recursive = TRUE)


library('Seurat')
library('dplyr')

illumina_ser <- readRDS("Final_fold/illumina_ser")

illumina_ser <- illumina_ser %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.3)

DimPlot(illumina_ser, reduction = "tsne", label = TRUE)

illumina_ser <- JoinLayers(object = illumina_ser)

cl_markers_illumina <- FindAllMarkers(illumina_ser, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

cl_markers_illumina %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC) %>% print(n =100)

new_ident <- setNames(c('imature Glutamatergic', 'cycling radial glia', 
                        'mature Glutamatergic','imature GABAergic','mature GABAergic',
                        'intermediate progenitor', 'radial glia', 'Cajal-Retzius'),
                      levels(illumina_ser))

illumina_ser <- RenameIdents(illumina_ser, new_ident)

plot1 <- FeaturePlot(illumina_ser, c("Lrfn5","Pcp4","Zfhx3",
                                     "Dgkg","Hs6st2","Atp1a2"),
                     ncol=3, reduction = "tsne")

DimPlot(illumina_ser, reduction = "tsne", label = TRUE) + NoLegend()

###################nanopore

nanopore_ser <- readRDS("Final_fold/nanopore_ser")

nanopore_ser <- nanopore_ser %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:9) %>%
  FindNeighbors(dims = 1:9) %>%
  FindClusters(resolution = 0.7)

DimPlot(nanopore_ser, reduction = "tsne", label = TRUE)

nanopore_ser <- JoinLayers(object = nanopore_ser)

cl_markers_nanopore <- FindAllMarkers(nanopore_ser, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

cl_markers_nanopore %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(n =100)


new_ident <- setNames(c('intermediate progenitor', 'imature GABAergic', 
                        'mature Glutamatergic','cycling radial glia','mature GABAergic',
                        'radial glia','Imature Glutamatergic', 'Cajal-Retzius'),
                      levels(nanopore_ser))
nanopore_ser <- RenameIdents(nanopore_ser, new_ident)

plot1 <- FeaturePlot(nanopore_ser, c("Lrfn5","Gad2","Sparc",
                                     "Spc25","Rpp25","Lhx5"),
                     ncol=3, reduction = "tsne")

DimPlot(nanopore_ser, reduction = "tsne", label = TRUE) + NoLegend()

write.csv(cl_markers_illumina, file = 'Final_fold/illumina_markers_table')

write.csv(cl_markers_nanopore, file = 'Final_fold/nanopore_markers_table')


nanopore_ser$orig.ident <- 'nanopore'

illumina_ser$orig.ident <- 'illumina'

nanollumina <- merge(illumina_ser, nanopore_ser)

nanollumina$orig.ident=='nanopore'

nanollumina <- nanollumina %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.9)

plot1 <- DimPlot(nanollumina, split.by = "orig.ident")
plot1


 <- JoinLayers(object = nanollumina)
cl_markers_nanollumina <- FindAllMarkers(, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

cl_markers_nanollumina %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(n =100)


new_ident <- setNames(c('intermediate progenitor', 'mature Glutamatergic', 
                        'cycling radial glia','cycling radial glia','radial glia',
                        'imature GABAergic','mature Glutamatergic', 'intermediate progenitor',
                        'mature GABAergic','imature Glutamatergic','mature Glutamatergic',
                        'intermediate progenitor','Cajal-Retzius'),
                      levels())
 <- RenameIdents(, new_ident)

plot1 <- DimPlot(, split.by = "orig.ident", reduction = "tsne", label = TRUE) + NoLegend()
plot1


head(@meta.data)

$CellType <- Idents()


levels($CellType)

######################################################starting 
nanollumina <- readRDS("Final_fold/nanollumina")

library(Seurat)

library(EnhancedVolcano)

library(SCPA)

########radial glia

radial_glia_i <- seurat_extract(
    nanollumina,
    meta1 = "orig.ident",
    value_meta1 = 'illumina',
    meta2 = "CellType",
    value_meta2 = "radial glia")


radial_glia_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "radial glia")


seurat1 <- CreateSeuratObject(radial_glia_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(radial_glia_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "Radial glia")

########intermediate progenitor

intermediate_progenitor_i <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "intermediate progenitor")


intermediate_progenitor_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "intermediate progenitor")


seurat1 <- CreateSeuratObject(intermediate_progenitor_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(intermediate_progenitor_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "intermediate progenitor")



########imature Glutamatergic

imature_Glutamatergic_i <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "imature Glutamatergic")


imature_Glutamatergic_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "imature Glutamatergic")


seurat1 <- CreateSeuratObject(imature_Glutamatergic_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(imature_Glutamatergic_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "imature Glutamatergic")


########mature Glutamatergic

mature_Glutamatergic_i <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "mature Glutamatergic")


mature_Glutamatergic_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "mature Glutamatergic")


seurat1 <- CreateSeuratObject(mature_Glutamatergic_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(mature_Glutamatergic_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "mature Glutamatergic")

########imature GABAergic

imature_GABAergic_i <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "imature GABAergic")


imature_GABAergic_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "imature GABAergic")


seurat1 <- CreateSeuratObject(imature_GABAergic_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(imature_GABAergic_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "imature GABAergic")


########Cajal-Retzius
nanollumina$CellType

Cajal-Retzius_i <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "Cajal-Retzius")


Cajal-Retzius_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "Cajal-Retzius")


seurat1 <- CreateSeuratObject(Cajal-Retzius_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(Cajal-Retzius_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "Cajal-Retzius")


########cycling radial glia

cycling_radial_glia_i <- seurat_extract(
  ,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "cycling radial glia")


cycling_radial_glia_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "cycling radial glia")


seurat1 <- CreateSeuratObject(cycling_radial_glia_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(cycling_radial_glia_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "cycling radial glia")

########mature GABAergic

mature_GABAergic_i <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'illumina',
  meta2 = "CellType",
  value_meta2 = "mature GABAergic")


mature_GABAergic_n <- seurat_extract(
  nanollumina,
  meta1 = "orig.ident",
  value_meta1 = 'nanopore',
  meta2 = "CellType",
  value_meta2 = "mature GABAergic")


seurat1 <- CreateSeuratObject(mature_GABAergic_i, project="DS1", min.features = 0, min.cells = 0)

seurat2 <- CreateSeuratObject(mature_GABAergic_n, project="DS2", min.features = 0, min.cells = 0)

seurat1$orig.ident <- "illumina"

seurat2$orig.ident <- "nanopore"


seurat3 <- merge(seurat1, seurat2)


aggregate_ifnb <- AggregateExpression(seurat3, group.by = "orig.ident", assays="RNA", return.seurat = TRUE)

aggregate_ifnb@meta.data

seurat3 <- NormalizeData(seurat3)

seurat3 <- JoinLayers(object = seurat3)

markers <- FindMarkers(seurat3, ident.1 = "illumina", group.by = "orig.ident", ident.2="nanopore")


EnhancedVolcano(markers,
                lab = rownames(markers),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = "mature GABAergic")

