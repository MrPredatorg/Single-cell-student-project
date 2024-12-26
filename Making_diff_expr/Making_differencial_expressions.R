######################################################starting 

library(Seurat)

library(EnhancedVolcano)

library(SCPA)

########radial glia
nanollumina <- readRDS("Final_fold/nanollumina")


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
  nanollumina,
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

