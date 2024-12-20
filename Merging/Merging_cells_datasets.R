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