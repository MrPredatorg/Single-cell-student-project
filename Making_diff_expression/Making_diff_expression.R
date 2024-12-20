library('Seurat')
library('dplyr')

illumina_ser <- readRDS("Final_fold/illumina_ser")

illumina_ser <- illumina_ser %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:8) %>%
  FindNeighbors(dims = 1:8) %>%
  FindClusters(resolution = 0.4)

DimPlot(illumina_ser, reduction = "tsne", label = TRUE)

illumina_ser <- JoinLayers(object = illumina_ser)

cl_markers_illumina <- FindAllMarkers(illumina_ser, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

cl_markers_illumina %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(n =100)

new_ident <- setNames(c('mature Glutamatergic', 'Cajal-Retzius', 
                        'radial glia','imature GABAergic','cycling radial glia',
                        'imature Glutamatergic', 'radial glia'),
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
  RunTSNE(dims = 1:8) %>%
  FindNeighbors(dims = 1:8) %>%
  FindClusters(resolution = 0.6)

DimPlot(nanopore_ser, reduction = "tsne", label = TRUE)

nanopore_ser <- JoinLayers(object = nanopore_ser)

cl_markers_nanopore <- FindAllMarkers(nanopore_ser, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

cl_markers_nanopore %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(n =100)


new_ident <- setNames(c('mature Glutamatergic', 'imature GABAergic', 
                        'radial glia','cycling radial glia','cycling radial glia',
                        'mature GABAergic', 'Cajal-Retzius'),
                      levels(nanopore_ser))
nanopore_ser <- RenameIdents(nanopore_ser, new_ident)

plot1 <- FeaturePlot(nanopore_ser, c("Lrfn5","Gad2","Sparc",
                                     "Spc25","Rpp25","Lhx5"),
                     ncol=3, reduction = "tsne")

DimPlot(nanopore_ser, reduction = "tsne", label = TRUE) + NoLegend()

write.csv(cl_markers_illumina, file = 'Final_fold/illumina_markers_table')

write.csv(cl_markers_nanopore, file = 'Final_fold/nanopore_markers_table')