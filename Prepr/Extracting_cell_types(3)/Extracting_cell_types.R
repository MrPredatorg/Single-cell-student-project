library('Seurat')
library('dplyr')

illumina_ser <- readRDS("Final_fold/illumina_ser")

nanopore_ser <- readRDS("Final_fold/nanopore_ser")


nanollumina_layer <- JoinLayers(object = nanollumina)

cl_markers_nanollumina <- FindAllMarkers(nanollumina_layer, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

cl_markers_nanollumina %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(n =100)


new_ident <- setNames(c('intermediate progenitor', 'mature Glutamatergic', 
                        'cycling radial glia','cycling radial glia','radial glia',
                        'imature GABAergic','mature Glutamatergic', 'intermediate progenitor',
                        'mature GABAergic','imature Glutamatergic','mature Glutamatergic',
                        'intermediate progenitor','Cajal-Retzius'),
                      levels(nanollumina_layer))
nanollumina_layer <- RenameIdents(nanollumina_layer, new_ident)

plot1 <- DimPlot(nanollumina_layer, split.by = "orig.ident", reduction = "tsne", label = TRUE) + NoLegend()
plot1


head(nanollumina_layer@meta.data)

nanollumina_layer$CellType <- Idents(nanollumina_layer)