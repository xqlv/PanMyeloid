library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)

seurat.obj <-
  readRDS("~/seurat.obj.rds")
table(seurat.obj$cell_type)
DimPlot(seurat.obj)

seurat.obj.sce <- as.SingleCellExperiment(seurat.obj)

sim <- slingshot(seurat.obj.sce, 
                 clusterLabels = 'cell_type', 
                 start.clus = 'cMO',
                 end.clus=c('MHCIIhigh Mac', 'TREM2+ Mac', 'CCR2+ Mac', 'ncMO', 'cDC1'),
                 reducedDim = 'UMAP')
sim
summary(sim$slingPseudotime_1)

plot(reducedDims(sim)$UMAP, col = c('#F59C80','#E84B35','#F89E1E','#92D2C3','#95CC5E','#00A489','#7E6148','#B09C85')[sim$cell_type], pch=16, asp = 1)
lines(SlingshotDataSet(sim), type = 'lineages', col = 'black', lwd = 2)