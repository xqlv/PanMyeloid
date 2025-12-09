library(destiny)
library(Biobase)
library(Seurat)

seurat_obj <- readRDS("~/merged_r_annotation.rds")
table(seurat_obj$disease)
seurat_obj <- subset(seurat_obj,
                     subset = disease %in% c('AS','AMI','ICM','ACM','DCM','HCM'))
table(seurat_obj$cell_type)
seurat_obj <- subset(seurat_obj,
                     subset = cell_type %in% c('MHCIIhi Mac', 'TREM2+ Mac', 'CCR2+ Mac', 'cMO', 'ncMO', 'cDC1', 'cDC2'))
seurat_obj <- NormalizeData(seurat_obj)

data <- seurat_obj@assays$RNA@data ##data <-GetAssayData(seurat_obj)
data <- data[VariableFeatures(seurat_obj),]
pd <- new('AnnotatedDataFrame', data = as.data.frame(seurat_obj@meta.data))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
destinyObj <- ExpressionSet(as.matrix(data), phenoData=pd, featureData=fd)

destinyObj$celltype<- seurat_obj$cell_type
destinyObj$cluster<- seurat_obj$cluster_num
destinyObj$disease<- seurat_obj$disease
destinyObj$sample<- seurat_obj$sample

dm <- DiffusionMap(destinyObj, k = 100)
dpt<-DPT(dm)






