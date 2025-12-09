library(Seurat)
library(tidyverse)
source("~/custom_function.R")

set.seed(717)
dir.create("~/deg_RNA_cell_type")

seurat.obj <- readRDS("~/merged_r_annotation.rds")

seurat.obj$celltype.disease <- paste(seurat.obj$cell_type, seurat.obj$disease, sep = "_")
table(seurat.obj$celltype.disease)

Idents(seurat.obj) <- seurat.obj$celltype.disease
table(Idents(seurat.obj))

DefaultAssay(seurat.obj) <- "RNA"
seurat.obj <- NormalizeData(seurat.obj)

#For cycle
outdir <- "~/deg_RNA_cell_type/DCM"
celltype <- c('TLF+ Mac', 'MHCIIhi Mac', 'TREM2+ Mac', 'CCR2+ Mac', 'cMO', 'ncMO', 'cDC1', 'cDC2', 'Mast')

for (celltype in celltype) {
  print(celltype)
  ident.1 <- paste(celltype, "DCM", sep = "_")
  ident.2 <- paste(celltype, "NC", sep = "_")
  temp_outdir <- file.path(outdir, celltype)
  dir.create(temp_outdir, recursive = T)
  set.seed(717)

  markers <- FindMarkers(seurat.obj,
                         ident.1 = ident.1,
                         ident.2 = ident.2,
                         #test.use = "MAST",
                         #latent.vars=c("study",'sample','gender'),
                         min.pct = 0,
                         logfc.threshold = 0)
  markers$celltype <- celltype
  markers$gene <- rownames(markers)
  write.csv(markers,file.path(temp_outdir, "markers.csv"))
  saveRDS(markers,file.path(temp_outdir, "markers.rds"))

  gsea.input <- markers

  #H
  gsea_res <- cat_gsea(gsea.input,
                       arrange_by = "avg_log2FC",
                       gene_name = "gene",
                       category = "H",
                       species = "human")
  #### 保存结果 ---
  # 分析结果
  saveRDS(gsea_res, file.path(temp_outdir, "gsea_H.rds"))

  # 导出表格
  write.csv(gsea_res, file.path(temp_outdir, "gsea_H_table.csv"),row.names = F)

  #C2
  gsea_res <- cat_gsea(gsea.input,
                       arrange_by = "avg_log2FC",
                       gene_name = "gene",
                       category = "C2",
                       species = "human")
  #### 保存结果 ---
  # 分析结果
  saveRDS(gsea_res, file.path(temp_outdir, "gsea_C2.rds"))

  # 导出表格
  write.csv(gsea_res, file.path(temp_outdir, "gsea_C2_table.csv"),row.names = F)

  #C5
  gsea_res <- cat_gsea(gsea.input,
                       arrange_by = "avg_log2FC",
                       gene_name = "gene",
                       category = "C5",
                       species = "human")
  #### 保存结果 ---
  # 分析结果
  saveRDS(gsea_res, file.path(temp_outdir, "gsea_C5.rds"))

  # 导出表格
  write.csv(gsea_res, file.path(temp_outdir, "gsea_C5_table.csv"),row.names = F)

}
