library(Seurat)
library(ggsci)
library(ggplot2)
library(lisi)
library(tidyverse)

#### input data ----
ACM.DanielReichart2022 <- readRDS("~/QC/ACM.DanielReichart2022/seurat_obj_QC.rds")
AMI.ChristophKuppe2022 <- readRDS("~/QC/AMI.ChristophKuppe2022/seurat_obj_QC.rds")
AS.LukaNicin2022 <- readRDS("~/QC/AS.LukaNicin2022/seurat_obj_QC.rds")
DCM.DanielReichart2022 <- readRDS("~/QC/DCM.DanielReichart2022/seurat_obj_QC.rds")
DCM.MarkChaffin2022 <- readRDS("~/QC/DCM.MarkChaffin2022/seurat_obj_QC.rds")
HCM.MarkChaffin2022 <- readRDS("~/QC/HCM.MarkChaffin2022/seurat_obj_QC.rds")
HFDonor.JunedhMAmrute2023 <- readRDS("~/QC/HFDonor.JunedhMAmrute2023/seurat_obj_QC.rds")
HFNRpost.JunedhMAmrute2023 <- readRDS("~/QC/HFNRpost.JunedhMAmrute2023/seurat_obj_QC.rds")
HFNRpre.JunedhMAmrute2023 <- readRDS("~/QC/HFNRpre.JunedhMAmrute2023/seurat_obj_QC.rds")
HFRpost.JunedhMAmrute2023 <- readRDS("~/QC/HFRpost.JunedhMAmrute2023/seurat_obj_QC.rds")
HFRpre.JunedhMAmrute2023 <- readRDS("~/QC/HFRpre.JunedhMAmrute2023/seurat_obj_QC.rds")
ICM.BridgetSimonson2023 <- readRDS("~/QC/ICM.BridgetSimonson2023/seurat_obj_QC.rds")
ICM.ChristophKuppe2022 <- readRDS("~/QC/ICM.ChristophKuppe2022/seurat_obj_QC.rds")
NC.BridgetSimonson2023 <- readRDS("~/QC/NC.BridgetSimonson2023/seurat_obj_QC.rds")
NC.ChristophKuppe2022 <- readRDS("~/QC/NC.ChristophKuppe2022/seurat_obj_QC.rds")
NC.DanielReichart2022 <- readRDS("~/QC/NC.DanielReichart2022/seurat_obj_QC.rds")
NC.MarkChaffin2022 <- readRDS("~/QC/NC.MarkChaffin2022/seurat_obj_QC.rds")
NC.MonikaLitviňuková2020 <- readRDS("~/QC/NC.MonikaLitviňuková2020/seurat_obj_QC.rds")

#### merge data ----
merged <- merge(ACM.DanielReichart2022,
                y=c(AMI.ChristophKuppe2022,
                    AS.LukaNicin2022,
                    DCM.DanielReichart2022,
                    DCM.MarkChaffin2022,
                    HCM.MarkChaffin2022,
                    HFDonor.JunedhMAmrute2023,
                    HFNRpost.JunedhMAmrute2023,
                    HFNRpre.JunedhMAmrute2023,
                    HFRpost.JunedhMAmrute2023,
                    HFRpre.JunedhMAmrute2023,
                    ICM.BridgetSimonson2023,
                    ICM.ChristophKuppe2022,
                    NC.BridgetSimonson2023,
                    NC.ChristophKuppe2022,
                    NC.DanielReichart2022,
                    NC.MarkChaffin2022,
                    NC.MonikaLitviňuková2020))
counts <- merged@assays$RNA@counts
metadata <- merged@meta.data
merged <- CreateSeuratObject(counts = counts,
                             meta.data = metadata)

#### black list ----
MT_gene <- read.csv("~/Mitochondria_gene.txt")
MT_gene <- MT_gene[,1]
RB_gene <- read.csv("~/Ribosome_gene.txt")
RB_gene <- RB_gene[,1]
HSP_gene <- read.csv("~/Heat_shock_protein_gene.txt")
HSP_gene <- HSP_gene[,1]

#### harmony ----
ACM.DanielReichart2022 <- NormalizeData(ACM.DanielReichart2022, normalization.method = "LogNormalize", scale.factor = 10000)
ACM.DanielReichart2022 <- FindVariableFeatures(ACM.DanielReichart2022, selection.method = "vst", nfeatures = 6000)
VF_ACM.DanielReichart2022 <- VariableFeatures(ACM.DanielReichart2022)

AMI.ChristophKuppe2022 <- NormalizeData(AMI.ChristophKuppe2022, normalization.method = "LogNormalize", scale.factor = 10000)
AMI.ChristophKuppe2022 <- FindVariableFeatures(AMI.ChristophKuppe2022, selection.method = "vst", nfeatures = 6000)
VF_AMI.ChristophKuppe2022 <- VariableFeatures(AMI.ChristophKuppe2022)

DCM.DanielReichart2022 <- NormalizeData(DCM.DanielReichart2022, normalization.method = "LogNormalize", scale.factor = 10000)
DCM.DanielReichart2022 <- FindVariableFeatures(DCM.DanielReichart2022, selection.method = "vst", nfeatures = 6000)
VF_DCM.DanielReichart2022 <- VariableFeatures(DCM.DanielReichart2022)

DCM.MarkChaffin2022 <- NormalizeData(DCM.MarkChaffin2022, normalization.method = "LogNormalize", scale.factor = 10000)
DCM.MarkChaffin2022 <- FindVariableFeatures(DCM.MarkChaffin2022, selection.method = "vst", nfeatures = 6000)
VF_DCM.MarkChaffin2022 <- VariableFeatures(DCM.MarkChaffin2022)

HCM.MarkChaffin2022 <- NormalizeData(HCM.MarkChaffin2022, normalization.method = "LogNormalize", scale.factor = 10000)
HCM.MarkChaffin2022 <- FindVariableFeatures(HCM.MarkChaffin2022, selection.method = "vst", nfeatures = 6000)
VF_HCM.MarkChaffin2022 <- VariableFeatures(HCM.MarkChaffin2022)

ICM.ChristophKuppe2022 <- NormalizeData(ICM.ChristophKuppe2022, normalization.method = "LogNormalize", scale.factor = 10000)
ICM.ChristophKuppe2022 <- FindVariableFeatures(ICM.ChristophKuppe2022, selection.method = "vst", nfeatures = 6000)
VF_ICM.ChristophKuppe2022 <- VariableFeatures(ICM.ChristophKuppe2022)

NC.MonikaLitviňuková2020 <- NormalizeData(NC.MonikaLitviňuková2020, normalization.method = "LogNormalize", scale.factor = 10000)
NC.MonikaLitviňuková2020 <- FindVariableFeatures(NC.MonikaLitviňuková2020, selection.method = "vst", nfeatures = 6000)
VF_NC.MonikaLitviňuková2020 <- VariableFeatures(NC.MonikaLitviňuková2020)

HFDonor.JunedhMAmrute2023 <- NormalizeData(HFDonor.JunedhMAmrute2023, normalization.method = "LogNormalize", scale.factor = 10000)
HFDonor.JunedhMAmrute2023 <- FindVariableFeatures(HFDonor.JunedhMAmrute2023, selection.method = "vst", nfeatures = 6000)
VF_HFDonor.JunedhMAmrute2023 <- VariableFeatures(HFDonor.JunedhMAmrute2023)

NC.DanielReichart2022 <- NormalizeData(NC.DanielReichart2022, normalization.method = "LogNormalize", scale.factor = 10000)
NC.DanielReichart2022 <- FindVariableFeatures(NC.DanielReichart2022, selection.method = "vst", nfeatures = 6000)
VF_NC.DanielReichart2022 <- VariableFeatures(NC.DanielReichart2022)

NC.MarkChaffin2022 <- NormalizeData(NC.MarkChaffin2022, normalization.method = "LogNormalize", scale.factor = 10000)
NC.MarkChaffin2022 <- FindVariableFeatures(NC.MarkChaffin2022, selection.method = "vst", nfeatures = 6000)
VF_NC.MarkChaffin2022 <- VariableFeatures(NC.MarkChaffin2022)

VF_dataset <- Reduce(intersect, list(VF_ACM.DanielReichart2022, 
                                     VF_AMI.ChristophKuppe2022, 
                                     VF_DCM.DanielReichart2022, 
                                     VF_DCM.MarkChaffin2022,
                                     VF_HCM.MarkChaffin2022,
                                     VF_ICM.ChristophKuppe2022,
                                     VF_NC.MonikaLitviňuková2020,
                                     VF_HFDonor.JunedhMAmrute2023,
                                     VF_NC.DanielReichart2022,
                                     VF_NC.MarkChaffin2022))

VF_r <- setdiff(x=VF_dataset,y=MT_gene)
VF_r <- setdiff(x=VF_r,y=RB_gene)
VF_r <- setdiff(x=VF_r,y=HSP_gene)

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
VariableFeatures(merged) <- VF_r
merged <- ScaleData(merged,
                     features = VF_r,
                     vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb",
                                         "S.Score","G2M.Score"))
merged <- RunPCA(merged, features = VF_r)
merged <- merged |>
  harmony::RunHarmony(
    assay.use = "RNA",
    reduction.use = "pca",
    group.by.vars = "sample",
    plot_convergence = TRUE
  )

ElbowPlot(merged, ndims = 50)
merged <- RunUMAP(merged, reduction = 'harmony', dims = 1:30)

umap <- data.frame(Embeddings(merged,reduction = "umap"))
meta_data <- data.frame(sample=merged$sample,donor=merged$donor,dataset=merged$dataset)
LISI <- compute_lisi(umap, meta_data, c("sample",'donor','dataset'))






