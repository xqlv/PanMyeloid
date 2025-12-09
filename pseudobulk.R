rm(list=ls())
library(Seurat)
library(tidyverse)
library(survminer)
library(GeneOverlap)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(sva)
library(DESeq2)
set.seed(717)

#### TLF+ Mac ----
scRNA = readRDS("~/merged_r_annotation.rds") 
dim(scRNA)
table(scRNA$study)
scRNA <- subset(scRNA, subset = study %in% c('Chaffin M et al.','Kuppe et al.',
                                             'Reichart et al.','Simonson et al.'))

table(scRNA$disease)
scRNA$condition <- ''
scRNA$condition[scRNA$disease %in% c('AMI','ICM','ACM','DCM','HCM')] <- 'HF'
scRNA$condition[scRNA$disease %in% c('NC')] <- 'NC'
table(scRNA$condition)
table(scRNA$cell_type)
scRNA <- subset(scRNA, subset = cell_type %in% c('TLF+ Mac'))
dim(scRNA)

scRNA <- SetIdent(scRNA, value = scRNA@meta.data$sample)

keep <- scRNA[["RNA"]]@counts[rowSums(scRNA@assays[["RNA"]]@counts > 1) >= 10, ] %>% rownames

bs = split(colnames(scRNA),scRNA$sample)  #sample
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(scRNA) %in% bs[[x]]
    rowSums(as.matrix(scRNA@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
head(ct)
rownames(ct) <- rownames(scRNA)
bulk_counts <- ct

bulk_counts <- bulk_counts[keep, ] 
coldata <- scRNA@meta.data[,c('sample','donor','region','condition','dataset','gender','age','cell_type','study')]
coldata <- unique(coldata)
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(bulk_counts), ]
stopifnot(all(colnames(bulk_counts) == rownames(coldata)))  

cts <- ComBat_seq(bulk_counts %>% as.matrix(), batch = coldata$study)

coldata$condition <- as.factor(coldata$condition)  
coldata$gender <- as.factor(coldata$gender)
coldata$age <- as.factor(coldata$age)

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ gender + age + condition  
)
dds <- estimateSizeFactors(dds)  
dds <- DESeq(dds)  

res <- results(dds, contrast = c("condition", "HF", "NC"))
res <- data.frame(res)
res$gene <- rownames(res)
res$Group <- 'HF'
res$cell_type <- 'TLF+ Mac'

DEG_deseq2 <- res %>%
  mutate(Type = if_else(padj > 0.05, "stable",
                        if_else(abs(log2FoldChange) < 0.25, "stable",
                                if_else(log2FoldChange >= 0.25, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Symbol")
