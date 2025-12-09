library(Seurat)
library(ggsci)
library(ggplot2)
library(lisi)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(tibble)

merged <- readRDS("~/mouse/clustering/sub_merged3_anno.rds")
merged <- NormalizeData(merged)
markers1 <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
table(markers1$cluster)
markers1$cluster <- factor(markers1$cluster,
                           levels = c('TLF_MP','MHCII_MP','Trem2_MP','Ccr2_MP',
                                      'IFNIC_MP','Prolif_MP',
                                      'cMO','ncMO','cDC1','cDC2','NP'))

merged2 <- readRDS("~/merged_r_annotation.rds")
table(merged2$cluster_num)
table(merged2$cell_type)
merged2$cell_type2 <- ''
merged2$cell_type2 <- merged2$cell_type
merged2$cell_type2 <- as.character(merged2$cell_type2)
merged2$cell_type2[merged2$cluster_num %in% c('c01','c02','c03','c04','c07')] <- 'TLF+ Mac'
merged2$cell_type2[merged2$cluster_num %in% c('c05')] <- 'IFNIC Mac'
merged2$cell_type2[merged2$cluster_num %in% c('c06')] <- 'Prolif Mac'
table(merged2$cell_type2)
merged2$cell_type2 <- factor(merged2$cell_type2,
                             levels = c('TLF+ Mac','MHCIIhi Mac','TREM2+ Mac','CCR2+ Mac',
                                        'IFNIC Mac','Prolif Mac',
                                        'cMO','ncMO','cDC1','cDC2','Mast'))
Idents(merged2) <- merged2$cell_type2
table(Idents(merged2))
merged2 <- NormalizeData(merged2)
markers2 <- FindAllMarkers(
  merged2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

df_gene1=markers1[markers1$p_val_adj <= 1e-5 & markers1$avg_log2FC >= 0.75,]
df_gene2=markers2[markers2$p_val_adj <= 1e-5 & markers2$avg_log2FC >= 0.75,]

cluster1=names(table(df_gene1$cluster))
cluster2=names(table(df_gene2$cluster))

library(homologene)
inTax  <-  10090
outTax  <-  9606

df_gene1_convert <-data.frame(cluster='cluster', gene='gene')
for (i in cluster1) {
  output <- homologene(df_gene1[df_gene1$cluster %in% c(i),]$gene, inTax = inTax, outTax = outTax)  
  df <- data.frame(cluster=i, gene=unique(output[,'9606']))
  df_gene1_convert <- rbind(df_gene1_convert, df)
}

df_gene1_convert <- df_gene1_convert[-1,]

## jaccard
df_ja=c()

for (i in cluster1) {
  ja=c()
  for (j in cluster2) {
    a=df_gene1_convert[df_gene1_convert$cluster==i,]
    a=a$gene
    b=df_gene2[df_gene2$cluster==j,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    ja=c(ja,jaccard)
  }
  df_ja=rbind(df_ja,ja)
}

rownames(df_ja)=paste0(cluster1)
colnames(df_ja)=paste0(cluster2)

library(pheatmap)
pheatmap(df_ja,
         cluster_cols = F,cluster_rows = F,scale = "none",
         show_rownames = T, show_colnames =T,
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         breaks = seq(0,0.18,length.out = 100),border_color = NA) 
