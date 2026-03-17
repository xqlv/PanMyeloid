#### GO analysis ----
TREM2_MP <- readRDS('~/MP/TLF_MP/clustering/markers.rds')

library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
set.seed(717)
deg <- markers
head(deg)
deg <- deg[deg$cluster=="c07",]
gene <- deg$gene
bp <-
  enrichGO(
    gene,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
    # ,
    # qvalueCutoff = 0.2
  )
