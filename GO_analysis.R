#### GO analysis ----
library(Seurat)

TREM2_MP <- readRDS('~/MP/TREM2_MP/clustering/TREM2_MP_annotation.rds')
TREM2_MP <- NormalizeData(TREM2_MP)
Idents(TREM2_MP) <- TREM2_MP$subclusters
markers <- FindAllMarkers(
  TREM2_MP,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
set.seed(717)
deg <- markers
head(deg)
deg <- deg[deg$cluster=="PLA2G7hiTREM2+_Mac",]
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
term <- bp@result

term_select <- term[term$Description %in% c('regulation of neuroinflammatory response',
                                            'plasma lipoprotein particle organization',
                                            'regulation of lipid metabolic process',
                                            'lipoprotein biosynthetic process',
                                            'macrophage derived foam cell differentiation'),]
df <- term_select
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
df$pvalue <- as.numeric(df$pvalue)
ggplot(data = df,
       aes(x = -log10(pvalue),
           y = reorder(Description,-log10(pvalue))))  +
  geom_bar(
    stat="identity",
    alpha=1,
    fill= '#CFC1A2',  #3b89e2 189f7f
    width = 0.6) +
  geom_text(aes(x=labelx,
                y=labely,
                label = df$Description),
            size=3.5,
            hjust =0)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12),
  )+
  xlab("-log10(pvalue)")+
  ggtitle("PLA2G7hiTREM2+_Mac")+
  scale_x_continuous(expand = c(0,0))