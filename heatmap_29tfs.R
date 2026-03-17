deg_TLF_MP <- read.csv('~/pseudobulk/DESeq2/celltype/markers_TLFMP_HF.csv',row.names = 1)
deg_MHCII_MP <- read.csv('~/pseudobulk/DESeq2/celltype/markers_MHCIIMP_HF.csv',row.names = 1)
deg_TREM2_MP <- read.csv('~/pseudobulk/DESeq2/celltype/markers_TREM2MP_HF.csv',row.names = 1)
deg_CCR2_MP <- read.csv('~/pseudobulk/DESeq2/celltype/markers_CCR2MP_HF.csv',row.names = 1)
deg_cMO <- read.csv('~/pseudobulk/DESeq2/celltype/markers_cMO_HF.csv',row.names = 1)
deg_ncMO <- read.csv('~/pseudobulk/DESeq2/celltype/markers_ncMO_HF.csv',row.names = 1)
deg_DC <- read.csv('~/pseudobulk/DESeq2/celltype/markers_DC_HF.csv',row.names = 1)

deg_TLF_MP$cluster <- 'TLF+ Mac'
deg_MHCII_MP$cluster <- 'MHCIIhi Mac'
deg_TREM2_MP$cluster <- 'TREM2+ Mac'
deg_CCR2_MP$cluster <- 'CCR2+ Mac'
deg_cMO$cluster <- 'cMO'
deg_ncMO$cluster <- 'ncMO'
deg_DC$cluster <- 'DC'

deg_celltype <- rbind(deg_TLF_MP, deg_MHCII_MP, deg_TREM2_MP, deg_CCR2_MP, deg_cMO, deg_ncMO, deg_DC)

tfs_selected <- unique(c("ELF1","NR3C1","ETV6","JDP2","ETV5","NFIL3","NFIB",
                         "HDX","PRDM1","IRF2","FOXO1","FOXO3","SREBF2","PKNOX2",
                         "BCL11A","ELK4","ATF4","ERG","BCL6","IRF8",
                         "ETV6","NR3C1","ELF1","ETV5","JDP2","XBP1","NFIL3",
                         "HDX","IRF2","ATF4","CEBPB","NFIB","ELK3","E2F7","NFKB1",
                         "ELF4","BCL6","E2F3","E2F1","PRDM1","RUNX3"))

deg_celltype_tfs <- deg_celltype[deg_celltype$gene %in% tfs_selected,]

library(tidyverse)
colnames(deg_celltype_tfs)
df1 <- deg_celltype_tfs[,c(8,12,3)]
data1 <- df1 %>% pivot_wider(names_from = "cluster", values_from = "log2FoldChange")
data1 <- as.matrix(data1)
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
data1 <- t(data1)

df2 <- deg_celltype_tfs[,c(8,12,7)]
data2 <- df2 %>% pivot_wider(names_from = "cluster", values_from = "padj")
data2 <- as.matrix(data2)
rownames(data2) <- data2[,1]
data2 <- data2[,-1]
data2 <- t(data2)

class(data1) <- "numeric"
data1[data1 >= 1] <- 1
data1[data1 <= -1] <- -1

class(data2) <- "numeric"
data2[data2 < 0.05] <- "*"
data2[data2 >= 0.05] <- ""

data3 <- data1
data3[data3 >= 0.25 | data3 <= -0.25] <- "*"
data3[data3 !=  "*"] <- ""

p_text <- ifelse(data3 == "*" & data2 == "*", "*", "")
colnames(p_text) <- colnames(data3)
rownames(p_text) <- rownames(data3)

library('pheatmap')
library(RColorBrewer)
pheatmap(as.matrix(data1), 
         scale = "none", cluster_rows=F, cluster_cols = T,
         display_numbers = p_text, 
         fontsize_number = 15, number_color = 'black',
         # color = colorRampPalette(c("white","#FFE4B5","#D55E00","darkred"))(256),
         color = colorRampPalette(rev(c(brewer.pal(9, "RdBu"))))(256),
         cellwidth = 20, cellheight = 20, 
         border_color = "white", fontsize = 15, 
         # # display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)),
         # # annotation_col = annotation_col,
         # annotation_colors = ann_colors,
         # annotation_col = annotation_col,
         # treeheight_row = 10, treeheight_col =10,
         # gaps_col=c(1),
         # gaps_col = c(5,9,13,18,23),
         angle_col = c("90"),
         # clustering_method = "average",
         # main = "DEGs",
         show_rownames = T,show_colnames = T
)