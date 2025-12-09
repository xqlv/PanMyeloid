library(ktplots)
library(tidyr)
library(Seurat)

count_net <- read.delim("~/cellphoneDB/cell_type_all/output/count_network.txt", check.names = FALSE)
count_matrix<-spread(count_net, TARGET, count)
rownames(count_matrix) <- count_matrix$SOURCE
count_matrix <- count_matrix[, -1]
count_matrix <- as.matrix(count_matrix)
count_matrix2 <- count_matrix[,-c(1,6,7,8)]

library(pheatmap)
pheatmap(count_matrix2, show_rownames = T, show_colnames = T, scale="none", cluster_cols = T,        
         border_color='white', cluster_rows = T, fontsize_row = 14, fontsize_col = 14,         
         main = "Cell-cell interaction", treeheight_row = 0, family = 'Arial',         
         color = colorRampPalette(c("dodgerblue4",'peachpuff','deeppink4' ))( 1000 ),         
         treeheight_col = 0,         
         display_numbers = T, number_color="white",fontsize_number=10, 
         number_format="%.0f", legend_labels = c(0,300))