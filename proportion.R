library(Seurat)
library(ggsci)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

merged <- readRDS("~/merged_r_annotation.rds")

col2 <- c('#F59C80','#E84B35','#F89E1E','#92D2C3','#95CC5E','#00A489','#7E6148','#B09C85','#8491B4')

#### proportion of cells in disease ----
cell_type_prop <- as.data.frame(prop.table(table(merged$cell_type, merged$disease), margin = 2) * 100)
colnames(cell_type_prop) <- c('cell_type',"disease","prop")
table(cell_type_prop$disease)
cell_type_prop$disease <- factor(cell_type_prop$disease,
                                 levels = rev(c('NC','AS','AMI','ICM','ACM','DCM','HCM','HF_preLVAD_NR','HF_postLVAD_NR','HF_preLVAD_R','HF_postLVAD_R')))
table(cell_type_prop$cell_type)
cell_type_prop$cell_type <- factor(cell_type_prop$cell_type,
                                   levels = rev(c('TLF+ Mac','MHCIIhi Mac','TREM2+ Mac','CCR2+ Mac','cMO','ncMO','cDC1','cDC2','Mast')))
p1 <- ggplot(cell_type_prop,aes(x=disease,y=prop,fill=cell_type))+
  geom_bar(stat = 'identity',width = 0.7, size = 0)+
  theme_bw()+ scale_fill_manual(values = rev(col2))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",angle = 0,hjust = 0.5,vjust = 1),
        axis.text.y = element_text(color = "black"))+
  labs(x = 'disease', y = 'proportion(%)', fill = 'celltype')+coord_flip()
p1