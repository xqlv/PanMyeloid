#### MI time point prop ----
library(Seurat)
library(ggsci)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(clusterProfiler)
library(ggpubr)

TREM2_MP <- readRDS('/home/lvxq/heart_immu/pan_heart_disease/result_new2/MP/TREM2_MP/clustering/TREM2_MP_annotation.rds')

TREM2_MP_MI <- subset(TREM2_MP, subset = study %in% c('Kuppe et al.'))
table(TREM2_MP_MI$donor)

TREM2_MP_MI$days_after_infarction <- ''
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$disease %in% c('NC')] <- '0'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P10')] <- '6'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P11')] <- '31'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P12')] <- '31'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P13')] <- '45'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P14')] <- '62'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P15')] <- '11'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P16')] <- '4'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P18')] <- '153'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P19')] <- '40'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P2')] <- '5'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P20')] <- '166'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P3')] <- '2'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P4')] <- '88'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P5')] <- '101'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P6')] <- '3'
TREM2_MP_MI$days_after_infarction[TREM2_MP_MI$donor %in% c('P9')] <- '2'

meta_TREM2_MP_MI <- unique(TREM2_MP_MI@meta.data[,c('sample','days_after_infarction')])
rownames(meta_TREM2_MP_MI) <- meta_TREM2_MP_MI$sample

AMI <- subset(TREM2_MP_MI,subset = disease %in% "AMI")
sub_cluster_prop_AMI <- as.data.frame(prop.table(table(AMI$subclusters, AMI$sample), margin = 2) * 100)
colnames(sub_cluster_prop_AMI) <- c('sub_cluster',"sample","prop")
sub_cluster_prop_AMI$disease <- "AMI"

ICM <- subset(TREM2_MP_MI,subset = disease %in% "ICM")
sub_cluster_prop_ICM <- as.data.frame(prop.table(table(ICM$subclusters, ICM$sample), margin = 2) * 100)
colnames(sub_cluster_prop_ICM) <- c('sub_cluster',"sample","prop")
sub_cluster_prop_ICM$disease <- "ICM"

NC <- subset(TREM2_MP_MI,subset = disease %in% "NC")
sub_cluster_prop_NC <- as.data.frame(prop.table(table(NC$subclusters, NC$sample), margin = 2) * 100)
colnames(sub_cluster_prop_NC) <- c('sub_cluster',"sample","prop")
sub_cluster_prop_NC$disease <- "NC"

sub_cluster_prop <- rbind(sub_cluster_prop_NC, sub_cluster_prop_AMI, sub_cluster_prop_ICM)

meta_combine <- merge(sub_cluster_prop, meta_TREM2_MP_MI, by='sample')

# meta_combine <- meta_combine %>% arrange(days_after_infarction)
table(meta_combine$days_after_infarction)
meta_combine$days_after_infarction <- factor(meta_combine$days_after_infarction,
                                             levels = c('0','2','3','4','5','6','11','31','40','45','62','88','101','153','166'))


library(ggplot2)
library(dplyr)

colnames(meta_combine)
medians <- meta_combine %>%
  group_by(days_after_infarction, sub_cluster) %>%
  summarise(median_value = median(prop), .groups = 'drop')

ggplot(meta_combine, aes(x = days_after_infarction, y = prop, fill = sub_cluster)) +
  # geom_boxplot(aes(fill = sub_cluster),alpha = 0.5, position = position_dodge(0)) +  # 绘制箱型图
  geom_line(data = medians, aes(x = days_after_infarction, y = median_value, group = sub_cluster, color = sub_cluster), 
            size = 1.5, position = position_dodge(0)) +  # 绘制中位数折线
  geom_point(data = medians, aes(x = days_after_infarction, y = median_value, color = sub_cluster), 
             size = 3, position = position_dodge(0)) +  # 在中位数位置添加点
  labs(title = "", x = "days_after_infarction", y = "Proportion (%)") +
  theme_minimal() +
  theme(legend.title = element_blank())+theme_classic()+
  scale_color_manual(values = c('#CFC1A2','#BD696D','#80B18E','#9B86AD'))+
  scale_fill_manual(values = c('#CFC1A2','#BD696D','#80B18E','#9B86AD'))
