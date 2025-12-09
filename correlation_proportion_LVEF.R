library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

seurat.obj <-
  readRDS("~/merged_r_annotation.rds")
seurat.obj_LV <- subset(seurat.obj,
                        subset = region %in% c('LV'))

#### proportion cor with LVEF ----
sub_seurat.obj_LV <- subset(seurat.obj_LV,
                            subset = disease %in% c('AMI','ICM','ACM','DCM','HCM','NC'))

cluster_num_prop <- as.data.frame(prop.table(table(sub_seurat.obj_LV$cluster_num, sub_seurat.obj_LV$sample), margin = 2) * 100)
colnames(cluster_num_prop) <- c('cluster_num',"sample","prop")
cluster_num_prop <- cluster_num_prop[cluster_num_prop$cluster_num == 'c06',]
row.names(cluster_num_prop)  <- cluster_num_prop$sample

data_LVEF <- read.csv("~/cellproportion_cor_LVEF.csv")
rownames(data_LVEF) <- data_LVEF$X

samples <- Reduce(intersect, list(rownames(data_LVEF), rownames(cluster_num_prop)))
sample_LVEF <- data_LVEF[samples,]
sample_prop <- cluster_num_prop[samples,]

sample_combine <- cbind(sample_LVEF, sample_prop)
colnames(sample_combine)

library(ggpubr)
library(ggplot2)
ggscatter(sample_combine, 
          x = "LVEF", y = "prop",
          color = "red3",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson",
          # cor.coef.coord = c(25, 40),
          cor.coef.size = 5)+
  labs(x='LVEF(%)',y='Proportion of c06')


