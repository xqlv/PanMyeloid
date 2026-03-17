### data ---
library(Seurat)
library(ggsci)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(clusterProfiler)
library(ggpubr)

source("~/custom_function.R")
set.seed(717)

deg_AMI <- readRDS("~/pseudobulk/DESeq2/disease/AMI/markers_HF.rds")
deg_ICM <- readRDS("~/pseudobulk/DESeq2/disease/ICM/markers_HF.rds")
deg_ACM <- readRDS("~/pseudobulk/DESeq2/disease/ACM/markers_HF.rds")
deg_DCM <- readRDS("~/pseudobulk/DESeq2/disease/DCM/markers_HF.rds")
deg_HCM <- readRDS("~/pseudobulk/DESeq2/disease/HCM/markers_HF.rds")

gsea.input <- deg_AMI
colnames(gsea.input)
gsea_res_AMI <- cat_gsea(gsea.input,
                         arrange_by = "log2FoldChange",
                         gene_name = "gene",
                         category = "C5",
                         species = "human")

gsea.input <- deg_ICM
colnames(gsea.input)
gsea_res_ICM <- cat_gsea(gsea.input,
                         arrange_by = "log2FoldChange",
                         gene_name = "gene",
                         category = "C5",
                         species = "human")

gsea.input <- deg_ACM
colnames(gsea.input)
gsea_res_ACM <- cat_gsea(gsea.input,
                         arrange_by = "log2FoldChange",
                         gene_name = "gene",
                         category = "C5",
                         species = "human")

gsea.input <- deg_DCM
colnames(gsea.input)
gsea_res_DCM <- cat_gsea(gsea.input,
                         arrange_by = "log2FoldChange",
                         gene_name = "gene",
                         category = "C5",
                         species = "human")
gsea.input <- deg_HCM
colnames(gsea.input)
gsea_res_HCM <- cat_gsea(gsea.input,
                         arrange_by = "log2FoldChange",
                         gene_name = "gene",
                         category = "C5",
                         species = "human")


df_gsea_res_AMI <- gsea_res_AMI@result
df_gsea_res_AMI$disease <- 'AMI'
df_gsea_res_ICM <- gsea_res_ICM@result
df_gsea_res_ICM$disease <- 'ICM'
df_gsea_res_ACM <- gsea_res_ACM@result
df_gsea_res_ACM$disease <- 'ACM'
df_gsea_res_DCM <- gsea_res_DCM@result
df_gsea_res_DCM$disease <- 'DCM'
df_gsea_res_HCM <- gsea_res_HCM@result
df_gsea_res_HCM$disease <- 'HCM'

df_gsea_res <- rbind(df_gsea_res_AMI,
                     df_gsea_res_ICM,
                     df_gsea_res_ACM,
                     df_gsea_res_DCM,
                     df_gsea_res_HCM) %>% as.data.frame()
df_gsea_res$change <- ""
df_gsea_res$change <- "no"
df_gsea_res[df_gsea_res$NES > 1 & df_gsea_res$pvalue < 0.05 & df_gsea_res$qvalue < 0.25,]$change <- "up"
df_gsea_res[df_gsea_res$NES < -1 & df_gsea_res$pvalue < 0.05 & df_gsea_res$qvalue < 0.25,]$change <- "down"

library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(vroom)
library(forcats)

df_gsea_res_select <- df_gsea_res[df_gsea_res$Description %in% c('GOBP_OXIDATIVE_PHOSPHORYLATION',
                                                                 'GOMF_FATTY_ACID_BINDING',
                                                                 'GOBP_TYPE_2_IMMUNE_RESPONSE',
                                                                 'GOBP_INTERLEUKIN_13_PRODUCTION',
                                                                 'GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_12_PRODUCTION',
                                                                 'GOBP_ACTIVATION_OF_PROTEIN_KINASE_B_ACTIVITY',
                                                                 'GOBP_MITOTIC_CELL_CYCLE_CHECKPOINT_SIGNALING',
                                                                 'GOBP_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION'),]
colnames(df_gsea_res_select)
library(Hmisc)
df_gsea_res_select$Description <- tolower(df_gsea_res_select$Description)
df_gsea_res_select$Description <- gsub("gobp_", "", df_gsea_res_select$Description)

df_gsea_res_select$disease <- factor(df_gsea_res_select$disease,
                                     levels = rev(c('AMI','ICM','ACM','DCM','HCM')))
df_gsea_res_select$Description <- factor(df_gsea_res_select$Description,
                                         levels = c('oxidative_phosphorylation',
                                                    'fatty_acid_binding',
                                                    'type_2_immune_response',
                                                    'interleukin_13_production',
                                                    'negative_regulation_of_interleukin_12_production',
                                                    'activation_of_protein_kinase_b_activity',
                                                    "mitotic_cell_cycle_checkpoint_signaling",
                                                    'regulation_of_endothelial_cell_migration'))
df_gsea_res_select %>% 
  ggplot(aes(x = NES, y = disease)) +
  geom_vline(xintercept = 0, color = 'grey', linewidth = 1) +
  geom_segment(aes(x = 0, xend = NES, y = disease, yend = disease), color = 'grey', linewidth = 1) +
  geom_point(aes(color = disease, size = -log10(p.adjust))) +
  # scale_color_gradient(low = 'blue', high = 'red', breaks = seq(3, 9, 2), limits = c(1, 9), name = 'Significance\n(-log10 FDR)') +
  scale_color_manual(values = rev(c('#AFD796','#61BAAF','#EC74B1','#FFC179','#00C4E4')), name = 'Disease')+
  scale_size_continuous(name = 'Significance\n(-log10 FDR)', 
                        # limits = c(30, 60), 
                        range = c(2, 5)) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-2, 3, 1), expand = c(0, 0)) +
  labs(x = 'Normalized Enrichment Score (NES)', y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1, color = 'grey'),
    panel.spacing.x = unit(0.1, units = 'in'),
    strip.background = element_rect(linewidth = 0.5, color = 'grey'),
    strip.text = element_text(size = 7),
    legend.position = "right",  # 图例位置（right/bottom）
    axis.ticks = element_line(color = 'black'),
    axis.text.x = element_text(colour = 'black', size = 8),
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.title.x = element_text(color = 'black', size = 8),
    legend.title = element_text(color = 'black', size = 8),
    legend.text = element_text(color = 'black', size = 8)
  ) +
  facet_wrap(~Description, ncol = 4)




