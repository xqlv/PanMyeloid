### data ---
gsea_TREM2_MP_AMI <- read.csv('~/MP/TREM2_MP/pseudobulk_deg/GSEA/gsea_C5_AMI.csv')
gsea_TREM2_MP_AMI$disease <- 'AMI'
gsea_TREM2_MP_ICM <- read.csv('~/MP/TREM2_MP/pseudobulk_deg/GSEA/gsea_C5_ICM.csv')
gsea_TREM2_MP_ICM$disease <- 'ICM'
gsea_TREM2_MP_ACM <- read.csv('~/MP/TREM2_MP/pseudobulk_deg/GSEA/gsea_C5_ACM.csv')
gsea_TREM2_MP_ACM$disease <- 'ACM'
gsea_TREM2_MP_DCM <- read.csv('~/MP/TREM2_MP/pseudobulk_deg/GSEA/gsea_C5_DCM.csv')
gsea_TREM2_MP_DCM$disease <- 'DCM'
gsea_TREM2_MP_HCM <- read.csv('~/MP/TREM2_MP/pseudobulk_deg/GSEA/gsea_C5_HCM.csv')
gsea_TREM2_MP_HCM$disease <- 'HCM'

gsea_TREM2_MP_combine <- rbind(gsea_TREM2_MP_AMI,gsea_TREM2_MP_ICM,gsea_TREM2_MP_ACM,gsea_TREM2_MP_DCM,gsea_TREM2_MP_HCM)
colnames(gsea_TREM2_MP_combine)
gsea_TREM2_MP_combine <- gsea_TREM2_MP_combine[,c(1,2,3,13,4,5,6,7,8,9,10,11,12)]

library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(vroom)
library(forcats)

gsea_TREM2_MP_combine_select <- gsea_TREM2_MP_combine[gsea_TREM2_MP_combine$Description %in% c("GOBP_WOUND_HEALING",
                                                                                               "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
                                                                                               'GOBP_INSULIN_RECEPTOR_SIGNALING_PATHWAY',
                                                                                               'GOBP_LIPOPROTEIN_BIOSYNTHETIC_PROCESS',
                                                                                               'GOMF_MHC_PROTEIN_COMPLEX_BINDING'),]
colnames(gsea_TREM2_MP_combine_select)
library(Hmisc)
gsea_TREM2_MP_combine_select$Description <- tolower(gsea_TREM2_MP_combine_select$Description)
gsea_TREM2_MP_combine_select$Description <- gsub("gobp_", "", gsea_TREM2_MP_combine_select$Description)
gsea_TREM2_MP_combine_select$Description <- gsub("gocc_", "", gsea_TREM2_MP_combine_select$Description)
gsea_TREM2_MP_combine_select$Description <- gsub("gomf_", "", gsea_TREM2_MP_combine_select$Description)

gsea_TREM2_MP_combine_select$disease <- factor(gsea_TREM2_MP_combine_select$disease,
                                               levels = rev(c('AMI','ICM','ACM','DCM','HCM')))
gsea_TREM2_MP_combine_select$Description <- factor(gsea_TREM2_MP_combine_select$Description,
                                                   levels = c("wound_healing",
                                                              "collagen_containing_extracellular_matrix",
                                                              'insulin_receptor_signaling_pathway',
                                                              'lipoprotein_biosynthetic_process',
                                                              'mhc_protein_complex_binding'))
gsea_TREM2_MP_combine_select %>% 
  ggplot(aes(x = NES, y = disease)) +
  geom_vline(xintercept = 0, color = 'grey', linewidth = 1) +
  geom_segment(aes(x = 0, xend = NES, y = disease, yend = disease), color = 'grey', linewidth = 1) +
  geom_point(aes(color = disease, size = -log10(p.adjust))) +
  # scale_color_gradient(low = 'blue', high = 'red', breaks = seq(3, 9, 2), limits = c(1, 9), name = 'Significance\n(-log10 FDR)') +
  scale_color_manual(values = rev(c('#AFD796','#61BAAF','#EC74B1','#FFC179','#00C4E4')), name = 'Disease')+
  scale_size_continuous(name = 'Significance\n(-log10 FDR)', 
                        # limits = c(30, 60), 
                        range = c(2, 5)) +
  scale_x_continuous(limits = c(-2, 3), breaks = seq(-2, 3, 1), expand = c(0, 0)) +
  labs(x = 'Normalized Enrichment Score (NES)', y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1, color = 'grey'),
    panel.spacing.x = unit(0.1, units = 'in'),
    strip.background = element_rect(linewidth = 0.5, color = 'grey'),
    strip.text = element_text(size = 7),
    axis.ticks = element_line(color = 'black'),
    axis.text.x = element_text(colour = 'black', size = 8),
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.title.x = element_text(color = 'black', size = 8),
    legend.title = element_text(color = 'black', size = 8),
    legend.text = element_text(color = 'black', size = 8)
  ) +
  facet_wrap(~Description, ncol = 3)

