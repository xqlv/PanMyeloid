library(Augur)
library(Seurat)
library(tidyverse)

### AMI_NC ----
Kuppe_et_al <- readRDS('~/Kuppe_et_al.rds')
table(Kuppe_et_al$disease)
AMI_NC <- subset(Kuppe_et_al,
                 subset = disease %in% c('AMI','NC'))
expr <- Seurat::GetAssayData(AMI_NC)
meta <- AMI_NC@meta.data %>% droplevels()
augur1 = calculate_auc(expr, meta, cell_type_col = "cell_type", label_col = "disease",n_threads = 8)
augur1$AUC
saveRDS(augur1,'~/Augur/augurAMI.rds')

### ICM_NC ----
Kuppe_et_al <- readRDS('~/Kuppe_et_al.rds')
table(Kuppe_et_al$disease)
ICM_NC <- subset(Kuppe_et_al,
                 subset = disease %in% c('ICM','NC'))
expr <- Seurat::GetAssayData(ICM_NC)
meta <- ICM_NC@meta.data %>% droplevels()
augur2 = calculate_auc(expr, meta, cell_type_col = "cell_type", label_col = "disease",n_threads = 8)
augur2$AUC
saveRDS(augur2,'~/Augur/augurICM_Kuppe_et_al.rds')

### ACM_NC ----
Reichart_et_al <- readRDS('~/Reichart_et_al.rds')
ACM_NC <- subset(Reichart_et_al,
                 subset = disease %in% c('ACM','NC'))
expr <- Seurat::GetAssayData(ACM_NC)
meta <- ACM_NC@meta.data %>% droplevels()
augur3 = calculate_auc(expr, meta, cell_type_col = "cell_type", label_col = "disease",n_threads = 8)
augur3$AUC
saveRDS(augur3,'~/Augur/augurACM.rds')

### DCM_NC ----
Reichart_et_al <- readRDS('~/Reichart_et_al.rds')
DCM_NC <- subset(Reichart_et_al,
                 subset = disease %in% c('DCM','NC'))
expr <- Seurat::GetAssayData(DCM_NC)
meta <- DCM_NC@meta.data %>% droplevels()
augur4 = calculate_auc(expr, meta, cell_type_col = "cell_type", label_col = "disease",n_threads = 8)
augur4$AUC
saveRDS(augur4,'~/Augur/augurDCM_Reichart_et_al.rds')

### HCM_NC ----
ChaffinM_et_al <- readRDS('~/ChaffinM_et_al.rds')
HCM_NC <- subset(ChaffinM_et_al,
                 subset = disease %in% c('HCM','NC'))
expr <- Seurat::GetAssayData(HCM_NC)
meta <- HCM_NC@meta.data %>% droplevels()
augur5 = calculate_auc(expr, meta, cell_type_col = "cell_type", label_col = "disease",n_threads = 8)
augur5$AUC
saveRDS(augur5,'~/Augur/augurHCM.rds')

#### vis ----
plot_lollipop(augur1) +
  geom_segment(aes(xend = cell_type, yend = 0.5), size = 0.5) +
  geom_point(size = 4, aes(color = cell_type)) +  
  scale_color_manual(values = c('#F59C80','#E84B35','#F89E1E','#92D2C3','#95CC5E','#00A489','#7E6148','#B09C85','#8491B4'))+
  theme(axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black'))
