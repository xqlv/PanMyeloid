library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

seurat.obj <-
  readRDS("~/merged_r_annotation.rds")

#### proportion by sample ----
AS.LukaNicin2022 <- subset(seurat.obj,subset = dataset %in% "AS.LukaNicin2022")
cluster_num_prop_AS.LukaNicin2022 <- as.data.frame(prop.table(table(AS.LukaNicin2022$cluster_num, AS.LukaNicin2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_AS.LukaNicin2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_AS.LukaNicin2022$dataset <- "AS.LukaNicin2022"

AMI.ChristophKuppe2022 <- subset(seurat.obj,subset = dataset %in% "AMI.ChristophKuppe2022")
cluster_num_prop_AMI.ChristophKuppe2022 <- as.data.frame(prop.table(table(AMI.ChristophKuppe2022$cluster_num, AMI.ChristophKuppe2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_AMI.ChristophKuppe2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_AMI.ChristophKuppe2022$dataset <- "AMI.ChristophKuppe2022"

ICM.BridgetSimonson2023 <- subset(seurat.obj,subset = dataset %in% "ICM.BridgetSimonson2023")
cluster_num_prop_ICM.BridgetSimonson2023 <- as.data.frame(prop.table(table(ICM.BridgetSimonson2023$cluster_num, ICM.BridgetSimonson2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_ICM.BridgetSimonson2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_ICM.BridgetSimonson2023$dataset <- "ICM.BridgetSimonson2023"

ICM.ChristophKuppe2022 <- subset(seurat.obj,subset = dataset %in% "ICM.ChristophKuppe2022")
cluster_num_prop_ICM.ChristophKuppe2022 <- as.data.frame(prop.table(table(ICM.ChristophKuppe2022$cluster_num, ICM.ChristophKuppe2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_ICM.ChristophKuppe2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_ICM.ChristophKuppe2022$dataset <- "ICM.ChristophKuppe2022"

ACM.DanielReichart2022 <- subset(seurat.obj,subset = dataset %in% "ACM.DanielReichart2022")
cluster_num_prop_ACM.DanielReichart2022 <- as.data.frame(prop.table(table(ACM.DanielReichart2022$cluster_num, ACM.DanielReichart2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_ACM.DanielReichart2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_ACM.DanielReichart2022$dataset <- "ACM.DanielReichart2022"

DCM.DanielReichart2022 <- subset(seurat.obj,subset = dataset %in% "DCM.DanielReichart2022")
cluster_num_prop_DCM.DanielReichart2022 <- as.data.frame(prop.table(table(DCM.DanielReichart2022$cluster_num, DCM.DanielReichart2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_DCM.DanielReichart2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_DCM.DanielReichart2022$dataset <- "DCM.DanielReichart2022"

DCM.MarkChaffin2022 <- subset(seurat.obj,subset = dataset %in% "DCM.MarkChaffin2022")
cluster_num_prop_DCM.MarkChaffin2022 <- as.data.frame(prop.table(table(DCM.MarkChaffin2022$cluster_num, DCM.MarkChaffin2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_DCM.MarkChaffin2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_DCM.MarkChaffin2022$dataset <- "DCM.MarkChaffin2022"

HCM.MarkChaffin2022 <- subset(seurat.obj,subset = dataset %in% "HCM.MarkChaffin2022")
cluster_num_prop_HCM.MarkChaffin2022 <- as.data.frame(prop.table(table(HCM.MarkChaffin2022$cluster_num, HCM.MarkChaffin2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_HCM.MarkChaffin2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_HCM.MarkChaffin2022$dataset <- "HCM.MarkChaffin2022"

NC.MonikaLitviňuková2020 <- subset(seurat.obj,subset = dataset %in% "NC.MonikaLitviňuková2020")
cluster_num_prop_NC.MonikaLitviňuková2020 <- as.data.frame(prop.table(table(NC.MonikaLitviňuková2020$cluster_num, NC.MonikaLitviňuková2020$sample), margin = 2) * 100)
colnames(cluster_num_prop_NC.MonikaLitviňuková2020) <- c('cluster_num',"sample","prop")
cluster_num_prop_NC.MonikaLitviňuková2020$dataset <- "NC.MonikaLitviňuková2020"

NC.BridgetSimonson2023 <- subset(seurat.obj,subset = dataset %in% "NC.BridgetSimonson2023")
cluster_num_prop_NC.BridgetSimonson2023 <- as.data.frame(prop.table(table(NC.BridgetSimonson2023$cluster_num, NC.BridgetSimonson2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_NC.BridgetSimonson2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_NC.BridgetSimonson2023$dataset <- "NC.BridgetSimonson2023"

NC.ChristophKuppe2022 <- subset(seurat.obj,subset = dataset %in% "NC.ChristophKuppe2022")
cluster_num_prop_NC.ChristophKuppe2022 <- as.data.frame(prop.table(table(NC.ChristophKuppe2022$cluster_num, NC.ChristophKuppe2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_NC.ChristophKuppe2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_NC.ChristophKuppe2022$dataset <- "NC.ChristophKuppe2022"

NC.DanielReichart2022 <- subset(seurat.obj,subset = dataset %in% "NC.DanielReichart2022")
cluster_num_prop_NC.DanielReichart2022 <- as.data.frame(prop.table(table(NC.DanielReichart2022$cluster_num, NC.DanielReichart2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_NC.DanielReichart2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_NC.DanielReichart2022$dataset <- "NC.DanielReichart2022"

NC.MarkChaffin2022 <- subset(seurat.obj,subset = dataset %in% "NC.MarkChaffin2022")
cluster_num_prop_NC.MarkChaffin2022 <- as.data.frame(prop.table(table(NC.MarkChaffin2022$cluster_num, NC.MarkChaffin2022$sample), margin = 2) * 100)
colnames(cluster_num_prop_NC.MarkChaffin2022) <- c('cluster_num',"sample","prop")
cluster_num_prop_NC.MarkChaffin2022$dataset <- "NC.MarkChaffin2022"

HFDonor.JunedhMAmrute2023 <- subset(seurat.obj,subset = dataset %in% "HFDonor.JunedhMAmrute2023")
cluster_num_prop_HFDonor.JunedhMAmrute2023 <- as.data.frame(prop.table(table(HFDonor.JunedhMAmrute2023$cluster_num, HFDonor.JunedhMAmrute2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_HFDonor.JunedhMAmrute2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_HFDonor.JunedhMAmrute2023$dataset <- "HFDonor.JunedhMAmrute2023"

HFNRpre.JunedhMAmrute2023 <- subset(seurat.obj,subset = dataset %in% "HFNRpre.JunedhMAmrute2023")
cluster_num_prop_HFNRpre.JunedhMAmrute2023 <- as.data.frame(prop.table(table(HFNRpre.JunedhMAmrute2023$cluster_num, HFNRpre.JunedhMAmrute2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_HFNRpre.JunedhMAmrute2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_HFNRpre.JunedhMAmrute2023$dataset <- "HFNRpre.JunedhMAmrute2023"

HFNRpost.JunedhMAmrute2023 <- subset(seurat.obj,subset = dataset %in% "HFNRpost.JunedhMAmrute2023")
cluster_num_prop_HFNRpost.JunedhMAmrute2023 <- as.data.frame(prop.table(table(HFNRpost.JunedhMAmrute2023$cluster_num, HFNRpost.JunedhMAmrute2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_HFNRpost.JunedhMAmrute2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_HFNRpost.JunedhMAmrute2023$dataset <- "HFNRpost.JunedhMAmrute2023"

HFRpre.JunedhMAmrute2023 <- subset(seurat.obj,subset = dataset %in% "HFRpre.JunedhMAmrute2023")
cluster_num_prop_HFRpre.JunedhMAmrute2023 <- as.data.frame(prop.table(table(HFRpre.JunedhMAmrute2023$cluster_num, HFRpre.JunedhMAmrute2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_HFRpre.JunedhMAmrute2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_HFRpre.JunedhMAmrute2023$dataset <- "HFRpre.JunedhMAmrute2023"

HFRpost.JunedhMAmrute2023 <- subset(seurat.obj,subset = dataset %in% "HFRpost.JunedhMAmrute2023")
cluster_num_prop_HFRpost.JunedhMAmrute2023 <- as.data.frame(prop.table(table(HFRpost.JunedhMAmrute2023$cluster_num, HFRpost.JunedhMAmrute2023$sample), margin = 2) * 100)
colnames(cluster_num_prop_HFRpost.JunedhMAmrute2023) <- c('cluster_num',"sample","prop")
cluster_num_prop_HFRpost.JunedhMAmrute2023$dataset <- "HFRpost.JunedhMAmrute2023"

cluster_num_prop_dataset <- rbind(cluster_num_prop_AS.LukaNicin2022,
                                  cluster_num_prop_AMI.ChristophKuppe2022,
                                  cluster_num_prop_ICM.BridgetSimonson2023,
                                  cluster_num_prop_ICM.ChristophKuppe2022,
                                  cluster_num_prop_ACM.DanielReichart2022,
                                  cluster_num_prop_DCM.DanielReichart2022,
                                  cluster_num_prop_DCM.MarkChaffin2022,
                                  cluster_num_prop_HCM.MarkChaffin2022,
                                  cluster_num_prop_NC.MonikaLitviňuková2020,
                                  cluster_num_prop_NC.BridgetSimonson2023,
                                  cluster_num_prop_NC.ChristophKuppe2022,
                                  cluster_num_prop_NC.DanielReichart2022,
                                  cluster_num_prop_NC.MarkChaffin2022,
                                  cluster_num_prop_HFDonor.JunedhMAmrute2023,
                                  cluster_num_prop_HFNRpre.JunedhMAmrute2023,
                                  cluster_num_prop_HFNRpost.JunedhMAmrute2023,
                                  cluster_num_prop_HFRpre.JunedhMAmrute2023,
                                  cluster_num_prop_HFRpost.JunedhMAmrute2023)
cluster_num_prop_dataset$dataset <- factor(cluster_num_prop_dataset$dataset,
                                           levels = c('NC.MonikaLitviňuková2020',
                                                      'AS.LukaNicin2022',
                                                      'NC.ChristophKuppe2022',
                                                      'AMI.ChristophKuppe2022',
                                                      'ICM.ChristophKuppe2022',
                                                      'NC.BridgetSimonson2023',
                                                      'ICM.BridgetSimonson2023',
                                                      'NC.DanielReichart2022',
                                                      'ACM.DanielReichart2022',
                                                      'DCM.DanielReichart2022',
                                                      'NC.MarkChaffin2022',
                                                      'DCM.MarkChaffin2022',
                                                      'HCM.MarkChaffin2022',
                                                      'HFDonor.JunedhMAmrute2023',
                                                      'HFNRpre.JunedhMAmrute2023',
                                                      'HFNRpost.JunedhMAmrute2023',
                                                      'HFRpre.JunedhMAmrute2023',
                                                      'HFRpost.JunedhMAmrute2023'))
cluster_num_prop_dataset$disease <- ''
cluster_num_prop_dataset$disease <- gsub("\\..*", "", cluster_num_prop_dataset$dataset)

cluster_num_prop_dataset$study <- ''
cluster_num_prop_dataset$study[cluster_num_prop_dataset$dataset %in% c('AS.LukaNicin2022','NC.MonikaLitviňuková2020')] <- 'Nicin et al., 2022'
cluster_num_prop_dataset$study[cluster_num_prop_dataset$dataset %in% c('AMI.ChristophKuppe2022','ICM.ChristophKuppe2022',
                                                                       'NC.ChristophKuppe2022')] <- 'Kuppe et al., 2022'
cluster_num_prop_dataset$study[cluster_num_prop_dataset$dataset %in% c('ICM.BridgetSimonson2023','NC.BridgetSimonson2023')] <- 'Simonson et al., 2023'
cluster_num_prop_dataset$study[cluster_num_prop_dataset$dataset %in% c('ACM.DanielReichart2022','DCM.DanielReichart2022',
                                                                       'NC.DanielReichart2022')] <- 'Reichart et al., 2022'
cluster_num_prop_dataset$study[cluster_num_prop_dataset$dataset %in% c('DCM.MarkChaffin2022','HCM.MarkChaffin2022',
                                                                       'NC.MarkChaffin2022')] <- 'Chaffin M et al., 2022'
cluster_num_prop_dataset$study[cluster_num_prop_dataset$dataset %in% c('HFDonor.JunedhMAmrute2023','HFNRpre.JunedhMAmrute2023',
                                                                       'HFNRpost.JunedhMAmrute2023','HFRpre.JunedhMAmrute2023',
                                                                       'HFRpost.JunedhMAmrute2023')] <- 'Amrute et al., 2023'
library(ggpubr)

c06_prop <- cluster_num_prop_dataset[cluster_num_prop_dataset$cluster_num == "c06",]
c06_prop <- c06_prop[c06_prop$study %in% c('Simonson et al., 2023'),]
c06_prop$disease <- factor(c06_prop$disease, levels = c('NC','ICM'))
ggplot(c06_prop,aes(disease,prop,color=disease))+
  geom_boxplot(width = 0.5)+
  geom_point(size=0.1)+
  geom_jitter(width = 0.25)+
  # facet_grid(~sub_cluster_num)+
  stat_compare_means(comparisons = list(
    c("NC","ICM")),
    method = "t.test",label = "p.signif")+
  scale_color_manual(values = c('#F9BCB8','#5569B1'))+
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,color = 'black',size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.title.y = element_text(size = 10),
        # legend.position = "none",
        # legend.direction = "vertical",
        # legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,color = 'black',size = 10))+
  labs(y="proportion of c06 (%)", title = 'Simonson et al., 2023')



