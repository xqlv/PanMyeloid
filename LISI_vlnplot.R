library(ggpubr)
library(ggplot2)
library(ggsci)

LISI1 <- readRDS("~/integration/LISI1.rds")
LISI2 <- readRDS("~/integration/LISI2.rds")
LISI3 <- readRDS("~/integration/LISI3.rds")
LISI4 <- readRDS("~/integration/LISI4.rds")

df1 <- data.frame(step='step1', LISI=LISI1$sample)
df2 <- data.frame(step='step2', LISI=LISI2$sample)
df3 <- data.frame(step='step3', LISI=LISI3$sample)
df4 <- data.frame(step='step4', LISI=LISI4$sample)
df <- rbind(df1,df2,df3,df4)

ggviolin(df, "step", "LISI",  
         color = "step",
         # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot")+scale_color_npg()

