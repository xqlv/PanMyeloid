#### corrplot ----
library(corrplot)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

data <- read.csv("~/cellproportion_cor/data_cell_type.csv")
rownames(data) <- data$X
data <- data[,-1]
data_cell <- data
adata <- cor (data_cell, method="pearson")
testRes = cor.mtest(data_cell, method="pearson",conf.level = 0.95)
corrplot(adata)
corrplot(adata, method = "circle", col = rev(colorRampPalette(c(brewer.pal(11,'RdBu')))(200)),
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt", #number.font = 0.25,cl.cex = 0.25,
         p.mat = testRes$p, diag = T, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')
corrplot(adata, method = "number", type = "lower", col = rev(colorRampPalette(c(brewer.pal(11,'RdBu')))(200)),
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE', #number.font = 0.25, cl.cex = 0.25,
         add = T)

p <- corrplot(adata, method = "number", type = "lower", col = rev(colorRampPalette(c(brewer.pal(11,'RdBu')))(200)),
              tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE', #number.font = 0.25, cl.cex = 0.25,
              add = T)