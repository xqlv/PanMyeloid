library(tidyverse)
library(reshape2)
library(rstatix)
library(ggpubr)
library(GGally)
library(ggsci)
source("~/custom_function.R")
source("~/custom_plot_function.R")

#读入数据
####合并所有差异基因----
readCSV <- function(dir_dta){
  file_list <- list.files(path=dir_dta,full.names=T)
  varSave_func <- function(x){
    table_x <- read.csv(file=x,sep=",",header = T)
  }
  a<-invisible(lapply(file_list,FUN=varSave_func))
  b<-as.data.frame(a[[1]])
  for (i in 2:length(a)){
    c<-rbind(b,a[[i]])
    b <- c
  }
  return(b)
}

##########调用函数
dir_dta <- "~/pseudobulk/DESeq2/disease/Waterfallplot/"  #目录下为Findmarker出的结果
subset <-readCSV(dir_dta)

# subset <- subset[which(subset$p_val_adj != "NA"),]
subset <- na.omit(subset)

colnames(subset)
log2FC <- 0.25 #log2FC最低标准0.25，最高1.5
subset <- subset[abs(subset$log2FoldChange) >= log2FC & subset$padj < 0.05, ]
# subset <- subset[subset$pct.1 >= 0.1 | subset$pct.2 >= 0.1,]
subset$condition <- "No change"
subset[subset$log2FoldChange >= log2FC, ]$condition  <- "Upregulated"
subset[subset$log2FoldChange <= -log2FC, ]$condition  <- "Downregulated"
table(subset$condition,subset$Type)
# subset$celltype <- substring(subset$celltype, 0, 3)
subset$condition=paste(subset$cluster,subset$condition,sep = "_")
subset <- subset[,-1]
# colnames(subset) <- c("PValue","logFC","pct.1","pct.2","PValue_adj","cluster","gene","condition")
marker_condition <- subset


##画图
group1_name = "Upregulated"
group2_name = "Downregulated"
plot1_d = 40
plot2_d = 20
text.size = 4

deg.group1=marker_condition[str_detect(marker_condition$condition,paste0("_",group1_name)),]
deg.group2=marker_condition[str_detect(marker_condition$condition,paste0("_",group2_name)),]

### high expression genes in group1 #####################################################
celltype.order.name=names(sort(table(deg.group1$cluster)))
celltype.order=1:length(celltype.order.name)
names(celltype.order)=celltype.order.name

gene.shared=sort(table(deg.group1$gene))[sort(table(deg.group1$gene)) >= 2]
gene.shared=names(gene.shared)
gene.unique=sort(table(deg.group1$gene))[sort(table(deg.group1$gene)) < 2]
gene.unique=names(gene.unique)

group1.gene.shared=deg.group1[deg.group1$gene %in% gene.shared,c("gene","condition","cluster")]
write.csv(group1.gene.shared, '~/pseudobulk/DESeq2/disease/group1_up.gene.shared.csv')

gene.shared.freq=as.data.frame(sort(table(group1.gene.shared$gene)))
colnames(gene.shared.freq)=c("gene","freq")
gene.shared.freq$gene=as.character(gene.shared.freq$gene)
gene.order=c()
for (freqi in rev(unique(gene.shared.freq$freq))) {
  one.freq.df=as.data.frame(gene.shared.freq %>% filter(freq == freqi))
  if (dim(one.freq.df)[1] > 1) {
    one.gene.left=c()
    for (genei in one.freq.df[,"gene"]) {
      one.gene.distribution=group1.gene.shared[group1.gene.shared$gene %in% genei,]
      one.gene.left=append(one.gene.left,min(celltype.order[one.gene.distribution$cluster]))
    }
    one.gene.leftcelltype.df=as.data.frame(one.gene.left)
    one.gene.leftcelltype.df$gene=one.freq.df[,"gene"]
    one.gene.leftcelltype.df=one.gene.leftcelltype.df%>%arrange(one.gene.left)
    gene.order=append(gene.order,one.gene.leftcelltype.df$gene)
  }else{
    gene.order=append(gene.order,one.freq.df[1,"gene"])
  }
}

group1.gene.unique=deg.group1[deg.group1$gene %in% gene.unique,c("gene","condition","cluster")]
write.csv(group1.gene.unique, '~/pseudobulk/DESeq2/disease/group1_up.gene.unique.csv')

group1.gene.unique$cluster=factor(group1.gene.unique$cluster,levels = celltype.order.name)
group1.gene.unique=group1.gene.unique%>%arrange(cluster)
gene.order=append(gene.order,group1.gene.unique$gene)

group1.plot.data=deg.group1[,c("gene","condition","cluster")]
group1.plot.data$gene=factor(group1.plot.data$gene,levels = rev(gene.order))
group1.plot.data$cluster=factor(group1.plot.data$cluster,levels = celltype.order.name)

p1=group1.plot.data%>%ggplot(aes(x=cluster,y=gene))+
  geom_stripped_cols(odd = "#d9d9d9",even ="white",alpha=0.5)+
  geom_tile(aes(width=0.9),color='#C30D23',fill='#C30D23')+
  geom_hline(yintercept = length(gene.unique)+0.5,linetype=5)+
  scale_y_discrete(paste0("shared: ",length(gene.shared)," genes; unique: ",length(gene.unique)," genes"))+
  scale_x_discrete(expand = c(0,0))+
  labs(title = paste0(group1_name," DEGs"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.y.left = element_text(size = 16),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    plot.title = element_text(size = 18,hjust = 0.5)
  )

plot1.lower.text.df=as.data.frame(table(group1.gene.unique$cluster))
group1.gene.shared$cluster=factor(group1.gene.shared$cluster,levels = celltype.order.name)
plot1.upper.text.df=as.data.frame(table(group1.gene.shared$cluster))
plot1_h = length(gene.unique)

p3 = p1+
  geom_text(data = plot1.lower.text.df,mapping = aes(x=Var1,y=plot1_d,label=Freq),size=text.size)+
  geom_text(data = plot1.upper.text.df,mapping = aes(x=Var1,y=plot1_d+plot1_h,label=Freq),size=text.size)
p3

### high expression genes in group2 ############################################
celltype.order.name=names(sort(table(deg.group2$cluster)))
celltype.order=1:length(celltype.order.name)
names(celltype.order)=celltype.order.name

gene.shared=sort(table(deg.group2$gene))[sort(table(deg.group2$gene)) >= 2]
gene.shared=names(gene.shared)
gene.unique=sort(table(deg.group2$gene))[sort(table(deg.group2$gene)) < 2]
gene.unique=names(gene.unique)

group2.gene.shared=deg.group2[deg.group2$gene %in% gene.shared,c("gene","condition","cluster")]
write.csv(group2.gene.shared, '~/pseudobulk/DESeq2/disease/group2_down.gene.shared.csv')

gene.shared.freq=as.data.frame(sort(table(group2.gene.shared$gene)))
colnames(gene.shared.freq)=c("gene","freq")
gene.shared.freq$gene=as.character(gene.shared.freq$gene)
gene.order=c()
for (freqi in rev(unique(gene.shared.freq$freq))) {
  one.freq.df=as.data.frame(gene.shared.freq %>% filter(freq == freqi))
  if (dim(one.freq.df)[1] > 1) {
    one.gene.left=c()
    for (genei in one.freq.df[,"gene"]) {
      one.gene.distribution=group2.gene.shared[group2.gene.shared$gene %in% genei,]
      one.gene.left=append(one.gene.left,min(celltype.order[one.gene.distribution$cluster]))
    }
    one.gene.leftcelltype.df=as.data.frame(one.gene.left)
    one.gene.leftcelltype.df$gene=one.freq.df[,"gene"]
    one.gene.leftcelltype.df=one.gene.leftcelltype.df%>%arrange(one.gene.left)
    gene.order=append(gene.order,one.gene.leftcelltype.df$gene)
  }else{
    gene.order=append(gene.order,one.freq.df[1,"gene"])
  }
}

group2.gene.unique=deg.group2[deg.group2$gene %in% gene.unique,c("gene","condition","cluster")]
write.csv(group2.gene.unique, '~/pseudobulk/DESeq2/disease/group2_down.gene.unique.csv')

group2.gene.unique$cluster=factor(group2.gene.unique$cluster,levels = celltype.order.name)
group2.gene.unique=group2.gene.unique%>%arrange(cluster)
gene.order=append(gene.order,group2.gene.unique$gene)

group2.plot.data=deg.group2[,c("gene","condition","cluster")]
group2.plot.data$gene=factor(group2.plot.data$gene,levels = rev(gene.order))
group2.plot.data$cluster=factor(group2.plot.data$cluster,levels = celltype.order.name)

p2=group2.plot.data%>%ggplot(aes(x=cluster,y=gene))+
  geom_stripped_cols(odd = "#d9d9d9",even ="white",alpha=0.5)+
  geom_tile(aes(width=0.9),color='#026EB7',fill='#026EB7')+
  geom_hline(yintercept = length(gene.unique)+0.5,linetype=5)+
  scale_y_discrete(paste0("shared: ",length(gene.shared)," genes; unique: ",length(gene.unique)," genes"))+
  scale_x_discrete(expand = c(0,0))+
  labs(title = paste0(group2_name," DEGs"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.y.left = element_text(size = 16),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    plot.title = element_text(size = 18,hjust = 0.5)
  )

plot2.lower.text.df=as.data.frame(table(group2.gene.unique$cluster))
group2.gene.shared$cluster=factor(group2.gene.shared$cluster,levels = celltype.order.name)
plot2.upper.text.df=as.data.frame(table(group2.gene.shared$cluster))
plot2_h = length(gene.unique)

p4 = p2+
  geom_text(data = plot2.lower.text.df,mapping = aes(x=Var1,y=plot2_d,label=Freq),size=text.size)+
  geom_text(data = plot2.upper.text.df,mapping = aes(x=Var1,y=plot2_d+plot2_h,label=Freq),size=text.size)
p4

p3 + p4




