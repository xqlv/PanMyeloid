library(Seurat)
library(SeuratData)
library(tidyverse)
library(RColorBrewer)
library(ggsci)
set.seed(717)

merged <- readRDS('~/cellphoneDB/merged_total.rds')
table(merged$disease)
table(merged$cluster_num)

seurat_obj <- subset(merged, subset = cluster_num %in% c('c03', 'Endothelial'))

library(nichenetr)

output <- "~/nichenetr/c03_Endothelial/"
dir.create(output, recursive = T)

table(seurat_obj$disease)
seurat_obj$group <- ''
seurat_obj$group[seurat_obj$disease %in% c('AMI','ICM','ACM','DCM','HCM')] <- 'HF'
seurat_obj$group[seurat_obj$disease %in% c('NC')] <- 'NC'
table(seurat_obj$group)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)

table(seurat_obj$cluster_num)
Idents(seurat_obj) <- seurat_obj$cluster_num
table(Idents(seurat_obj))

#### Prepare NicheNet analysis ----
ligand_target_matrix = readRDS("~/nichenetr/ligand_target_matrix.rds")
ligand_target_matrix[1:5, 1:5] # target genes in rows, ligands in columns
lr_network = readRDS("~/nichenetr/lr_network.rds")
head(lr_network)
weighted_networks = readRDS("~/nichenetr/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

#### Perform the NicheNet analysis ----
## receiver
receiver = "Endothelial"
expressed_genes_receiver = get_expressed_genes(receiver, 
                                               seurat_obj, 
                                               pct = 0.10,
                                               assay_oi = "RNA")
saveRDS(expressed_genes_receiver,
        file.path(output, "expressed_genes_receiver.rds"))
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c('c03')

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.10, "RNA") # lapply to get the expressed genes of every sender cell type separately here
saveRDS(list_expressed_genes_sender, file.path(output, "list_expressed_genes_sender.rds"))
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


#### Define a gene set of interest ----
seurat_obj_receiver = subset(seurat_obj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

condition_oi = "HF"
condition_reference = "NC"

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, 
                                ident.1 = condition_oi, 
                                ident.2 = condition_reference, 
                                min.pct = 0.10) %>% rownames_to_column("gene")
saveRDS(DE_table_receiver, file.path(output, "DE_table_receiver.rds"))

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 &
                                            abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#### Define a set of potential ligands ----
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors = intersect(receptors, expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands &
                                            to %in% expressed_receptors) %>% pull(from) %>% unique()

#### Perform NicheNet ligand activity analysis ----
ligand_activities = predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# best_upstream_ligands <- c(best_upstream_ligands, "CCL14", "SPP1")
DotPlot(seurat_obj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

## Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(
  get_weighted_ligand_target_links,
  geneset = geneset_oi,
  ligand_target_matrix = ligand_target_matrix,
  n = 200
) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33
)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets, order_ligands] %>% t()

# order_targets enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
set.seed(717)
dir.create("~/nichenetr/c03_Endothelial/function")
bp <-
  enrichGO(
    order_targets,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
    # ,
    # qvalueCutoff = 0.2
  )
term <- bp@result
write.table(
  term,
  file = "~/nichenetr/c03_Endothelial/function/order_targets_gobp.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(bp, "~/nichenetr/c03_Endothelial/function/order_targets_gobp.rds")

pathways <- c("epithelial cell proliferation",
              "wound healing",
              "nitric oxide biosynthetic process",
              "vasculogenesis",
              "nitric oxide metabolic process",
              "positive regulation of vasculogenesis",
              "endothelial cell migration",
              'regulation of vasculogenesis')

df <- subset(term,
             subset = Description %in% pathways)
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
df$pvalue <- as.numeric(df$pvalue)
p1 <- ggplot(data = df,
             aes(x = -log10(pvalue),
                 y = reorder(Description,-log10(pvalue))))  +
  geom_bar(
    stat="identity",
    alpha=1,
    fill= '#ee6470',  #3b89e2 189f7f
    width = 0.8) +
  geom_text(aes(x=labelx,
                y=labely,
                label = df$Description),
            size=3.5,
            hjust =0)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12),
  )+
  xlab("-log10(pvalue)")+
  ggtitle("")+
  scale_x_continuous(expand = c(0,0))
p1





