rm(list = ls())
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(Seurat)
library(tidyverse)
library(patchwork)
setwd("/Users/mardream_75/Desktop/Rdata")

tdm <- readRDS("tdm-original-RNA.rds")

DimPlot(tdm, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "orig.ident") + NoLegend()
DimPlot(tdm, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident") + NoLegend()


?FindMarkers

Idents(tdm)

#定义组别
tdm$celltype <- Idents(tdm)
tdm$celltype.group <- paste(tdm$celltype, tdm$orig.ident, sep = "_")
Idents(tdm) <- "celltype.group"

#差异分析,ident.1分子，ident.2分母,wilcox秩和检验,最低表达比例达多少才纳入

tdeg <- FindMarkers(tdm,ident.1 = 'T_P-TDM',ident.2 = 'T_M-TDM', verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)

cd4deg <- FindMarkers(tdm,ident.1 = 'CD4+ T_P-TDM',ident.2 = 'CD4+ T_M-TDM', verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)
head(cd4deg)

#循环查询每组别的差异基因
cellfordeg <- levels(tdm$celltype)
for(i in 1:length(cellfordeg))
{
  CELLDEG <- FindMarkers(tdm,ident.1 = paste0(cellfordeg[i],'_P-TDM'),ident.2 = paste0(cellfordeg[i],'_M-TDM'), verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)
  write.csv(CELLDEG,paste0("./241220差异基因/",cellfordeg[i],".CSV")) 
}

#筛选10个最大差异基因
library(dplyr)
top10 <- CELLDEG %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10

#差异基因可视化
tdm <- ScaleData(tdm, features = rownames(tdm))
DoHeatmap(tdm, features = top10, size = 3)

#总组间差异
VlnPlot(tdm, features = top10, idents = 'M-TDM')
FeaturePlot(tdm,features = top10[1:5],split.by = 'orig.ident')
DotPlot(tdm,features = top10,split.by ='orig.ident',cols = c('blue','yellow','pink'))

#火山图
EnhancedVolcano(tdeg, lab = rownames(tdeg), x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3.0, labSize = 6.0,
                title = 'diff_T')


