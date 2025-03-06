rm(list = ls())
BiocManager::install("EnhancedVolcano")
library(ggrepel)
library(EnhancedVolcano)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
setwd("/Users/mardream_75/Desktop/Rdata")
setwd("C:/Users/Mardream/OneDrive/Desktop/Rdata")
tdm <- readRDS("tdm-25-treg.rds")

DimPlot(tdm, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "orig.ident") + NoLegend()
DimPlot(tdm, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident") + NoLegend()

#制图
tdm <- JoinLayers(tdm)
tdm.markers <- FindAllMarkers(tdm, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25) 
library(dplyr)
top5 <- tdm.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
tdm <- c(object = tdm, " " = "Unknowns")
DoHeatmap(tdm, features = top5$gene, size = 3, angle = -50, hjust = 0.8) + scale_fill_gradientn(colors = c("white","grey","firebrick3"))

?FindMarkers

Idents(tdm)

#定义组别
tdm$celltype <- Idents(tdm)
tdm$celltype.group <- paste(tdm$celltype, tdm$orig.ident, sep = "_")
Idents(tdm) <- "celltype.group"

#差异分析,ident.1分子，ident.2分母,wilcox秩和检验,最低表达比例达多少才纳入


cd4deg <- FindMarkers(tdm,ident.1 = 'CD4+ T_P-TDM',ident.2 = 'CD4+ T_M-TDM', verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)
head(cd4deg)

cd8deg <- FindMarkers(tdm,ident.1 = 'CD8+ T_P-TDM',ident.2 = 'CD8+ T_M-TDM', verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)
head(cd8deg)

tregdeg <- FindMarkers(tdm,ident.1 = 'Treg_P-TDM',ident.2 = 'Treg_M-TDM', verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)
head(tregdeg)


#循环查询每组别的差异基因
cellfordeg <- levels(tdm$celltype)
for(i in 1:length(cellfordeg))
{
  if (cellfordeg[i] == "Lymphocyte") {next} #Lymphocyte 在MTDM中没有
  CELLDEG <- FindMarkers(tdm,ident.1 = paste0(cellfordeg[i],'_P-TDM'),ident.2 = paste0(cellfordeg[i],'_M-TDM'), verbose = TRUE, test.use = 'wilcox', min.pct = 0.1)
  write.csv(CELLDEG,paste0("./25年/差异基因/",cellfordeg[i],".CSV")) 
}

#筛选10个最大差异基因
library(dplyr)
top10 <- CELLDEG %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10

#手动筛选每组
cd4top10 <- cd4deg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
cd4top10

cd8top10 <- cd48deg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
cd8top10

treg <- tregdeg %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
tregtop10


#火山图
EnhancedVolcano(cd4deg, lab = rownames(cd4deg), x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3.0, labSize = 6.0,
                title = 'diff_CD4+ T')

EnhancedVolcano(cd8deg, lab = rownames(cd8deg), x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3.0, labSize = 6.0,
                title = 'diff_CD8+ T')

EnhancedVolcano(tregdeg, lab = rownames(tregdeg), x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3.0, labSize = 6.0,
                title = 'diff_Treg')


#手动读入
cd4 <- read_csv("./25年/差异基因/CD4+ T.CSV")
cd8 <- read_csv("./25年/差异基因/CD8+ T.CSV")
treg <- read_csv("./25年/差异基因/Treg.CSV")
th <- read_csv("./25年/差异基因/Th.CSV")
cd4top15 <- cd4 %>% top_n(n = 15, wt = avg_log2FC) %>% pull(`...1`) #读入csv会多一列...1，里面存着基因
cd4top15

pdfout <- function(data1,custom_title) {
  p <- EnhancedVolcano(data1, lab = data1$`...1`, x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0,
                selectLab = data1$`...1`[data1$p_val_adj < 0.05 & abs(data1$avg_log2FC) > 1],
                max.overlaps = Inf,
                col = c("grey20","grey20","grey20","red"),
                labCol = 'black', labFace = 'bold', boxedLabels = TRUE, drawConnectors = TRUE,
                endsConnectors = "last", widthConnectors = 0.8,
                title = paste0('Differential gene_',custom_title), subtitle = NULL) + theme(legend.position = "none")
  ggsave(paste0('./25年/差异基因/图/',custom_title,'-volcano.pdf'), device = "pdf", width = 15, height = 7, dpi = 600, limitsize = FALSE)
}

pdfout(cd4,"CD4+ T")
pdfout(cd8,"CD8+ T")
pdfout(treg,"Treg")
pdfout(th,"Th")

#EnhancedVolcano(cd4, lab = cd4$`...1`, x = 'avg_log2FC', y = 'p_val_adj',
#                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 2.0, labSize = 3.0,
#                legendPosition = 'right', legendLabSize = 5, legendIconSize = 5,
#                selectLab = c('Igkc'),
#                title = 'diff_CD4+ T') #同理，以...1的数据
#差异基因可视化
#tdm <- ScaleData(tdm, features = rownames(tdm))
#DoHeatmap(tdm, features = top10, size = 3)


#总组间差异
#VlnPlot(tdm, features = top10, idents = 'M-TDM')
#FeaturePlot(tdm,features = top10[1:5],split.by = 'orig.ident')
#DotPlot(tdm,features = top10,split.by ='orig.ident',cols = c('blue','yellow','pink'))


