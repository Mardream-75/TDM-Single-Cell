rm(list = ls())

BiocManager::install("RColorBrewer")
BiocManager::install("magrittr")
BiocManager::install("patchwork")
BiocManager::install("tidydr")
library(devtools)
install.packages('NMF')
BiocManager::install('BiocNeighbors')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github('immunogenomics/presto')
devtools::install_github("jinworks/CellChat",type="binary")

library(Seurat)
# library(SeuratData)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
library(tidydr)


setwd("/Users/mardream_75/Desktop/Rdata")
setwd("C:/Users/Mardream/OneDrive/Desktop/Rdata")
tdm <- readRDS("tdm-241219-cd4cd8.rds")

#新增cell_type后续用
tdm$seurat_clusters <- as.factor(tdm$seurat_clusters)
new.cluster.ids <- c("B", "B","Neutrophil","B","CD4+ T","CD4+ T","B",
                     "Red Blood","CD8+ T","B","Red Blood","Nature Killer","Macrophage",
                     "Plasma","Macrophage","Unknown","T",
                     "Unknown","Neutrophil","Red Blood","Unknown")
cluster_mapping <- setNames(new.cluster.ids, seq_along(new.cluster.ids) - 1)
cell_type_vector <- cluster_mapping[as.character(tdm$seurat_clusters)]
cell_type_df <- data.frame(cell_type = cell_type_vector, row.names = colnames(tdm))
tdm <- AddMetaData(tdm, metadata = cell_type_df)

View(tdm@meta.data)

mtdm.object <- subset(tdm, subset = orig.ident == "M-TDM")
ptdm.object <- subset(tdm, subset = orig.ident == "P-TDM")

mtdm.data.input <- GetAssayData(mtdm.object, assay = "RNA", slot = 'data')
mtdm.meta <- mtdm.object@meta.data[,c("cell_type","orig.ident")]
mtdm.meta$CellType %<>% as.vector(.)

#细胞互作
mtdm.cellchat <- createCellChat(object = mtdm.data.input)
mtdm.cellchat <- addMeta(mtdm.cellchat, meta = mtdm.meta)
mtdm.cellchat <- setIdent(mtdm.cellchat, ident.use = "cell_type")

levels(mtdm.cellchat@idents)
mtdmgroupSize <- as.numeric(table(mtdm.cellchat@idents))
mtdmgroupSize

#选择数据库，人 CellChatDB.human
mtdm.cellchat@DB <- CellChatDB.mouse
#showDatabaseCategory(CellChatDB.mouse)
#dplyr::glimpse(CellChatDB.mouse$interaction)

#Secreted Signaling 旁分泌，自分泌； ECM-Receptor 外泌体； Cell-Cell Contact 直接接触； Non-Protein Signaling
#如果要选择
#CellChatDB.use <- subsetDB(CellChatDB.use, search= "Secreted Signaling")
#mtdm.cellchat@DB <- CellChatDB.use

mtdm.cellchat <- subsetData(mtdm.cellchat, features = NULL)
future::plan("multisession", workers = 8) #处理CPU线程
mtdm.cellchat <- identifyOverExpressedGenes(mtdm.cellchat)
mtdm.cellchat <- identifyOverExpressedInteractions(mtdm.cellchat)
mtdm.cellchat <- smoothData(mtdm.cellchat, adj = PPI.mouse)

#cellchat分析
options(future.globals.maxSize = 3 * 1024^3)  # 内存设置为 3GB，不然后面报错
mtdm.cellchat <- computeCommunProb(mtdm.cellchat,raw.use = T)
#raw.use=F 使用预处理后的数据，可能会引进非生物学因素

mtdm.cellchat <- filterCommunication(mtdm.cellchat, min.cells = 10) #通讯细胞过少没意义
mtdm.cellchat <- computeCommunProbPathway(mtdm.cellchat)
mtdm.cellchat <- aggregateNet(mtdm.cellchat)#计算聚合网络
mtdm.cellchat <- netAnalysis_computeCentrality(mtdm.cellchat, slot.name = "netP")

group1.net <- subsetCommunication(mtdm.cellchat)
write.csv(group1.net, file = "./25年/互作/MTDM_net_inter_raw.useT.csv")
saveRDS(mtdm.cellchat,"mtdm.cellchat.rds")

#PTDM
ptdm.data.input <- GetAssayData(ptdm.object, assay = "RNA", slot = 'data')
ptdm.meta <- ptdm.object@meta.data[,c("cell_type","orig.ident")]
ptdm.meta$CellType %<>% as.vector(.)

#细胞互作
ptdm.cellchat <- createCellChat(object = ptdm.data.input)
ptdm.cellchat <- addMeta(ptdm.cellchat, meta = ptdm.meta)
ptdm.cellchat <- setIdent(ptdm.cellchat, ident.use = "cell_type")

levels(ptdm.cellchat@idents)
ptdmgroupSize <- as.numeric(table(ptdm.cellchat@idents))
ptdmgroupSize
ptdm.cellchat@DB <- CellChatDB.mouse

ptdm.cellchat <- subsetData(ptdm.cellchat, features = NULL)
future::plan("multisession", workers = 8) #处理CPU线程
ptdm.cellchat <- identifyOverExpressedGenes(ptdm.cellchat)
ptdm.cellchat <- identifyOverExpressedInteractions(ptdm.cellchat)
ptdm.cellchat <- smoothData(ptdm.cellchat, adj = PPI.mouse)

#cellchat分析
options(future.globals.maxSize = 3 * 1024^3)  # 内存设置为 3GB，不然后面报错
ptdm.cellchat <- computeCommunProb(ptdm.cellchat,raw.use = T)

ptdm.cellchat <- filterCommunication(ptdm.cellchat, min.cells = 10) #通讯细胞过少没意义
ptdm.cellchat <- computeCommunProbPathway(ptdm.cellchat)
ptdm.cellchat <- aggregateNet(ptdm.cellchat)#计算聚合网络
ptdm.cellchat <- netAnalysis_computeCentrality(ptdm.cellchat, slot.name = "netP")

group2.net <- subsetCommunication(ptdm.cellchat)
write.csv(group2.net, file = "./25年/互作/PTDM_net_inter_raw.useT.csv")
saveRDS(ptdm.cellchat,"ptdm.cellchat.rds")


#可视化
pdf("./25年/互作/mtdm.pdf",width = 13, height = 6)
par(mfrow = c(1,2))
netVisual_circle(mtdm.cellchat@net$count,
                 vertex.weight = mtdmgroupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Number of interactions")

netVisual_circle(mtdm.cellchat@net$weight,
                 vertex.weight = mtdmgroupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Interaction strength")

dev.off()

pdf("./25年/互作/ptdm.pdf",width = 13, height = 6)
par(mfrow = c(1,2))
netVisual_circle(ptdm.cellchat@net$count,
                 vertex.weight = ptdmgroupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Number of interactions")

netVisual_circle(ptdm.cellchat@net$weight,
                 vertex.weight = ptdmgroupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Interaction strength")
dev.off()


netVisual_bubble(mtdm.cellchat,
                 targets.use = "CD4+ T",
                 remove.isolate = FALSE,
                 font.size = 14)

netVisual_bubble(ptdm.cellchat,
                 targets.use = "CD4+ T",
                 remove.isolate = FALSE,
                 font.size = 14)


#合并
object.list <- list("M-TDM" = mtdm.cellchat, "P-TDM" = ptdm.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#比较 c(2,1)换顺序，这样红色才是实验组
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(2,1))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(2,1), measure = "weight")
gg1 + gg2
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(2,1))
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(2,1), measure = "weight")
dev.off()

#差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat, comparison = c(2,1))
h2 <- netVisual_heatmap(cellchat, comparison = c(2,1), measure = "weight")
h1+h2

#保守和特异性信号通路
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, comparison = c(2,1))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T, comparison = c(2,1))
gg1 + gg2

#弦图改热图，直观
library(pheatmap)
diff.count <- cellchat@net$`M-TDM`$count - cellchat@net$`P-TDM`$count
pheatmap(diff.count,
         treeheight_col = "0", treeheight_row = "0",
         cluster_rows = T,
         cluster_cols = T)
