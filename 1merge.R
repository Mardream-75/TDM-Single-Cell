BiocManager::install("multtest")
BiocManager::install("Seurat")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("patchwork")
BiocManager::install("tidyverse")

library("multtest")
library("Seurat")
library("dplyr")
library("ggplot2")
library("patchwork")
library("tidyverse")

rm(list = ls())
setwd("C:/Users/Mardream/OneDrive/Desktop/Rdata")
setwd("/Users/mardream_75/Desktop/Rdata")

mtdm.data <- Read10X(data.dir = "./M-TDM/filtered_feature_bc_matrix")
ptdm.data <- Read10X(data.dir = "./P-TDM/filtered_feature_bc_matrix")

#10X的三文件读取
#至少3个细胞有该基因，每个细胞至少有200个基因
mtdm <- CreateSeuratObject(counts = mtdm.data, project = "M-TDM", min.cells =3, min.features = 200)
ptdm <- CreateSeuratObject(counts = ptdm.data, project = "P-TDM", min.cells =3, min.features = 200)

table(mtdm@meta.data$orig.ident)
table(ptdm@meta.data$orig.ident)

#计算线粒体RNA含量，人源换MT
mtdm[["percent.mt"]] <- PercentageFeatureSet(mtdm, pattern = "^mt-")
ptdm[["percent.mt"]] <- PercentageFeatureSet(ptdm, pattern = "^mt-")

#质控降维标准化
mtdm <- subset(mtdm,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
ptdm <- subset(ptdm,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#查看剩余细胞量（across xxx）
mtdm
ptdm

mtdm <- NormalizeData(mtdm, normalization.method = "LogNormalize", scale.factor = 10000)
ptdm <- NormalizeData(ptdm, normalization.method = "LogNormalize", scale.factor = 10000)
mtdm <- FindVariableFeatures(mtdm, selection.method = "vst", nfeatures = 2000)
ptdm <- FindVariableFeatures(ptdm, selection.method = "vst", nfeatures = 2000)

#merge
tdm.anchors <- FindIntegrationAnchors(object.list = list(mtdm,ptdm), dims = 1:20)
tdm.intergrated <- IntegrateData(anchorset = tdm.anchors, dims = 1:20)

#rm(mtdm,mtdm.data,ptdm,ptdm.data)
#pca
tdm <- tdm.intergrated

#DefaultAssay(tdm) <- "integrated"
#应该用RNA
DefaultAssay(tdm) <- "RNA"
tdm <- ScaleData(tdm, features = rownames(tdm))
tdm <- RunPCA(tdm, npcs = 30, verbose = FALSE)
print(tdm[["pca"]], dims = 1:5, nfeatures = 5) 

#Heatmap
VizDimLoadings(tdm, dims = 1:2, reduction = "pca") 
DimHeatmap(tdm, dims = 1, cells = 500, balanced = TRUE) 
DimHeatmap(tdm, dims = 1:15, cells = 500, balanced = TRUE) 

tdm <- JackStraw(tdm, num.replicate = 100)
tdm <- ScoreJackStraw(tdm, dims = 1:20)
JackStrawPlot(tdm, dims = 1:15)

tdm <- FindNeighbors(tdm, reduction = "pca", dims = 1:30)
tdm <- FindClusters(tdm, resolution = 0.5)
tdm <- RunUMAP(tdm, reduction = "pca",dims= 1:30)
tdm <- RunTSNE(tdm, dims= 1:30)

#UMAP
p1 <- DimPlot(tdm, reduction = "umap",label = T, split.by = "orig.ident")
p2 <- DimPlot(tdm, reduction = "umap", group.by = "orig.ident")
p1|p2

#TSNE
p3 <- DimPlot(tdm, reduction = "tsne",label = T, split.by = "orig.ident")
p4 <- DimPlot(tdm, reduction = "tsne", group.by = "orig.ident")
p3|p4

#手动注释
#注释热图
tdm <- JoinLayers(tdm)
tdm.markers <- FindAllMarkers(tdm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
library(dplyr)
top10 <- tdm.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(tdm, features = top10$gene) + NoLegend()

#注释小提琴图
VlnPlot(tdm, features = top10$gene[1:20],pt.size=0)

DimPlot(tdm,label = T,reduction = "tsne")
#手动注释
#http://117.50.127.228/CellMarker/
tdmrename <- tdm
levels(tdmrename)
new.cluster.ids <- c("B", "B","Neutrophil","B","CD4+ T","CD4+ T","B",
                     "Red Blood","CD8+ T","B","Red Blood","Nature Killer","Macrophage",
                     "Plasma","Macrophage"," ","T",
                     " ","Neutrophil","Red Blood"," ")

new.cluster.ids <- c("B", "B","Neutrophil","B","T","T","B",
                     "Red Blood","T","B","Red Blood","Nature Killer","Macrophage",
                     "Plasma","Macrophage"," ","T",
                     " ","Neutrophil","Red Blood"," ")

names(new.cluster.ids) <- levels(tdmrename)

tdmrename <- RenameIdents(tdmrename, new.cluster.ids)
DimPlot(tdmrename, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(tdmrename, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(tdmrename,file="tdm-241219-T.rds")
saveRDS(tdm,file="tdm-241219-original-RNA.rds")

#再细分用
new.cluster.ids <- c("Nature Killer", "B","Neutrophil",
                     "T","T","Red Blood",
                     "UNKNOWN","T","CD8 T",
                     "Macrophage","T","Treg",
                     "T","Plasma","T helper",
                     "Th17")
