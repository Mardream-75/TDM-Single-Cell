#GSVA

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(dplyr)
library(enrichplot)

rm(list = ls())

setwd("/Users/mardream_75/Desktop/Rdata")
#读入文件
t.marker <- read.csv('./241220差异基因/T.csv',row.names = 1)
#head(treg.marker)

#自定义可视化函数
erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[5],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[6],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p <- p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") +
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
         x="Richfactor",y="term.description",title="Enrichment Process")
  
  pr <-pr + theme_bw()
  pr
}


#倒数100个基因
mykegg <- function(input){
  inputgene <- input %>% top_n(n = -100, wt = avg_log2FC) %>% rownames()
  gene.df <- bitr(inputgene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
  #organism mmu小鼠，hsa人
  ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = 'mmu',
                      pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,
                      minGSSize = 10,maxGSSize = 500,use_internal_data = F)
  ekegg <- setReadable(ekegg,'org.Mm.eg.db','ENTREZID')
  write.csv(ekegg@result,'./富集/kegg.result.csv')
  erich2plot(ekegg@result)
}

mygo <- function(input){
  inputgene <- input %>% top_n(n = -100, wt = avg_log2FC) %>% rownames()
  gene.df <- bitr(inputgene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
  #GO
  ggoMF <- enrichGO(gene.df$SYMBOL, org.Mm.eg.db, keyType = "SYMBOL", ont = "MF",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                    minGSSize = 10,maxGSSize = 500,readable = FALSE, pool = FALSE)
  write.csv(ggoMF@result,'./富集/GO.MF.result.csv')
  erich2plot(ggoMF@result)
}

#原始代码

tgene <- t.marker %>% top_n(n = -100, wt = avg_log2FC) %>% rownames()
gene.df <- bitr(tgene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)

#View(gene.df)



#KEGG organism mmu小鼠，hsa人
ekegg <- enrichKEGG(gene.df$ENTREZID, organism = 'mmu', keyType = 'kegg',
                    pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,
                    minGSSize = 10,maxGSSize = 500,use_internal_data = F)
#ekegg <- setReadable(ekegg,'org.Mm.eg.db','ENTREZID')
write.csv(ekegg@result,'./241220/富集/T-kegg.result.csv')
#pathway <- ekegg@result 看一下数据，行列与自定义函数要一致
erich2plot(ekegg@result)
pdf(file = "./241219/富集/kegg.pdf", width = 10, height = 10)
dotplot(ekegg,showCategory=20)
barplot(ekegg,showCategory=20,drop=T)
cnetplot(ekegg,categorySize= "pvalue", foldChange=NULL)

#KEGG亚通路图
write.table(ekegg@result$ID, file = "./241220/富集/T-KEGG_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
ekegg@result$ID
browseKEGG(ekegg@result,"mmu04657")

#GO
ggoMF <- enrichGO(gene.df$SYMBOL, org.Mm.eg.db, keyType = "SYMBOL", ont = "MF",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                  minGSSize = 10,maxGSSize = 500,readable = FALSE, pool = FALSE)
write.csv(ggoMF@result,'./241220/富集/T-GO.MF.result.csv')
erich2plot(ggoMF@result)
dotplot(ggoMF,showCategory=20)
barplot(ggoMF,showCategory=20,drop=T)
cnetplot(ggoMF,categorySize= "pvalue", foldChange=NULL)
