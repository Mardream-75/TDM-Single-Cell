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
setwd("C:/Users/Mardream/OneDrive/Desktop/Rdata")

erich2plot <- function(data1, data2, custom_title) {
  library(ggplot2)
  # 处理前一个 data（上调基因）
  data1 <- data1[order(data1$qvalue, decreasing = F)[1:10], ]
  data1$type <- "Upregulated"  # 标识为上调基因
  # 处理后一个 data（下调基因）
  data2 <- data2[order(data2$qvalue, decreasing = F)[1:10], ]
  data2$RichFactor <- -data2$RichFactor
  data2$type <- "Downregulated"  # 标识为下调基因
  data4plot <- rbind(data1, data2)
  # 按 type 和 RichFactor 排序
  data4plot <- data4plot[order(data4plot$type, data4plot$RichFactor, decreasing = FALSE), ]
  # 转换 Description 为因子，以确保 Y 轴顺序正确
  data4plot$Description <- factor(data4plot$Description, levels = unique(data4plot$Description))
  # 绘图
  p <- ggplot(data4plot, aes(RichFactor, Description, color = type))
  p <- p + geom_point()
  pbubble <- p + geom_point(aes(size = Count, color = type))
  pr <- pbubble + 
    scale_size(range = c(1, 7.5)) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "darkgreen"),
      breaks = c("Upregulated","Downregulated")) +
    labs(color = "Gene Type",
         size = "Observed.gene.count",
         x = "Richfactor",
         title = custom_title,
         subtitle = "Enrichment Process")
  pr <- pr + theme_bw()
  pr
}

go2plot <- function(data1, data2, custom_title) {
  library(ggplot2)
  # 处理前一个 data（上调基因）
  data1 <- data1[order(data1$qvalue, decreasing = F), ]
  data1 <- do.call(rbind, by(data1, data1$ONTOLOGY, function(df) {
    df[1:min(5, nrow(df)), ]
  }))
  data1$GeneRatio <- apply(data1, 1, function(x) {
    as.numeric(strsplit(x[4], '/')[[1]][1]) /
      as.numeric(strsplit(x[4], '/')[[1]][2])
  })
  data1$type <- "Upregulated"  # 标识为上调基因
  # 处理后一个 data（下调基因）
  data2 <- data2[order(data2$qvalue, decreasing = F), ] 
  data2 <- do.call(rbind, by(data2, data2$ONTOLOGY, function(df) {
    df[1:min(5, nrow(df)), ]
  }))
  data2$GeneRatio <- apply(data2, 1, function(x) {
    as.numeric(strsplit(x[4], '/')[[1]][1]) /
      as.numeric(strsplit(x[4], '/')[[1]][2])
  })
  data2$GeneRatio <- -data2$GeneRatio
  data2$type <- "Downregulated"  # 标识为下调基因
  data4plot <- rbind(data1, data2)
  #柱状图
  p <- ggplot(data4plot, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = type)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(
      values = c("Upregulated" = "red", "Downregulated" = "darkgreen"),
      breaks = c("Upregulated","Downregulated")) +  # 上调红色，下调蓝色
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +            # 按 ONTOLOGY 分面
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.text.y = element_text(angle = 0),                              
      axis.text.y = element_text(size = 8),                               
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = "grey")
    ) +
    labs(
      x = "Gene Ratio",                                                   
      y = "GO Term Description",                                           
      fill = "Direction",
      title = custom_title
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  #气泡图
  p <- ggplot(data4plot, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = type)) +
    geom_point(alpha = 0.7) +  # 使用geom_point()绘制气泡图
    scale_size(range = c(1, 7.5)) +  # 设置气泡大小范围
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "darkgreen"),
      breaks = c("Upregulated", "Downregulated")
    ) +  
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +            # 按 ONTOLOGY 分面
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.text.y = element_text(angle = 0),                              
      axis.text.y = element_text(size = 8),                               
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = "grey")
    ) +
    labs(
      x = "Gene Ratio",                                                   
      y = "GO Term Description",                                           
      size = "Gene Count",
      color = "Direction",
      title = custom_title
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
  ggsave(paste0('./25年/富集基因/GO/',custom_title,'-dot.pdf'), plot = p, device = "pdf", width = 10, height = 10, dpi = 600, limitsize = FALSE)
}

#上下调各100个基因
mykegg <- function(marker,title){
  marker$foldChange <- 2^marker$avg_log2FC
  upgene <- marker %>% top_n(n = 100, wt = avg_log2FC) %>% rownames()
  dngene <- marker %>% top_n(n = -100, wt = avg_log2FC) %>% rownames()
  foldChange_up <- marker$foldChange[rownames(marker) %in% upgene]
  foldChange_dn <- marker$foldChange[rownames(marker) %in% dngene]
  foldChange_vector <- c(foldChange_up, foldChange_dn)
  names(foldChange_vector) <- c(upgene, dngene)
  upgene.df <- bitr(upgene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
  dngene.df <- bitr(dngene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
  #organism mmu小鼠，hsa人
  #分别富集
  upekegg <- enrichKEGG(unique(upgene.df$ENTREZID), organism = 'mmu',
                      pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,
                      minGSSize = 10,maxGSSize = 500,use_internal_data = F)
  dnekegg <- enrichKEGG(unique(dngene.df$ENTREZID), organism = 'mmu',
                        pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,
                        minGSSize = 10,maxGSSize = 500,use_internal_data = F)
  #ID转基因名
  upekegg <- setReadable(upekegg,'org.Mm.eg.db','ENTREZID')
  dnekegg <- setReadable(dnekegg,'org.Mm.eg.db','ENTREZID')
  #删除物种名
  upekegg@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", upekegg@result$Description)
  dnekegg@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", dnekegg@result$Description)
  write.csv(upekegg@result, paste0('./25年/富集基因/',title,'-up-kegg.result.csv'))
  write.csv(dnekegg@result, paste0('./25年/富集基因/',title,'-down-kegg.result.csv'))
  print(erich2plot(upekegg@result,dnekegg@result,title))
  # 合并上调和下调基因的富集结果
  ekegg_result <- rbind(upekegg@result, dnekegg@result)
  ekegg <- new("enrichResult", result = ekegg_result)
  ekegg@gene <- unique(c(upekegg@gene, dnekegg@gene))
  # 合并基因数据框
  gene.df <- rbind(upgene.df, dngene.df)
  # 提取富集分析中的基因 ID
  common_genes <- intersect(ekegg@gene, gene.df$ENTREZID)
  # 对数处理
  foldChange_vector <- log2(foldChange_vector)
  #自行修改最大值以美化图表
  foldChange_vector <- pmin(foldChange_vector, 2.5)
  cnetplot(ekegg,
           categorySize = "pvalue",
           foldChange = foldChange_vector,
           showCategory = 20,
           max.overlaps = 50) +
    scale_color_gradient2(name = "log2(Fold Change)",
                          low = "green", mid = "grey", high = "red",
                          midpoint = 1, limits = c(0,2.5), oob = scales::oob_keep)
  # 自定义颜色渐变，超出范围的颜色设置为红色
}

mygo <- function(marker,title){
  marker$foldChange <- 2^marker$avg_log2FC
  upgene <- marker %>% top_n(n = 100, wt = avg_log2FC) %>% rownames()
  dngene <- marker %>% top_n(n = -100, wt = avg_log2FC) %>% rownames()
  upgene.df <- bitr(upgene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
  dngene.df <- bitr(dngene,fromType = "SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
  foldChange_up <- marker$foldChange[rownames(marker) %in% upgene]
  foldChange_dn <- marker$foldChange[rownames(marker) %in% dngene]
  foldChange_vector <- c(foldChange_up, foldChange_dn)
  names(foldChange_vector) <- c(upgene, dngene)
  #GO ont: BP生物学过程，MF分子功能，CC细胞组分
  upggoMF <- enrichGO(upgene.df$SYMBOL, org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL",
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.2,
                    minGSSize = 10,maxGSSize = 500,readable = TRUE, pool = FALSE)
  dnggoMF <- enrichGO(dngene.df$SYMBOL, org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL",
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.2,
                    minGSSize = 10,maxGSSize = 500,readable = TRUE, pool = FALSE)
  write.csv(upggoMF@result, paste0('./25年/富集基因/',title,'-up-go.result.csv'))
  write.csv(dnggoMF@result, paste0('./25年/富集基因/',title,'-down-go.result.csv'))
  go2plot(upggoMF@result,dnggoMF@result,title)
  # 合并上调和下调基因的富集结果
  ggoMF_result <- rbind(upggoMF@result, dnggoMF@result)
  ggoMF <- new("enrichResult", result = ggoMF_result)
  ggoMF@gene <- unique(c(upggoMF@gene, dnggoMF@gene))
  # 合并基因数据框
  gene.df <- rbind(upgene.df, dngene.df)
  # 提取富集分析中的基因 ID
  common_genes <- intersect(ggoMF@gene, gene.df$ENTREZID)
  # 对数处理
  #foldChange_vector <- log2(foldChange_vector)
  #自行修改最大值以美化图表
  foldChange_vector <- pmin(foldChange_vector, 2.5)
  p <- cnetplot(ggoMF,
           categorySize = "pvalue",
           foldChange = foldChange_vector,
           showCategory = 15,
           max.overlaps = 50,
           layout = "kk") +
    scale_color_gradient2(name = "Fold Change",
                          low = "green", mid = "grey", high = "red",
                          midpoint = 1, limits = c(0,2.5), oob = scales::oob_keep)
  ggsave(paste0('./25年/富集基因/GO/',title,'-cnetplot.pdf'), plot = p, device = "pdf", width = 10, height = 10, dpi = 600, limitsize = FALSE)
}

#读入文件
cd4.marker <- read.csv('./25年/差异基因/CD4+ T.CSV',row.names = 1)
cd8.marker <- read.csv('./25年/差异基因/CD8+ T.CSV',row.names = 1)
mac.marker <- read.csv('./25年/差异基因/Macrophage.CSV',row.names = 1)
b.marker <- read.csv('./25年/差异基因/B.CSV',row.names = 1)
mykegg(cd4.marker,"CD4+ T")
mykegg(cd8.marker,"CD8+ T")
mykegg(mac.marker,"Macrophage")
mykegg(b.marker,"B")
mygo(cd4.marker,"CD4+ T")
mygo(cd8.marker,"CD8+ T")
mygo(mac.marker,"Macrophage")
mygo(b.marker,"B")
