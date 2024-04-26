library(readr)
counts <- read_csv("GSE158359-counts.csv", 
                   col_types = cols(Length = col_skip()))

contracts <- read_csv("contracts.csv")

View(counts)
library(tidyverse)
write_csv(counts, "counts.csv")

library(readr)
## ----setup2, message=FALSE, eval=TRUE---------------------------------------------------
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(htmlwidgets)

# #加载parallel包
# library(parallel)
# 
# #detectCores函数可以告诉你你的CPU可使用的核数
# clnum<-detectCores() 
# 
# #设置参与并行的CPU核数目，这里我们使用了所有的CPU核，也就是我们刚才得到的clnum，具体到这个案例，clnum=4
# cl <- makeCluster(getOption("cl.cores", clnum));
# 
# # #关闭并行计算
# # stopCluster(cl);

countdata <- read.csv("counts.csv",header = T,row.names = 1)
# head(countdata)
# ent<- select(Mus.musculus, keys=rownames(countdata), columns=c("SYMBOL"),
#              keytype="ENTREZID")
# merged_table <- merge(ent, countdata, by.x = "ENSEMBL", by.y = "row.names", all.x = TRUE)
# 
# # Find rows with NA values in ENTREZID
# missing_entrezid <- is.na(merged_table$ENTREZID)
# 
# # Remove rows with missing ENTREZID
# merged_table <- merged_table[!missing_entrezid, ]
# 
# # Find duplicate rows based on ENTREZID
# duplicates <- duplicated(merged_table$ENTREZID)
# 
# # Remove duplicate rows
# merged_table <- merged_table[!duplicates, ]
# 
# # Set ENTREZID as row names and remove the ENTREZID column
# rownames(merged_table) <- merged_table$ENTREZID
# merged_table$ENTREZID <- NULL
# merged_table$ENSEMBL <- NULL
# countdata <- merged_table


# # 过滤表达矩阵
# nrow(countdata)
# countdata <- countdata[rowSums(countdata)>10,]
# nrow(countdata)
# countdata <- countdata[rowSums(countdata)>50,]
# nrow(countdata)
#读取样品信息
coldata <- read.csv("contracts.csv",header = T,row.names = 1)
head(coldata)

# 删减样本
# coldata
# str(coldata)
# coldata <- coldata[-c(1),]
# coldata
group = as.factor(coldata$group)

# group = as.factor(coldata)
group
# countdata <- countdata[,-c(1)]

#生成DGElist对象
# # x <- DGEList(counts = countdata,group = as.factor(coldata$group))
# x <- DGEList(counts = countdata,group = group)
# x
# samplenames <- colnames(x)
# geneid <- rownames(x)
# genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM","ENSEMBL"),
#                 keytype="SYMBOL")
# # head(genes)
# ## ----removedups-------------------------------------------------------------------------
# # genes <- genes[!duplicated(genes$ENSEMBL),]
# # # genes <- genes[!duplicated(genes$ENTREZID),]
# summary(genes)
# ## ----assigngeneanno---------------------------------------------------------------------
# x$genes <- genes
# x
# # rownames(x) <- x$genes$ENTREZID
# library(org.Mm.eg.db)
# # View(org.Mm.eg.db)

## ----cpm--------------------------------------------------------------------------------
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE, prior.count=2)
# rpkm <- rpkm(x)
## ----lcpm-------------------------------------------------------------------------------
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

## ----zeroes-----------------------------------------------------------------------------
table(rowSums(x$counts==0)==9)

# ## ----filter-----------------------------------------------------------------------------
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----filterplot1, fig.height=4, fig.width=8, fig.cap="每个样本过滤前的原始数据
# （A）和过滤后（B）的数据的log-CPM值密度。竖直虚线标出了过滤步骤中所用阈值（相当于CPM值为约0.2）。"----
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

## ----normalize--------------------------------------------------------------------------
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
lcpm <- cpm(x, log=TRUE)

boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

## ----MDS1, fig.height=4, fig.width=8, fig.cap="以样品分组上色并标记的log-CPM值在维度1和2的MDS图（A）和以测序泳道上色并标记的维度3和4的MDS图（B）。图中的距离对应于领先倍数变化（leading fold-change），默认情况下也就是前500个在每对样品之间差异最大的基因的平均（均方根）log2倍数变化。"----
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
# col.lane <- lane
# levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
# col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
# plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
# title(main="B. Sequencing lanes")

## ----GlimmaMDSplot----------------------------------------------------------------------
glMDSPlot(lcpm, labels=colnames(x), 
          groups=x$samples[,c(1)], launch=TRUE)

## ----design-----------------------------------------------------------------------------
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

## ----contrasts--------------------------------------------------------------------------

contr.matrix <- makeContrasts(
  treatedvsuntreated = FGF19-Ctrl, 
  levels = colnames(design))
contr.matrix

## ----voom, fig.height=4, fig.width=8, fig.cap="图中绘制了每个基因的均值（x轴）和方差（y轴），显示了在该数据上使用`voom`前它们之间的相关性（左），以及当运用`voom`的权重后这种趋势是如何消除的（右）。左侧的图是使用`voom`函数绘制的，它为log-CPM转换后的数据拟合线性模型并提取残差方差。然后，对方差取四次方根（或对标准差取平方根），并相对每个基因的平均表达作图。均值通过平均计数加上2再进行log2转换计算得到。右侧的图使用`plotSA`绘制了log2残差标准差与log-CPM均值的关系。在这两幅图中，每个黑点表示一个基因。左侧图中，红色曲线展示了用于计算voom权重的估算所得的均值-方差趋势。右侧图中，由经验贝叶斯算法得到的平均log2残差标准差由水平蓝线标出。 "----
par(mfrow=c(1,2))
# v <- voom(x, design, plot=TRUE)
v <- voomWithQualityWeights(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


## ----decidetests------------------------------------------------------------------------
summary(decideTests(efit))


dt <- decideTests(efit)
summary(dt)

# write.fit(efit, dt, file="results-eft.csv",sep =",")
## ----toptables--------------------------------------------------------------------------
treatedvsuntreated <- topTable(efit, coef=1, n=Inf)

write.csv(treatedvsuntreated,"treatedvsuntreated.csv")
# basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(treatedvsuntreated)
# head(basal.vs.ml)

## ----MDplot, fig.keep='none'------------------------------------------------------------
plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], 
       xlim=c(-8,13))

## ----GlimmaMDplot-----------------------------------------------------------------------
glMDPlot(efit, coef=1, status=dt, main=colnames(efit)[1],
         side.main="ENSEMBL", counts=lcpm, groups=group, launch=FALSE)




## ----softwareinfo-----------------------------------------------------------------------
sessionInfo()

treatedvsuntreated$SYMBOL <- rownames(treatedvsuntreated)
write.csv(treatedvsuntreated,"treatedvsuntreated.csv")
#==============================================
#               九十二 基因功能注释           #
#==============================================
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Mm.eg.db)
library(GSEABase)
library(clusterProfiler)
set.seed(1503)
# dta_HFD.Fructose_Chow <- read.csv("HFD.Fructose_Chow.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Glucose_Chow <- read.csv("HFD.Glucose_Chow.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Fructose_HFD.Glucose <- read.csv("HFD.Fructose_HFD.Glucose.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD_Chow <- read.csv("HFD_Chow.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Fructose_HFD <- read.csv("HFD.Fructose_HFD.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子 
# dta_HFD.Glucose_HFD <- read.csv("HFD.Glucose_HFD.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# #==============================================
# 
# ## ------------------------------------------------------------------------
# siggenename <- rownames(dta)
# length(siggenename)
# keytypes(org.Mm.eg.db)
# 
# 
# ids <- bitr(siggenename, fromType="ENSEMBL", toType=c( "SYMBOL"), OrgDb="org.Mm.eg.db")
# head(ids)
# go <-  bitr(siggenename,fromType = "ENSEMBL",toType = c("SYMBOL","GO","ONTOLOGY"),OrgDb = "org.Mm.eg.db")
# head(go)
# 
# eg <-  bitr(siggenename, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
# head(eg)
# #====================================
# # GO功能富集分析，输入Gene ID 为GI号 #
# #====================================
# gene <- eg$ENTREZID
# gene.df <- bitr(gene, fromType = "ENTREZID",
#                 toType = c("ENSEMBL", "SYMBOL"),
#                 OrgDb = org.Mm.eg.db)
# head(gene.df)
# ggo <- groupGO(gene     = gene,
#                OrgDb    = org.Mm.eg.db,
#                ont      = "MF",
#                level    = 3,
#                readable = TRUE)
# head(ggo)
# 
# # View(as.data.frame(ggo))
# ego <- enrichGO(gene          = gene,
#                 OrgDb         = org.Mm.eg.db,
#                 ont           = "CC",
#                 pAdjustMethod = "BH",
#                 readable      = TRUE)
# head(ego)
# View(as.data.frame(ego))
# ## GO功能富集可视化
# barplot(ggo, drop=TRUE, showCategory=12)
# ## ----fig.height=5, fig.width=8-------------------------------------------
# 
# ego
# #====================================
# #           KEGG富集分析            #
# #====================================
# search_kegg_organism('mmu', by='kegg_code')
# 
# # KEGG pathway over-representation analysis over-representation analysis (ORA) 
# kk <- enrichKEGG(gene         = gene,
#                  organism     = "mmu",
#                  pvalueCutoff = 0.05)
# head(kk)
# kk@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kk@result$Description)
# 
# View(kk@result)
# barplot(kk, showCategory=20)
# dotplot(kk,showCategory=20) #这里可以画出金葡菌感染来
# write_csv(kk@result,"HFD.Fructose_Chow_kegg_result_ORA.csv")
# 



# KEGG pathway gene set enrichment analysis (GSEA)
#==============================================
#               GSEA  分析                     #
#==============================================
#====================================
#             方法1  GSEAbase       #
#====================================
#BiocManager::install("GSEABase")
library(GSEABase)
library(clusterProfiler)
library(org.Mm.eg.db)
#calculate ranks
# #根据上面的几个数据分别进行gsea替换下面的res
# dta_HFD.Fructose_Chow <- read.csv("HFD.Fructose_Chow.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Glucose_Chow <- read.csv("HFD.Glucose_Chow.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Fructose_HFD.Glucose <- read.csv("HFD.Fructose_HFD.Glucose.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD_Chow <- read.csv("HFD_Chow.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Fructose_HFD <- read.csv("HFD.Fructose_HFD.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子
# dta_HFD.Glucose_HFD <- read.csv("HFD.Glucose_HFD.csv",header = T,row.names = 1,stringsAsFactors = F)#字符串不要转化为因子

gsea_analysis <- function(data_file) {
  dta <- read.csv(data_file, header = T, row.names = 1, stringsAsFactors = FALSE)
  res <- dta
  
  ranks_RNAseq <- sign(res$logFC) * -log10(res$P.Value)
  
  genenames <- res$SYMBOL
  
  ranks_RNAseq <- setNames(object = ranks_RNAseq, nm = genenames)
  
  ranks_RNAseq <- sort(ranks_RNAseq, decreasing = TRUE)
  
  library(homologene)
  genelist <- names(ranks_RNAseq)
  homologues <- mouse2human(genelist)
  
  ranks_RNAseq <- ranks_RNAseq[names(ranks_RNAseq) %in% homologues$mouseGene]
  names(ranks_RNAseq) <- homologues$humanGene[match(names(ranks_RNAseq), homologues$mouseGene)]
  
  output_file <- paste0("RNASeq_ranks_", basename(data_file), ".rnk")
  write.table(data.frame(GeneName = names(ranks_RNAseq), rank = ranks_RNAseq),
              file.path(output_file), col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
  
  assign(paste0("ranks_RNAseq_", tools::file_path_sans_ext(basename(data_file))), ranks_RNAseq, envir = .GlobalEnv)
}
gsea_analysis("treatedvsuntreated.csv")



gsea_analysis_FC <- function(data_file) {
  dta <- read.csv(data_file, header = T, row.names = 1, stringsAsFactors = FALSE)
  res <- dta
  
  ranks_RNAseq <- res$logFC
  
  genenames <- res$SYMBOL
  
  ranks_RNAseq <- setNames(object = ranks_RNAseq, nm = genenames)
  
  ranks_RNAseq <- sort(ranks_RNAseq, decreasing = TRUE)
  
  library(homologene)
  genelist <- names(ranks_RNAseq)
  homologues <- mouse2human(genelist)
  
  ranks_RNAseq <- ranks_RNAseq[names(ranks_RNAseq) %in% homologues$mouseGene]
  names(ranks_RNAseq) <- homologues$humanGene[match(names(ranks_RNAseq), homologues$mouseGene)]
  
  output_file <- paste0("RNASeq_ranks_FC_", basename(data_file), ".rnk")
  write.table(data.frame(GeneName = names(ranks_RNAseq), rank = ranks_RNAseq),
              file.path(output_file), col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
  
  assign(paste0("ranks_RNAseq_FC_", tools::file_path_sans_ext(basename(data_file))), ranks_RNAseq, envir = .GlobalEnv)
}

gsea_analysis_FC("treatedvsuntreated.csv")

# ranks_RNAseq_FC_HFD_Chow


# # dta_HFD.Fructose_Chow
# res <- dta_HFD.Fructose_Chow #这里替换不同的数据 看看把下面包装为函数
# ranks_RNAseq_HFD.Fructose_Chow = sign(res$logFC) * -log10(res$P.Value)
# 
# # ranks_RNAseq = res$logFC
# #gene names from the TCGA set contain gene name and entrez gene ids separated by ‘|’
# # for all subsequent enrichment analysis we need to have just one id.  Separate the names
# # into their two ids.
# genenames <- res$SYMBOL
# 
# # Create ranks vector
# ranks_RNAseq <- setNames(object = ranks_RNAseq, nm = genenames)
# 
# #sort ranks in decreasing order
# ranks_RNAseq <- sort(ranks_RNAseq, decreasing = TRUE)
# 
# # 安装和加载 homologene 包
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# BiocManager::install("homologene")
# library(homologene)
# genelist <- names(ranks_RNAseq)
# homologues <- mouse2human(genelist)
# 
# # 将同源基因信息整合到 ranks_RNAseq 中
# ranks_RNAseq <- ranks_RNAseq[names(ranks_RNAseq) %in% homologues$mouseGene]
# names(ranks_RNAseq) <- homologues$humanGene[match(names(ranks_RNAseq), homologues$mouseGene)]
# 
# # 将 ranks_RNAseq 保存到文件中
# write.table(data.frame(GeneName = names(ranks_RNAseq), rank = ranks_RNAseq), file.path("RNASeq_ranks.rnk"),
#             col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
# 
# 
# # data(geneList, package="DOSE")
# 
# ## ----eval = FALSE--------------------------------------------------------
# 
# ranks_RNAseq_ID = sign(res$logFC) * -log10(res$PValue)
# # ranks_RNAseq = res$logFC
# #gene names from the TCGA set contain gene name and entrez gene ids separated by ‘|’
# # for all subsequent enrichment analysis we need to have just one id.  Separate the names 
# # into their two ids.
# geneID <- res$ENTREZID
# 
# # Create ranks vector
# ranks_RNAseq_ID <- setNames(object = ranks_RNAseq_ID, nm = geneID)
# 
# #sort ranks in decreasing order
# ranks_RNAseq_ID <- sort(ranks_RNAseq_ID, decreasing = TRUE)
# 
# library(clusterProfiler)
# kk2 <- gseKEGG(geneList     = ranks_RNAseq_ID,
#                organism     = 'mmu',
#                minGSSize    = 1,#调小一点才有jakstat3
#                pvalueCutoff = 0.05)
# ??gseKEGG
# 
# 
# head(kk2)
# 

#读取gmt文件
# gmtfile <- "~/bioinfo/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt"
# 
# 
# gmtfile <- "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM.v2023.1.Hs.gmt"
# geneset <- read.gmt(gmtfile)

#GSEA分析


# # c2kegg
# gmtfile <- "~/work/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt"
# geneset <- read.gmt(gmtfile)
# egmt_c2kegg_HFD_Chow <- GSEA(ranks_RNAseq_HFD_Chow, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.05, verbose = F)
# egmt_c2kegg_HFD_Fructose_Chow <- GSEA(ranks_RNAseq_HFD.Fructose_Chow, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.05, verbose = F)
# egmt_c2kegg_HFD_Glucose_Chow <- GSEA(ranks_RNAseq_HFD.Glucose_Chow, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.05, verbose = F)
# egmt_c2kegg_HFD_Fructose_HFD.Glucose <- GSEA(ranks_RNAseq_HFD.Fructose_HFD.Glucose, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.05, verbose = F)
# egmt_c2kegg_HFD_Fructose_HFD <- GSEA(ranks_RNAseq_HFD.Fructose_HFD, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.05, verbose = F)
# egmt_c2kegg_HFD_Glucose_HFD <- GSEA(ranks_RNAseq_HFD.Glucose_HFD, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.05, verbose = F)
# 
# 
# # c5
# gmtfile <- "~/work/msigdb_v2023.1.Hs_GMTs/c5.go.cc.v2023.1.Hs.symbols.gmt"
# geneset <- read.gmt(gmtfile)
# 
# egmt_c5cc_HFD_Chow <- GSEA(ranks_RNAseq_HFD_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5cc_HFD_Fructose_Chow <- GSEA(ranks_RNAseq_HFD.Fructose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5cc_HFD_Glucose_Chow <- GSEA(ranks_RNAseq_HFD.Glucose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5cc_HFD_Fructose_HFD.Glucose <- GSEA(ranks_RNAseq_HFD.Fructose_HFD.Glucose, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5cc_HFD_Fructose_HFD <- GSEA(ranks_RNAseq_HFD.Fructose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5cc_HFD_Glucose_HFD <- GSEA(ranks_RNAseq_HFD.Glucose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# 
# # C5gomf
# gmtfile <- "~/work/msigdb_v2023.1.Hs_GMTs/c5.go.mf.v2023.1.Hs.symbols.gmt"
# geneset <- read.gmt(gmtfile)
# egmt_c5mf_HFD_Chow <- GSEA(ranks_RNAseq_HFD_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5mf_HFD_Fructose_Chow <- GSEA(ranks_RNAseq_HFD.Fructose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5mf_HFD_Glucose_Chow <- GSEA(ranks_RNAseq_HFD.Glucose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5mf_HFD_Fructose_HFD.Glucose <- GSEA(ranks_RNAseq_HFD.Fructose_HFD.Glucose, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5mf_HFD_Fructose_HFD <- GSEA(ranks_RNAseq_HFD.Fructose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5mf_HFD_Glucose_HFD <- GSEA(ranks_RNAseq_HFD.Glucose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# 
# # c5bp
# gmtfile <- "~/work/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt"
# geneset <- read.gmt(gmtfile)
# egmt_c5bp_HFD_Chow <- GSEA(ranks_RNAseq_HFD_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5bp_HFD_Fructose_Chow <- GSEA(ranks_RNAseq_HFD.Fructose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5bp_HFD_Glucose_Chow <- GSEA(ranks_RNAseq_HFD.Glucose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5bp_HFD_Fructose_HFD.Glucose <- GSEA(ranks_RNAseq_HFD.Fructose_HFD.Glucose, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5bp_HFD_Fructose_HFD <- GSEA(ranks_RNAseq_HFD.Fructose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c5bp_HFD_Glucose_HFD <- GSEA(ranks_RNAseq_HFD.Glucose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# 
# # c7
# gmtfile <- "~/work/msigdb_v2023.1.Hs_GMTs/c7.immunesigdb.v2023.1.Hs.symbols.gmt"
# geneset <- read.gmt(gmtfile)
# egmt_c7immu_HFD_Chow <- GSEA(ranks_RNAseq_HFD_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c7immu_HFD_Fructose_Chow <- GSEA(ranks_RNAseq_HFD.Fructose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c7immu_HFD_Glucose_Chow <- GSEA(ranks_RNAseq_HFD.Glucose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c7immu_HFD_Fructose_HFD.Glucose <- GSEA(ranks_RNAseq_HFD.Fructose_HFD.Glucose, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c7immu_HFD_Fructose_HFD <- GSEA(ranks_RNAseq_HFD.Fructose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_c7immu_HFD_Glucose_HFD <- GSEA(ranks_RNAseq_HFD.Glucose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# 
# # all
# gmtfile <- "~/work/msigdb_v2023.1.Hs_GMTs/msigdb.v2023.1.Hs.symbols.gmt"
# geneset <- read.gmt(gmtfile)
# egmt_all_HFD_Chow <- GSEA(ranks_RNAseq_HFD_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_all_HFD_Fructose_Chow <- GSEA(ranks_RNAseq_HFD.Fructose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_all_HFD_Glucose_Chow <- GSEA(ranks_RNAseq_HFD.Glucose_Chow, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_all_HFD_Fructose_HFD.Glucose <- GSEA(ranks_RNAseq_HFD.Fructose_HFD.Glucose, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_all_HFD_Fructose_HFD <- GSEA(ranks_RNAseq_HFD.Fructose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# egmt_all_HFD_Glucose_HFD <- GSEA(ranks_RNAseq_HFD.Glucose_HFD, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05, verbose = F)
# 


# gmtfile <- "genesets.v2023.1.Hs-M2-KEYWORD-human.gmt"
# geneset <- read.gmt(gmtfile)
# 
# #GSEA分析
# egmt_M2 <- GSEA(ranks_RNAseq, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.05)
# # 整合包画STAT3集合=====

# 重构================

# 定义一个读取GMT文件并进行GSEA分析的函数
perform_GSEA <- function(gmt_path, prefix, min_size, pval_cutoff, ranks_list) {
  geneset <- read.gmt(gmt_path)
  results <- list()
  
  for (ranks_name in ranks_list) {
    set.seed(1503)
    ranks <- get(ranks_name)
    result_name <- paste0("egmt_", prefix, "_", sub("ranks_RNAseq_", "", ranks_name))
    results[[result_name]] <- GSEA(ranks, TERM2GENE = geneset, minGSSize = min_size, pvalueCutoff = pval_cutoff, verbose = F,eps = 0.0,seed = TRUE)
  }
  
  return(results)
}

??GSEA


# 定义一个将结果保存为ridgeplot的函数
# save_ridgeplots <- function(results_list, category_count = 15, output_dir = "~/ridgeplots") {
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir)
#   }
#   library(enrichplot)
#   library(ggplot2)
#   # for (result_name in names(results_list)) {
#     plot <- ridgeplot(results_list[[result_name]], showCategory = category_count) + ggtitle(paste("ridgeplot for", result_name))
#     ggsave(file.path(output_dir, paste(result_name, "pdf", sep = ".")), plot, width = 10, height = 10, units = "in", dpi = 1200)
#   }
# }

# 设置GMT文件和分析参数
gmt_files <- list(
  c2kegg = "~/work/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
  c5cc = "~/work/msigdb_v2023.1.Hs_GMTs/c5.go.cc.v2023.1.Hs.symbols.gmt",
  c5mf = "~/work/msigdb_v2023.1.Hs_GMTs/c5.go.mf.v2023.1.Hs.symbols.gmt",
  c5bp = "~/work/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt",
  c7immu = "~/work/msigdb_v2023.1.Hs_GMTs/c7.immunesigdb.v2023.1.Hs.symbols.gmt",
  all = "~/work/msigdb_v2023.1.Hs_GMTs/msigdb.v2023.1.Hs.symbols.gmt"
)

# 设置GSEA分析参数
ranks_names <- ls(pattern = "^ranks_RNAseq_")

# 进行GSEA分析并保存结果
all_results <- list()
for (prefix in names(gmt_files)) {
  gmt_file <- gmt_files[[prefix]]
  min_size <- if (prefix == "c2kegg") 1 else 10
  results <- perform_GSEA(gmt_file, prefix, min_size, 0.05, ranks_names)
  all_results <- c(all_results, results)
}
all_results

lapply(all_results, function(x) {
  if (is.null(x)) {
    "NULL"
  } else {
    length(x)
  }
})
# save_ridgeplots <- function(results_list, category_count = 15, output_dir = "~/ridgeplots") {
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir)
#   }
#   library(enrichplot)
#   library(ggplot2)
#   
#   for (result_name in names(results_list)) {
#     if (!is.null(results_list[[result_name]]) && length(results_list[[result_name]]) > 0) {
#       plot <- ridgeplot(results_list[[result_name]], showCategory = category_count) + ggtitle(paste("ridgeplot for", result_name))
#       ggsave(file.path(output_dir, paste(result_name, "pdf", sep = ".")), plot, width = 10, height = 10, units = "in", dpi = 1200)
#     } else {
#       message(paste("No data available for", result_name, "- skipping..."))
#     }
#   }
# }
# 
# # 保存ridgeplots
# save_ridgeplots(all_results)



# ===============原始画图代码================
#绘制GSEA图
library(enrichplot)
library(ggplot2)
# 
# # HFD_Fructose_HFD
# ridgeplot(egmt_c2kegg_HFD_Fructose_HFD, showCategory=15) + ggtitle("ridgeplot for c2kegg_HFD_Fructose_HFD")
# ggsave("ridgeplot_c2kegg_HFD_Fructose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5mf_HFD_Fructose_HFD, showCategory=15) + ggtitle("ridgeplot for c5mf_HFD_Fructose_HFD")
# ggsave("ridgeplot_c5mf_HFD_Fructose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5bp_HFD_Fructose_HFD, showCategory=15) + ggtitle("ridgeplot for c5bp_HFD_Fructose_HFD")
# ggsave("ridgeplot_c5bp_HFD_Fructose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5cc_HFD_Fructose_HFD, showCategory=15) + ggtitle("ridgeplot for c5cc_HFD_Fructose_HFD")
# ggsave("ridgeplot_c5cc_HFD_Fructose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c7immu_HFD_Fructose_HFD, showCategory=15) + ggtitle("ridgeplot for c7immu_HFD_Fructose_HFD")
# ggsave("ridgeplot_c7immu_HFD_Fructose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# 
# ridgeplot(egmt_c2kegg_HFD_Chow, showCategory=15) + ggtitle("ridgeplot for c2kegg_HFD_Chow")
# ggsave("ridgeplot_c2kegg_HFD_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5mf_HFD_Chow, showCategory=15) + ggtitle("ridgeplot for c5mf_HFD_Chow")
# ggsave("ridgeplot_c5mf_HFD_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5bp_HFD_Chow, showCategory=15) + ggtitle("ridgeplot for c5bp_HFD_Chow")
# ggsave("ridgeplot_c5bp_HFD_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5cc_HFD_Chow, showCategory=15) + ggtitle("ridgeplot for c5cc_HFD_Chow")
# ggsave("ridgeplot_c5cc_HFD_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c7immu_HFD_Chow, showCategory=15) + ggtitle("ridgeplot for c7immu_HFD_Chow")
# ggsave("ridgeplot_c7immu_HFD_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# 
# # HFD_Glucose_HFD
# ridgeplot(egmt_c2kegg_HFD_Glucose_HFD, showCategory=15) + ggtitle("ridgeplot for c2kegg_HFD_Glucose_HFD")
# ggsave("ridgeplot_c2kegg_HFD_Glucose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5mf_HFD_Glucose_HFD, showCategory=15) + ggtitle("ridgeplot for c5mf_HFD_Glucose_HFD")
# ggsave("ridgeplot_c5mf_HFD_Glucose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5bp_HFD_Glucose_HFD, showCategory=15) + ggtitle("ridgeplot for c5bp_HFD_Glucose_HFD")
# ggsave("ridgeplot_c5bp_HFD_Glucose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5cc_HFD_Glucose_HFD, showCategory=15) + ggtitle("ridgeplot for c5cc_HFD_Glucose_HFD")
# ggsave("ridgeplot_c5cc_HFD_Glucose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c7immu_HFD_Glucose_HFD, showCategory=15) + ggtitle("ridgeplot for c7immu_HFD_Glucose_HFD")
# ggsave("ridgeplot_c7immu_HFD_Glucose_HFD.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# 
# # HFD_Glucose_Chow
# ridgeplot(egmt_c2kegg_HFD_Glucose_Chow, showCategory=15) + ggtitle("ridgeplot for c2kegg_HFD_Glucose_Chow")
# ggsave("ridgeplot_c2kegg_HFD_Glucose_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5mf_HFD_Glucose_Chow, showCategory=15) + ggtitle("ridgeplot for c5mf_HFD_Glucose_Chow")
# ggsave("ridgeplot_c5mf_HFD_Glucose_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5bp_HFD_Glucose_Chow, showCategory=15) + ggtitle("ridgeplot for c5bp_HFD_Glucose_Chow")
# ggsave("ridgeplot_c5bp_HFD_Glucose_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5cc_HFD_Glucose_Chow, showCategory=15) + ggtitle("ridgeplot for c5cc_HFD_Glucose_Chow")
# ggsave("ridgeplot_c5cc_HFD_Glucose_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c7immu_HFD_Glucose_Chow, showCategory=15) + ggtitle("ridgeplot for c7immu_HFD_Glucose_Chow")
# ggsave("ridgeplot_c7immu_HFD_Glucose_Chow.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# 
# # HFD_Fructose_HFD_Glucose
# ridgeplot(egmt_c2kegg_HFD_Fructose_HFD.Glucose, showCategory=15) + ggtitle("ridgeplot for c2kegg_HFD_Fructose_HFD.Glucose")
# ggsave("ridgeplot_c2kegg_HFD_Fructose_HFD.Glucose.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5mf_HFD_Fructose_HFD.Glucose, showCategory=15) + ggtitle("ridgeplot for c5mf_HFD_Fructose_HFD.Glucose")
# ggsave("ridgeplot_c5mf_HFD_Fructose_HFD.Glucose.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5bp_HFD_Fructose_HFD.Glucose, showCategory=15) + ggtitle("ridgeplot for c5bp_HFD_Fructose_HFD.Glucose")
# ggsave("ridgeplot_c5bp_HFD_Fructose_HFD.Glucose.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c5cc_HFD_Fructose_HFD.Glucose, showCategory=15) + ggtitle("ridgeplot for c5cc_HFD_Fructose_HFD.Glucose")
# ggsave("ridgeplot_c5cc_HFD_Fructose_HFD.Glucose.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# ridgeplot(egmt_c7immu_HFD_Fructose_HFD.Glucose, showCategory=15) + ggtitle("ridgeplot for c7immu_HFD_Fructose_HFD.Glucose")
# ggsave("ridgeplot_c7immu_HFD_Fructose_HFD.Glucose.pdf", width = 10, height = 10, units = "in", dpi = 1200)
# 
# 



# 导出gsearesult===================================================================================================
# Assign eac  h element of all_results to the global environment
list2env(all_results, envir = .GlobalEnv)



library(enrichplot)
library(ggplot2)

# 定义保存 ridgeplot 图表的函数
save_ridgeplots <- function(results_list, category_count = 15, output_dir = "ridgeplots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  for (result_name in names(results_list)) {
    result <- results_list[[result_name]]
    # 检查结果是否为空或NULL，并且是否有行
    if (!is.null(result)  && nrow(result@result) > 0) {
      plot <- ridgeplot(result, showCategory = category_count) + ggtitle(paste("ridgeplot for", result_name))
      ggsave(file.path(output_dir, paste(result_name, "pdf", sep = ".")), plot, width = 10, height = 10, units = "in", dpi = 300)
    } else {
      message(paste("Skipping", result_name, "because it's NULL, not a GSEAresult, or has no rows."))
    }
  }
}
# 调用 save_ridgeplots 函数来保存 all_results 中的图表
save_ridgeplots(all_results)


# Now all elements of all_results are available as individual variables in the global environment
# 保存csv和rds文件===================================================================================================
# # Get the variables starting with "egmt"
# egmt_variables <- ls(pattern = "^egmt")
# 
# # Define the output directory
# output_directory <- "gsearesult"
# # Check if the output directory exists, if not, create it
# if (!dir.exists(output_directory)) {
#   dir.create(output_directory)
# }
# # Loop through the egmt variables and write CSV and RDS files
# for (variable in egmt_variables) {
#   # Access the 'result' slot of the S4 object
#   result <- get(variable)@result
# 
#   # Build the output filenames
#   output_filename_csv <- paste0(output_directory, "/", variable, ".csv")
#   output_filename_rds <- paste0(output_directory, "/", variable, ".rds")
# 
#   # Write the CSV file
#   write.csv(result, file = output_filename_csv)
# 
#   # Save the object as an RDS file
#   saveRDS(result, file = output_filename_rds)
# }
# 

# 定义输出目录
output_directory <- "gsearesult"
# 检查输出目录是否存在，如果不存在则创建
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# 循环遍历 all_results 中的每个 GSEAresult 对象
for (result_name in names(all_results)) {
  result <- all_results[[result_name]]
  
  # 检查结果是否有数据
  if (!is.null(result) && class(result)[1] == "gseaResult" && nrow(result@result) > 0) {
    # 构建输出文件名
    output_filename_csv <- paste0(output_directory, "/", result_name, ".csv")
    output_filename_rds <- paste0(output_directory, "/", result_name, ".rds")
    
    # 写入 CSV 文件
    write.csv(result@result, file = output_filename_csv)
    
    # 保存为 RDS 文件
    saveRDS(result, file = output_filename_rds)
    
  }
}


# egmt_all_HFD_Chow@result===================================================================================================
# write.csv(egmt_all_HFD_Chow@result, file = "egmt_all_HFD_Chow.csv")
# write.csv(egmt_all_HFD_Glucose_Chow@result, file = "egmt_all_HFD_Glucose_Chow.csv")
# write.csv(egmt_all_HFD_Fructose_HFD@result, file = "egmt_all_HFD_Fructose_HFD.csv")
# write.csv(egmt_all_HFD_Glucose_HFD@result, file = "egmt_all_HFD_Glucose_HFD.csv")
# write.csv(egmt_all_HFD_Fructose_HFD.Glucose@result, file = "egmt_all_HFD_Fructose_HFD.Glucose.csv")
# write.csv(egmt_all_HFD_Fructose_Chow@result, file = "egmt_all_HFD_Fructose_Chow.csv")
# 
# write.csv(egmt_c2kegg_HFD_Chow@result, file = "egmt_c2kegg_HFD_Chow.csv")
# write.csv(egmt_c2kegg_HFD_Glucose_Chow@result, file = "egmt_c2kegg_HFD_Glucose_Chow.csv")
# write.csv(egmt_c2kegg_HFD_Fructose_HFD@result, file = "egmt_c2kegg_HFD_Fructose_HFD.csv")
# write.csv(egmt_c2kegg_HFD_Glucose_HFD@result, file = "egmt_c2kegg_HFD_Glucose_HFD.csv")
# write.csv(egmt_c2kegg_HFD_Fructose_HFD.Glucose@result, file = "egmt_c2kegg_HFD_Fructose_HFD.Glucose.csv")
# write.csv(egmt_c2kegg_HFD_Fructose_Chow@result, file = "egmt_c2kegg_HFD_Fructose_Chow.csv")
# 
# write.csv(egmt_c5bp_HFD_Chow@result, file = "egmt_c5bp_HFD_Chow.csv")
# write.csv(egmt_c5bp_HFD_Glucose_Chow@result, file = "egmt_c5bp_HFD_Glucose_Chow.csv")
# write.csv(egmt_c5bp_HFD_Fructose_HFD@result, file = "egmt_c5bp_HFD_Fructose_HFD.csv")
# write.csv(egmt_c5bp_HFD_Glucose_HFD@result, file = "egmt_c5bp_HFD_Glucose_HFD.csv")
# write.csv(egmt_c5bp_HFD_Fructose_HFD.Glucose@result, file = "egmt_c5bp_HFD_Fructose_HFD.Glucose.csv")
# write.csv(egmt_c5bp_HFD_Fructose_Chow@result, file = "egmt_c5bp_HFD_Fructose_Chow.csv")
# 
# write.csv(egmt_c5cc_HFD_Chow@result, file = "egmt_c5cc_HFD_Chow.csv")
# write.csv(egmt_c5cc_HFD_Glucose_Chow@result, file = "egmt_c5cc_HFD_Glucose_Chow.csv")
# write.csv(egmt_c5cc_HFD_Fructose_HFD@result, file = "egmt_c5cc_HFD_Fructose_HFD.csv")
# write.csv(egmt_c5cc_HFD_Glucose_HFD@result, file = "egmt_c5cc_HFD_Glucose_HFD.csv")
# write.csv(egmt_c5cc_HFD_Fructose_HFD.Glucose@result, file = "egmt_c5cc_HFD_Fructose_HFD.Glucose.csv")
# wrrote.csv(eg,t_c5cc_HFD_Fructose_Chow@result, file = "egmt_c5cc_HFD_Fructose_Chow.csv")
# 
# write.csv(egmt_c5mf_HFD_Chow@result, file = "egmt_c5mf_HFD_Chow.csv")
# write.csv(egmt_c5mf_HFD_Glucose_Chow@result, file = "egmt_c5mf_HFD_Glucose_Chow.csv")
# write.csv(egmt_c5mf_HFD_Fructose_HFD@result, file = "egmt_c5mf_HFD_Fructose_HFD.csv")
# write.csv(egmt_c5mf_HFD_Glucose_HFD@result, file = "egmt_c5mf_HFD_Glucose_HFD.csv")
# write.csv(egmt_c5mf_HFD_Fructose_HFD.Glucose@result, file = "egmt_c5mf_HFD_Fructose_HFD.Glucose.csv")
# write.csv(egmt_c5mf_HFD_Fructose_Chow@result, file = "egmt_c5mf_HFD_Fructose_Chow.csv")
# 
# write.csv(egmt_c7immu_HFD_Chow@result, file = "egmt_c7immu_HFD_Chow.csv")
# write.csv(egmt_c7immu_HFD_Glucose_Chow@result, file = "egmt_c7immu_HFD_Glucose_Chow.csv")
# write.csv(egmt_c7immu_HFD_Fructose_HFD@result, file = "egmt_c7immu_HFD_Fructose_HFD.csv")
# write.csv(egmt_c7immu_HFD_Glucose_HFD@result, file = "egmt_c7immu_HFD_Glucose_HFD.csv")
# write.csv(egmt_c7immu_HFD_Fructose_HFD.Glucose@result, file = "egmt_c7immu_HFD_Fructose_HFD.Glucose.csv")
# write.csv(egmt_c7immu_HFD_Fructose_Chow@result, file = "egmt_c7immu_HFD_Fructose_Chow.csv")


# 点图=======================
dotplot(egmt_c2kegg_HFD_Fructose_HFD, showCategory=30) + ggtitle("dotplot for GSEA")
## convert gene ID to Symbol

p1 <- cnetplot(egmt_c2kegg_HFD_Fructose_HFD, color.params = list(foldChange= ranks_RNAseq_FC_HFD.Fructose_HFD))

HFD.Fructose_HFD
# ## convert gene ID to Symbol
# edox <- setReadable(egmt_c2kegg_HFD_Fructose_HFD, 'org.Mm.eg.db',keyType = "SYMBOL")
# 
# p1 <- cnetplot(edox, foldChange=)
# ## categorySize can be scaled by 'pvalue' or 'geneNum'
# p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
# p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

p2 <- cnetplot(egmt_c2kegg_HFD_Fructose_HFD, categorySize="pvalue", color.params = list(foldChange= ranks_RNAseq_FC_HFD.Fructose_HFD))
p2
p3 <- cnetplot(egmt_c2kegg_HFD_Fructose_HFD, foldChange=ranks_RNAseq_FC_HFD.Fructose_HFD, circular = TRUE, colorEdge = TRUE)
p3
# cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))



p1 <- heatplot(egmt_c2kegg_HFD_Fructose_HFD, showCategory=5)
p2 <- heatplot(egmt_c2kegg_HFD_Fructose_HFD, foldChange=ranks_RNAseq_FC_HFD.Fructose_HFD, showCategory=5)
heatplot(egmt_c2kegg_HFD_Fructose_HFD, foldChange=ranks_RNAseq_FC_HFD.Fructose_HFD, showCategory=20)

cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

edox2 <- pairwise_termsim(egmt_c2kegg_HFD_Fructose_HFD)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

options(repr.plot.width=6,repr.plot.height=4)

# ===============画图gesa plot2=================

??gseaNb()

library(GseaVis)
library(ggthemes)
library(ggplot2)
View(egmt@result)
# devtools::install_github("junjunlab/GseaVis")
??GseaVis
load("./rdata/2024年4月26日17点06分.RData")
library(httpgd)
hgd()
hgd_view()
hgd_browse()
dev.list()
dev.off()

gseaNb(object = egmt_c5bp_treatedvsuntreated,geneSetID = 'GOBP_MACROPHAGE_ACTIVATION',subPlot = 2,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.8)
HALLMARK_IL6_JAK_STAT3_SIGNALING

gseaNb(object = egmt_all_treatedvsuntreated,geneSetID = 'HALLMARK_IL6_JAK_STAT3_SIGNALING',subPlot = 2,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.8)

gseaplot2(egmt_all_HFD_Fructose_Chow, geneSetID = "GOBP_RESPONSE_TO_BACTERIUM", pvalue_table = F, base_size = 14, color = "black", 
          title = "GOBP_RESPONSE_TO_BACTERIUM", subplots = c(1,2))

gseaNb(object = all_results$egmt_c2kegg_HFD.Fructose_HFD,
       geneSetID = 'KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION',subPlot = 2,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.8)
??gseaNb
gseaplot2(all_results$egmt_c2kegg_HFD.Fructose_HFD,
          geneSetID = 'KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION',, pvalue_table = F, base_size = 14, color = "black", 
          title = "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION", subplots = c(1,2))

# GOBP_RESPONSE_TO_BACTERIUM
gseaNb(object = egmt_all_HFD_Fructose_Chow,
       geneSetID = 'GOBP_RESPONSE_TO_BACTERIUM',subPlot = 2,curveCol = colors_transparent,
       addPval = TRUE,pDigit = 3)

colors_transparent <- c("#0000FF99", "#FF000099","#32323299")

# 打印所选颜色的 HEX 值
print(colors_transparent)

gseaNb(object = all_results$egmt_all_Fru_Ctrl_table, geneSetID = c("HALLMARK_IL6_JAK_STAT3_SIGNALING","GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN",
                                                                   "GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"),subPlot = 2,newGsea = FALSE,
       curveCol = colors_transparent,
       addPval = TRUE,pDigit = 3,legend.position = "top")
??gseaNb
ggsave("stat3.pdf",width = 4,height = 4,units = "in",dpi = 1200)

gseaNb(object = all_results$egmt_all_Fru_Ctrl_table, geneSetID = c("GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"
                                                                   ,"GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN","KEGG_JAK_STAT_SIGNALING_PATHWAY",
                                                                   "HALLMARK_IL6_JAK_STAT3_SIGNALING"),subPlot = 2)
gseaNb(object = all_results$egmt_all_Fru_Ctrl_table, geneSetID = c("GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"
                                                                   ,"GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN","KEGG_JAK_STAT_SIGNALING_PATHWAY",
                                                                   "HALLMARK_IL6_JAK_STAT3_SIGNALING"),subPlot = 2)
KEGG_JAK_STAT_SIGNALING_PATHWAY


# GOBP_NEUTROPHIL_MIGRATION
# GOBP_GRANULOCYTE_MIGRATION
# GOBP_DENDRITIC_CELL_MIGRATION
# GOBP_POSITIVE_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION
# GOBP_MONONUCLEAR_CELL_MIGRATION
# GOBP_MYELOID_LEUKOCYTE_MIGRATION
# GOBP_LEUKOCYTE_MIGRATION


gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_GRANULOCYTE_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7) 
# 保存pdf,4.5 3 英寸
ggsave("GOBP_GRANULOCYTE_MIGRATION.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)

gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_NEUTROPHIL_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)
gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_DENDRITIC_CELL_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)
ggsave("GOBP_DENDRITIC_CELL_MIGRATION.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)

gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_POSITIVE_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)
ggsave("GOBP_POSITIVE_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)
gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_MONONUCLEAR_CELL_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)
ggsave("GOBP_MONONUCLEAR_CELL_MIGRATION.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)
gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_MYELOID_LEUKOCYTE_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)
ggsave("GOBP_MYELOID_LEUKOCYTE_MIGRATION.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)
gseaNb(all_results$egmt_c5bp_Fru_Ctrl_table, geneSetID = c("GOBP_LEUKOCYTE_MIGRATION"),subPlot = 2,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)
ggsave("GOBP_LEUKOCYTE_MIGRATION.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)
"GOBP_NEUTROPHIL_MIGRATION","GOBP_GRANULOCYTE_MIGRATION",
"GOBP_DENDRITIC_CELL_MIGRATION",
,"GOBP_MONONUCLEAR_CELL_MIGRATION","GOBP_MYELOID_LEUKOCYTE_MIGRATION","GOBP_LEUKOCYTE_MIGRATION"

# GSE9988_LPS_VS_LOW_LPS_MONOCYTE_UP
gseaNb(object = all_results$egmt_all_Fru_Ctrl_table,
       geneSetID = 'GSE9988_LPS_VS_LOW_LPS_MONOCYTE_UP',subPlot = 2,,newGsea = FALSE,addPval = TRUE,pDigit = 3,pHjust = 1,pvalY = 0.7)

ggsave("GSE9988_LPS_VS_LOW_LPS_MONOCYTE_UP.pdf",width = 4.5,height = 3,units = "in",dpi = 1200)

# 挑选颜色=====

# 如果还没有安装 RColorBrewer 包, 需要先安装
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer")
}

# 加载 RColorBrewer 包
library(RColorBrewer)

# 选择一组色盲友好的颜色组合
color_palette <- brewer.pal(3, "Dark2")

# 打印所选颜色的 HEX 值
print(color_palette)
barplot(rep(1, length(color_palette)), col = color_palette, border = NA)

"#1B9E77" "#D95F02" "#7570B3"
# 创建一个从蓝色到红色的颜色组合
color_palette <- colorRampPalette(c("blue", "red"))(3)

# 打印所选颜色的 HEX 值
print(color_palette)

gseaNb(egmt_all, geneSetID = c("HALLMARK_IL6_JAK_STAT3_SIGNALING","GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN",
                               "GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"),subPlot = 2,newGsea = FALSE,
       curveCol = color_palette,
       addPval = TRUE,pDigit = 3)


# 原始颜色
# 加载 grDevices 包，它提供了 adjustcolor() 函数
library(grDevices)

# 原始颜色
colors <- c("#1B9E77", "#D95F02", "#7570B3")
colors <- brewer.pal(8, "Dark2")
"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
colors <- colorRampPalette(c("blue", "grey","red"))(5)
"#0000FF" "#FF0000" "#BEBEBE"
"#0000FF" "#5F5FDE" "#BEBEBE" "#DE5F5F" "#FF0000"


colors <- c("blue", "red","grey")

# 设置 alpha 值
alpha <- 0.6 # 50% 的透明度
#1B9E77     #D95F02     #7570B3     #E7298A     #66A61E     #E6AB02     #A6761D     #666666 
"#1B9E7780" "#D95F0280" "#7570B380" "#E7298A80" "#66A61E80" "#E6AB0280" "#A6761D80" "#66666680"
# 使用 adjustcolor() 函数转换颜色
colors_transparent <- sapply(colors, adjustcolor, alpha.f = alpha)
colors_transparent
"#0000FF80" "#FF000080" "#BEBEBE80" 
#0000FF     #5F5FDE     #BEBEBE     #DE5F5F     #FF0000 
"#0000FF80" "#5F5FDE80" "#BEBEBE80" "#DE5F5F80" "#FF000080" 
#0000FF     #7F007F     #FF0000 
"#0000FF80" "#7F007F80" "#FF000080"


blue         red        grey 
"#0000FF99" "#FF000099" "#BEBEBE99" 

colors_transparent <- c("#0000FF99", "#FF000099","#32323299")

# 打印所选颜色的 HEX 值
print(colors_transparent)

gseaNb(egmt_all, geneSetID = c("HALLMARK_IL6_JAK_STAT3_SIGNALING","GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN",
                               "GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"),subPlot = 2,newGsea = FALSE,
       curveCol = colors_transparent,
       addPval = TRUE,pDigit = 3)


# REACTOME_MITOCHONDRIAL_IRON_SULFUR_CLUSTER_BIOGENESIS
# WP_MITOCHONDRIAL_COMPLEX_I_ASSEMBLY_MODEL_OXPHOS_SYSTEM
# GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY
gseaNb(egmt_all, geneSetID = c("REACTOME_MITOCHONDRIAL_IRON_SULFUR_CLUSTER_BIOGENESIS","WP_MITOCHONDRIAL_COMPLEX_I_ASSEMBLY_MODEL_OXPHOS_SYSTEM",
                               "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY"),subPlot = 2,newGsea = FALSE,
       curveCol = colors_transparent,
       addPval = TRUE,pDigit = 3)


# GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_PERMEABILITY_INVOLVED_IN_APOPTOTIC_PROCESS
# GOBP_ESTABLISHMENT_OF_MITOCHONDRION_LOCALIZATION
gseaNb(egmt_all, geneSetID = c("GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_PERMEABILITY_INVOLVED_IN_APOPTOTIC_PROCESS",
                               "GOBP_ESTABLISHMENT_OF_MITOCHONDRION_LOCALIZATION"),subPlot = 2,newGsea = FALSE,
       curveCol = colors_transparent,
       addPval = TRUE,pDigit = 3)

# 3.6.5
# 加载 grDevices 包，它提供了 col2rgb() 函数
library(grDevices)

# 原始颜色
color <- "#BEBEBE99"

# 将颜色转换为 RGB 格式
rgb_values <- col2rgb(color)

# 调整 RGB 值来使颜色更深
# 数值越低，颜色越深
rgb_values <- pmax(rgb_values - 80, 0)

# 将 RGB 值转回 HEX 格式
darker_color <- rgb(rgb_values[1,], rgb_values[2,], rgb_values[3,], maxColorValue = 255)

# 添加原始颜色的 alpha 值
darker_color <- paste0(darker_color, substr(color, 8, 9))
"#96969699"
"#82828299"
"#32323299"
# 打印新的颜色值
print(darker_color)
# # 小鼠基因集合
# gmtfile <- "~/bioinfo/msigdb_v2023.1.Mm_GMTs/m2.all.v2023.1.Mm.entrez.gmt"
# geneset <- read.gmt(gmtfile)
# egmt <- GSEA(ranks_RNAseq_ID, TERM2GENE = geneset, minGSSize = 10, pvalueCutoff = 0.99, verbose = F)
# # 
# ??GSEA

#查看结果
gsea.out.df <- egmt@result
View(gsea.out.df)
write_csv(gsea.out.df,"c5livergseahm.csv")
write.csv(egmt_all@result, file = "egmt_all_result.csv")

library(knitr)
BiocManager::install("kableExtra")
library(kableExtra)

egmt_all@result %>%
  kable("html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

gsea.out.df$ID
# 
# Accessing and modifying nested slots in S4 objects
kk2@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kk2@result$Description)
View(kk2@result)
# ====================================================================================================


# 6 4.17
# c5.go.bp.v2023.1.Hs.symbols.gmt
gseaplot2(egmt_c5bp, geneSetID = c("GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"
                                   ,"GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"), pvalue_table = F, base_size = 14, color = "black"
          , subplots = c(1,2))


gseaplot2(egmt_all, geneSetID = c("GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN"
                                  ,"GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN",
                                  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                                  "KEGG_JAK_STAT_SIGNALING_PATHWAY"
                                  # ,"GOBP_DEFENSE_RESPONSE_TO_BACTERIUM","GOBP_FRUCTOSE_METABOLIC_PROCESS"
                                  # ,"GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS"
), pvalue_table = F, base_size = 14, color = "black"
, subplots = c(1,2))
# 大小 高3，宽4.5 用AI调整  

# 6 4.17
# GSE9988_LPS_VS_LOW_LPS_MONOCYTE_UP
gseaplot2(egmt_all, geneSetID = "GSE9988_LPS_VS_LOW_LPS_MONOCYTE_UP", pvalue_table = F, base_size = 14, color = "black", 
          title = "LPS_VS_LOW_LPS_MONOCYTE_UP", subplots = c(1,2))

HP_ABNORMAL_MACROPHAGE_MORPHOLOGY
GOBP_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN
GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN

# GOBP_DEFENSE_RESPONSE_TO_BACTERIUM
gseaplot2(egmt, geneSetID = "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM", pvalue_table = F, base_size = 14, color = "black", 
          title = "DEFENSE RESPONSE TO BACTERIUM", subplots = c(1,2))
# GOBP_FRUCTOSE_METABOLIC_PROCESS
gseaplot2(egmt, geneSetID = "GOBP_FRUCTOSE_METABOLIC_PROCESS", pvalue_table = F, base_size = 14, color = "black", 
          title = "Fructose Metabolic Process", subplots = c(1,2))
# GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS
gseaplot2(egmt, geneSetID = "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS", pvalue_table = F, base_size = 14, color = "black", 
          title = "NEGATIVE REGULATION OF 
IMMUNE EFFECTOR PROCESS", subplots = c(1,2))


# HALLMARK_IL6_JAK_STAT3_SIGNALING
gseaplot2(egmt_all, geneSetID = "HALLMARK_IL6_JAK_STAT3_SIGNALING", pvalue_table = F, base_size = 14, color = "black", 
          title = "IL6 JAK STAT3 SIGNALING", subplots = c(1,2))

gseaplot2(kk2, geneSetID = c("mmu05100","mmu05132","mmu05034","mmu04530"),pvalue_table = T ,
          subplots = c(1,2,3),rel_heights = c(2, 0.5, 0.5)
)

gseaplot2(kk2, geneSetID = c("mmu05150","mmu04936","mmu04630"),pvalue_table = F ,
          subplots = c(1,2,3)
          # ,rel_heights = c(2, 0.5, 0.5)
)
# JAK,金葡菌，alcoholic 

gseaplot2(kk2, geneSetID = c("mmu05150"),pvalue_table = T ,
          subplots = c(1,2),title = 
            # ,rel_heights = c(2, 0.5, 0.5)
)

# HP_ABNORMAL_MACROPHAGE_MORPHOLOGY
gseaplot2(egmt_all, geneSetID = "HP_ABNORMAL_MACROPHAGE_MORPHOLOGY", pvalue_table = F, base_size = 14, color = "black", 
          title = "HP ABNORMAL MACROPHAGE MORPHOLOGY", subplots = c(1,2))
# GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_SUPERFAMILY_BINDING
gseaplot2(egmt_c5mf, geneSetID = "GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_SUPERFAMILY_BINDING", pvalue_table = F, base_size = 14, color = "black", 
          title = "GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_SUPERFAMILY_BINDING", subplots = c(1,2))
# GOMF_HYDROLASE_ACTIVITY_HYDROLYZING_N_GLYCOSYL_COMPOUNDS
gseaplot2(egmt_c5mf, geneSetID = "GOMF_HYDROLASE_ACTIVITY_HYDROLYZING_N_GLYCOSYL_COMPOUNDS", pvalue_table = F, base_size = 14, color = "black", 
          title = "GOMF_HYDROLASE_ACTIVITY_HYDROLYZING_N_GLYCOSYL_COMPOUNDS", subplots = c(1,2))
# GOMF_HYDROLASE_ACTIVITY_ACTING_ON_GLYCOSYL_BONDS
gseaplot2(egmt_c5mf, geneSetID = "GOMF_HYDROLASE_ACTIVITY_ACTING_ON_GLYCOSYL_BONDS", pvalue_table = F, base_size = 14, color = "black", 
          title = "GOMF_HYDROLASE_ACTIVITY_ACTING_ON_GLYCOSYL_BONDS", subplots = c(1,2))

# 加标题=====

# Define the geneSetID
geneSetID <- c("mmu04630")

# Get the description for the specified geneSetID
description <- kk2@result[kk2@result$ID == geneSetID, "Description"]

# Create a filename with the description, date, and size
filename <- paste0(description, "_", Sys.Date(), "_4.5x3.pdf")

# Plot the GSEA plot with the description as the title
p <- gseaplot2(kk2, geneSetID = geneSetID, pvalue_table = F, base_size = 14, color = "black", title = description, subplots = c(1,2))

# Save the plot with ggsave
ggsave(filename, plot = p, width = 4.5, height = 3)



gseaplot2(kk2, geneSetID = c("mmu05340","mmu01240"),pvalue_table = T ,
          subplots = c(1,2,3)
          # ,rel_heights = c(2, 0.5, 0.5)
)

c("mmu05100","mmu05340","mmu04530","mmu01240")

# GseaVis=====



# mmu04510
# Focal adhesion
# Define function to plot GSEA with description as title
# mmu05132 Salmonella infection
# mmu05034 Alcoholism
# mmu05340 Primary immunodeficiency
# mmu05100 Bacterial invasion of epithelial cells
# 
# mmu05340
# mmu04530
# mmu05034 Alcoholism
# mmu01240
# Biosynthesis of cofactors


??gseaplot2
# ====
library(EGSEA)
names(v)
v

## ----collections--------------------------------------------------------------
library(EGSEAdata)
## ----collectionlist-----------------------------------------------------------
info = egsea.data("mouse", returnInfo = TRUE)
names(info)
info$msigdb$info$collections

## ----loadegsea, message=FALSE, warning=FALSE----------------------------------

library(htmlwidgets)
library(DT)
library(plotly)
## ----indexing-----------------------------------------------------------------
# gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="mouse", 
#                      msigdb.gsets=c("c2", "c3", "c4", "c5", "c6","c7"), go.part = TRUE)
# gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="mouse", 
# msigdb.gsets=c("c2", "c3", "c4", "c5", "c6","c7"), go.part = TRUE,gsdb.gsets = "all",kegg.updated = TRUE)
gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="mouse", 
                     msigdb.gsets=c("c2", "c5"), go.part = TRUE,gsdb.gsets = "none",kegg.updated = TRUE)


# ??buildIdx
names(gs.annots)

# ## ----exploresets--------------------------------------------------------------
# class(gs.annots$c2)
# summary(gs.annots$c2)
# show(gs.annots$c2)
# s = getSetByName(gs.annots$c2, "SMID_BREAST_CANCER_LUMINAL_A_DN")
# class(s)
# names(s)
# names(s$SMID_BREAST_CANCER_LUMINAL_A_DN)

# ## ----indexclass---------------------------------------------------------------
# slotNames(gs.annots$c2)

## ----symbolmap----------------------------------------------------------------
colnames(v$genes)
# symbolsMap = v$genes[, c(1, 2)]
symbolsMap = v$genes[, c(1, 3)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

## ----base---------------------------------------------------------------------
egsea.base()

## ----selectbasemethods--------------------------------------------------------
baseMethods = egsea.base()[-2]
baseMethods

## ----combine------------------------------------------------------------------
egsea.combine()

## ----sort---------------------------------------------------------------------
egsea.sort()

## ----egseatest----------------------------------------------------------------
gsa = egsea(voom.results=v, contrasts=contr.matrix,  
            gs.annots=gs.annots, symbolsMap=symbolsMap,
            baseGSEAs=baseMethods, sort.by="med.rank",
            num.threads = 32, report = FALSE,sum.plot.axis = "med.rank",sum.plot.cutoff = 110,keep.set.scores = TRUE,keep.base = TRUE,keep.limma = TRUE
            ,display.top = 25)

## ----showegsea----------------------------------------------------------------
show(gsa)

## ----summariseegsea-----------------------------------------------------------
summary(gsa)

# ORA

top.Table = topTable(efit, coef = 1, number = Inf, p.value = 0.05,
                     lfc = 1)
deGenes = as.character(top.Table$ENTREZID)
logFC = top.Table$logFC
names(logFC) = deGenes
# build the gene set collection index
gs.annots_ORA = buildIdx(entrezIDs = deGenes, species = "mouse",
                         msigdb.gsets = "none", kegg.exclude = c("Metabolism"))
## [1] "Building KEGG pathways annotation object ... "
# perform the ORA analysis set report = TRUE to generate
# the EGSEA interactive report
gsa_ORA = egsea.ora(geneIDs = deGenes, universe = as.character(v$genes[,1]), logFC = logFC, title = "X24IL13-X24", gs.annots = gs.annots_ORA,
                    symbolsMap = top.Table[, c(1, 2)], display.top = 5, report.dir = "./il13-egsea-ora-report",
                    num.threads = 4, report = FALSE)

show(gsa_ORA)
summary(gsa_ORA)

## ----topsets------------------------------------------------------------------
names(gs.annots)
topSets(gsa, gs.label="c2", contrast = "treatedvsuntreated", names.only=TRUE)

View(topSets(gsa, gs.label="kegg", contrast = "treatedvsuntreated", names.only=FALSE,number = Inf) )

## ----topsetslim---------------------------------------------------------------
t = topSets(gsa, contrast = "treatedvsuntreated",
            names.only=FALSE, number = Inf, verbose = FALSE)

## ----topsets2-----------------------------------------------------------------
topSets(gsa, gs.label="kegg", contrast="treatedvsuntreated", sort.by="med.rank"
        ,names.only=FALSE, number = Inf, verbose = FALSE)

View(gsa@results[["kegg"]][["test.results"]][["treatedvsuntreated"]])

View(gsa@results[["kegg"]][["test.results"]][["treatedvsuntreated"]])
write.csv(gsa@results[["kegg"]][["test.results"]][["treatedvsuntreated"]],"kegg.csv")

topSets(gsa, gs.label="kegg", contrast="treatedvsuntreated", sort.by="med.rank",names.only = TRUE,number = 10,verbose = FALSE)

topSets(gsa, gs.label="c5CC", contrast="treatedvsuntreated", sort.by="med.rank",names.only = TRUE,number = 10)
topSets(gsa, gs.label="c5BP", contrast="treatedvsuntreated", sort.by="med.rank",names.only = TRUE,number = 10)
topSets(gsa, gs.label="c5MF", contrast="treatedvsuntreated", sort.by="med.rank",names.only = TRUE,number = 10)

## ----heatmaps-----------------------------------------------------------------
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_UP", gs.label="c5MF",
            contrast = "treatedvsuntreated", file.name = "GO_GLUTATHIONE_TRANSFERASE_ACTIVITY", format="png")
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_DN", gs.label="c2",
            contrast = "treatedvsuntreated", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_DN", format="png")


## ----pathwayplot1, eval=FALSE-------------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "BasalvsLP", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction")

## ----pathwayplot2, eval=FALSE-------------------------------------------------
plotPathway(gsa, gene.set = "ErbB signaling pathway",
            contrast = "treatedvsuntreated", gs.label = "kegg",
            file.name = "12ErbB signaling pathway")


# 批量画图===============================

top10_pathways <- topSets(gsa, gs.label="kegg", contrast="treatedvsuntreated", sort.by="med.rank",names.only = TRUE,number = 10)

for (i in 1:10) {
  pathway <- top10_pathways[i]
  file_name <- paste0(i,"_",pathway)
  plotPathway(gsa, gene.set = pathway, contrast = "treatedvsuntreated", gs.label = "kegg", file.name = file_name)
}
# ==============================
topSets(gsa, gs.label="kegg", contrast="treatedvsuntreated", sort.by="med.rank")
gsa@results[["kegg"]][["test.results"]][["treatedvsuntreated"]]

## ----mdsplot------------------------------------------------------------------
plotMethods(gsa, gs.label = "c2", contrast = "treatedvsuntreated", 
            file.name = "mds_c2_treatedvsuntreated", format="pdf")
plotMethods(gsa, gs.label = "c5BP", contrast = "BasalvsLP", 
            file.name = "mds_c5_treatedvsuntreated", format="pdf")

## ----keggsummaryplot1---------------------------------------------------------
plotSummary(gsa, gs.label = "kegg", contrast = "treatedvsuntreated", 
            file.name = "summary_kegg_treatedvsuntreated", format="pdf",x.cutoff = 50,use.names = TRUE,x.axis = "med.rank")
??plotSummary
## ----c2summaryplot2-----------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3, 
            file.name = "summary_c2_LPvsML", 
            x.axis = "med.rank", format="png")

## ----c2summaryplot3-----------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3, 
            file.name = "summary_sig_c2_LPvsML", 
            x.axis = "med.rank", x.cutoff=300, format="png")

## ----summaryplotkegg1and2-----------------------------------------------------
plotSummary(gsa, gs.label = "kegg", contrast = c(1,2), 
            file.name = "summary_kegg_1vs2", format="png")

## ----gographs-----------------------------------------------------------------
plotGOGraph(gsa, gs.label="c5BP", contrast = 1, file.name="BasalvsLP-c5BP-top-", format="png")
plotGOGraph(gsa, gs.label="c5CC", contrast = 1, file.name="BasalvsLP-c5CC-top-", format="png")

## ----summarybarplot-----------------------------------------------------------
plotBars(gsa, gs.label = "c2", contrast="treatedvsuntreated", file.name="treatedvsuntreated-c2-bars", format="png")
plotBars(gsa, gs.label = "kegg", contrast="treatedvsuntreated", file.name="treatedvsuntreated-kegg", format="pdf")

## ----summaryheatmap-----------------------------------------------------------
plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir",
                   file.name="summary_heatmaps_c2", format="png")
plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir",
                   file.name="summary_heatmaps_kegg", format="pdf")




#========探索S4对象

dds = as.DESeqDataSet(x)
rse = as(dds, "RangedSummarizedExperiment")
rse
# BiocManager::install("iSEE")

suppressPackageStartupMessages({
  library(iSEE)
})
app <- iSEE(rse)
shiny::runApp(app)



# ======



