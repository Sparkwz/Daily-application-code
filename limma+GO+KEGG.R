
#引用包
library(limma)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
logFC_cutoff <- 1         #logFC过滤条件(logFC=1: Fold change=2     logFC=0.585: Fold change=1.5)
adj.P.Val.CutOff <- 0.05         #fdr临界值

#读取文件,并对输入文件整理
rt=read.table("symbol.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=log2(data+1)
table(colnames(data))

#提取对照组和实验组数据
control_cols <- grep("^C-", colnames(data), value = TRUE)
experimental_cols <- grep("^F-", colnames(data), value = TRUE)
data_control <- data[,which(colnames(data) %in% control_cols)] ##对照组数据
data_case <- data[,which(colnames(data) %in% experimental_cols)]##实验组数据

#设置control和case分组矩阵
group_list <- c(rep("control",dim(data_control)[2]),rep("case",dim(data_case)[2]))
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(data)

#构建比较矩阵
contrast.matrix <- makeContrasts(case-control,levels = design)

#构建线性拟合模型
fit <- lmFit(data,design) # 非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)# 用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = "case - control",n = Inf)

# 定义差异基因
dif <- DEG %>%
  mutate(Significant = case_when(
  logFC > logFC_cutoff & P.Value < adj.P.Val.CutOff ~ "Up",
  logFC < -logFC_cutoff & P.Value < adj.P.Val.CutOff ~ "Down",
  TRUE~"Not-Change"))

# 查看统计的差异基因数量
table(dif$Significant)

# 将结果保存为Txt文件
write.table(dif, file="dif.txt", sep="\t", row.names=T, quote=F)

#GO+KEGG分析
#引用包
pvalueFilter=0.05           #p值过滤条件
qvalueFilter=0.05           #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

rt=read.table("dif.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=as.vector(rownames(rt))
head(genes)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="GO-barplot.pdf",width = 10,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GO-bubble.pdf",width = 10,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#KEGG analysis
#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="kegg-barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

#气泡图
pdf(file="kegg-bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

