library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tinyarray)

#选择项目
project <- "TCGA-BLCA"
setwd("~/Spark")
if (!dir.exists(paste0("./", project))) {
  dir.create(paste0("./", project))
}
#TCGAbiolinks下载 project 数据
query <- GDCquery(project = project,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts"
)
#多线程下载
GDCdownload(query, files.per.chunk = 100)
data <- GDCprepare(query = query)

#提取基因表达量
mrna <- data[data@rowRanges$gene_type=="protein_coding",]

#提取tpm_unstrand数据（不应该使用基因表达counts数据或log(counts)数据，否则会造成样本之间的表达量不可比）
Matrix <- assay(mrna,"tpm_unstrand")

# 将矩阵转换为数据框，并保留原始列名
tpm <- data.frame(Matrix, check.names = FALSE)

#提取基因名(新版TCGAbiolinks下载的数据集自带gene name，无需手动转换)
symbol <- rowData(mrna)$gene_name
tpm$gene_name=mrna@rowRanges$gene_name
tpm=aggregate(tpm,.~gene_name,"sum") #对相同基因的表达量求和

#查看转换后的TPM数据
head(tpm[,1:5])
rownames(tpm) <-tpm$gene_name
tpm <- tpm[,-1]
tpm[1:4,1:4]

#基因过滤
#仅保留在一半以上样本里表达的基因
nrow(tpm)
tpm <- tpm[apply(tpm, 1, function(x) sum(x > 0) > 0.5*ncol(tpm)), ]
nrow(tpm)

#使用logCPM或logTPM数据(不同样本间进行基因表达量比较使用统一的归一化手段才有可比性)
#tpm <- log2(tpm+1)

#去除正常样本
Group <- make_tcga_group(tpm)
table(Group)
tpm <- tpm[, Group=='tumor']
tpm <- tpm[apply(tpm,1, function(x){sum(x>0)>0.5*ncol(tpm)}),]#去除样本后再次过滤
tpm[1:4,1:4]


#临床信息
clinical <- colData(data)
meta <- as.data.frame(clinical)
nrow(meta)
length(unique(meta$sample))
meta <- distinct(meta,sample,.keep_all = T)

#样本处理
meta <- meta %>%
  mutate(overall_survival = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)) %>%
  mutate(event = vital_status)

#去掉生存信息不全或者生存时间小于30天的样本
k1 <- meta$overall_survival >= 30;table(k1) #筛选小于30days样本
k2 <- !(is.na(meta$overall_survival)|is.na(meta$overall_survival));table(k2) #筛选生存信息不全样本
meta <- meta[k1 & k2,]

#规范结局事件为0/1
table(meta$event)
meta$event <- ifelse(meta$event == "Dead", 1, 0)
table(meta$event)

#规范生存时间按月份
range(meta$overall_survival)
meta$overall_survival <- meta$overall_survival/30.4375
range(meta$overall_survival)

#匹配表达数据及临床信息
head(rownames(meta))
head(colnames(tpm))
s <- intersect(rownames(meta),colnames(tpm));length(s) #获取两者交集
tpm <- tpm[,s]
meta <- meta[s,]
dim(tpm)
dim(meta)
identical(rownames(meta),colnames(tpm))

#保存表达数据
write.table(tpm, paste0("./", project,"/symbol_unlog.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#保存临床信息
table(colnames(meta))
write.csv(meta_df, paste0("./", project,"/clinical.csv"), row.names = T) #临床注释信息
