library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(msigdbr)

#选择项目
class <-"TCGA-BLCA"

#数据下载
query <- GDCquery(
  project = class,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query = query)
data <- GDCprepare(query = query)
if (!dir.exists(paste0("./", class))) {
  dir.create(paste0("./", class))
}
Exp <- assay(data) %>% as.data.frame() # 提取数据表达
ann <- rowRanges(data) #提取基因注释
ann <- as.data.frame(ann)
rownames(ann) <- ann$gene_id
ann <- ann[rownames(Exp),]
write.table(ann, file = paste0("./", class, "/ann.txt"), row.names = FALSE)
Exp <- cbind(data.frame(Gene = ann$gene_name), Exp)
write.table(Exp, file = paste0("./", class, "/exp.txt"), row.names = FALSE) #表达矩阵
clinical <- GDCquery_clinic(project= class, type = "clinical") #提取临床信息
write.csv(clinical, paste0("./", class,"/clinical.csv"), row.names = F) #临床注释信息