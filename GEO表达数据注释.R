library(devtools)
install_github("jmzeng1314/idmap1") ##bioconductor包
library(idmap1)
install_github("jmzeng1314/idmap2") #全部表达芯片的soft文件
library(idmap2)
library(GEOmirror)
library(Biobase)
library(base)
setwd(...)
#设置下载时限
options(timeout = 600)
geoChina('GSE13507')
load('GSE13507_eSet.Rdata')
gset
a=exprs(gset[[1]])
a[1:4,1:4]
gset[[1]]@annotation
#将平台号输入下方（若idmap1无法获取注释数据，改用idmap2获取）
b=idmap1::getIDs("GPL6102")
#b=idmap2::get_soft_IDs("GPL6102")
head(b)
#使用match函数找到对应的symbol
  symbol_annotations <- b$symbol[match(rownames(a), b$ID)]
#将符号注释添加到a的行名中
rownames(a) <- symbol_annotations
#删去无法注释基因
a <- a[!is.na(rownames(a)), ]
write.table(a, file = "GSE13507_bulk.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
