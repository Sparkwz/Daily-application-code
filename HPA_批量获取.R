library(BiocStyle)
library(HPAanalyze)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
rm(list = ls())
options(stringsAsFactors = F)

dir.create("..\\HPA-output")
gene="NEK8"
#获得HPA网站中该基因的xml文件
hpa_target_gene<-hpaXmlGet(gene)
#根据标本选定组织名称
tissue="Urinary bladder"
#将xml中组织染色的信息提取出来
hpa_target_gene_fig_url<-hpaXmlTissueExpr(hpa_target_gene)
hpa_target_gene_fig_url_1<-as.data.frame(hpa_target_gene_fig_url[[1]])
hpa_target_gene_fig_url_1[1:6,1:18]
hpa_target_gene_fig_url_2<-as.data.frame(hpa_target_gene_fig_url[[2]])
hpa_target_gene_fig_url_2[1:6,1:18]
#选择自己感兴趣的组织
hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_1[hpa_target_gene_fig_url_1$tissueDescription2==tissue,]
hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_2[hpa_target_gene_fig_url_2$tissueDescription2==tissue,]
#为该组织该基因单独建个文件夹储存
picDir <- paste('..\\HPA-output\\',gene, tissue,"IHC-2\\", sep = "_")
if (!dir.exists(picDir)) {  
     dir.create(picDir)
  }

for (i in 1:nrow(hpa_target_gene_fig_url_tissue)) {  
  file_url<-hpa_target_gene_fig_url_tissue$imageUrl[i]  
  file_dir<-paste(picDir,gene,tissue,hpa_target_gene_fig_url_tissue$patientId[i],hpa_target_gene_fig_url_tissue$tissueDescription1[i],hpa_target_gene_fig_url_tissue$tissueDescription2[i],".tiff",sep = "_")  
  download.file(url = file_url,destfile = file_dir,mode = "wb")}
write.csv(hpa_target_gene_fig_url_tissue,paste(picDir,gene,"IHC-2_result_tab.csv",sep = "_"))
