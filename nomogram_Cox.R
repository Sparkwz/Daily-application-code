# 加载包
library(rms)
library(survival)
library(pROC)
library(timeROC)
library(ggplot2)

# 设置工作目录
setwd("//")

################ 1. 读取训练集数据 ################
#读取数据
training_dataset <- read.table("nomogram.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(training_dataset)

################ 2. 单因素Cox回归 ################
#单因素分析(单个分别读取)
vars_to_test <- c("Tumor_distance", "POST_NLR", "AMP") #纳入变量名称
univ_formulas <- sapply(vars_to_test, function(x) as.formula(paste0("Surv(OS, Censor) ~ ", x)))
univ_models <- lapply(univ_formulas, function(x) coxph(x, data = training_dataset))
univ_results <- lapply(univ_models, function(x) {
  x <- summary(x)
  HR <- round(exp(x$coef[1]), 2)
  lower <- round(exp(x$conf.int[,"lower .95"]), 2)
  upper <- round(exp(x$conf.int[,"upper .95"]), 2)
  p.value <- signif(x$wald["pvalue"], 3)
  data.frame(HR = HR, lower = lower, upper = upper, p.value = p.value)
})
univ_results_df <- do.call(rbind, univ_results)
rownames(univ_results_df) <- vars_to_test
print(univ_results_df)

################ 3. 多因素Cox回归 ################
# 选择纳入多因素的变量
multiv_model <- coxph(Surv(OS, Censor) ~ Tumor_distance + POST_NLR + AMP, data = training_dataset)
summary(multiv_model)

################ 4. 构建列线图 ################
# 设置数据分布
dd <- datadist(training_dataset)
options(datadist = "dd")

# 使用rms包构建模型
fit <- cph(Surv(OS, Censor) ~ Tumor_distance + POST_NLR + AMP, data = training_dataset, surv = TRUE)

# 绘制列线图（1年、3年、5年生存概率）
surv <- Survival(fit)
nom <- nomogram(fit, fun = list(
  function(x) surv(12, x),
  function(x) surv(36, x),
  function(x) surv(60, x)
), funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"))

pdf("nomogram_training.pdf", width = 8, height = 6)
plot(nom)
dev.off()

################ 5. 绘制ROC曲线 ################
#生成分析需要的输入格式，三列数据分别为生存时间、生存状态、Cox预测值
time_roc <- timeROC(  
  T=training_dataset$OS,  
  delta=training_dataset$Censor,  
  marker=fit$linear.predictors,  
  cause = 1,  
  weighting = "marginal",  
  times = c(12,36,60),  
  ROC=TRUE,  
  iid=TRUE
)
time_ROC_df <- data.frame(  
  TP_1year = time_roc$TP[, 1],  
  FP_1year = time_roc$FP[, 1],  
  TP_3year = time_roc$TP[, 2],  
  FP_3year = time_roc$FP[, 2],  
  TP_5year = time_roc$TP[, 3],  
  FP_5year = time_roc$FP[, 3]
  )
###提取AUC的百分之95CI
#1年
ci95low <- sprintf("%.3f", time_roc$AUC[[1]]-(1.96*time_roc[["inference"]][["vect_sd_1"]][["t=12"]]))
ci95up <- sprintf("%.3f", time_roc$AUC[[1]]+(1.96*time_roc[["inference"]][["vect_sd_1"]][["t=12"]]))
ci951year <- paste("(",ci95low,"-",ci95up,")")
#3年
ci95low <- sprintf("%.3f", time_roc$AUC[[2]]-(1.96*time_roc[["inference"]][["vect_sd_1"]][["t=36"]]))
ci95up <- sprintf("%.3f", time_roc$AUC[[2]]+(1.96*time_roc[["inference"]][["vect_sd_1"]][["t=36"]]))
ci953year <- paste("(",ci95low,"-",ci95up,")")
#5年
ci95low <- sprintf("%.3f", time_roc$AUC[[3]]-(1.96*time_roc[["inference"]][["vect_sd_1"]][["t=60"]]))
ci95up <- sprintf("%.3f", time_roc$AUC[[3]]+(1.96*time_roc[["inference"]][["vect_sd_1"]][["t=60"]]))
ci955year <- paste("(",ci95low,"-",ci95up,")")

p1 <- ggplot(data = time_ROC_df) +  
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +  
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +  
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +  
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +  
  theme_classic() +  
  annotate("text",           
           x = 0.65, y = 0.25, size = 4.5,           
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc$AUC[[1]]),ci951year), color = "#BC3C29FF"  
           ) +  
  annotate("text",           
           x = 0.65, y = 0.15, size = 4.5,           
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc$AUC[[2]]),ci953year), color = "#0072B5FF"  ) +  
  annotate("text",           
           x = 0.65, y = 0.05, size = 4.5,           
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc$AUC[[3]]),ci955year), color = "#E18727FF"  ) +  
  labs(x = "False positive rate", y = "True positive rate") +  
  theme(    
    axis.text = element_text(face = "plain", size = 15, color = "black"),    
    axis.title.x = element_text(face = "bold", size = 15, color = "black", margin = margin(c(15, 0, 0, 0))),    
    axis.title.y = element_text(face = "bold", size = 15, color = "black", margin = margin(c(0, 15, 0, 0)))  
    )
p1
ggsave("./ROC.pdf",plot = p1,height = 8,width = 8)
