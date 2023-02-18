# 机器学习模型绘制校正曲线（Calibration Curve），多分类onevsrest 
# Generate multiple calibration plots with colored/shaded 95 percent confidence intervals
rm(list = ls())
library(runway)
library(multiROC)

dt <- read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\bls\\mripredictions\\outputsave\\test-cal-f3_f4.csv")
head(dt)

#### 三分类哑变量处理
# Merge true labels and predicted values
library(fastDummies)
true_label <- fastDummies::dummy_cols(dt$outcomes)
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label[-1], dt[-1])
write.csv(final_df,'C:\\Users\\wane1\\Documents\\file\\sci\\bls\\mripredictions\\outputsave\\test-cal-f4.csv')

# 4.1 multi_roc and multi_pr function
roc_res <- multi_roc(dt, force_diag=T)
pr_res <- multi_pr(dt, force_diag=T)

# 4.2 Confidence Intervals
# 4.2.1 List of AUC results
unlist(roc_res$AUC)

# 4.2.2 Bootstrap
roc_ci_res <- roc_ci(dt, conf= 0.95, type='basic', R = 100, index = 4)
pr_ci_res <- pr_ci(dt, conf= 0.95, type='basic', R = 100, index = 4)

# 4.2.3 Output All Results
roc_auc_with_ci_res <- roc_auc_with_ci(dt, conf= 0.95, type='bca', R = 100)
roc_auc_with_ci_res
pr_auc_with_ci_res <- pr_auc_with_ci(dt, conf= 0.95, type='bca', R = 100)
pr_auc_with_ci_res

# 5 Plots
plot_roc_df <- plot_roc_data(roc_res)
plot_pr_df <- plot_pr_data(pr_res)

# 5.2 Plot
# 5.2.1 ROC Plot
library(ggplot2)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) + geom_path(aes(color = Group, linetype=Method)) + geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') + theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.justification=c(1, 0), legend.position=c(.95, .05), legend.title=element_blank(), legend.background = element_rect(fill=NULL, size=0.5, linetype="solid", colour ="black"))

# 5.2.2 PR Plot
ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + geom_path(aes(color = Group, linetype=Method), size=1.5) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.justification=c(1, 0), legend.position=c(.95, .05), legend.title=element_blank(), legend.background = element_rect(fill=NULL, size=0.5, linetype="solid", colour ="black"))


data(multi_model_dataset)
df <- read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\bls\\mripredictions\\outputsave\\test-cal-f4_1.csv")
head(df)

library(reshape2) #  首先加载一下reshape2包
aql <- melt(df,  id.vars = c("outcomes"), variable.name = "classes", value.name = "predictions") 
head(aql)  # 查看数据前6列
tail(aql)   # 查看数据后6列
df <- df0[,complete.cases(t(df0))]  # 提取不含空值的列
aql <- aql[complete.cases(aql),]  # 提取不含空值的行
write.csv(aql,'C:\\Users\\wane1\\Documents\\file\\sci\\bls\\mripredictions\\outputsave\\test-cal-f4_2.csv')

cal_plot_multi(aql, outcome = 'outcomes',model = 'classes', prediction = 'predictions', n_bins = 10000)
threshperf_plot_multi(aql, outcome = 'outcomes', prediction = 'predictions', model = 'classes')

library(riskRegression)



