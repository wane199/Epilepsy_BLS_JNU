# [生存分析之landmark](https://mp.weixin.qq.com/s?__biz=Mzg2OTY4NzAxNA==&mid=2247483674&idx=1&sn=b71b104e19dc5ee688b2d4c9733b9695&chksm=ce980a93f9ef83859827e1c4fa89c8f2c9a8ec08819060d5a04a356d69535a43bdde9395f74e&mpshare=1&scene=1&srcid=0113tqJAW7uPgIqaH9kLSHrh&sharer_sharetime=1673612760201&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd)
# [生存分析5：限制平均生存时间Restricted Mean Survival Time,RMST](https://mp.weixin.qq.com/s?__biz=MzkyNTIyNDYwNQ==&mid=2247485393&idx=1&sn=695f8d9a7d2bb12b394cc4c1b8167f82&chksm=c1c89c23f6bf1535764ca7fa18e888d1a57d7b9c98b47eb7359e10aafa5fa82eae17b483fe19&mpshare=1&scene=1&srcid=01132ByPutJJcv2Hv5rmwgvV&sharer_sharetime=1673612698209&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd)
rm(list = ls())
library(survival)

# 数据集
# data(colon)
colon <- colon
str(colon)

# 将分层变量rx由3个水平变为2个水平
oldvals <- c("Obs", "Lev", "Lev+5FU")
newvals <- factor(c("观察组", "治疗组", "治疗组"))
colon$newrx <- newvals[match(colon$rx, oldvals)]

# 生存率估计，按newrx进行分层
fit <- survfit(Surv(time, status) ~ newrx, data = colon)
summary(fit)

library(jskm) # 进行landmark分析
# 简单设置
jskm(fit, marks = T, cut.landmark = 720)
# 图片格式修改调整
jskm(fit,
  marks = F, pval = T, table = T, label.nrisk = "No. at risk",
  size.label.nrisk = 8, xlabs = "days", ylabs = "Survival",
  ystrataname = "rx", ystratalabs = c("观察组", "治疗组"),
  legendposition = c(0.85, 0.9), timeby = 360, ylines = c(0.25, 1),
  cut.landmark = 720
)

################################################
# 程序包的安装及数据的加载
# install.packages("survRM2")
# install.packages("tidyverse")
# install.packages("survminer")
library(survRM2)
library(tidyverse)
library(survminer)

# data<- read_csv("Checkmate057.csv")
# fkdt <- tibble(subject = as.factor(1:10),
#                time = sample(4:20, 10, replace = T),
#                status = sample(c("Censor", rep("Event", 2)), 10,
#                                replace = T))

dt <- read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\cph\\TLE234group.csv", sep = ",", header = TRUE)
# 设置因子的水平标签
dt$Rel._in_5yrs <- factor(dt$Rel._in_5yrs,
                          levels = c(0, 1),
                          labels = c("未复发", "复发"))
fkdt <- transform(dt, subject = as.factor(1:nrow(dt)))
colnames(fkdt)
# 游泳图（Swimmer Plots）
ggplot(fkdt, aes(subject, Follow_up_timemon)) +
  geom_bar(stat = "identity", width = 0.25) +
  geom_point(
    data = fkdt,
    aes(subject, Follow_up_timemon, color = factor(Rel._in_5yrs), shape = factor(Rel._in_5yrs)),
    size = 2
  ) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "grey20", size = 4, angle = 30, hjust = 1, vjust = 0, face = "plain"), 
    legend.title = element_blank(),
    legend.position = "bottom"
  )
# 可视化数据框热图
library(funkyheatmap)
library(kableExtra)
data <- dt[,6:24]
data <- transform(data, ID = as.factor(1:nrow(data)))
data <- data %>% 
  rownames_to_column("id") %>%
  # arrange(desc()) %>%
  head(20)
funky_heatmap(dt[,6:24])




