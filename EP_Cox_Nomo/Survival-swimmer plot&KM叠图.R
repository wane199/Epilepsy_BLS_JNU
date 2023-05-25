# Creating a survival-swimmer plot in R(http://rstudio-pubs-static.s3.amazonaws.com/389060_a0ea812d8f8d490393bd1c65a6dcddef.html)
# Importing the required libraries
rm(list = ls())
library(magrittr)
library(stringi)
library(readr)   # Reading in the dataset
library(ggplot2) # Viewing the dataset
library(forcats) # Sorting factors
library(RColorBrewer) # Plot colours
library(dplyr, warn.conflicts=FALSE)   # Manipulating the dataframes
library(purrr, warn.conflicts=FALSE)   # Manipulating dataframe metadata
library(zoo, warn.conflicts=FALSE)     # Filling in  NA values
library(reshape2) # Reformmating dataframes 

# Importing the dataset
swimmer_file = "https://blogs.sas.com/content/graphicallyspeaking/files/2014/06/Swimmer_93.txt"
col.names = c("subjectID", "stage", "startTime", "endTime", 
              "isContinued", "responseType", "responseStartTime", "responseEndTime", "Durable")
df <- readr::read_lines(swimmer_file) %>%
  # Split by line recursion (\r\n)
  stringi::stri_split(fixed="\r\n", simplify=TRUE) %>%
  # Take only lines starting with a number (sample id)
  .[grepl("^[0-9]+", .)] %>%
  # Remove spaces from response column
  gsub(pattern="\\sresponse", replacement="_response") %>%
  # Remove spaces from stage column
  gsub(pattern="Stage\\s",  replacement="Stage_") %>%
  # Some lines missing 'Stage' and 'isContinued' column. 
  # Replace any set of 8 or more spaces with ' . '
  gsub(pattern="\\s{8,}", replacement=' . ') %>%
  # Split strings by spaces, do not include empty strings as columns
  stringi::stri_split(fixed=" ", simplify=TRUE, omit_empty=TRUE) %>%
  # Convert to dataframe
  as.data.frame(stringsAsFactors=FALSE) %>%
  # Set the column names
  purrr::set_names(col.names) %>%
  # We need to do some more cleaning up of the dataframe
  # Convert all . to NAs
  dplyr::na_if(".") %>%
  # Fill NAs in Stage column
  dplyr::mutate(stage=zoo::na.locf(stage)) %>%
  # Turn isContinued into boolean
  dplyr::mutate(isContinued=dplyr::if_else(isContinued=="FilledArrow", TRUE, FALSE, missing=FALSE)) %>%
  # Convert stage variable to factor, remove underscore
  dplyr::mutate(stage = as.factor(gsub(pattern="_", replacement=" ", x=stage))) %>%
  # Remove underscore from response types 
  dplyr::mutate(responseType = gsub("_", " ", responseType)) %>%
  # Change Durable from character to numeric
  dplyr::mutate(Durable = as.numeric(Durable)) %>%
  # Change Time variables from character to numeric
  dplyr::mutate_at(vars(dplyr::ends_with("Time")), as.numeric)

# Viewing the data.
df %>% dplyr::glimpse()



# 可视化临床试验中常用的游泳图（swimmer plot）(https://mp.weixin.qq.com/s?__biz=MzkzNzMxNjgxMA==&mid=2247487555&idx=1&sn=bd76a790ef77e6e6ad2c19ada8d320ac&chksm=c29008e6f5e781f08f0e4ef892b4e483a477586dcc2fd231240a3fcc8716e70df7de606c6653&scene=21#wechat_redirect)
# 游泳者图/时间轨迹图（swimmer plots）(https://zhuanlan.zhihu.com/p/578958546)
# install.packages("swimplot")
library(swimplot)
library(ggplot2)

## 不同id患者的治疗方式及其治疗结束时间
knitr::kable(head(ClinicalTrial.Arm,10))

## 不同id患者事件发生时间
knitr::kable(head(ClinicalTrial.AE,10))

## 不同id患者治疗反应开始和结束时间以及反应类型
knitr::kable(head(ClinicalTrial.Response,10))

# 基础绘图
swimmer_plot(df=ClinicalTrial.Arm,
             id='id',
             end='End_trt',
             fill='lightblue',
             width=.85)

# 添加分组信息
arm_plot <- swimmer_plot(df=ClinicalTrial.Arm,
                         id='id',
                         end='End_trt',
                         name_fill='Arm',
                         id_order='Arm',
                         col="black",alpha=0.75,width=.8)

arm_plot

# 分面
swim_plot_stratify <-swimmer_plot(df=ClinicalTrial.Arm,
                                  id='id',end='End_trt',
                                  name_fill='Arm',
                                  id_order ='increasing',
                                  col="black",alpha=0.75,width=.8,
                                  base_size = 14,
                                  stratify= c('Age','Sex'))

swim_plot_stratify

# 添加时间点信息
AE_plot <- arm_plot + 
  swimmer_points(df_points=ClinicalTrial.AE,
                 id='id',time='time',
                 name_shape = 'event',
                 size=2.5,fill='white',col='black')
AE_plot

# 反应时间分组
arm_plot + swimmer_points(df_points=ClinicalTrial.AE,
                          id='id',time='time',
                          name_shape = 'event',
                          size=2.5,fill='white',name_col = 'Related')

# 添加时间线信息
Response_plot <- arm_plot +
  swimmer_lines(df_lines=ClinicalTrial.Response,id='id',
                start ='Response_start',end='Response_end',
                name_col='Response',size=1)

Response_plot

# 同时添加时间点和时间线信息
Response_plot_with_points <- Response_plot+
  swimmer_points_from_lines(df_lines=ClinicalTrial.Response,id='id',
                            start = 'Response_start',end = 'Response_end', 
                            cont = 'Continued_response',
                            name_col='Response',size=2)

Response_plot_with_points

# 添加箭头
AE_plot+
  swimmer_arrows(df_arrows=ClinicalTrial.Arm,
                 id='id',arrow_start='End_trt',
                 cont = 'Continued_treatment',
                 name_col='Arm',type ="open",cex=1)

# 美化
AE_plot <-  AE_plot +
  scale_fill_manual(name="Treatment",values=c("#e41a1c", "#377eb8",'#4daf4a'))+
  scale_color_manual(name="Treatment",values=c("#e41a1c", "#377eb8",'#4daf4a')) +
  scale_shape_manual(name="Adverse event",values=c(21,24,17),breaks=c('AE','SAE','Death'))

AE_plot

AE_plot <-  AE_plot +
  scale_fill_manual(name="Treatment",values=c("#e41a1c", "#377eb8",'#4daf4a'))+
  scale_color_manual(name="Treatment",values=c("#e41a1c", "#377eb8",'#4daf4a')) +
  scale_shape_manual(name="Adverse event",values=c(21,24,17),breaks=c('AE','SAE','Death'))

AE_plot

Response_plot_with_points <- Response_plot_with_points+guides(fill = guide_legend(override.aes = list(shape = NA)))
Response_plot_with_points

Response_plot_with_points <- Response_plot_with_points+
  annotate("text", x=3.5, y=20.45, label="Continued response",size=3.25)+
  annotate("text",x=2.5, y=20.25, label=sprintf('\u2192'),size=8.25)+
  coord_flip(clip = 'off', ylim = c(0, 17))
Response_plot_with_points

Response_plot_with_points <- Response_plot_with_points + theme(axis.title.y=element_blank(),
                                                               axis.text.y=element_blank(),
                                                               axis.ticks.y=element_blank()) +labs(y="Time since enrollment (months)") 
Response_plot_with_points

#Overriding legends to have colours for the events and no points in the lines
p1 <- arm_plot + swimmer_points(df_points=ClinicalTrial.AE,id='id',time='time',name_shape =
                                  'event',size=2.5,col='black',name_fill = 'event') +
  scale_shape_manual(values=c(21,22,23),breaks=c('AE','SAE','Death'))


p1 +scale_fill_manual(name="Treatment",values=c('grey40',"#e41a1c", "#377eb8",1,'#4daf4a','grey90'))
#####################################
# 叠中叠，一行代码实现带有风险表的生存曲线插入子图(https://mp.weixin.qq.com/s?__biz=MzU4OTc0OTg2MA==&mid=2247500015&idx=1&sn=21902096c48ac75b40fa6af7195e289a&chksm=fdca4be4cabdc2f20e8e67760f503a4b3660ed8acfab96ae965b65ca4569295e9102c803e74f&mpshare=1&scene=1&srcid=0513p68Pd85YCXbx50UjJlB7&sharer_sharetime=1684423067245&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd)
# 1. 加载数据
rm(list = ls())
dt <- read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\cph\\cph2\\TLE220group.csv",sep = ";")
str(dt)
library(showtext)
showtext_auto()
# 2. 拟合曲线
# 使用survfit2()函数来拟合曲线。
library(ggsurvfit)
fit <- survfit2(Surv(Follow_up_timemon, Rel._in_5yrs == 1) ~ 1, data = dt)
fit # 可输出部分信息

# 3. 简单生存曲线
fit %>%
  ggsurvfit()

# 4. 生存曲线排版
# 可以在主图上面添加带有风险表的子图。
p1 <- survfit(Surv(Follow_up_timemon, Rel._in_5yrs == 1) ~ 1, data = dt) %>% 
  ggsurvfit(linewidth = 0.8, linetype_aes = TRUE, color = "cyan") + 
  add_risktable(risktable_stats = "n.risk",
                stats_label = list(n.risk = "复发风险病例数")) +
  add_censor_mark(color = "cyan") +
  add_confidence_interval() +
  scale_x_continuous(breaks = seq(0,60,12)) + 
  labs(y = "无复发概率",
    x = "术后随访时间（月）") +
  # add_quantile() +
  theme_classic(base_size = 14) +
  theme(legend.position=c(0.2,0.2)) +
  ggsci::scale_fill_nejm()
p1
p2 <- survfit(Surv(Follow_up_timemon, Rel._in_5yrs == 1) ~ 1, data = dt) %>% 
  ggsurvfit(linewidth = 0.8, type = "cumhaz", linetype_aes = TRUE, color = "red") +
  add_risktable(risktable_stats = "n.risk",
                stats_label = list(n.risk = "复发风险病例数"),
                risktable_height = 0.20) +
  add_censor_mark(color = "red") +
  add_confidence_interval() +
  scale_x_continuous(breaks = seq(0,60,12)) + 
  labs(y = "复发概率",
       x = "术后随访时间（月）") +
  theme_classic(base_size = 10) +
  theme(legend.position='none') +
  ggsci::scale_fill_nejm()

library(patchwork)
wrap_plots(ggsurvfit_build(p1)) +
  inset_element(ggsurvfit_build(p2), 
                left = 0.55, bottom = 0.60, right = 1.0, top = 1.0)

# https://www.danieldsjoberg.com/ggsurvfit/
p <- survfit2(Surv(time, status) ~ surg, data = df_colon) |>
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable() +
  add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.75)
p +
  # limit plot to show 8 years and less
  coord_cartesian(xlim = c(0, 8)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    title = "Recurrence by Time From Surgery to Randomization",
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.2),
                     expand = c(0.015, 0)) +
  scale_x_continuous(breaks = 0:10, 
                     expand = c(0.02, 0))

str(dt)

# 基础绘图
swimmer_plot(df=dt,
             id='ID',
             end='Follow_up_timemon',
             fill='lightblue',
             width=.85)

# 添加分组信息
arm_plot <- swimmer_plot(df=dt,
                         id='ID',
                         end='Follow_up_timemon',name_fill='Rel._in_5yrs',
                         # name_fill='Rel._in_5yrs',
                         # id_order='Rel._in_5yrs',
                         col="black",alpha=0.75,width=.8)

arm_plot
AE_plot <- arm_plot + 
  swimmer_points(df_points=dt,
                 id='ID',time='Follow_up_timemon',
                 name_shape = 'Rel._in_5yrs',
                 size=2.5,fill='white',col='black')
AE_plot










