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
# Let’s have a look at the data using the glimpse tool. We can see that the data frame is ‘tidy’. There is one row for every observation. A tidy dataframe may mean there are multiple entries for a given subject.
df %>% dplyr::glimpse()



# 可视化临床试验中常用的游泳图（swimmer plot）(https://mp.weixin.qq.com/s?__biz=MzkzNzMxNjgxMA==&mid=2247487555&idx=1&sn=bd76a790ef77e6e6ad2c19ada8d320ac&chksm=c29008e6f5e781f08f0e4ef892b4e483a477586dcc2fd231240a3fcc8716e70df7de606c6653&scene=21#wechat_redirect)





install.packages("swimplot")
install.packages("ggplot2")

library(swimplot)
library(ggplot2)











#####################################
# 叠中叠，一行代码实现带有风险表的生存曲线插入子图(https://mp.weixin.qq.com/s?__biz=MzU4OTc0OTg2MA==&mid=2247500015&idx=1&sn=21902096c48ac75b40fa6af7195e289a&chksm=fdca4be4cabdc2f20e8e67760f503a4b3660ed8acfab96ae965b65ca4569295e9102c803e74f&mpshare=1&scene=1&srcid=0513p68Pd85YCXbx50UjJlB7&sharer_sharetime=1684423067245&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd)
# 1. 加载数据
rm(list = ls())
dt <- read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\cph\\cph2\\TLE220group.csv",sep = ";")
str(dt)

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
  add_risktable(risktable_stats = "n.risk") +
  add_censor_mark() +
  add_confidence_interval() +
  # add_quantile() +
  theme_classic(base_size = 14) +
  theme(legend.position=c(0.2,0.2)) +
  ggsci::scale_fill_nejm()
p1
p2 <- survfit(Surv(Follow_up_timemon, Rel._in_5yrs == 1) ~ 1, data = dt) %>% 
  ggsurvfit(linewidth = 0.8, type = "cumhaz", linetype_aes = TRUE, color = "cyan") +
  add_risktable(risktable_stats = "n.risk",
                risktable_height = 0.20) +
  add_censor_mark() +
  add_confidence_interval() +
  theme_classic(base_size = 10) +
  theme(legend.position='none') +
  ggsci::scale_fill_nejm()

library(patchwork)
wrap_plots(ggsurvfit_build(p1)) +
  inset_element(ggsurvfit_build(p2), 
                left = 0.50, bottom = 0.55, right = 1.0, top = 1.0)
















