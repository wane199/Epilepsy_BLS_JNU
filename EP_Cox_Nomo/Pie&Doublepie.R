# (多层饼图)[https://www.jianshu.com/p/7410c023df7b]
rm(list = ls())
library(ggplot2) # 绘图
library(ggsci) # 配色
library(showtext)
showtext_auto(enable = TRUE)
font_add(family = "YaHei", regular = "msyh.ttc")

# 构建测试数据
dat <- data.frame(
  x = rep("b", 21),
  y = rep("a", 21),
  z = rep("c", 21),
  cat1 = paste("c", 1:21, sep = "_"),
  # cat1 = c('进展','进展','进展','机遇','机遇','机遇','挑战'),
  cat2 = c(
    "进展", "进展", "进展", "进展", "进展", "进展", "进展",
    "机遇", "机遇", "机遇", "机遇", "机遇", "机遇", "机遇", "机遇",
    "挑战", "挑战", "挑战", "挑战", "挑战", "挑战"
  ),
  value1 = rep(1, 21),
  value2 = rep(1, 21)
)

# 分别求所占百分比
dat1 <- aggregate(dat$value1, by = list(dat$cat1), FUN = sum)
dat1$per1 <- dat1$x / sum(dat1$x)

# for循环构建标签的相对位置
for (i in seq(nrow(dat1), 1)) {
  if (i == nrow(dat1)) {
    dat1$per.y1[i] <- dat1$per1[i] / 2
  } else {
    dat1$per.y1[i] <- sum(dat1$per1[(i + 1):nrow(dat1)]) + dat1$per1[i] / 2
  }
}

# 构建标签后合并数据
dat1$label1 <- paste(dat1$Group.1, "(", round(dat1$per1 * 100, 2), "%", ")", sep = "")
dat <- merge(dat, dat1[, c(1, 3, 4, 5)], by.x = "cat1", by.y = "Group.1")

# 重复操作
dat2 <- aggregate(dat$value2, by = list(dat$cat2), FUN = sum)
dat2$per2 <- dat2$x / sum(dat2$x)

for (i in seq(nrow(dat2), 1)) {
  if (i == nrow(dat2)) {
    dat2$per.y2[i] <- dat2$per2[i] / 2
  } else {
    dat2$per.y2[i] <- sum(dat2$per2[(i + 1):nrow(dat2)]) + dat2$per2[i] / 2
  }
}

dat2$label2 <- paste(dat2$Group.1, "(", round(dat2$per2 * 100, 2), "%", ")", sep = "")
dat <- merge(dat, dat2[, c(1, 3, 4, 5)], by.x = "cat2", by.y = "Group.1")

# 绘图
ggplot(dat) +
  # 绘制柱状图
  geom_bar(
    aes(y,
      ifelse(cat2 == "a3", per2, per2 / 2),
      fill = cat2
    ),
    stat = "identity", width = 1.3
  ) +
  # 添加标签
  geom_text(
    aes(1.25, as.numeric(per.y2),
      label = label2
    ),
    size = 2.5, color = "black"
  ) +
  # 绘制柱状图
  geom_bar(aes(x, per1, fill = cat1),
    stat = "identity", width = .8, color = "white"
  ) +
  # 添加标签
  geom_text(aes(2, as.numeric(per.y1), label = label1),
    size = 2.5, color = "black"
  ) +
  # 设置Y轴刻度
  scale_y_continuous(labels = scales::percent) +
  coord_polar(theta = "y") + # 转换坐标轴
  theme_void() +
  scale_fill_igv() + # 设置填充色
  theme(legend.position = "none") # 隐藏图例


##################
browsers <- structure(list(browser = structure(c(
  3L, 3L, 3L, 3L, 2L, 2L,
  2L, 1L, 5L, 5L, 4L
), .Label = c(
  "Chrome", "Firefox", "MSIE",
  "Opera", "Safari"
), class = "factor"), version = structure(c(
  5L,
  6L, 7L, 8L, 2L, 3L, 4L, 1L, 10L, 11L, 9L
), .Label = c(
  "Chrome 10.0",
  "Firefox 3.5", "Firefox 3.6", "Firefox 4.0", "MSIE 6.0", "MSIE 7.0",
  "MSIE 8.0", "MSIE 9.0", "Opera 11.x", "Safari 4.0", "Safari 5.0"
), class = "factor"), share = c(
  10.85, 7.35, 33.06, 2.81, 1.58,
  13.12, 5.43, 9.91, 1.42, 4.55, 1.65
), ymax = c(
  10.85, 18.2, 51.26,
  54.07, 55.65, 68.77, 74.2, 84.11, 85.53, 90.08, 91.73
), ymin = c(
  0, 10.85, 18.2, 51.26, 54.07, 55.65, 68.77, 74.2, 84.11, 85.53,
  90.08
)), .Names = c("browser", "version", "share", "ymax", "ymin"), row.names = c(NA, -11L), class = "data.frame")

theme_set(theme_classic())
ggplot(browsers) +
  geom_rect(aes(fill = version, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  coord_polar(theta = "y") +
  xlim(c(0, 4))

ggplot(browsers) +
  geom_bar(aes(x = factor(1), fill = browser), width = 1) +
  coord_polar(theta = "y")

ggplot(browsers) +
  geom_rect(aes(fill = version, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(aes(fill = browser, ymax = ymax, ymin = ymin, xmax = 3, xmin = 0)) +
  xlim(c(0, 4)) +
  theme(aspect.ratio = 1)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)
ggplot(browsers) +
  geom_rect(aes(fill = version, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(aes(fill = browser, ymax = ymax, ymin = ymin, xmax = 3, xmin = 0)) +
  xlim(c(0, 4)) + scale_fill_manual(values = mycolors) +
  theme(aspect.ratio = 1) +
  coord_polar(theta = "y")


#######################
# (R语言绘制旭日图（嵌套多层的饼图或圆环图）)[https://mp.weixin.qq.com/s?__biz=MzIxNzc1Mzk3NQ==&mid=2247484687&idx=1&sn=48820a0c5f03106849681168450cfde6&chksm=97f5b517a0823c01af54228ae6cfcbadf32920b63f7b1267275d1b83d5d5ff9ded436d8d898d&mpshare=1&scene=1&srcid=0403BSAGMnLCRzESzK7nzbEB&sharer_sharetime=1680488674175&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd]
#ggsunburst 包的旭日图，参考自
#https://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot/37015211

#rPython 仅在 Linux 或 Mac 下可用 install.packages() 安装
#Windows 下可通过该方法安装：https://github.com/cjgb/rPython-win
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('rPython')) install.packages('rPython')
install.packages('http://genome.crg.es/~didac/ggsunburst/ggsunburst_0.0.9.tar.gz', repos = NULL, type = 'source')
library(ggsunburst)

#模拟数据
browsers<-structure(list(
  browser = structure(c(3L, 3L, 3L, 3L, 2L, 2L, 2L, 1L, 5L, 5L, 4L), .Label = c('Chrome', 'Firefox', 'MSIE', 'Opera', 'Safari'),),
  version = structure(c(5L, 6L, 7L, 8L, 2L, 3L, 4L, 1L, 10L, 11L, 9L), .Label = c('Chrome 10.0', 'Firefox 3.5', 'Firefox 3.6', 'Firefox 4.0', 'MSIE 6.0', 'MSIE 7.0', 'MSIE 8.0', 'MSIE 9.0', 'Opera 11.x', 'Safari 4.0', 'Safari 5.0'),),
  share = c(10.85, 7.35, 33.06, 2.81, 1.58, 13.12, 5.43, 9.91, 1.42, 4.55, 1.65),
  ymax = c(10.85, 18.2, 51.26, 54.07, 55.65, 68.77, 74.2, 84.11, 85.53, 90.08, 91.73),
  ymin = c(0, 10.85, 18.2, 51.26, 54.07, 55.65, 68.77, 74.2, 84.11, 85.53, 90.08)),
  .Names = c('browser', 'version', 'share', 'ymax', 'ymin'), row.names = c(NA, -11L), class = 'data.frame')

browsers$browser <- browsers$parent
write.table(browsers, file = 'browsers.csv', row.names = FALSE, sep = ',')

sb <- sunburst_data('browsers.csv', type = 'node_parent', sep = ',', node_attributes = c('browser', 'size'))

#为内部节点添加名称以着色
sb$rects[!sb$rects$leaf,]$browser <- sb$rects[!sb$rects$leaf,]$name

#作图
p <- sunburst(sb, rects.fill.aes = 'browser', node_labels = TRUE, node_labels.min = 15)
p + geom_text(data = sb$leaf_labels,
              aes(x = x, y = 0.1, label=paste(size, '%'), angle = angle, hjust = hjust), size = 2)

#继续展示 plotly 包交互式图形
library(plotly)
d <- data.frame(
  ids = c(
    '进展', '机遇', '挑战', '进展 - Football', 'Soccer',
    '进展 - Rugby', '机遇 - Football', 'Rugby',
    '机遇 - American Football','挑战 - 公平','执行方面的挑战','责任'),
  labels = c(
    '进展', '机遇', '挑战', 'Football', 'Soccer', 'Rugby',
    'Football', 'Rugby', 'American<br>Football', '公平','执行方面的挑战','责任'),
  parents = c(
    '', '', '', '进展', '进展', '进展', '机遇',
    '机遇', '机遇','挑战','挑战','挑战'),
  stringsAsFactors = FALSE
)
write.table(d, file = 'C:\\Users\\wane1\\Documents\\file\\sci\\aiep\\AI展望.csv', row.names = FALSE, sep = ',')

p <- plot_ly(d, ids = ~ids, labels = ~labels, parents = ~parents, type = 'sunburst')
p


