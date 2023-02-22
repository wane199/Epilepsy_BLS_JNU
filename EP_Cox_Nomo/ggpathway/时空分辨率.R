# 时空分辨率示意图
# 跨尺度(https://blog.51cto.com/u_15535576/5116197)
df <- data.frame(
  x = c("10^-4", "10^-3", "10^-2", "10^-1"),
  y = c(0.0001, 0.001, 0.01, 0.1)
)
df

lwd_pt <- .pt*72.27/96
theme_set(theme_classic(base_size = 22, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))

ggplot(df, aes(x = x, y = y)) +
  geom_col(aes(fill = x)) +
  scale_y_continuous(
    labels = c(
      expression(italic(0)),
      expression(10^-4),
      expression(10^-3),
      expression(10^-2),
      expression(10^-1) # 1 %*% 
    ),
    position = "right",
    expand = c(0, 0),
    # breaks = c(0, 0.001, 0.002, 0.003, 0.004),
    # limits = c(0, 0.005)
  ) +
  coord_flip()+ xlab("时间分辨率(秒)")+ylab("空间分辨率(米)")+  theme(legend.position = "none") + 
  geom_text(x = 1, y = 7, label = "fontsize = 10",  size = 10/.pt) 
  # scale_y_reverse(expand=c(0,0),
  #                 position="right")+

ggplot(data = rate,aes(x=reorder(地区,地区生产总值.)))+
  geom_bar(aes(y=地区生产总值.,fill=地区),stat = "identity")+
  labs(title = "2015年山东各地市生产总值及占比情况")+
  geom_text(aes(label=rate$地区生产总值.,y=地区生产总值.-500))+
  geom_text(aes(label=b,y=地区生产总值.+600))+
  coord_flip()+

  scale_y_continuous(limits = c(0,10000),breaks = seq(0,10000,2500))+
  theme_minimal()+theme(legend.position = "none")


a <- 10
i <- seq(-4, -1, 1)
x <- a * 10^i
x <- seq(10^(-4), 10^(-1), by = 10000)

# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

# The dataset is provided in the gapminder library
library(gapminder)
data <- gapminder %>% filter(year=="2007") %>% dplyr::select(-year)

# Most basic bubble plot
data %>%
  arrange(desc(pop)) %>%
  mutate(country = factor(country, country)) %>%
  ggplot(aes(x=gdpPercap, y=lifeExp, size=pop, fill=continent)) +
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(range = c(.1, 24), name="Population (M)") +
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  # theme_ipsum() +
  scale_x_reverse(labels = c(expression(10^-4),
                                expression(10^-3),
                                expression(10^-2),
                                expression(10^-1),
                                expression(10^0),
                                expression(10^1),
                                expression(10^2),
                                expression(10^3),
                                expression(10^4)),
                             breaks = seq(0, 50000, 6250)) + 
  scale_y_continuous(labels = c(expression(10^-4),
                                expression(10^-3),
                                expression(10^-2),
                                expression(10^-1)),
                     breaks = c(32,48,64,80), limits = c(30,90),expand = c(0, 0),
                     position = "right") + coord_flip() +
  theme(axis.ticks.length.x = unit(-0.20, 'cm'),axis.ticks.length.y = unit(-0.20, 'cm'),
            axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))) +  #修改y轴刻度朝内
  theme(axis.line = element_line(arrow = arrow(length = unit(0.5, 'cm'))))  + #坐标轴尾端为箭头
  geom_text(x = 1, y = 7, label = "fontsize = 10",  size = 10/.pt) +
  theme(legend.position="bottom") +
  ylab("空间分辨率(米)") +
  xlab("时间分辨率(秒)") + 
  theme(legend.position = "none")




