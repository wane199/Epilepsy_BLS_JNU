# https://github.com/cxli233/ggpathway
# https://mp.weixin.qq.com/s?__biz=Mzg3MzQzNTYzMw==&mid=2247500255&idx=1&sn=a7e14dafa645ccf5a9ed11bf5e2c1d12&chksm=cee29941f99510577ef8f1223108c8df337dd08760fc6ffc02e9ad322cddcfc2baf066849ce0&mpshare=1&scene=1&srcid=12052bMoDDkCPnt1WoZd7Yj3&sharer_sharetime=1670383300989&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd

# 加载R包
rm(list = ls())
library(tidyverse)
library(igraph)
library(ggraph)
library(readxl)
library(viridis)
library(RColorBrewer)
library(rcartocolor)

# 构建边数据
example1_edge_table <- tribble(
  ~from, ~to,  ~label,
  "Glc6P", "6P-gluconolactone",  "Glc6PHD",
  "6P-gluconolactone", "6P-glucoconate",  "6P-gluconolactonase",
  "6P-glucoconate", "Ru5P", "6P-gluconateDH"
)
# 构建点数据
example1_nodes_table <- tribble(
  ~name, ~x,  ~y,
  "Glc6P", 1, 0,
  "6P-gluconolactone", 2, 0,  
  "6P-glucoconate", 3, 0,
  "Ru5P", 4, 0
)

# 文件整合
example1_network <- graph_from_data_frame(
  d = example1_edge_table,
  vertices = example1_nodes_table,
  directed = T
)

# 绘制基础通路图
ggraph(example1_network, layout = "manual", 
       x = x, y = y) +
  geom_node_text(aes(label = name), hjust = 0.5) +
  geom_edge_link(aes(label = example1_edge_table$label), 
                 angle_calc = 'along',
                 label_dodge = unit(2, 'lines'),
                 arrow = arrow(length = unit(0.5, 'lines')), 
                 start_cap = circle(4, 'lines'),
                 end_cap = circle(4, 'lines')) +
  theme_void()  
ggsave("./EP_Cox_Nomo/Pentose_1.svg", height = 2, width = 6.5, bg = "white")
ggsave("./EP_Cox_Nomo/Pentose_1.png", height = 2, width = 6.5, bg = "white")

###################
# 绘制PPP代谢途径
# 导入数据
example2_edges <- read_excel("./EP_Cox_Nomo/OPPP_edges.xlsx")
example2_nodes <- read_excel("./EP_Cox_Nomo/OPPP_nodes.xlsx")

# 构建新标签
example2_nodes <- example2_nodes %>% 
  mutate(label = str_remove(name, "_\\d"))

# 文件整合
example2_network <- graph_from_data_frame(
  d = example2_edges,
  vertices = example2_nodes,
  directed = T)
# 绘制PPP代谢途径图
ggraph(example2_network, layout = "kk") +
  geom_node_point(size = 3, aes(fill = as.factor(carbons)), 
                  alpha = 0.8, shape = 21, color = "grey20") +
  geom_node_text(aes(label = label), hjust = 0.5, repel = T) +
  geom_edge_link(label_dodge = unit(2, 'lines'),
                 arrow = arrow(length = unit(0.4, 'lines')), 
                 start_cap = circle(1, 'lines'),
                 end_cap = circle(2, 'lines')) +
  scale_fill_manual(values = carto_pal(7, "Vivid")) +
  labs(fill = "Carbons") +
  theme_void()  

###################
# TCA途径
# 导入数据
example3_edges <- read_excel("./EP_Cox_Nomo/ggpathway/TCA_cycle_edges.xlsx")
example3_nodes <- read_excel("./EP_Cox_Nomo/ggpathway/TCA_cycle_nodes.xlsx")

# 构建标签文本
example3_nodes <- example3_nodes %>% 
  mutate(label = str_remove(name, "_\\d"))

# 整合数据
example3_network <- graph_from_data_frame(
  d = example3_edges,
  vertices = example3_nodes,
  directed = T)

# 绘制TCA途径图
ggraph(example3_network, layout = "manual",
       x = x, y = y) +
  geom_node_point(size = 3, aes(fill = as.factor(carbons)), 
                  alpha = 0.8, shape = 21, color = "black") +
  geom_edge_link(arrow = arrow(length = unit(0.4, 'lines')), 
                 start_cap = circle(0.5, 'lines'),
                 end_cap = circle(0.5, 'lines'), 
                 width = 1.1, alpha = 0.5) +
  geom_node_text(aes(label = label), hjust = 0.5, repel = T) +
  annotate(geom = "text", label = "TCA Cycle", 
           x = 0, y = 0, size = 5, fontface = "bold") +
  scale_fill_manual(values = carto_pal(7, "Vivid")) +
  labs(fill = "Carbons") +
  theme_void() +
  coord_fixed()
ggsave("./EP_Cox_Nomo/ggpathway/TCA_2.0.png", height = 4, width = 5, bg = "white")

# Subsetting nodes and edges
# We can simplify this by removing the cofactors.
example3_nodes_trim <- example3_nodes %>% 
  filter(carbons != "cofactor")

example3_edges_trim <- example3_edges %>% 
  filter(from %in% example3_nodes_trim$name &
           to %in% example3_nodes_trim$name)

example3_network_trim <- graph_from_data_frame(
  d = example3_edges_trim,
  vertices = example3_nodes_trim,
  directed = T
)

ggraph(example3_network_trim, layout = "manual",
       x = x, y = y) +
  geom_node_point(size = 3, aes(fill = as.factor(carbons)), 
                  alpha = 0.8, shape = 21, color = "grey20") +
  geom_edge_link(arrow = arrow(length = unit(0.4, 'lines')), 
                 start_cap = circle(0.5, 'lines'),
                 end_cap = circle(1, 'lines'), 
                 width = 1.1, alpha = 0.5) +
  geom_node_text(aes(label = label), hjust = 0.5, repel = T) +
  annotate(geom = "text", label = "TCA Cycle", 
           x = 0, y = 0, size = 5, fontface = "bold") +
  scale_fill_manual(values = carto_pal(7, "Vivid")) +
  labs(fill = "Carbons") +
  theme_void() +
  coord_fixed()

ggsave("./EP_Cox_Nomo/ggpathway/TCA_2.svg", height = 4, width = 5, bg = "white")
ggsave("./EP_Cox_Nomo/ggpathway/TCA_2.png", height = 4, width = 5, bg = "white")






