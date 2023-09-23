# https://mp.weixin.qq.com/s?__biz=MzkxOTM5MzQwNQ==&mid=2247489440&idx=1&sn=a371882fe1a885f4946810c93755a3b2&chksm=c1a392b2f6d41ba4ac5fd39fee3c4a258b29297f6d34972634b8ee00e4ef00a2b6dbfa7d4d69&mpshare=1&scene=1&srcid=04224cSkGsEVh4lUG7l3rtN2&sharer_sharetime=1682151773725&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("xmc811/Scillus", ref = "development")
rm(list = ls())
library(tidyverse)
library(Scillus) # scRNA-seq数据的处理和可视化
library(Seurat)
library(magrittr)
library(purrr)

# 4.1 文件路径
a <- list.files("/media/wane/KINGSTON1/GSE128531_RAW", full.names = TRUE)
m <- tibble(
  file = a,
  sample = stringr::str_remove(basename(a), ".csv.gz"),
  group = rep(c("CTCL", "Normal"), each = 3)
)

DT::datatable(m)
# 4.2 配色
pal <- tibble(
  var = c("sample", "group", "seurat_clusters"),
  pal = c("Set2", "Set1", "Paired")
)

# 4.3 读入数据
scRNA <- load_scfile(m)
map(scRNA, print)

# 4.4 数据长度
# 这里的长度等于metadata的行数，即m的行数。🤨
length(scRNA)

# 5QC可视化
# 5.1 线粒体基因可视化
plot_qc(scRNA, metrics = "percent.mt")
plot_qc(scRNA,
  metrics = "percent.mt",
  plot_type = "box" # "combined", "box" or "violin"
)

plot_qc(scRNA,
  metrics = "percent.mt",
  group_by = "group" # "sample", 其他在metadata中的列
)

# 5.2 nFeature可视化
plot_qc(scRNA, metrics = "nFeature_RNA")

plot_qc(scRNA,
  metrics = "nFeature_RNA",
  group_by = "group",
  pal_setup = "Accent"
)

# 5.3 nCount可视化
plot_qc(scRNA, metrics = "nCount_RNA")
plot_qc(scRNA,
  metrics = "nCount_RNA",
  plot_type = "density"
) +
  scale_x_log10()

# 6过滤与整合
# 6.1 过滤
# subset参数的语法与Seurat对象的subset()函数是一样的。🤒
# 用的时候，会自动绘制barplot以显示过滤前后的细胞数。😉
scRNA_f <- filter_scdata(scRNA, subset = nFeature_RNA > 500 & percent.mt < 10)
# 6.2 标准化处理
# 接着就是做一下Normalize，FindVariableFeatures，CellCycleScoring等标准化处理了。🧐
scRNA_f %<>%
  purrr::map(.f = NormalizeData) %>%
  purrr::map(.f = FindVariableFeatures) %>%
  purrr::map(
    .f = CellCycleScoring,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes
  )

# 6.3 整合数据
scRNA_int <- IntegrateData(
  anchorset = FindIntegrationAnchors(
    object.list = scRNA_f,
    dims = 1:30, k.filter = 50
  ),
  dims = 1:30
)
scRNA_int %<>%
  ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

scRNA_int %<>%
  RunPCA(npcs = 50, verbose = T)

scRNA_int %<>%
  RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 30) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

scRNA_int

# 7因子化处理（可选步骤）
# 主要是处理一下metadata的数据，这样作图会更好看一些，如果你没有metadata，可以不做这一步。😘
m %<>%
  mutate(group = factor(group, levels = c("Normal", "CTCL")))

scRNA_int %<>%
  refactor_seurat(metadata = m)

scRNA_int

# 4降维及其可视化
# 4.1 初步绘图
plot_scdata(scRNA_int, pal_setup = pal)
plot_scdata(scRNA_int, pal_setup = "Dark2")

# 4.3 分组可视化
plot_scdata(scRNA_int, color_by = "group", pal_setup = pal)
# 4.4 分面可视化
plot_scdata(scRNA_int, split_by = "sample", pal_setup = pal)
# 4.5 手动配色
plot_scdata(scRNA_int, color_by = "sample", 
            pal_setup = c("red","orange","yellow","green","blue","purple"))

# 5统计及其可视化
# 5.1 按sample统计
plot_stat(scRNA_int, 
          plot_type = "group_count"
          ## 三种，"group_count", "prop_fill", and "prop_multi"
)
# 5.2 按cluster统计
plot_stat(scRNA_int, "group_count", group_by = "seurat_clusters", pal_setup = pal)
# 5.3 堆叠柱形图
plot_stat(scRNA_int, 
          plot_type = "prop_fill", 
          pal_setup = c("grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20"))
# 5.4 按cluster和sample统计
plot_stat(scRNA_int, plot_type = "prop_multi", pal_setup = "Set3")
# 5.5 按cluster和group统计
plot_stat(scRNA_int, plot_type = "prop_fill", group_by = "group")
# 5.6 换个配色
plot_stat(scRNA_int, plot_type = "prop_multi", 
          group_by = "group", pal_setup = c("sienna","bisque3"))

# 6热图及其可视化
# 6.1 寻找marker
# 我们首先要用Seurat包的FindAllMarkers确定一下marker，分析方法也比较多，包括wilcox, roc, t, poisson, DESeq2等。😷
markers <- FindAllMarkers(scRNA_int, logfc.threshold = 0.25, min.pct = 0.1, only.pos = F)
# 6.2 热图可视化
# 在热图中，每一行代表一个基因，每一列代表一种细胞，默认绘制每种细胞的前8个基因。🤓
# 细胞可以通过sort_var进行排序，默认设置为c("seurat_clusters")，即细胞按cluster进行排序。🙃
# anno_var用来指定注释数据或者metadata中的相关数据。🥰
plot_heatmap(dataset = scRNA_int, 
             markers = markers,
             sort_var = c("seurat_clusters","sample"),
             anno_var = c("seurat_clusters","sample","percent.mt","S.Score","G2M.Score"),
             anno_colors = list("Set2",                                             
                                # RColorBrewer palette
                                c("red","orange","yellow","purple","blue","green"), 
                                # color vector
                                "Reds",
                                c("blue","white","red"),                            
                                # Three-color gradient
                                "Greens"))
plot_heatmap(dataset = scRNA_int,
             n = 6,
             markers = markers,
             sort_var = c("seurat_clusters","sample"),
             anno_var = c("seurat_clusters","sample","percent.mt"),
             anno_colors = list("Set2",
                                c("red","orange","yellow","purple","blue","green"),
                                "Reds"),
             hm_limit = c(-1,0,1),
             hm_colors = c("purple","black","yellow"))

# 7富集分析
# 7.1 制定cluster的GO分析
plot_cluster_go(markers, cluster_name = "1", org = "human", ont = "CC")
# 7.2 所有cluster的GO分析
plot_all_cluster_go(markers, org = "human", ont = "BP")

# 8GSEA分析
# 8.1 差异分析
# 首先我们要用find_diff_genes做一下差异分析，寻找差异基因。😜
de <- find_diff_genes(dataset = scRNA_int, 
                      clusters = as.character(0:7),
                      comparison = c("group", "CTCL", "Normal"),
                      logfc.threshold = 0,   # threshold of 0 is used for GSEA
                      min.cells.group = 1)   # To include clusters with only 1 cell

gsea_res <- test_GSEA(de, 
                      pathway = pathways.hallmark)

# 8.2 GSEA结果可视化
plot_GSEA(gsea_res, p_cutoff = 0.1, colors = c("#0570b0", "grey", "#d7301f"))

# 9.1 小试牛刀
plot_measure(dataset = scRNA_int, 
             measures = c("KRT14","percent.mt"), 
             group_by = "seurat_clusters", 
             pal_setup = pal)

plot_measure_dim(dataset = scRNA_int, 
                 measures = c("nFeature_RNA","nCount_RNA","percent.mt","KRT14"),
                 split_by = "sample")
















