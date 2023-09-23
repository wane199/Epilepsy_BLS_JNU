# https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247492295&idx=1&sn=d4e81588c0ac2906849c8bc44d079209&chksm=9b3c9b97ac4b1281ab462f866538f096ba2473a3d1426a27d8d97ccc636bfd081e1ba8934dfd&scene=21#wechat_redirect
library(multtest)
if (!require(multtest)) install.packages("multtest")
if (!require(Seurat)) install.packages("Seurat")
if (!require(dplyr)) install.packages("dplyr")
if (!require(mindr)) install.packages("mindr")
if (!require(mindr)) install.packages("tidyverse")
##### 自动读取cellranger(LINUX)输出的feature barcode matric
rm(list = ls())
pbmc.data <- Read10X(data.dir = "C:\\Users\\wane1\\Downloads\\第二课各类数据结构及读取方法\\filtered_gene_bc_matrices\\hg19\\")
# 自动读取10X的数据，是一些tsv与mtx文件
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "biomamba")

#### 仅有一个稀疏矩阵时的读取方法#####
matrix_data <- read.table("C:\\Users\\wane1\\Downloads\\第二课各类数据结构及读取方法\\single_cell_datamatrix.txt", sep = "\t", header = T, row.names = 1)
dim(matrix_data)
seurat_obj <- CreateSeuratObject(counts = matrix_data)

###### 读取RDS文件############
# rm(list = ls())
pbmc <- readRDS("C:\\Users\\wane1\\Downloads\\第二课各类数据结构及读取方法\\panc8.rds")
# saveRDS(pbmc,"pbmc.rds")

#############
library(tidyverse)
str(pbmc)
library(mindr)
(out <- capture.output(str(pbmc)))
out2 <- paste(out, collapse = "\n")
mm(gsub("\\.\\.@", "# ", gsub("\\.\\. ", "#", out2)), type = "markdown", root = "Seurat")


# 单样本分析
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
rm(list = ls())
# 下载并解压数据
getwd()
setwd('C:\\Users\\wane1\\Downloads\\第三讲单样本分析\\')
download.file('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz','my.pbmc.gz')
untar(gunzip("my.pbmc.gz"))
# 读入数据并创建Seurat分析对象
pbmc.data <- Read10X(data.dir = "C:\\Users\\wane1\\Downloads\\第三讲单样本分析\\filtered_gene_bc_matrices\\hg19\\")
# Load the PBMC dataset
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

ncol(pbmc)
ncol(pbmc.data)
lalala <- as.data.frame(pbmc[["RNA"]]@counts)
write.table(lalala,'mycount.txt',sep = '\t')#表达矩阵可以这么存出来

# 质控数据及可视化
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
#鼠源的需要换成mt
head(pbmc@meta.data,5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
CombinePlots(plots = list(plot1, plot2)) 

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
CombinePlots(plots = list(plot1, plot2)) 

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
ncol(as.data.frame(pbmc[["RNA"]]@counts))

# 计算、分群
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#CLR、RC
#normalizes the feature expression measurements for each cell by the total expression
#pbmc[["RNA"]]@data

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#for PCA DoHeatmap
top10 <- head(VariableFeatures(pbmc), 10)# Identify the 10 most highly variable genes


plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
## When using repel, set xnudge and ynudge to 0 for optimal results
plot1 + plot2# plot variable features with and without labels

pbmc <- ScaleData(pbmc, features = rownames(pbmc))
## Centering and scaling data matrix
#pbmc <- ScaleData(pbmc) ##only4VariableFeatures

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)

ElbowPlot(pbmc) 

pbmc <- FindNeighbors(pbmc, dims = 1:10)
## Computing nearest neighbor graph
## Computing SNN
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 分群后的可视化：tsne+umap
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

# 寻找marker基因并对cluster进行重命名
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

if(!require(dplyr))install.packages("dplyr")
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#http://biocc.hrbmu.edu.cn/CellMarker/
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

sessionInfo()
saveRDS(pbmc,'pbmc.rds')
pbmc<- readRDS('pbmc.rds')
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

###########单纯的merge#################
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
##########准备用于拆分的数据集#########并
#pbmc <subset(pbmc,downsample 50)
ifnb <readRDS('pbmcrenamed.rds')
ifnb.list <Splitobject(ifnb,split.by "group")
C57 <ifnb.list$c57
AS1 <ifnb.list$AS1
######简单merge########
#不具有去批次效应功能
pbmc <-merge(C57,y c(AS1),add.cell.ids c("C57","AS1"),