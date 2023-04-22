# https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499086&idx=1&sn=85f1a292ab82fc3c5980c82069057a79&chksm=9b3c841eac4b0d08bb8b821f67b2740a677b040d310e414c28e1fec7bd2a5cf4e0cef4ba0b0f&mpshare=1&scene=1&srcid=04225hR6eUh2LxTU03RZyBTz&sharer_sharetime=1682147335292&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#.rd
rm(list = ls())
getwd()
# 预处理
# 1.1 读取数据
if (!require(WGCNA)) BiocManager::install("WGCNA")
enableWGCNAThreads() # 查看WGCNA可以调用多少线程工作
## Allowing parallel execution with up to 7 working processes.
options(stringsAsFactors = FALSE)
# Read in the female liver data set
femData <- read.csv("C:\\Users\\wane1\\Downloads\\testdata\\LiverFemale3600.csv", sep = ";")
# Take a quick look at what is in the data set:
dim(femData)
names(femData)[1:10]
femData[1:10, 1:10] # 这个数据前面几列是包含基因注释信息的，需要整理成行名是样本，列名是基因名的矩阵
datExpr0 <- as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) <- femData$substanceBXH
rownames(datExpr0) <- names(femData)[-c(1:8)]
datExpr0[1:5, 1:5]

# 1.2 检查是否包含缺失值
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # 返回TRUE则不包含缺失值
# 反之需要删除有问题的样品与基因
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  }
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr0[1:5, 1:5]

# 1.3 聚类 观察是否有离群值
sampleTree <- hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12, 9)
# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree,
  main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2
)
# 1.4 去除聚类高度高于15的样本
# Plot a line to show the cut
plot(sampleTree,
  main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2
) + abline(h = 15, col = "red")

## integer(0)
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
datExpr[1:5, 1:5]

# 1.5 加载注释信息
traitData <- read.csv("C:\\Users\\wane1\\Downloads\\testdata\\ClinicalTraits.csv", sep = ";")
dim(traitData) # 有三十多列注释
traitData[1:5, 1:5]
names(traitData)
# remove columns that hold information we do not need.
allTraits <- traitData[, -c(31, 16)]
allTraits <- allTraits[, c(2, 11:36)]
dim(allTraits)
names(allTraits) # 取出需要的
# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples <- rownames(datExpr) # 取出样本名
traitRows <- match(femaleSamples, allTraits$Mice)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]
datTraits[1:5, 1:5] # 行名是样本名，每列是一种注释
collectGarbage()

sampleTree2 <- hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = FALSE)
# 给注释中的数字变量赋予颜色，感觉与我的hcluster热图无异
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
  groupLabels = names(datTraits),
  main = "Sample dendrogram and trait heatmap"
)

# 构建network、识别module
# 2.1 自动化分块式
# 2.1.1 确定power（β）
# Choosing the soft-thresholding power: analysis of network topology
# Constructing a weighted gene network entails the choice of the soft thresholding power
# to which co-expression similarity is raised to calculate adjacency
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2)) # 给定一些值用于测试
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
  main = paste("Scale independence")
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers, cex = cex1, col = "red"
)
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
  main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
# 要选择lowest power for which the scale-free topology index curve  attens out upon reaching a high value
sft$powerEstimate

# 2.1.2 一步法网络构建
net <- blockwiseModules(datExpr,
  power = sft$powerEstimate,
  TOMType = "unsigned",
  minModuleSize = 30, # 最小模块的基因数
  reassignThreshold = 0,
  mergeCutHeight = 0.25, # threshold for merging of modules
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "femaleMouseTOM",
  verbose = 3,
  maxBlockSize = ncol(datExpr) # 可以处理最大基因的数量，数据集大于这一参数会break
  # 对应内存需求如下
  # 16GB workstation should handle up to 20000 probes;
  # 32GB workstation should handle perhaps 30000.
  # 4GB standard  8000-10000 probes,
) # 这一步参数众多，需要的话可以?一下

# 返回值中net$colors 包含了模块的分配情况
# net$MEs 包含模块的特征基因.
table(net$colors) # 看看这些基因被分为了哪些模块
# 注意，0是不属于任何模块的基因，也就是后面画出灰色的部分，其余模块颜色按模块基因从大到小排列
net$dendrograms[[1]] # hierarchical clustering dendrogram存放位置

mergedColors <- labels2colors(net$colors) # 把数字返回成颜色
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)

# 2.1.3 修改、存储
# recutBlockwiseTrees可以用于重新计算，比从头计算节省时间
# 存储文件用于下游分析
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
  file = "C:\\Users\\wane1\\Downloads\\testdata\\FemaleLiver-02-networkConstruction-auto.RData"
)

# 2.2 步进法 有助于自定义各类参数
# 2.2.1 还是一样要确定β
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))+
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+abline(h=0.90,col="red")

# this line corresponds to using an R^2 cut-off of h
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))+
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 2.2.2 adjacency
adjacency = adjacency(datExpr, power = sft$powerEstimate)
adjacency[1:5,1:5]

# 2.2.3 Topological Overlap Matrix (TOM)
#去噪、防止假关联的出现将adjacency转换为Topological Overlap Matrix（反应similarity）, 并用1-TOM反应dissimilarity:
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
TOM[1:5,1:5]

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")#利用dissTOM聚类更快
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)#这种聚类树每个leaf代表一个基因

# 2.2.4 色块聚类树
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)#每个模块最小的基因数量
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)# Convert numeric lables into colors
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 2.2.5 合并相似表达的模块
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEs[1:5,1:5]

MEDiss = 1-cor(MEs)# Calculate dissimilarity of module eigengenes,注意计算的是不详细性
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")+abline(h=0.25, col = "red")

####开始merge
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)

#merge为返回的新的WGCNA对象
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)#查看合并前后的模块聚类图

# 2.2.6 保存核心文件用于下游分析
moduleColors = mergedColors# Rename to moduleColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "C:\\Users\\wane1\\Downloads\\testdata\\FemaleLiver-02-networkConstruction-stepByStep.RData")
