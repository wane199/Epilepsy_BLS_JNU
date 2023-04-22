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
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# 1.4 去除聚类高度高于15的样本
# Plot a line to show the cut
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)+abline(h = 15, col = "red")

## integer(0)
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datExpr[1:5,1:5]

# 1.5 加载注释信息
traitData = read.csv("C:\\Users\\wane1\\Downloads\\testdata\\ClinicalTraits.csv", sep = ";")
dim(traitData)#有三十多列注释
traitData[1:5,1:5]
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36) ]
dim(allTraits)
names(allTraits)#取出需要的
# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr)#取出样本名
traitRows = match(femaleSamples, allTraits$Mice)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
datTraits[1:5,1:5]#行名是样本名，每列是一种注释
collectGarbage()

sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
#给注释中的数字变量赋予颜色，感觉与我的hcluster热图无异
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
