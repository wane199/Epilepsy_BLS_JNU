#### Rcpi: R/Bioconductor Package as an Integrated Informatics Platform for Drug Discovery ####
# https://www.bioconductor.org/packages/devel/bioc/vignettes/Rcpi/inst/doc/Rcpi.html#Applications_in_chemogenomics
##### Installation #####
install.packages("BiocManager")
BiocManager::install("Rcpi")
BiocManager::install("Rcpi", dependencies = c("Imports", "Enhances"))

##### Predicting protein subcellular localization #####
library("Rcpi")

# Load FASTA files
extracell <- readFASTA(system.file(
  "vignettedata/extracell.fasta",
  package = "Rcpi"
))
mitonchon <- readFASTA(system.file(
  "vignettedata/mitochondrion.fasta",
  package = "Rcpi"
))

length(extracell)

length(mitonchon)

extracell <- extracell[(sapply(extracell, checkProt))]
mitonchon <- mitonchon[(sapply(mitonchon, checkProt))]
length(extracell)

length(mitonchon)

# Calculate APAAC descriptors
x1 <- t(sapply(extracell, extractProtAPAAC))
x2 <- t(sapply(mitonchon, extractProtAPAAC))
x <- rbind(x1, x2)

# Make class labels
labels <- as.factor(c(rep(0, length(extracell)), rep(1, length(mitonchon))))

set.seed(123)

# Split training and test set
tr.idx <- c(
  sample(1:nrow(x1), round(nrow(x1) * 0.75)),
  sample(nrow(x1) + 1:nrow(x2), round(nrow(x2) * 0.75))
)
te.idx <- setdiff(1:nrow(x), tr.idx)

x.tr <- x[tr.idx, ]
x.te <- x[te.idx, ]
y.tr <- labels[tr.idx]
y.te <- labels[te.idx]

library("randomForest")
rf.fit <- randomForest(x.tr, y.tr, cv.fold = 5)
print(rf.fit)

# Predict on test set
rf.pred <- predict(rf.fit, newdata = x.te, type = "prob")[, 1]

# Plot ROC curve
library("pROC")
plot.roc(y.te, rf.pred, grid = TRUE, print.auc = TRUE)

#### Applications in chemoinformatics ####
##### Regression modeling in QSRR study of retention indices #####
library("Rcpi")

RI.smi <- system.file(
  "vignettedata/RI.smi",
  package = "Rcpi"
)
RI.csv <- system.file(
  "vignettedata/RI.csv",
  package = "Rcpi"
)

x.mol <- readMolFromSmi(RI.smi, type = "mol")
x.tab <- read.table(RI.csv, sep = "\t", header = TRUE)
y <- x.tab$RI

# Calculate selected molecular descriptors
x <- suppressWarnings(cbind(
  extractDrugALOGP(x.mol),
  extractDrugApol(x.mol),
  extractDrugECI(x.mol),
  extractDrugTPSA(x.mol),
  extractDrugWeight(x.mol),
  extractDrugWienerNumbers(x.mol),
  extractDrugZagrebIndex(x.mol)
))

# Run regression on training set
library("caret")
library("pls")

# Cross-validation settings
ctrl <- trainControl(
  method = "repeatedcv", number = 5, repeats = 10,
  summaryFunction = defaultSummary
)

# Train a PLS model
set.seed(123)
pls.fit <- train(
  x, y,
  method = "pls", tuneLength = 10, trControl = ctrl,
  metric = "RMSE", preProc = c("center", "scale")
)

# Print cross-validation results
print(pls.fit)

# Number of components vs. RMSE
print(plot(pls.fit, asp = 0.5))

# Plot experimental RIs vs predicted RIs
plot(y, predict(pls.fit, x),
     xlim = range(y), ylim = range(y),
     xlab = "Experimental RIs", ylab = "Predicted RIs"
)
abline(a = 0, b = 1)

##### In silico toxicity classification for drug discovery #####
library("Rcpi")

fdamdd.smi <- system.file("vignettedata/FDAMDD.smi", package = "Rcpi")
fdamdd.csv <- system.file("vignettedata/FDAMDD.csv", package = "Rcpi")

x.mol <- readMolFromSmi(fdamdd.smi, type = "mol")
x.smi <- readMolFromSmi(fdamdd.smi, type = "text")
y <- as.factor(paste0("class", scan(fdamdd.csv)))

# Calculate molecular fingerprints
x1 <- extractDrugEstateComplete(x.mol)
x2 <- extractDrugMACCSComplete(x.mol)
x3 <- extractDrugOBFP4(x.smi, type = "smile")

library("caret")

# Remove near-zero variance variables
x1 <- x1[, -nearZeroVar(x1)]
x2 <- x2[, -nearZeroVar(x2)]
x3 <- x3[, -nearZeroVar(x3)]

# Split training and test set
set.seed(123)
tr.idx <- sample(1:nrow(x1), round(nrow(x1) * 0.75))
te.idx <- setdiff(1:nrow(x1), tr.idx)
x1.tr <- x1[tr.idx, ]
x1.te <- x1[te.idx, ]
x2.tr <- x2[tr.idx, ]
x2.te <- x2[te.idx, ]
x3.tr <- x3[tr.idx, ]
x3.te <- x3[te.idx, ]
y.tr <- y[tr.idx]
y.te <- y[te.idx]

# SVM classification on training sets
library("kernlab")

# Cross-validation settings
ctrl <- trainControl(
  method = "cv", number = 5, # repeats = 100,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# SVM with RBF kernel
svm.fit1 <- train(
  x1.tr, y.tr,
  method = "svmRadial", trControl = ctrl,
  metric = "ROC", preProc = c("center", "scale")
)
svm.fit2 <- train(
  x2.tr, y.tr,
  method = "svmRadial", trControl = ctrl,
  metric = "ROC", preProc = c("center", "scale")
)
svm.fit3 <- train(
  x3.tr, y.tr,
  method = "svmRadial", trControl = ctrl,
  metric = "ROC", preProc = c("center", "scale")
)

# Print cross-validation results
print(svm.fit1)
print(svm.fit2)
print(svm.fit3)

# Predict on test set
svm.pred1 <- predict(svm.fit1, newdata = x1.te, type = "prob")[, 1]
svm.pred2 <- predict(svm.fit2, newdata = x2.te, type = "prob")[, 1]
svm.pred3 <- predict(svm.fit3, newdata = x3.te, type = "prob")[, 1]

# Generate colors
library("RColorBrewer")
pal <- brewer.pal(3, "Set1")

# ROC curves of different fingerprints
library("pROC")
plot(smooth(roc(y.te, svm.pred1)), col = pal[1], grid = TRUE)
plot(smooth(roc(y.te, svm.pred2)), col = pal[2], grid = TRUE, add = TRUE)
plot(smooth(roc(y.te, svm.pred3)), col = pal[3], grid = TRUE, add = TRUE)

##### Clustering of molecules based on structural similarities #####
library("Rcpi")
mols <- readMolFromSDF(system.file(
  "compseq/tyrphostin.sdf",
  package = "Rcpi"
))

simmat <- diag(length(mols))

for (i in 1:length(mols)) {
  for (j in i:length(mols)) {
    fp1 <- extractDrugEstate(mols[[i]])
    fp2 <- extractDrugEstate(mols[[j]])
    tmp <- calcDrugFPSim(fp1, fp2, fptype = "compact", metric = "tanimoto")
    simmat[i, j] <- tmp
    simmat[j, i] <- tmp
  }
}

mol.hc <- hclust(as.dist(1 - simmat), method = "ward.D")

library("ape") # Tree visualization of clusters
clus5 <- cutree(mol.hc, 5) # Cut dendrogram into 5 clusters

# Generate colors
library("RColorBrewer")
pal5 <- brewer.pal(5, "Set1")
plot(as.phylo(mol.hc),
     type = "fan",
     tip.color = pal5[clus5],
     label.offset = 0.1, cex = 0.7
)

##### Structure-based chemical similarity searching #####
library("Rcpi")

mol <- system.file("compseq/DB00530.sdf", package = "Rcpi")
moldb <- system.file("compseq/tyrphostin.sdf", package = "Rcpi")

rank1 <- searchDrug(
  mol, moldb,
  cores = 4, method = "fp",
  fptype = "maccs", fpsim = "tanimoto"
)
rank2 <- searchDrug(
  mol, moldb,
  cores = 4, method = "fp",
  fptype = "fp2", fpsim = "cosine"
)
rank3 <- searchDrug(
  mol, moldb,
  cores = 4, method = "mcs",
  mcssim = "tanimoto"
)

head(rank1)

head(rank2)

head(rank3)

# Convert SDF format to SMILES format
convMolFormat(
  infile = mol, outfile = "DB00530.smi", from = "sdf", to = "smiles"
)
convMolFormat(
  infile = moldb, outfile = "tyrphostin.smi", from = "sdf", to = "smiles"
)

smi1 <- readLines("DB00530.smi")
smi2 <- readLines("tyrphostin.smi")[92] # Select the #92 molecule
calcDrugMCSSim(smi1, smi2, type = "smile", plot = TRUE)

#### Applications in chemogenomics ####
##### Predicting drug-target interaction by integrating chemical and genomic spaces #####
library("Rcpi")

gpcr <- read.table(system.file(
  "vignettedata/GPCR.csv",
  package = "Rcpi"
),
header = FALSE, as.is = TRUE
)

head(gpcr)

library("igraph")
library("reshape")
# remotes::install_github("gastonstat/arcdiagram")
library("arcdiagram")

g <- graph.data.frame(gpcr[1:(nrow(gpcr) / 2), ], directed = FALSE)
edgelist <- get.edgelist(g)
vlabels <- V(g)$name
vgroups <- c(rep(0, 95), rep(1, 223))
vfill <- c(rep("#8B91D4", 95), rep("#B2C771", 223))
vborders <- c(rep("#6F74A9", 95), rep("#8E9F5A", 223))
degrees <- degree(g)

xx <- data.frame(vgroups, degrees, vlabels, ind = 1:vcount(g))
yy <- arrange(xx, desc(vgroups), desc(degrees))
new.ord <- yy$ind

arcplot(
  edgelist,
  ordering = new.ord, labels = vlabels,
  cex.labels = 0.1, show.nodes = TRUE,
  col.nodes = vborders, bg.nodes = vfill,
  cex.nodes = log10(degrees) + 0.1,
  pch.nodes = 21, line = -0.5, col.arcs = hsv(0, 0, 0.2, 0.25)
)

library("Rcpi")

gpcr <- read.table(system.file(
  "vignettedata/GPCR.csv",
  package = "Rcpi"
),
header = FALSE, as.is = TRUE
)

protid <- unique(gpcr[, 1])
drugid <- unique(gpcr[, 2])

protseq <- getSeqFromKEGG(protid, parallel = 5)
drugseq <- getSmiFromKEGG(drugid, parallel = 50)

x0.prot <- cbind(
  t(sapply(unlist(protseq), extractProtAPAAC)),
  t(sapply(unlist(protseq), extractProtCTriad))
)

x0.drug <- cbind(
  extractDrugEstateComplete(readMolFromSmi(textConnection(drugseq))),
  extractDrugMACCSComplete(readMolFromSmi(textConnection(drugseq))),
  extractDrugOBFP4(drugseq, type = "smile")
)

# Generate drug x / protein x / y
x.prot <- matrix(NA, nrow = nrow(gpcr), ncol = ncol(x0.prot))
x.drug <- matrix(NA, nrow = nrow(gpcr), ncol = ncol(x0.drug))
for (i in 1:nrow(gpcr)) x.prot[i, ] <- x0.prot[which(gpcr[, 1][i] == protid), ]
for (i in 1:nrow(gpcr)) x.drug[i, ] <- x0.drug[which(gpcr[, 2][i] == drugid), ]

y <- as.factor(c(rep("pos", nrow(gpcr) / 2), rep("neg", nrow(gpcr) / 2)))

x <- getCPI(x.prot, x.drug, type = "combine")

##### Compound-protein interaction (CPI) descriptors #####
# Protein-protein interaction (PPI) descriptors
library("caret")
x <- x[, -nearZeroVar(x)]

# Cross-validation settings
ctrl <- trainControl(
  method = "cv", number = 5, repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Train a random forest classifier
library("randomForest")

set.seed(1006)
rf.fit <- train(
  x, y,
  method = "rf", trControl = ctrl,
  metric = "ROC", preProc = c("center", "scale")
)

print(rf.fit)

rf.pred <- predict(rf.fit$finalModel, x, type = "prob")[, 1]

library("pROC")
plot(smooth(roc(y, rf.pred)), grid = TRUE, print.auc = TRUE)












































