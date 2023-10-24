#### MR之后生信分析 ####
## MR之后根据目标遗传变异SNP寻找HUB基因，利用基因疾病关联数据库，进行疾病相关基因的研究。
## 做生物信息分析，然后研究功能通路的改变，再用实验验证，本节首先绘制基因-疾病网络图、进行富集分析等，
# https://www.disgenet.org/static/disgenet2r/disgenet2r.html
##### Installation and first run #####
# library(devtools)
# install_bitbucket("ibi_group/disgenet2r")
# 
rm(list = ls())
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
  email = "wane199@outlook.com",
  password = "wu19940101"
)
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

##### Retrieving Gene-Disease Associations from DisGeNET #####
###### Searching by gene ######
data1 <- gene2disease(
  gene = 6334, vocabulary = "ENTREZ",
  database = "ALL"
)
class(data1)
data1

results <- extract(data1)
head(results, 3)
###### Visualizing the diseases associated to a single gene ######
# Figure 1: The Gene-Disease Network for the SCN8A gene 
plot(data1,
  class = "Network",
  prop = 20
)
# Figure 2: The Disease Class Network for the SCN8A Gene
plot(data1,
  class = "DiseaseClass",
  prop = 3
)

###### Exploring the evidences associated to a gene ######
data1 <- gene2evidence(
  gene = "SCN8A",
  vocabulary = "HGNC",
  disease = "C0014544",
  database = "ALL",
  score = c(0.3, 1)
)
data1
results <- extract(data1)
head(results)

myListOfGenes <- c("SCN8A", "GRIN2B", "ABCB1", "PCDH19", "SLC12A5")
data2 <- gene2disease(
  gene = myListOfGenes,
  score = c(0.2, 1),
  verbose = TRUE
)
data2

# Visualizing the diseases associated to multiple genes
# Figure 3: The Gene-Disease Network for a list of genes belonging to the voltage-gated potassium channel family
plot(data2,
  class = "Network",
  prop = 10
)

# Figure 4: The Gene-Disease Heatmap for a list of genes belonging to the voltage-gated potassium channel family
plot(data2,
  class = "Heatmap",
  limit = 100, nchars = 50
)

# Figure 5: The Gene-Disease Class Heatmap for a list of genes belonging to the voltage-gated potassium channel family
plot(data2,
  class = "DiseaseClass", nchars = 60
)

# Searching by disease
data3 <- disease2gene(
  disease = "C0014544",
  database = "CURATED", # CURATED/ALL
  score = c(0.4, 1)
)
data3

data4 <- disease2gene(
  disease = "181500", vocabulary = "OMIM",
  database = "CURATED",
  score = c(0.4, 1)
)
data4 <- disease2gene(
  disease = "345.9", vocabulary = "ICD9CM",
  database = "CURATED",
  score = c(0.4, 1)
)
data4 <- disease2gene(
  disease = "C3020", vocabulary = "NCI",
  database = "CURATED",
  score = c(0.4, 1)
)
data4
# Visualizing the genes associated to a single disease
# Figure 6: The Gene-Disease Network for genes associated to Epilepsy
plot(data3,
  prop = 10
)
# Figure 7: The Protein Class-Disease Network for genes associated to Epilepsy
plot(data3,
  class = "ProteinClass"
)

##### Exploring the evidences associated to a disease #####
data3 <- disease2evidence(
  disease = "C0014544",
  type = "GDA",
  database = "CURATED",
  score = c(0.4, 1)
)
data3
data3 <- disease2evidence(
  disease = "C0014544",
  gene = c("SCN8A", "GRIN2B"),
  type = "GDA",
  database = "CURATED",
  score = c(0.4, 1)
)
data3
results <- extract(data3)

###### Searching multiple diseases ######
diseasesOfInterest <- c("C0036341", "C0002395", "C0030567", "C0005586", "C0014544")
data5 <- disease2gene(
  disease = diseasesOfInterest,
  database = "CURATED",
  score = c(0.4, 1),
  verbose = TRUE
)
###### Visualizing the genes associated to multiple diseases ######
plot(data5,
  class = "Network",
  prop = 10
)
plot(data5,
  class = "Heatmap",
  limit = 30,
  cutoff = 0.2
)
# Figure 10: The Protein Class-Disease Heatmap for genes associated to a list of diseases
plot(data5,
  class = "ProteinClass"
)

##### Retrieving Variant-Disease Associations from DisGeNET #####
###### Searching by variant ######
data6 <- variant2disease( variant= "rs2844697", # rs2844697、rs9273368
                          database = "ALL")
data6

###### Visualizing the diseases associated to a single variant ######
# Figure 11: The Variant-Disease Network for the variant rs121913279
plot( data6, 
      class = "Network",
      prop  = 10)

# Figure 12: The Variant- Disease Class Network for the variant rs121913279
plot( data6,
      class = "DiseaseClass",
      prop = 3)

##### Exploring the evidences associated to a variant #####
###### Searching multiple variants ######
data7 <- variant2disease(
  variant = c("rs121913013", "rs1060500621",
              "rs199472709", "rs72552293",
              "rs74315445", "rs199472795"),
  database = "ALL")

###### Visualizing the diseases associated to multiple variants ######
plot( data7,
      class = "Network",
      prop = 10)














##### Get DisGeNET data version #####
get_disgenet_version()








