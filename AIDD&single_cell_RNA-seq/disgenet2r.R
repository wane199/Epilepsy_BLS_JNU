# https://www.disgenet.org/static/disgenet2r/disgenet2r.html
# Installation and first run
library(devtools)
install_bitbucket("ibi_group/disgenet2r")

library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
  email = "wane199@outlook.com",
  password = "wu19940101"
)
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

# Retrieving Gene-Disease Associations from DisGeNET
# Searching by gene
data1 <- gene2disease(
  gene = 6334, vocabulary = "ENTREZ",
  database = "ALL"
)
class(data1)
data1

results <- extract(data1)
head(results, 3)
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

# Exploring the evidences associated to a gene
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

# Exploring the evidences associated to a disease
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

# Searching multiple diseases
diseasesOfInterest <- c("C0036341", "C0002395", "C0030567", "C0005586", "C0014544")
data5 <- disease2gene(
  disease = diseasesOfInterest,
  database = "CURATED",
  score = c(0.4, 1),
  verbose = TRUE
)
# Visualizing the genes associated to multiple diseases
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

# Retrieving Variant-Disease Associations from DisGeNET
# Searching by variant
data6 <- variant2disease( variant= "rs121913279",
                          database = "CURATED")
data6










