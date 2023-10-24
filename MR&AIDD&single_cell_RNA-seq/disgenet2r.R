#### MR之后生信分析 ####
## MR之后根据目标遗传变异SNP寻找HUB基因，利用基因疾病关联数据库，进行疾病相关基因的研究。
## 做生物信息分析，然后研究功能通路的改变，再用实验验证，本节首先绘制基因-疾病网络图、进行富集分析等，
# https://www.disgenet.org/static/disgenet2r/disgenet2r.html
##### Installation and first run #####
# library(devtools)
# install_bitbucket("ibi_group/disgenet2r")
rm(list = ls())
library(disgenet2r) # [bitbucket::ibi_group/disgenet2r] v0.99.3
disgenet_api_key <- get_disgenet_api_key(
  email = "wane199@outlook.com",
  password = "wu19940101"
)
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)

#### Retrieving Gene-Disease Associations from DisGeNET ####
##### Searching by gene #####
# https://www.ncbi.nlm.nih.gov/gene/
data1 <- gene2disease(
  gene = 6323, vocabulary = "ENTREZ", # SCN1A
  database = "ALL"
)
class(data1)
data1

results <- extract(data1)
head(results, 3)
###### Visualizing the diseases associated to a single gene ######
# Figure 1: The Gene-Disease Network for the SCN1A gene
plot(data1,
  class = "Network",
  prop = 20
)
# Figure 2: The Disease Class Network for the SCN1A Gene
plot(data1,
  class = "DiseaseClass",
  prop = 3
)

###### Exploring the evidences associated to a gene ######
data1 <- gene2evidence(
  gene = "SCN1A",
  vocabulary = "HGNC",
  disease = "C0014544", # MeSH
  database = "ALL",
  score = c(0.3, 1)
)
data1
results <- extract(data1)
head(results)

##### Searching multiple genes #####
myListOfGenes <- c("SCN1A", "TBC1D24", "El1", "SAMD12", "STARD7")
data2 <- gene2disease(
  gene = myListOfGenes,
  score = c(0.2, 1),
  verbose = TRUE
)
data2

###### Visualizing the diseases associated to multiple genes ######
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
  class = "DiseaseClass", limit = 10, nchars = 60
)

##### Searching by disease #####
# https://www.ncbi.nlm.nih.gov/mesh/
data3 <- disease2gene(
  disease = "D004827",
  vocabulary = "MESH",
  database = "ALL", # CURATED/ALL
  score = c(0.4, 1)
)
data3

data4 <- disease2gene(
  disease = "2538", vocabulary = "DO",
  database = "ALL",
  score = c(0.4, 1)
)
data4 <- disease2gene(
  disease = "345.1", vocabulary = "ICD10", # http://icd9.chrisendres.com/index.php
  database = "ALL",
  score = c(0.4, 1)
)
data4 <- disease2gene(
  disease = "G40.8", vocabulary = "ICD10", # https://www.medsci.cn/sci/icd-10.do
  database = "ALL",
  score = c(0.4, 1)
)
data4

###### Visualizing the genes associated to a single disease ######
# Figure 6: The Gene-Disease Network for genes associated to Epilepsy
plot(data3,
  prop = 10
)
# Figure 7: The Protein Class-Disease Network for genes associated to Epilepsy
plot(data3,
  class = "ProteinClass"
)

###### Exploring the evidences associated to a disease ######
data3 <- disease2evidence(
  disease = "2538",
  type = "GDA",
  database = "CURATED",
  score = c(0.4, 1)
)
data3
data3 <- disease2evidence(
  disease = "D004827",
  gene = c("SCN1A", "TBC1D24"),
  type = "GDA",
  database = "CURATED",
  score = c(0.4, 1)
)
data3
results <- extract(data3)

##### Searching multiple diseases #####
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

#### Retrieving Variant-Disease Associations from DisGeNET ####
##### Searching by variant #####
data6 <- variant2disease(
  variant = "rs2844697", # rs2844697、rs9273368
  database = "ALL"
)
data6

###### Visualizing the diseases associated to a single variant ######
# Figure 11: The Variant-Disease Network for the variant rs121913279
plot(data6,
  class = "Network",
  prop  = 10
)

# Figure 12: The Variant- Disease Class Network for the variant rs121913279
plot(data6,
  class = "DiseaseClass",
  prop = 3
)

###### Exploring the evidences associated to a variant ######
# variant2evidence()

##### Searching multiple variants #####
data7 <- variant2disease(
  variant = c(
    "rs121913013", "rs1060500621",
    "rs199472709", "rs72552293",
    "rs74315445", "rs199472795"
  ),
  database = "ALL"
)

###### Visualizing the diseases associated to multiple variants ######
# Figure 13: The Variant-Disease Network for a list of variants
plot(data7,
  class = "Network",
  prop = 10
)
# Figure 14: The Variant-Gene-Disease Network for a list of variants
plot(data7,
  class = "Network",
  showGenes = T
)
# Figure 15: The Variant-Disease Heatmap for a list of variants
plot(data7,
  class = "Heatmap",
  prop = 10
)
# Figure 16: The Variant-Disease Class Heatmap for a list of variants
plot(data7,
  class = "DiseaseClass"
)

##### Searching by disease #####
data8 <- disease2variant(
  disease = c("C0752166"),
  database = "CLINVAR",
  score = c(0, 1)
)
data8
###### Visualizing the variants associated to a single disease ######
# Figure 17: The Variant-Disease Network for a single disease
plot(data8,
  class = "Network"
)

###### Explore the evidences associated to a single disease ######
data8 <- disease2evidence(
  disease = "C0036341",
  variant = "rs1344706",
  type = "VDA",
  year = c(2019, 2020),
  database = "ALL",
  score = c(0.1, 1)
)
data8

data8 <- disease2evidence(
  disease = "C0036341",
  gene = "CACNA1C",
  type = "VDA",
  database = "BEFREE",
  score = c(0.1, 1)
)
data8

##### Searching multiple diseases #####
data8 <- disease2variant(
  disease = c("C3150943", "C1859062", "C2678485", "C4015695"),
  database = "CURATED",
  score = c(0.75, 1)
)
data8

###### Visualizing the variants associated to multiple diseases ######
# Figure 18: The Variant-Disease Network for a list of diseases
plot(data8,
  class = "Network"
)

# Figure 19: The Variant-Gene-Disease Network for a list of diseases
plot(data8,
  class = "Network",
  showGenes = T
)

# Figure 20: The Variant-Disease Heatmap for a list of diseases
plot(data8,
  class = "Heatmap",
  prop = 10
)

#### Retrieving Disease-Disease Associations from DisGeNET ####
##### Searching DDAs by genes for a single disease #####
data9 <- disease2disease_by_gene(
  disease = "C0010674",
  database = "CURATED", ndiseases = 5
)
data9

qr <- extract(data9)
head(qr[c("disease1_name", "disease2_name", "jaccard_genes", "ngenes", "pvalue")])

# Figure 21: The Disease-Disease Network by shared genes for MELAS Syndrome
plot(data9,
  class = "Network",
  layout = "layout.lgl",
  prop = 0.4
)

# Figure 22: The Diseasome plot for MELAS Syndrome
plot(data9,
  class = "Diseasome",
  layout = "layout.lgl",
  prop = 5
)

# Figure 23: Disease-Disease Barplot visualizing the top 5 diseases associated to MELAS Syndrome
plot(data9,
  class = "Barplot"
)

##### Searching DDAs via genes for multiple diseases #####
###### Visualizing the diseases associated to multiple diseases ######
diseasesOfInterest <- c("C0013182", "C0013221", "C3658290", "C0860207", "C1274933")
data10 <- disease2disease_by_gene(
  disease = diseasesOfInterest,
  database = "CURATED", ndiseases = 5
)

# Figure 24: The Disease-Disease Network by shared genes for a list of diseases
plot(data10,
  class = "Network",
  prop  = 0.1
)

# Figure 25: The Diseasome plot for a list of diseases
plot(data10,
  class = "Diseasome",
  prop  = 0.1
)

# Figure 26: The Venn Diagram for 3 diseases
plot(data10,
  class = "Venn"
)

# Figure 27: The Venn Diagram for a list of 4 diseases
data11 <- disease2disease_by_gene(
  disease = c("C0018801", "C0028754", "C0004153", "C0011849"),
  database = "CURATED",
  ndiseases = 5
)
plot(data11,
  class = "Venn"
)

# Figure 28: The Venn Diagram for a list of 5 diseases
data11 <- disease2disease_by_gene(
  disease = c("C0018801", "C0028754", "C0004153", "C0011849", "C0021368"),
  database = "CURATED",
  ndiseases = 5
)
plot(data11,
  class = "Venn"
)

###### Searching DDAs via shared variants for a single disease ######
data12 <- disease2disease_by_variant(
  disease = "C0011860",
  database = "CURATED", ndiseases = 3
)
data12

# Figure 29: An example of Disease-Disease Network
plot(data12,
  layout = "layout.lgl",
  prop = 0.2
)

data13 <- disease2disease_by_variant(
  disease = c("C0018801", "C0020538"),
  database = "CURATED", ndiseases = 5
)

data13
plot(data13,
  layout = "layout.lgl",
  class = "Diseasome",
  prop = 0.1
)

##### Performing a disease enrichment #####
list_of_genes <- disease2gene(disease = "C0004352", database = "GENOMICS_ENGLAND")

list_of_genes <- list_of_genes@qresult$gene_symbol

res_enrich <- disease_enrichment(
  entities = list_of_genes, vocabulary = "HGNC",
  database = "CURATED"
)

table1 <- res_enrich@qresult[1:10, c("Description", "FDR", "Ratio", "BgRatio")]

# Figure : The Enrichment plot for a list of genes
plot(res_enrich, class = "Enrichment", count = 3, cutoff = 0.05, nchars = 70)

res_enrich <- disease_enrichment(
  entities = list_of_genes, vocabulary = "HGNC",
  database = "ALL"
)

table1 <- res_enrich@qresult[1:10, c("Description", "FDR", "Ratio", "BgRatio")]

list_of_variants <- disease2variant(disease = "C0004352", database = "CLINVAR")

list_of_variants <- as.character(list_of_variants@qresult$variantid)

res_enrich <- disease_enrichment(
  entities = list_of_variants, vocabulary = "DBSNP",
  database = "CURATED"
)

table1 <- res_enrich@qresult[1:10, c("Description", "FDR", "Ratio", "BgRatio")]

# Figure : The Enrichment plot for a list of genes
plot(res_enrich, class = "Enrichment", count = 6, cutoff = 0.05)

##### Expanding DisGeNET data with other resources available in the LOD cloud #####
##### Retrieving the pathways associated to a disease #####
dis2path <- disease2pathway(
  disease = "C0018801",
  database = "ALL", score = c(0.5, 1)
)
dis2path
# head(dis2path@qresult)
# qr <- extract(dis2path)
# head(qr, 3)

##### Retrieving the diseases associated to a given pathway #####
path2dis <- pathway2disease(
  pathway = "WP1591",
  database = "CURATED",
  score = c(0.9, 1)
)
path2dis
head(path2dis@qresult)
qr <- extract(path2dis)
head(qr, 3)

##### Retrieve the drug targets for a disease #####
library("SPARQL") # SPARQL_1.16.tar.gz
dis2com <- disease2compound(
  disease = c("C0020538"),
  database = "CURATED"
)
dis2com
qr <- extract(dis2com)
head(qr, 3)

##### Retrieving Gene Ontology data for a disease #####
##### Retrieving the biological processes associated to a disease #####
disease2bp <- disease2biologicalprocess(
  disease = "C0036341",
  database = "CURATED",
  score = c(0.6, 1)
)
disease2bp
disease2bp <- extract(disease2bp)
head(disease2bp[c("gene_product_label", "go_label")])

##### Retrieving the molecular functions of the genes associated to a disease #####
disease2mf <- disease2molecularfunction(
  disease = "C0002395",
  database = "UNIPROT"
)

disease2mf
disease2mf <- extract(disease2mf)
head(disease2mf[c("gene_product_label", "go_label")])

##### Retrieving the molecular functions of the genes associated to a disease #####
disease2cc <- disease2cellcomponent(
  disease = "C0002395",
  database = "MGD"
)

disease2cc <- extract(disease2cc)
head(disease2cc[c("gene_product_label", "go_label")])

##### Get DisGeNET data version #####
get_disgenet_version()
