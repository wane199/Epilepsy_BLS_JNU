##### TwoSample MR Practical #####
# https://zhuanlan.zhihu.com/p/605966100
# 判断是否已经安装了“pacman”包I如果没有就安装它
if (!require("pacman")) install.packages("pacman", update = F, ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/")
# 加载“pacman”包，用于方便加载其他的R包
library("pacman") # Package Management Tool
# 使用“p_load”函数加载所需的R包
p_load(
  ggplot2, # CRAN v3.4.3
  readr # CRAN v2.1.4
)
p_load_gh("MRCIEU/MRInstruments")


# 在加载这些包之前，“pacman“会先检查是否已经安装，如果没有会自动安装
install.packages(c("devtools", "knitr", "rmarkdown"))
library(devtools) # Tools to Make Developing R Packages Easier
install_github(c("MRCIEU/TwoSampleMR", "MRCIEU/MRInstruments"))
install_github("WSpiller/MRPracticals", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
install_github("WSpiller/MRPracticals", build_opts = c("--no-resave-data", "--no-manual"))
library(TwoSampleMR) # Two Sample MR Functions and Interface to MR Base Database
library(MRInstruments) # [github::MRCIEU/MRInstruments] v0.3.2
library(MRPracticals) # [github::WSpiller/MRPracticals] v0.0.1

## Step 1: 获取暴露摘要估计
### GWAS 目录
data(gwas_catalog)
exposure_gwas <- subset(gwas_catalog, grepl("Blood", Phenotype_simple))
exposure_gwas <- subset(gwas_catalog, grepl("Neale", Author))

exposure_gwas <- subset(gwas_catalog, grepl("Locke", Author) &
  Phenotype_simple == "Body mass index")

head(exposure_gwas[, c(7:12, 18:21)])

exposure_gwas <- subset(exposure_gwas, grepl("EA", Phenotype_info))
exposure_gwas <- exposure_gwas[!grepl("women", exposure_gwas$Phenotype_info), ]
exposure_gwas <- exposure_gwas[!grepl("men", exposure_gwas$Phenotype_info), ]
head(exposure_gwas[, c(7:12, 18:21)])

exposure_gwas <- exposure_gwas[exposure_gwas$pval < 5 * 10^-8, ]

exposure_data <- format_data(exposure_gwas)

exposure_data <- clump_data(exposure_data, clump_r2 = 0.001)

exposure_data <- clump_data(exposure_data, clump_r2 = 0.001)

quick_extract <- extract_instruments(2)

## Step 2 获取结局摘要估计
ao <- available_outcomes()
head(ao)[, c(3, 4, 6, 8, 9, 20)]

outcome_gwas <- subset(ao, grepl("Systolic", trait))
head(outcome_gwas)[, c(3, 4, 6, 8, 9, 11, 16, 20)]

outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "UKB-a:360"
)

### 代理LD （LD proxys）的说明
H_data <- harmonise_data(
  exposure_dat = exposure_data,
  outcome_dat = outcome_data
)

### 重复条目
H_data <- power_prune(H_data)

## Step 3 进行MR分析
### 获取效应估计值
mr_results <- mr(H_data)
mr_results

mr(H_data, method_list = c("mr_egger_regression", "mr_ivw"))

mr_method_list()

head(mr_method_list())[, 1:2]

### 生成具有 95% 置信区间的OR值
generate_odds_ratios(mr_results) # 请注意，本例中的分析使用的是连续结局变量

### 敏感性分析
mr_pleiotropy_test(H_data)

mr_heterogeneity(H_data, method_list = c("mr_egger_regression", "mr_ivw"))

### 生成MR结果的图示
plot1 <- mr_scatter_plot(mr_results, H_data)
plot1

res_single <- mr_singlesnp(H_data)
plot2 <- mr_forest_plot(res_single)
plot2

res_loo <- mr_leaveoneout(H_data)
plot3 <- mr_leaveoneout_plot(res_loo)
plot3

plot4 <- mr_funnel_plot(res_single)
plot4



##### 药物靶基因孟德尔随机化 #####
# https://www.bilibili.com/video/BV1BN411p7ZD/?p=6&spm_id_from=pageDriver&vd_source=23f183f0c5968777e138f31842bde0a0
# 判断是否已经安装了“pacman”包I如果没有就安装它
if (!require("pacman")) install.packages("pacman", update = F, ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/")
# 加载“pacman”包，用于方便加载其他的R包
library("pacman") # Package Management Tool
# 使用“p_load”函数加载所需的R包
p_load(
  ggplot2, # CRAN v3.4.3
  readr # CRAN v2.1.4
)
p_load_gh("MRCIEU/MRInstruments")

library(TwoSampleMR) # Two Sample MR Functions and Interface to MR Base Database
library(Matrix)
library(MendelianRandomization)

# 暴露数据
# ebi-a-GCST90018926	2021	Type 2 diabetes	NA	490,089	24,167,560
expofile <- "ebi-a-GCST90018926"
# 靶基因位点
# TBC1D24: GRCh38.p14 (GCF_000001405.40)	16	NC_000016.10 (2475127..2505730)
# HNF4A: GRCh38.p14 (GCF_000001405.40)	20	NC_000020.11 (44355699..44434596)
chr_pos <- 20 # 染色体位置
pos_start <- 44355699 # 开始位置
pos_end <- 44434596 # 结束位置

# 结局数据
# ebi-a-GCST90018840	2021	Epilepsy	NA	458,310	24,186,492
# ukb-b-20124	2018	Heel bone mineral density (BMD) T-score, automated	MRC-IEU	265,753	9,851,867
outcfile <- "ebi-a-GCST90018840"

# 在线读取gwas数据
# 提取检测的SNP,存到biof_exposure_dat中
biof_exposure_dat <- extract_instruments(
  outcome = expofile,
  clump = FALSE
)

# 进行连锁不平衡刷除
biof_exposure_dat <- clump_data(biof_exposure_dat,
  clump_kb = 100, # 定义连锁不平衡窗口大小
  clump_r2 = 0.3, # 定义连锁不平衡的R平方阈值
  clump_p1 = 5e-08, # 保留p值最小的SNP
  clump_p2 = 5e-08 # 删除p值大于阈值的SNP
)

# 获取去除连锁不平衡后的行数和列数
dim(biof_exposure_dat)
view(biof_exposure_dat)

# 提取药物靶点周围的SNP
# subset函数从biof_exposure_dat筛选SNP
# chr.exposure==chr_pos:chromosome与药物靶点chromosome相等
# pos.exposure>pos_start-100000:SNP位置大于靶点start位置左侧100k!
# pos . exposure < pos_end + 100000: SNP位置小于靶点end位置右侧100kb
# 这样提取出药物靶点周围±100kb的SNP
Drug_Target_SNP <- subset(biof_exposure_dat,
                          chr.exposure == chr_pos &
                                   pos.exposure > pos_start - 100000 &
                                   pos.exposure < pos_end + 100000)
                          
# 将结果写入csv文件
# Drug_Target_SNP:提取的药物靶点周围SNP
# "Drug_Target_SNP.csv": 输出的csv文件名
write.csv(Drug_Target_SNP, file = "Drug_Target_SNP,csv")

# load outcome data找到和结局相关性的snp
# 从Outcome数据中提取与药物靶点SNP相关的表型数据
biof_Outcome_dat <- extract_outcome_data(
  snps = Drug_Target_SNP$SNP, # 药物靶点SNP
  outcomes = outcfile)  # 表型数据文件
  
# harmonize and merge 数据
# 确保SNP对暴露和结果的效应基于同一等位基因
harmonise_dat <- harmonise_data(
    exposure_dat = Drug_Target_SNP, # 药物靶点SNP数据
    outcome_dat = biof_Outcome_dat, # 相关表型数据
    action = 1)   # 调和方法

# 查看harmonise后的数据
View(harmonise_dat)

# 将结果写入文件
write.table(harmonise_dat,
            "expo_outc_harmonise.csv",
            row.names = F,
            sep = "\t",
            quote = F)

# 去除混杂因素
p_load(MendelianRandomization,purrr,readr)
# 设置变量grp_size,用于设置每个分组的SNP数量
grp_size <- 100

# 将data中的SNP分成多个大小为grp_size的分组，保存到变量grps中
grps <- split(harmonise_dat$SNP,ceiling(seq_along(harmonise_dat$SNP)/grp_size))

# 通过phenoscanner API对每个分组进行关联分析
# map_dfr将结果合并为一个数据框
# 对grps中的每个分组运行phenoscanner API进行关联分析
results <- map_dfr(grps, ~phenoscanner(snpquery=.×, #,x表示传入的分组
                                       catalogue="GwAS",#在GNAS目录中援索
                                       pvalue=1e-05,#p值阔值设置为1e-05
                                       proxies="None",#不使用代理SNP
                                       r2=0.8,#相关性阀值设置为0.8
                                       build=38)$resu1ts)#基因组版本为b38

# 将关联结果写入文件confounder07.csy
write_csv(results,"confounder.csv")
View(results)

# 设置变量remove_snps,用于保存要移除的混杂因素snp
remove_snps <- c("rs766466","rs2456973","rs9739070")

# 移除含有混杂SNP的行
harmonise_dat <- harmonise_dat[!harmonise_dat$SNP %in% remove_snps,]

# 将移除混杂因素后的结果写入clean_confounder07.csv
write_csv(harmonise_dat,"clean_confounder.csv")

res = mr(harmonise_dat)
res
write.csv(res, file="MRres.csv")

# mendelian randomization(MR)analysis
MR_result <- mr(harmonise_dat)
# View(MR_result)
result_OR=generate_odds_ratios(MR_result)
# 靶基因对疾病的关系，药物是抑制剂还是促进剂，
result_ORSor1=1/result_OR$or 

result_OR$or_lci951 = result_OR$or1-1.96*result_OR$se
result_ORSor_uci951 = result_OR$or1+1.96*result_OR$se
write.table(result_OR[,5:ncol(result_OR)],"MR_OR.xls",
            row.names = F, sep = "\t", quote = F)

mr_heterogeneity(harmonise_dat, method_list=c("mr_egger_regression","mr_ivw")) # 异质性检验
# outlier test
run_mr_presso(harmonise_dat, NbDistribution=1000) # 偏倚检验
pleio <-mr_pleiotropy_test(harmonise_dat) # 多效性检验--MR egger
# View(pleio)
write.csv(pleio,file="pleio.csv") # 多效性结果

single <- mr_leaveoneout(harmonise_dat)
mr_leaveoneout_plot(single) # 留一法检验敏感性
write.csv(single, file = "single.csv")

p1 <-mr_scatter_plot(res,harmonise_dat) # 散点图
ggsave(p1[[1]], file = "scatter.pdf", width = 8, height = 8)



