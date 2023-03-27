# (TwoSample MR Practical)[https://zhuanlan.zhihu.com/p/605966100]
install.packages(c("devtools", "knitr", "rmarkdown"))
library(devtools)
install_github(c("MRCIEU/TwoSampleMR", "MRCIEU/MRInstruments"))
install_github("WSpiller/MRPracticals", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
install_github("WSpiller/MRPracticals", build_opts = c("--no-resave-data", "--no-manual"))
library(TwoSampleMR)
library(MRInstruments)
library(MRPracticals)

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
mr_results<-mr(H_data)
mr_results

mr(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

mr_method_list()

head(mr_method_list())[,1:2]

### 生成具有 95% 置信区间的OR值
generate_odds_ratios(mr_results) # 请注意，本例中的分析使用的是连续结局变量

### 敏感性分析
mr_pleiotropy_test(H_data)

mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

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



