#### SVMR&MVMR ####
library(TwoSampleMR)
library(openxlsx)
#### 单变量MR（暴露-结局）####
## 批量MR分析
# 自动一次提取三个变量的工具变量，只需要把id填入
Multiple_exp_dat <- extract_instruments(
  outcomes = c(
    "finn-b-RHEUMA_ENDPOINTS",
    "ieu-a-832",
    "finn-b-M13_SLE",
    "finn-b-M13_ANKYLOSPON"
  )
) # 风湿学终点,RA,SLE,AS
PSS <- extract_instruments(outcomes = "finn-b-L12_LOCALSCLERODERMA", p1 = 5e-05) # 硬皮病
vasculitis <- extract_instruments(outcomes = "finn-b-L12_VASCULITISNAS", p1 = 5e-05) # 血管炎
Multiple_exp_dat_all <- rbind(Multiple_exp_dat, PSS, vasculitis)

# 提取三个变量合并SNP的结局
# 当然这边也可以填入多个研究结局 outcomes=c("ieu-a-2","ieu-a-7")

outcome <- extract_outcome_data(
  snps = Multiple_exp_dat_all$SNP,
  outcomes = c(
    "ebi-a-GCST90000514",
    "ukb-b-10144",
    "ukb-b-19886"
  )
) # GERD,有炎症,无炎症
# harmonise # 【以这个函数,查看原理】
Multiple_har <- harmonise_data(Multiple_exp_dat_all, outcome)
# 进行mr分析
Multiple_res <- mr(Multiple_har)
Multiple_res <- generate_odds_ratios(Multiple_res)
# 由于结局是二分类，所以生成OR值
OR <- generate_odds_ratios(Multiple_res)
# 异质性检验
heterogeneity <- mr_heterogeneity(Multiple_har)
res_MRPRESSO <- run_mr_presso(Multiple_har)
# 多效性检验
pleiotropy <- mr_pleiotropy_test(Multiple_har)
pleiotropy
# 散点图
mr_scatter_plot(Multiple_res, Multiple_har)
# 逐个剔除检验，做留一图
leaveoneout <- mr_leaveoneout(Multiple_har)
mr_leaveoneout_plot(leaveoneout)
# 森林图
Multiple_res_single <- mr_singlesnp(Multiple_har)
mr_forest_plot(Multiple_res_single)
# 漏斗图
mr_funnel_plot(Multiple_res_single)
# 文件输出
write.xlsx(Multiple_res, file = "单变量MR结果(XY).xlsx")
write.xlsx(heterogeneity, file = "反向heterogeneity(XY).xlsx")
write.xlsx(res_MRPRESSO, file = "res_MRPRESSO(XY).xlsx")
write.xlsx(pleiotropy, file = "反向pleiotropy(XY).xlsx")
write.xlsx(Multiple_har, file = "反向Multiple_har(XY).xlsx")


#### 多变量MR####
# 混杂因素：体脂率(ebi-a-GCST003435)，BMI(ieu-a-2)，抑郁症(ieu-b-102)，久坐行为(ukb-a-5)
## 两两多变量
list_id <- c("ebi-a-GCST003435", "ieu-a-94", "ieu-b-102", "ukb-a-5")
exp_dat_MV <- NULL
for (i in 1:length(list_id)) {
  id_exposure <- c("ieu-a-832", list_id[i]) # 三个暴露：克罗恩病，溃疡性结肠炎，BMI
  id_outcome <- "ebi-a-GCST90000514" # 定义结局变量ID
  exposure_dat <- mv_extract_exposures(id_exposure) # 显著相关的SNP提出
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) # 结局数据
  mvdat_MV1 <- mv_harmonise_data(exposure_dat, outcome_dat) # mvdat
  res <- mv_multiple(mvdat_MV1)
  res_OR_MV1 <- generate_odds_ratios(res$result)
  exp_dat_MV <- rbind(exp_dat_MV, res_OR_MV1) # exp_dat_MV
}
write.xlsx(exp_dat_MV, file = "MVMR结果(两两混杂).xlsx")
## 全部多变量
id_exposure <- c("ieu-a-832", "ebi-a-GCST003435", "ieu-a-2", "ieu-b-102", "ukb-a-5") #
id_outcome <- "ebi-a-GCST90000514" # 定义结局变量ID
exposure_dat <- mv_extract_exposures(id_exposure) # 三个显著相关的SNP提出
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) # 结局数据
mvdat_MV2 <- mv_harmonise_data(exposure_dat, outcome_dat) # mvdat
res <- mv_multiple(mvdat_MV2)
res_OR_MV2 <- generate_odds_ratios(res$result) # exp_dat_MV
write.xlsx(res_OR_MV2, file = "MVMR结果(全部混杂).xlsx")

#### 潜在中介分析（暴露-中介）####
# 自动一次提取三个变量的工具变量，只需要把id填入
Multiple_exp_dat_mediator <- extract_instruments(
  outcomes = "ieu-a-832"
) # RA
# 提取三个变量合并SNP的结局
outcome_mediator <- extract_outcome_data(
  snps = Multiple_exp_dat_mediator$SNP,
  outcomes = c(
    "ieu-a-2",
    "ieu-b-102",
    "ukb-a-5",
    "finn-b-J10_ASTHMA"
  )
) # BMI,抑郁症,久坐,哮喘
# harmonise # 【以这个函数,查看原理】
Multiple_har_mediator <- harmonise_data(Multiple_exp_dat_mediator, outcome_mediator)
# 进行mr分析
Multiple_res_mediator <- mr(Multiple_har_mediator)
Multiple_res_mediator <- generate_odds_ratios(Multiple_res_mediator)
write.xlsx(Multiple_res_mediator, file = "潜在中介MR结果(XY).xlsx")


#### MVMR ####
# 设置工作目录
setwd("D:/【科研学习】/【胃外科】/【孟德尔随机化】")
# 读入数据包
library(TwoSampleMR)
######################## IEU在线提取数据MVMR###########################
id_exposure <- c("ieu-a-30", "ukb-a-104", "ukb-a-225", "ebi-a-GCST90012044") # 三个暴露：克罗恩病，溃疡性结肠炎，吸烟
id_outcome <- "bbj-a-119" # 定义结局变量ID
exposure_dat <- mv_extract_exposures(id_exposure) # 三个显著相关的SNP提出
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) # 结局数据
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
res <- mv_multiple(mvdat)
res_OR <- generate_odds_ratios(res$result)
res_OR
write.csv(res_OR, file = "MVMR结果(克-溃-烟-TNF).csv")

################## MVMR所用的SNP是根据什么规则来的######################
exp1 <- extract_instruments(outcomes = "ieu-a-755")
exp2 <- extract_instruments(outcomes = "ukb-b-16878")
exp3 <- extract_instruments(outcomes = "ieu-a-974")

sum <- dim(exp1)[1] + dim(exp2)[1] + dim(exp3)[1]

sum

exp <- rbind(exp1, exp2, exp3) # 按列合并
exp <- exp["SNP"] # 提出素有SNP
exp <- unique(exp) # 去非重复SNP
dim(exp) # 计算最后SNP

#################### LASSO去除高度共线MR###############
mv_lasso_feature_selection(mvdat) # LASSO

mv_lasso <- mv_subset(
  mvdat,
  features = mv_lasso_feature_selection(mvdat), # LASSO后进行MR
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
mv_lasso
mv_lasso_OR <- generate_odds_ratios(mv_lasso$result)
mv_lasso_OR
#################### 单个SNP的MVMR################
mv_residual(
  mvdat,
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
####################### PRESSO####################
# mr_presso可以完成多变量
if (!require("devtools")) {
  install.packages("devtools")
} else {}
devtools::install_github("rondolab/MR-PRESSO")

SummaryStats <- cbind(
  mvdat[["outcome_beta"]],
  mvdat[["exposure_beta"]][, 1],
  mvdat[["exposure_beta"]][, 2],
  mvdat[["exposure_beta"]][, 3],
  mvdat[["exposure_se"]][, 1],
  mvdat[["exposure_se"]][, 2],
  mvdat[["exposure_se"]][, 3],
  mvdat[["outcome_se"]]
)
SummaryStats <- data.frame(SummaryStats)

library(MRPRESSO)
mr_presso(
  BetaOutcome = "X1",
  BetaExposure = c("X2", "X3", "X4"),
  SdOutcome = "X8",
  SdExposure = c("X5", "X6", "X7"),
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = SummaryStats,
  NbDistribution = 1000,
  SignifThreshold = 0.05
)
################# MendelianRandomization包的几种方法############

if (!requireNamespace("MendelianRandomization")) {
  install.packages("MendelianRandomization")
}
library(MendelianRandomization)
# 读取实例数据
# 准备数据
MRMVInputObject <- mr_mvinput(
  bx = cbind(ldlc, hdlc, trig),
  bxse = cbind(ldlcse, hdlcse, trigse),
  by = chdlodds,
  byse = chdloddsse
)

MRMVInputObject

MRMVInputObject_1 <- mr_mvinput(
  bx = cbind(SummaryStats$X2, SummaryStats$X3, SummaryStats$X4),
  bxse = cbind(SummaryStats$X5, SummaryStats$X6, SummaryStats$X7),
  by = SummaryStats$X1,
  byse = SummaryStats$X8
)
MRMVInputObject_1
# IVW方法
MRMVObject <- mr_mvivw(MRMVInputObject,
  model = "default",
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)

MRMVObject

MRMVObject <- mr_mvivw(MRMVInputObject_1,
  model = "default",
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
MRMVObject

# egger方法
MRMVObject <- mr_mvegger(
  MRMVInputObject,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
MRMVObject

MRMVObject <- mr_mvegger(
  MRMVInputObject_1,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
MRMVObject

# LASSO
MRMVObject <- mr_mvlasso(
  MRMVInputObject,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject

MRMVObject <- mr_mvlasso(
  MRMVInputObject_1,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject

# median
MRMVObject <- mr_mvmedian(
  MRMVInputObject,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject

MRMVObject <- mr_mvmedian(
  MRMVInputObject_1,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject

#################### RMVMR###########################
remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
# rawdat_mvmr<-rawdat_rmvmr
head(rawdat_mvmr)
F.data <- format_mvmr(
  BXGs = rawdat_mvmr[, c(1, 2, 3)],
  BYG = rawdat_mvmr[, 7],
  seBXGs = rawdat_mvmr[, c(4, 5, 6)],
  seBYG = rawdat_mvmr[, 8],
  RSID = rawdat_mvmr[, 9]
)

head(F.data)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)












# 下边给大家放了TwoSampleMR本地的代码，个人感觉很难用，不灵活，有兴趣可以研究
################### 本地暴露读取#################
exp_l <- mv_extract_exposures_local(
  c("exp_a.txt", "exp_b.txt"),
  sep = " ",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  id_col = "id.exposure"
)

outcome_dat1 <- format_data(
  a,
  type = "outcome",
  snps = exp$"SNP",
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA_INSOMNIA",
  se_col = "SE_INSOMNIA",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_INSOMNIA",
)
