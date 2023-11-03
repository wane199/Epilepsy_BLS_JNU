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
library(doParallel) # 并行

## Step 1: 获取暴露摘要估计
### GWAS 目录
rm(list = ls())
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


##### 中介孟德尔随机化 #####
# https://mp.weixin.qq.com/s?__biz=Mzg5NjgwMzAzNQ==&mid=2247484020&idx=1&sn=71941a4d0d049a579da7324f09dad7c7&chksm=c07acbc0f70d42d60ad9251e4152feaadd7c908dce5b287e3ddaf2db01b12d4048ff36d100b2&mpshare=1&scene=1&srcid=0926XRKstADEkmt4HGh9iy7x&sharer_shareinfo=8d92a27f4cb892329551633b3390bf53&sharer_shareinfo_first=8d92a27f4cb892329551633b3390bf53#rd
# install packages
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
# install.packages("ggplot2")

# library
library(TwoSampleMR)
library(kableExtra)
library(ggplot2)
library(cowplot)

# set up your work directory
setwd("")

# two step 
###### step 1 X to Y TSMR(得到总效应，beta_all) ######
expo_rt<-extract_instruments(outcome="ukb-b-20124",clump = FALSE)
expo_rt<-clump_data(expo_rt,clump_kb = 10000,
                    clump_r2 = 0.001,
                    clump_p1 = 5e-08,
                    clump_p2 = 5e-08,
)

outc_rt <- extract_outcome_data(expo_rt$SNP, outcomes = "finn-b-M13_ATLOAXSUBLUX")


harm_rt <- harmonise_data(
  exposure_dat =  expo_rt, 
  outcome_dat = outc_rt,action=2)

mr_result<- mr(harm_rt)
View(mr_result)
result_or=generate_odds_ratios(mr_result)
write.table(result_or[,5:ncol(result_or)],"OR1.txt",row.names = F,sep = "\t",quote = F)

###### step 2 Y to X TSMR(确定可以做中介) ######
expo_rt2<-extract_instruments(outcome="finn-b-M13_ATLOAXSUBLUX",clump = FALSE)
expo_rt2<-clump_data(expo_rt,clump_kb = 10000,
                     clump_r2 = 0.001,
                     clump_p1 = 5e-08,
                     clump_p2 = 5e-08,
)

outc_rt2 <- extract_outcome_data(expo_rt2$SNP, outcomes = "ukb-b-20124")

harm_rt2 <- harmonise_data(
  exposure_dat =  expo_rt2, 
  outcome_dat = outc_rt2,action=2)

mr_result2<- mr(harm_rt2)
View(mr_result2)
result_or2=generate_odds_ratios(mr_result2)
write.table(result_or2[,5:ncol(result_or2)],"OR2.txt",row.names = F,sep = "\t",quote = F)

###### step 3 X to X' TSMR(得到beta1) ######
expo_rt3<-expo_rt
outc_rt3 <- extract_outcome_data(expo_rt3$SNP, outcomes = "ieu-b-4764")

harm_rt3 <- harmonise_data(
  exposure_dat =  expo_rt3, 
  outcome_dat = outc_rt3,action=2)
write.table(harm_rt3, "harmonise3.txt",row.names = F,sep = "\t",quote = F)

mr_result3<- mr(harm_rt3)
View(mr_result3)
result_or3=generate_odds_ratios(mr_result3)
write.table(result_or3[,5:ncol(result_or3)],"OR3.txt",row.names = F,sep = "\t",quote = F)

###### step 4 X' to Y TSMR(得到beta2) ######
expo_rt4<-extract_instruments(outcome="ieu-b-4764",clump = FALSE)
expo_rt4<-clump_data(expo_rt4,clump_kb = 10000,
                     clump_r2 = 0.001,
                     clump_p1 = 5e-08,
                     clump_p2 = 5e-08,
)

outc_rt4 <- extract_outcome_data(expo_rt4$SNP, outcomes = "finn-b-M13_ATLOAXSUBLUX")

harm_rt4 <- harmonise_data(
  exposure_dat =  expo_rt4, 
  outcome_dat = outc_rt4,action=2)
write.table(harm_rt4, "harmonise4.txt",row.names = F,sep = "\t",quote = F)

mr_result4<- mr(harm_rt4)
View(mr_result4)
result_or4=generate_odds_ratios(mr_result4)
write.table(result_or4[,5:ncol(result_or4)],"OR4.txt",row.names = F,sep = "\t",quote = F)

# 中介效应(乘积法) beta12=beta1*beta2
# https://www.bilibili.com/video/BV1Th4y1C7Je/?spm_id_from=333.788&vd_source=23f183f0c5968777e138f31842bde0a0
# 第一种中介效应 = 两样本（暴露 -> 中介）× 两样本（中介 -> 暴露）
product_method_PoE(EM_beta, EM_se, Mo_beta_total, Mo_se_total)

# 直接效应
beta_dir=beta_all-beta12


##### 药物靶基因孟德尔随机化 #####
# https://www.bilibili.com/video/BV1BN411p7ZD/?p=6&spm_id_from=pageDriver&vd_source=23f183f0c5968777e138f31842bde0a0
# 判断是否已经安装了“pacman”包I如果没有就安装它
rm(list = ls())
getwd()
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
# ebi-a-GCST90026412	2021	Severe autoimmune type 2 diabetes	NA	3,196	5,376,535
expofile <- "ebi-a-GCST90026412"

# 靶基因位点
# TBC1D24、HNF4A
# GABRG2：105.20220307	previous assembly	GRCh37.p13 (GCF_000001405.25)	5	NC_000005.9 (161494471..161582545)
chr_pos <- 5 # 染色体位置
pos_start <- 161494471 # 开始位置
pos_end <- 161582545 # 结束位置

# 结局数据
# ebi-a-GCST90018840	2021	Epilepsy	NA	458,310	24,186,492
# ukb-b-20124	2018	Heel bone mineral density (BMD) T-score, automated	MRC-IEU	265,753	9,851,867
# finn-b-FE	2021	Focal epilepsy	NA	—	16,380,452
outcfile <- "finn-b-FE"

# 在线读取gwas数据
# 提取检测的SNP,存到biof_exposure_dat中
biof_exposure_dat <- extract_instruments(
                    outcome = expofile,
                    clump = FALSE
)

# 读取本地暴露文件
{
biof_exposure_dat <- read_exposure_data(
            filename = "myexprosure.txt", # SNP
            sep = "\t",
            snp_Col = "SNP",
            beta_col = "Beta",
            se_col = "SE",
            effect_allele_col = "Effect_allele",
            other_allele_col = "Other_allele",
            pval_col = "P",
            eaf_col = "EAF",
            chr_col = "CHR",
            pos_col = "BP"
)

# data filter
biof_exposure_dat <- biof_exposure_dat[biof_exposure_dat$pval.exposure < 1e-8,]
biof_exposure_dat <- clump_data(biof_exposure_dat,clump_kb = 100,
                                 c1ump_r2 = 0.3)
}

# 进行连锁不平衡剔除
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
# chr.exposure == chr_pos:chromosome与药物靶点chromosome相等
# pos.exposure > pos_start - 100000:SNP位置大于靶点start位置左侧100kb
# pos.exposure < pos_end + 100000: SNP位置小于靶点end位置右侧100kb
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

# 读取本地结局数据
{
  biof_Outcome_dat <- read_outcome_data(
    snps = Drug_Target_SNP$SNP,
    filename = "CAD_OUTCOME.txt",
    sep = "\t",
    snp_col = "markername",
    beta_col ="beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "noneffect_allele", # other_allele
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value")
  
  # View(biof_Outcome_dat)
  write.csv(biof_Outcome_dat, file="biof_Outcome_dat.csv")
}

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

# 将移除混杂因素后的结果写入clean_confounder.csv
write_csv(harmonise_dat,"clean_confounder.csv")

res = mr(harmonise_dat)
res
write.csv(res, file="MRres.csv")

# mendelian randomization(MR)analysis
MR_result <- mr(harmonise_dat)
# View(MR_result)
result_OR=generate_odds_ratios(MR_result)
# 靶基因对疾病的关系:正/负相关，药物是抑制剂还是促进剂，
result_OR$or1=1/result_OR$or 

result_OR$or_lci951 = result_OR$or1-1.96*result_OR$se
result_OR$or_uci951 = result_OR$or1+1.96*result_OR$se
write.table(result_OR[,5:ncol(result_OR)],"MR_OR.xls",
            row.names = F, sep = "\t", quote = F)

mr_heterogeneity(harmonise_dat, method_list=c("mr_egger_regression","mr_ivw")) # 异质性检验
# outlier test
run_mr_presso(harmonise_dat, NbDistribution=1000) # 偏倚检验
pleio <- mr_pleiotropy_test(harmonise_dat) # 多效性检验--MR egger
# View(pleio)
write.csv(pleio,file="pleio.csv") # 多效性结果

single <- mr_leaveoneout(harmonise_dat)
mr_leaveoneout_plot(single) # 留一法检验敏感性
write.csv(single, file = "single.csv")

p1 <- mr_scatter_plot(res,harmonise_dat) # 散点图
p1
ggsave(p1[[1]], file = "scatter.pdf", width = 8, height = 8)


# 单SNP分析
res_single = mr_singlesnp(harmonise_dat)
# 绘制森林图，森林图可以表示每个SNP与结果的关联效应。
p2 <- mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]],file="forest.pdf",width=8,height=8)
# 绘制漏斗图
p3 <- mr_funnel_plot(singlesnp_results = res_single)
p3[[1]]
ggsave(p3[[1]],file="funnel_plot.pdf",width=8,height=8)

# 留一法敏感性分析
p4 <- mr_leaveoneout_plot(leaveoneout_results=mr_leaveoneout(harmonise_dat))
p4[[1]]
ggsave(p4[[1]],file="leaveoneout.pdf",width=8,height=8)

# 合并MR结果，绘制森林图
# 使用"p_load"函数加载所需的R包
p_load(forestploter,grid,ggplot2)

# 定义变量mrresultsfile,用于存储MR结果文件的文件名
mrresultsfile = "MRresults.txt"
# 使用read.table()函数读取MR结果文件
# 文件名为mrresultsfile变量的值"MRresults.txt"
# header = T表示文件第一行作为变量名
# sep = "\t”表示列之间使用制表符分隔
biofsci = read.table(mrresultsfile,header = T,sep = "\t")

# 处理P值的显示格式
# 对biofsci数据集的P.value字段进行处理
# 使用ifelse()进行条件判断
# 如果P.va1ue小于0.001显示"<0.001"
# 否则使用sprintf()格式化显示P值到小数点后4位
biofsci$P.value = ifelse(biofsci$P.value < 0.001, "<0.001",
                         sprintf("%.4f",biofsci$P.value))

# 处理NA值，使用ifelse()对多个字段进行处理
# 如果值为NA,则替换为""
biofsci$Target = ifelse(is.na(biofsci$Target),"",biofsci$Target)
biofsci$Method = ifelse(is.na(biofsci$Method),"",biofsci$Method)
biofsci$NSNP = ifelse(is.na(biofsci$NSNP),"",biofsci$NSNP)
biofsci$P.value = ifelse(is.na(biofsci$P.value),"",biofsci$p.value)

# 添加空格用于格式调整
# 使用rep()函数重复生成20个空格字符
# 使用paste()用空格拼接成一个字符串
# 赋值给新生成的变量biofsci''
biofscis$' ' <- paste(rep(" ", 20), collapse = " ")

# 处理OR值的显示格式
# 使用ifelse()对0R字段进行处理
# 如果0R是NA,显示""，否则显示格式化的0R值和95%可信区间
biofsci$'OR (95% CI)' <- ifelse(is.na(biofsci$OR),"",
                             sprintf("%.4f (%.4f to %.4f)",
                                     biofsci$OR,biofsci$OR_lci95,
                                     biofsci$OR_uci95))

# 设置森林图的主题，定义图形元素的是示格式
# 使用forest_theme()函数
# 定义基础字体大小，参考线颜色，脚注格式等参数
tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

# 绘制森林图
# 使用forest()函数，传递数据及参数
# biofsci[,c(1:5,9:10)]:选择需要的列数据
# est/lower/upper:OR值和置信区间
# arrowlab:曲线两个方向的标签文本
# sizes:点的大小
# ci_column:置信区间所在列
# ref_line:添加参考线
# xlim:x轴范围
# footnote:脚注
# theme:使用设置的主题对象
p=forest(biofsci[,c(1:5,9:10)],
         est = biofsci$OR,
         lower = biofsci$OR_lci95,
         upper = biofsci$OR_uci95,
         # arrow_lab = c("Placebo Better","Treatment Better").
         sizes = 0.4,
         ci_column = 6,
         ref_line = 1,
         xlim = c(0.05,3),
         footnote = "",
         theme = tm)
P

ggsave(p,file="forest.pdf",width=8,height=28)


# ..\\Downloads\\药物靶向MR\\药靶课程更新\\15分药靶（贝叶斯共定位分析）
library(TwoSampleMR)
library(data.table)
library(MRInstruments)
library(MRPRESSO)


ldsc_gwas<-extract_instruments(outcome="ieu-a-300",clump = FALSE)

ldsc_gwas<-clump_data(ldsc_gwas,clump_kb = 100,
                      clump_r2 = 0.3,
                      clump_p1 = 5e-08,
                      clump_p2 = 5e-08,
                      pop = "EUR")

HMGCR_gwas<-subset(ldsc_gwas,chr.exposure==5 & pos.exposure>74632993-100000 & pos.exposure<74657941+100000)
PCSK9_gwas<-subset(ldsc_gwas,chr.exposure==1 & pos.exposure>55505221-100000 & pos.exposure<55530525+100000)
NPC1L1_gwas<-subset(ldsc_gwas,chr.exposure==7 & pos.exposure>44552134-100000 & pos.exposure<44580929+100000)

HMGCR_gwas<-subset(HMGCR_gwas,eaf.exposure>0.01)
PCSK9_gwas<-subset(PCSK9_gwas,eaf.exposure>0.01)
NPC1L1_gwas<-subset(NPC1L1_gwas,eaf.exposure>0.01)

#ieu-a-7	Coronary heart disease
HMGCR_CHD_gwas <- extract_outcome_data(snps = HMGCR_gwas$SNP,
                                       outcomes = 'ieu-a-7')

dat <- harmonise_data(HMGCR_gwas,HMGCR_CHD_gwas)
res3 <- mr(dat)
res3
