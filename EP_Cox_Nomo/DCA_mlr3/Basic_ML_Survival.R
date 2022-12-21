# https://www.nature.com/articles/s41467-022-28421-6#Sec24
# Lasso, Ridge, Enet, StepCox, SurvivalSVM, CoxBoost, SuperPC, plsRcox, RSF, GBM(https://mp.weixin.qq.com/s?__biz=MzAwMjY4MDE2Mg==&mid=2247619851&idx=1&sn=2df8ba1c261dab389f630eeea91d1040&chksm=9ac5e106adb26810faecd8195976fc3910886079c926c3569865255936b2d4dc27b1bc8ed3d2&mpshare=1&scene=1&srcid=1123lWmWuQAxJ26m48HGkQQz&sharer_sharetime=1669216079860&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd)
library(glmnet)
library(mlr3)
library(mlr3verse) # lrn()函数生存分析机器学习算法汇总

# 1.Lasso
# 代码解析
#### 1.glmnet包
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1) # 核心函数(执行lasso回归)
#### 2.mlr3包
lrn <- lrn("classif.cv_glmnet", alpha = 1) # 核心函数(执行lasso回归)

# 2.Ridge
#代码解析
####1.glmnet包
cvfit <- cv.glmnet(x, y, family = "binomial",alpha = 0)#核心函数(执行Ridge回归)
####2.mlr3包
lrn <- lrn("classif.cv_glmnet",alpha = 0)#核心函数(执行Ridge回归)

# 3.Enet（弹性网络）
#代码解析
####1.glmnet包
cvfit <- cv.glmnet(x, y, family = "binomial",alpha = 0.5)#核心函数(执行Enet回归)
####2.mlr3包
lrn <- lrn("classif.cv_glmnet",alpha = 0.5)#核心函数(执行Enet回归)

# 4.StepCox
#代码解析
####1.survival包
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)}##核心函数
####2.mlr3包
#这个目前在mlr3中并不可以得到实现，所以我们只有采用基本函数来实现这个应用

# 5.survivalSVM
#代码解析
survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)##核心函数

# 6.CoxBoost
#代码解析
#github上有现成的R包可以完成这项操作
library(devtools)
install_github("binderh/CoxBoost")
RIF <- resample.CoxBoost(time=obs.time,status=obs.status,x=x,rep=100, maxstepno=200,multicore=FALSE,
                         mix.list=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99), 
                         stratum=group,stratnotinfocus=0,penalty=sum(obs.status)*(1/0.02-1),
                         criterion="hscore",unpen.index=NULL) 
#当然作者针对这里使用的是下面这套代码
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)


# 7.SuperPC(jedazard/superpc: Supervised Principal Components (github.com))
#代码解析
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)


# 8.plsRcox:高维数据中如果带有生存预后信息，那么可以通过plsRcox这个算法实现在高维数据中拟合cox模型
#代码实现
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))

# 9.RSF
#代码解析
library(randomForestSRC)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)

# 10.GBM
#代码解析
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)


###########################################
# [Survival Analysis in mlr3proba](https://mran.microsoft.com/snapshot/2020-01-15/web/packages/mlr3proba/vignettes/survival.html)
# Survival Tasks
# Unlike TaskClassif and TaskRegr which have a single ‘target’ argument, TaskSurv mimics the survival::Surv object and has three-four target arguments (dependent on censoring type)
rm(list = ls())
library(mlr3proba); library(mlr3); library(survival)
# type = "right" is default
TaskSurv$new(id = "right_censored", backend = survival::rats,
             time = "time", event = "status", type = "right")

task = TaskSurv$new(id = "interval_censored", backend = survival::bladder2[,-c(1, 7)],
                    time = "start", time2 = "stop", type = "interval2")
task
task$truth()[1:10]

# Train and Predict
# create task and learner
# veteran = mlr3misc::load_dataset("veteran", package = "survival")
veteran = survival::veteran
task_veteran = TaskSurv$new(id = "veteran", backend = veteran, time = "time", event = "status")
learner = lrn("surv.coxph")

# train/test split 
train_set = sample(task_veteran$nrow, 0.7 * task_veteran$nrow)
test_set = setdiff(seq_len(task_veteran$nrow), train_set)

# fit Cox PH and inspect model
learner$train(task_veteran, row_ids = train_set)
learner$model

prediction = learner$predict(task_veteran, row_ids = test_set)
prediction

# Evaluate - crank, lp, and distr
# In the previous example, Cox model predicts `lp` so `crank` is identical
all(prediction$lp == prediction$crank)
prediction$lp[1:10]

# These are evaluated with measures of discrimination and calibration.
# As all PredictionSurv objects will return crank, Harrell's C is the default measure.
prediction$score()

# distr is evaluated with probabilistic scoring rules.
measure = lapply(c("surv.graf", "surv.rmse"), msr)
prediction$score(measure)

# Often measures can be integrated over mutliple time-points, or return
# predictions for single time-points
measure = msr("surv.graf", times = 60)
prediction$score(measure)





 










# All Together Now
# Putting all of this together we can perform a (overly simplified) benchmark experiment to find the best learner for making predictions on a simulated dataset.
library(mlr3pipelines); library(mlr3); library(mlr3tuning); library(paradox)
set.seed(123)

task = tgen("simsurv")$generate(50)
composed_lrn_gbm = distrcompositor(lrn("surv.gbm", bag.fraction = 1, n.trees = 50L),
                                   "kaplan", "ph")

lrns = lapply(paste0("surv.", c("kaplan", "coxph", "parametric")), lrn)
lrns[[3]]$param_set$values = list(dist = "weibull", type = "ph")

design = benchmark_grid(tasks = task, learners = c(lrns, list(composed_lrn_gbm)),
                        resamplings = rsmp("cv", folds = 2))

bm = benchmark(design)
bm$aggregate(lapply(c("surv.harrellC","surv.graf","surv.grafSE"), msr))[,c(4, 7:9)]



library(mlr3); library(mlr3proba)
library(mlr3learners); library(mlr3pipelines)
library(mlr3verse)
kaplan = lrn("surv.kaplan")
cox = lrn("surv.coxph")
xgb = ppl("distrcompositor", learner = lrn("surv.xgboost"), estimator = "kaplan", form = "ph")
learners = list(cox, kaplan, xgb)
task = TaskSurv$new(id = "rats", backend = survival::rats[,1:4], time = "time", event = "status")
resample = rsmp("cv", folds = 5)
design = benchmark_grid(task, learners, resample)
bm = benchmark(design)
bm$aggregate(msr("surv.intlogloss"))

autoplot(res)








