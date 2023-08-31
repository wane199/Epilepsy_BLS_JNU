# 生存分析机器学习模型解释的一揽子解决方案-survex包(https://mp.weixin.qq.com/s?__biz=Mzg5MTg2ODA1MA==&mid=2247486431&idx=1&sn=82b895fc703ffaea487ec0c5539dc365&chksm=cfc78d22f8b0043438174d989db2d70c38ebea1fd3736140b1358e212c4897a1fb24a5e882f2&mpshare=1&scene=24&srcid=0707bP93ICgcyb1QDqzaixOp&sharer_sharetime=1688699505382&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd)
# survex: Explainable Machine Learning in Survival Analysis(https://github.com/ModelOriented/survex)
# 模型和解释器创建
library(survex)
library(survival)
set.seed(123)
vet <- survival::veteran
cph <- coxph(Surv(time, status)~., data = vet, model = TRUE, x = TRUE)
cph_exp <- explain(cph)
rsf <- randomForestSRC::rfsrc(Surv(time, status)~., data = vet)
rsf_exp <- explain(rsf)

# 进行预测
# 可以直接利用解释器预测风险
predict(cph_exp, veteran[1:2,], output_type="risk")
predict(rsf_exp, veteran[1:2,], output_type="risk")
# 可以直接预测时间点生存概率
predict(cph_exp, veteran[1:2,], output_type="survival", times=seq(1, 600, 100))
predict(rsf_exp, veteran[1:2,], output_type="survival", times=seq(1, 600, 100))
# 可以预测时间点累计风险情况
predict(cph_exp, veteran[1:2,], output_type="chf", times=seq(1, 600, 100))
predict(rsf_exp, veteran[1:2,], output_type="chf", times=seq(1, 600, 100))

# 测量性能
mp_cph <- model_performance(cph_exp)
mp_rsf <- model_performance(rsf_exp)
# 整体上模型评价指标相对于时间变化折线图
plot(mp_cph, mp_rsf)
# 条形图展示模型总体评价
plot(mp_cph, mp_rsf, metrics_type="scalar")

# 全局解释
# 变量重要性
model_parts_rsf <- model_parts(rsf_exp)
model_parts_cph <- model_parts(cph_exp)
# 变量重要性排序图-时间依赖折线图
plot(model_parts_cph,model_parts_rsf)
model_parts_rsf_auc <- model_parts(rsf_exp, loss_function=loss_one_minus_cd_auc, type="difference")
model_parts_cph_auc <- model_parts(cph_exp, loss_function=loss_one_minus_cd_auc, type="difference")
# NOTE: this may take a long time, so a progress bar is available. To enable it, use:
# progressr::with_progress(model_parts(rsf_exp, loss_function=loss_one_minus_cd_auc, type="difference"))
# 变量重要性排序-损失函数为动态AUC- y 轴上的值仅表示每个变量排列后损失函数的变化
plot(model_parts_cph_auc,model_parts_rsf_auc)

# 偏依赖
model_profile_cph <- model_profile(cph_exp, categorical_variables=c("trt", "prior"))
# 特征变量对于最终生存结果的影响图-分布生存曲线图
plot(model_profile_cph, facet_ncol = 1)

model_profile_rsf <- model_profile(rsf_exp, categorical_variables=c("trt", "prior"))
plot(model_profile_rsf, facet_ncol = 1, numerical_plot_type = "contour")

# 局部解释
# 局部变量归因
# SurvSHAP(t) SurvSHAP(t)对生存模型的扩展。它们将预测分解为各个变量。
predict_parts_cph_32 <- predict_parts(cph_exp, veteran[32,])
predict_parts_rsf_32 <- predict_parts(rsf_exp, veteran[32,])
# 32观测的shap值分布情况
plot(predict_parts_cph_32, predict_parts_rsf_32)

predict_parts_cph_12 <- predict_parts(cph_exp, veteran[12,])
predict_parts_rsf_12 <- predict_parts(rsf_exp, veteran[12,])
# 12观测的shap值分布情况
plot(predict_parts_cph_12, predict_parts_rsf_12)

# SurvLIME
predict_parts_cph_12_lime <- predict_parts(cph_exp, veteran[12,], type="survlime")
predict_parts_rsf_12_lime <- predict_parts(rsf_exp, veteran[12,], type="survlime")
# coxph模型survlime解释结果
plot(predict_parts_cph_12_lime, type="local_importance")
# 随机生存森林的survlime解释结果
plot(predict_parts_rsf_12_lime, type="local_importance")

# Ceteris paribus
predict_profile_cph_32 <- predict_profile(cph_exp, veteran[32,], categorical_variables=c("trt", "prior"))
# Ceteris paribus结果和偏依赖很类似
plot(predict_profile_cph_32, facet_ncol=1)


##############################
# survex包在mlr3流程下的用法
if (!require("ooplah")) install.packages("ooplah")
if (!require("dictionar6")) install.packages("dictionar6")
if (!require("set6")) install.packages("set6",
                                       repos = "https://mlr-org.r-universe.dev")
if (!require("param6")) install.packages("param6",
                                         repos = "https://mlr-org.r-universe.dev")
if (!require("distr6")) install.packages("distr6",
                                         repos = "https://mlr-org.r-universe.dev")
if (!require("mlr3")) install.packages("mlr3")
if (!require("mlr3proba")) install.packages("mlr3proba",
                                            repos = "https://mlr-org.r-universe.dev")
if (!require("mlr3extralearners")) install.packages("mlr3extralearners",
                                                    repos = "https://mlr-org.r-universe.dev")
if (!require("mlr3pipelines")) install.packages("mlr3pipelines")
library(mlr3proba)
library(mlr3extralearners)
library(mlr3pipelines)
library(survex)
library(survival)

# 为 mlr3proba 学习器创建解释器
# survex 允许以以下方式为 mlr3proba::LearnerSurv() 对象创建解释器：
veteran_task <- as_task_surv(veteran,
                             time = "time",
                             event = "status",
                             type = "right")
ranger_learner <- lrn("surv.ranger") 
ranger_learner$train(veteran_task)
ranger_learner_explainer <- explain(ranger_learner, 
                                    data = veteran[, -c(3,4)],
                                    y = Surv(veteran$time, veteran$status),
                                    label = "Ranger model")

# Ceteris Paribus参数的变量重要性图-具体行
ranger_learner_explainer |> predict_profile(veteran[1,]) |> plot(numerical_plot_type = "contours",
                                                                 variables = c("karno", "celltype"),
                                                                 facet_ncol = 2,
                                                                 subtitle = NULL)


# 为具有复合分布预测的学习器创建解释器
gbm_composite_learner <- as_learner(ppl(
  "distrcompositor",
  learner = lrn("surv.gbm"),
  estimator = "kaplan",
  form = "ph"
))
gbm_composite_learner$train(veteran_task)
# important!
class(gbm_composite_learner) <- c(class(gbm_composite_learner), "LearnerSurv") 
gbm_composite_learner_explainer <- explain(gbm_composite_learner,
                                           data = veteran[, -c(3,4)],
                                           y = Surv(veteran$time, veteran$status),
                                           label = "Composite GBM model")
# 复合分布预测的学习器也可以创建解释器
gbm_composite_learner_explainer |> model_parts(type = "ratio") |> plot(subtitle = NULL)

# 使用 mlr3proba 度量
# mlr3proba 包提供了许多额外的性能度量。它们可以使用 loss_adapt_mlr3proba() 函数适用于 survex。这些适应的度量可以与 model_parts() 和 model_performance() 函数一起使用。
loss_schmid <- loss_adapt_mlr3proba(msr("surv.schmid"))
gbm_composite_learner_explainer |> model_parts(loss = loss_schmid) -> model_part
# 自定义损失函数计算获得的重要性排名
model_part |> plot(subtitle = NULL)

loss_rcll <- loss_adapt_mlr3proba(msr("surv.rcll"))
custom_metrics <- c("Integrated Schmid Score" = loss_schmid, 
                    "Right-Censored Log loss" = loss_rcll)
perf_ranger <- model_performance(ranger_learner_explainer, metrics=custom_metrics)
perf_comp_gbm <- model_performance(gbm_composite_learner_explainer, metrics=custom_metrics) 
# 可以一次性定义多个指标进行分析，可以分析综合指标
plot(perf_ranger, perf_comp_gbm, metrics_type="scalar", subtitle = NULL)


################################
# survex包自定义解释器
# 自动解释器创建
# 在最好的情况下，已经实现了为您想要的模型创建解释器。这意味着一切都可以从模型对象中提取（有时需要设置额外参数）。例如，这是 survival 包中比例风险模型的情况：
library(survex)
library(survival)
cph <- coxph(Surv(time, status) ~ ., data = veteran, model = TRUE, x = TRUE)
auto_cph_explainer <- explain(cph)

# 手动创建解释器
# 下一个基本情况是使用所需模型可以进行所有类型的预测。然后可以使用explain_survival()函数手动创建解释器。让我们看一个例子-如何手动设置coxph解释器的所有参数：
cph <- coxph(Surv(time, status) ~ ., data=veteran)
# data must not include the target columns
veteran_data <- veteran[, -c(3,4)]
veteran_y <- Surv(veteran$time, veteran$status)
# set prediction functions of the required format
risk_pred <- function(model, newdata) predict(model, newdata, type = "risk")
surv_pred <- function(model, newdata, times) pec::predictSurvProb(model, newdata, times)
chf_pred <- function(model, newdata, times) -log(surv_pred(model, newdata, times))
# 手动创建的解释器
manual_cph_explainer <- explain_survival(model = cph,
                                         data = veteran_data,
                                         y = veteran_y,
                                         predict_function = risk_pred,
                                         predict_survival_function = surv_pred,
                                         predict_cumulative_hazard_function = chf_pred,
                                         label="manual coxph")

surv_pred_rsf <- transform_to_stepfunction(predict,
                                           type="survival",
                                           prediction_element = "survival",
                                           times_element = "time.interest")

# would also work 
# chf_pred_rsf <- transform_to_stepfunction(predict,
#                                           type="chf",
#                                           prediction_element = "chf",
#                                           times_element = "time.interest")
chf_pred_rsf <- function(model, newdata, times) {
  survival_to_cumulative_hazard(surv_pred_rsf(model, newdata, times))
}
times <- unique(veteran$times)
risk_pred_rsf <- risk_from_chf(chf_pred_rsf, times)
times <- unique(veteran$times)
risk_pred_rsf <- risk_from_chf(chf_pred_rsf, times)

library(randomForestSRC)
rsf <- rfsrc(Surv(time, status) ~ ., data = veteran)
manual_rsf_explainer <- explain_survival(model = rsf,
                                         data = veteran_data,
                                         y = veteran_y,
                                         predict_function = risk_pred_rsf,
                                         predict_survival_function = surv_pred_rsf,
                                         predict_cumulative_hazard_function = chf_pred_rsf,
                                         label = "manual rsf")


##### Breakdown #####
# 探索可解释性机器学习：Breakdown带你了解心脏病随机森林的预测关键(https://mp.weixin.qq.com/s?__biz=MzU4MTk3NzAxNA==&mid=2247485896&idx=1&sn=29f1eb56900bc5427901d769ae191654&chksm=fdbe1f11cac99607d9b69f698e2f51af72f1f07100e041f88a754e7fd68fb522250cadf6a2d0&scene=178&cur_album_id=3027341943300669440#rd)
library("DALEX")
library("iBreakDown")
set.seed(123)
model_titanic_glm <- glm(survived ~ gender + age + fare,
                         data = titanic_imputed, family = "binomial")
explain_titanic_glm <- explain(model_titanic_glm,
                               data = titanic_imputed,
                               y = titanic_imputed$survived,
                               label = "glm")

bd_glm <- break_down(explain_titanic_glm, titanic_imputed[1, ])
bd_glm
plot(bd_glm, max_features = 3)


library("randomForest")
set.seed(123)
# example with interaction
# classification for HR data
model <- randomForest(status ~ . , data = HR)
new_observation <- HR_test[1,]

explainer_rf <- explain(model,
                        data = HR[1:1000,1:5])

bd_rf <- break_down(explainer_rf,
                    new_observation)
head(bd_rf)
plot(bd_rf)


# 描述使用的疾病预测数据集
#「加载数据集和特征选择」
library(survival)
data(heart)
data <- heart[,c("age","surgery","transplant","event")]
data$age <- abs(data$age)
head(data)

# 运用随机森林模型进行心脏病预测
# 加载依赖库」
library(DALEX)
library(randomForest)
#「模型拟合」
rf <- randomForest(event~., data=data)
#「构建模型解释器」
# 构建解释器
rf_exp <- explain(rf,
                  data = data[, -4],
                  y = data$event,
                  label = "randomForest")
rf_exp

#「breakdown解释」
bd_rf <- break_down(rf_exp,new_observation = data[1, ])
head(bd_rf)
describe(bd_rf)
plot(bd_rf)
































































