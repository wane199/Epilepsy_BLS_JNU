#########################################
## https://mp.weixin.qq.com/s?__biz=Mzg3ODg5MzU5NA==&mid=2247485627&idx=1&sn=ea4c37193d22ebef6de4f9913ff4e369&chksm=cf0d87cef87a0ed85627aa18245e6ff13f045fa865d82c932baffe2c6a148399965b2f91b1d3&mpshare=1&scene=1&srcid=0508mw4LEN8NGw5BjMmrAaNA&sharer_sharetime=1683513791848&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd
## 自动机器学习（AutoML）(https://docs.h2o.ai/h2o/latest-stable/h2o-docs/explain.html?highlight=h2o%20pd_multi_plot)
# 基于h2o机器学习框架的自动机器学习方法，并对自动机器学习模型进行解释。
# 1、加载包和数据集
rm(list = ls())
library(h2o)
h2o.init() # 初始化

# 导入数据
dt <- h2o.importFile("C:\\Users\\wane1\\Documents\\file\\sci\\cph\\XML\\TLE234group_2019.csv")
# 分类变量转变为因子
dt <- na.omit(dt)
dt <- dt[c(7:24)] # 获取数据
vfactor <- c("oneyr", "side", "Sex")
dt[vfactor] <- lapply(dt[vfactor], factor)
# 批量数值转因子
for (i in names(dt)[c(-2, -3, -6:-8)]) {
  dt[, i] <- as.factor(dt[, i])
}
str(dt)
table(dt$oneyr) # 查看阳性结局

# 2、构建自动机器学习模型
# 2.1模型构建
am <- h2o.automl(
  y = "oneyr",
  training_frame = dt,
  max_models = 10
)

# 2.2  模型列表
b <- h2o.get_leaderboard(am) # 默认为排名前6的模型
b

# 2.3 选取最优模型
best <- h2o.get_best_model(am)

# 3、模型表现
perf <- h2o.performance(best)
perf

# 3.1 绘制ROC曲线
# AUC
h2o.auc(perf)

# 绘制ROC曲线
plot(perf, type = "roc")

# 3.2 变量重要性
h2o.varimp_plot(best)

# 3.3 基于排列的变量重要性
h2o.permutation_importance_plot(best, dt)

# 4、模型解释
# 4.1 SHAP summary plot
h2o.shap_summary_plot(best, dt)

# 4.2 单个样本的SHAP
h2o.shap_explain_row_plot(best, dt, row_index = 1)

# 4.3 部分依赖图（PDP）
h2o.pd_plot(best, dt, "AI_radscore")
# 多模型的PDP图
h2o.pd_multi_plot(am@leaderboard, dt, "AI_radscore")

# 4.4 个体条件期望图（ICE）
h2o.ice_plot(best,
  dt,
  show_pdp = TRUE,
  "side"
)

# 4.5 学习曲线
h2o.learning_curve_plot(best)


###########################################
# 基于Tidymodels的自动机器学习(https://mp.weixin.qq.com/s?__biz=Mzg3ODg5MzU5NA==&mid=2247485642&idx=1&sn=dcf03c1fddcfc3f6fbf4a24c3d6a21df&chksm=cf0d87bff87a0ea99924fdbecf737e9dbc5beb1294f8a32abb36f524f6ed2908b9c7ff5bafc4&cur_album_id=2792146949775884289&scene=190#rd)
# 1、加载R包及数据处理
library(tidymodels)
library(agua)
library(ggplot2)
library(h2o)
tidymodels_prefer()
theme_set(theme_classic())
# 启动h2o ，需要先安装h2o包
h2o_start()
# 数据加载及划分
dt <- read.csv("C:\\Users\\wane1\\Documents\\file\\sci\\cph\\XML\\TLE234group_2019.csv")
dt <- na.omit(dt)
dt <- dt[c(7:24)] # 获取数据
vfactor <- c("oneyr", "side", "Sex")
dt[vfactor] <- lapply(dt[vfactor], factor)
# 批量数值转因子
for (i in names(dt)[c(-2, -3, -6:-8)]) {
  dt[, i] <- as.factor(dt[, i])
}
str(dt)
table(dt$oneyr) # 查看阳性结局

set.seed(123)
split <- initial_split(dt, strata = oneyr)
train <- training(split)
test <- testing(split)

# 2、设置自动机器学习
# 2.1设置自动机器学习模型
auto_spec <-
  auto_ml() %>%
  set_engine("h2o", max_runtime_secs = 60, seed = 1) %>%
  set_mode("classification")

# 2.2 设置recipes
rec <- recipe(oneyr ~ ., data = train)

# 2.3 设置工作流
auto_wflow <-
  workflow() %>%
  add_model(auto_spec) %>%
  add_recipe(rec)

# 2.4 模型拟合
auto_fit <- fit(auto_wflow, data = train)

# 3、模型信息
# 3.1 提取模型拟合结果
extract_fit_parsnip(auto_fit) # 提取模型信息

# 3.2  预测
predict(auto_fit, new_data = test)

# 3.3 经交叉验证的模型根据按照不同指标进行排序
agua::rank_results(auto_fit) %>%
  filter(.metric == "auc") %>%
  arrange(rank)

# 3.4 每个模型的统计信息
collect_metrics(auto_fit, summarize = FALSE)

# 3.5 输出每个模型的预测值
tidy(auto_fit) %>%
  mutate(
    .predictions = map(.model, predict, new_data = head(test))
  )

# 3.6 基模型在stack集成模型中的权重
auto_fit %>%
  extract_fit_parsnip() %>%
  member_weights() %>%
  unnest(importance) %>%
  filter(type == "scaled_importance") %>%
  ggplot() +
  geom_boxplot(aes(value, algorithm)) +
  scale_x_sqrt() +
  labs(y = NULL, x = "scaled importance", title = "Member importance in stacked ensembles")

# 3.7 绘制性能比较图
autoplot(auto_fit, type = "rank", metric = c("auc", "accuracy")) +
  theme(legend.position = "none")



###################################
# 基于mlr3的全自动机器学习-mlr3automl代码示例(https://mp.weixin.qq.com/s?__biz=Mzg5MTg2ODA1MA==&mid=2247486673&idx=2&sn=ddfd4c42bc2104f4bfef160deced2e4f&chksm=cfc78a2cf8b0033af94511084f47a6edc27fc35731bb385017dbf644105b5beeae61693bd745&cur_album_id=2789575730761400322&scene=189#wechat_redirect)
# 创建模型
devtools::install_github('https://github.com/a-hanf/mlr3automl', dependencies = TRUE)
library(mlr3)
library(mlr3automl)

iris_task = tsk("iris")
iris_model = AutoML(iris_task)


# 训练，预测和重采样
# 在上面的代码中，我们的模型尚未训练。模型训练与mlr3中的训练和预测非常相似。我们使用train()方法，它需要一个训练索引向量作为附加参数。
train_indices = sample(1:iris_task$nrow, 2/3*iris_task$nrow)
iris_model$train(row_ids = train_indices)

predict_indices = setdiff(1:iris_task$nrow, train_indices)
predictions = iris_model$predict(row_ids = predict_indices)

resampling_result = iris_model$resample()

# 可解释模型
dalex_explainer = iris_model$explain(iml_package = "DALEX")
iml_explainer = iris_model$explain(iml_package = "iml")

# compute and plot feature permutation importance using DALEX
dalex_importance = DALEX::model_parts(dalex_explainer)
# 变量重要性
plot(dalex_importance)

# partial dependency plot using iml package
iml_pdp = iml::FeatureEffect$new(iml_explainer, feature="Sepal.Width", method="pdp")
# 变量预测价值图
plot(iml_pdp)

# 创建一个具有自定义学习器和固定时间预算的回归模型。每次超参数评估都会在10秒后停止。经过300秒的培训后，调优将停止，并返回迄今为止获得的最佳结果。
automl_model = AutoML(
  task=tsk("mtcars"),
  learner_list=c("regr.ranger", "regr.lm"),
  learner_timeout=10,
  runtime=300)

library(mlr3pipelines)
imbalanced_preproc = po("imputemean") %>>%
  po("smote") %>>%
  po("classweights", minor_weight=2)

automl_model = AutoML(task=tsk("pima"),
                      preprocessing = imbalanced_preproc) 


library(paradox)
new_params = ParamSet$new(list(
  ParamInt$new("classif.kknn.k",
               lower = 1, upper = 5, default = 3, tags = "kknn")))

my_trafo = function(x, param_set) {
  if ("classif.kknn.k" %in% names(x)) {
    x[["classif.kknn.k"]] = 2^x[["classif.kknn.k"]]
  }
  return(x)
}

automl_model = AutoML(
  task=tsk("iris"), 
  learner_list="classif.kknn",
  additional_params=new_params,
  custom_trafo=my_trafo)



























