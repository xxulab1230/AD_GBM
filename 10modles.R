# Load the mlr3verse library that provides a framework for machine learning in R
library(future)
library(mlr3verse)
# Load the tidyverse collection of packages for data manipulation and visualization
library(tidyverse)

# # Set the seed for random number generation to ensure reproducibility
{set.seed(3333) 
  # Create a classification task from the expr2_imp dataset with the target variable "Group"
  task<-as_task_classif(expr2_imp,
                        target = "Group")
  
  split<-partition(task)
  # Define the rows in the test set as the "holdout" group
  task$set_row_roles(split$test,"holdout")
  
  # View the feature types in the task
  task$feature_types
  
  # Convert the task data to a tibble for easier manipulation
  data<- task$data() |> as_tibble()
  
  # Convert character columns to numeric, assuming they represent integer features
  int_cols <- task$feature_types[type == "character"]$id
  data[, int_cols] <- lapply(data[, int_cols], as.numeric)
  # Add the converted data back to the task
  task$cbind(data[, int_cols])
  task
  # View the updated feature types in the task
  task$feature_types
  
  
  # Define the evaluation metric for the models (Area Under the ROC Curve)
  measure <- msr("classif.auc")
  inner_resample <- rsmp("cv", folds = 10)
  outer_resample <- rsmp("cv", folds = 10)
  # Define the hyperparameter tuning strategy (grid search)
  tuner <- tnr("grid_search",resolution=5)
  # Define the final set of performance metrics for model evaluation
  final_score <- msrs(c("classif.auc",
                        "classif.ce",
                        "classif.acc",
                        "classif.precision",
                        "classif.recall",
                        "classif.sensitivity", 
                        "classif.specificity"))
  
  
  # Check available cores for parallel processing
  availableCores()
  # Set up parallel processing with the multisession plan and 20 workers
  future::plan("multisession", workers=20)
  
  # 1.Define and tune a decision tree model, using a probability prediction type
  
  rr_classif.rpart<- 
    tune_nested(learner = lrn("classif.rpart",
                              predict_type = "prob",
                              id ="Classification_Tree"),
                search_space = ps(cp = p_dbl(lower = 0.001, upper = 0.1)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
  # Extract and format the results for the decision tree model
  res_classif.rpart<- 
    rr_classif.rpart$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  # Similar blocks of code are used to define and tune additional models 
  #(KNN, SVM, XGBoost, Random Forest, Glmnet, LDA, Logistic Regression, Naive Bayes, Neural Network)
  ##2.KNN
  rr_classif.kknn <- 
    tune_nested(learner = lrn("classif.kknn",
                              predict_type = "prob",
                              id="KNN"),
                search_space = ps(k = p_int(lower = 1, upper = 10)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
 
  
  res_classif.kknn<- 
    rr_classif.kknn$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  ##3.SVM
  rr_classif.svm<- 
    tune_nested(learner = lrn("classif.svm",
                              type = "C-classification",
                              kernel = "radial", 
                              predict_type = "prob",
                              id = "SVM"),
                search_space = ps(cost = p_dbl(log(0.1), log(10),
                                               trafo = function(x) exp(x)),
                                  gamma = p_dbl(log(0.1), log(10),
                                                trafo = function(x) exp(x))),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
  
  res_classif.svm<- 
    rr_classif.svm$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  ##4.xgboost
  
  rr_classif.xgboost<- 
    tune_nested(learner = lrn("classif.xgboost", 
                              predict_type = "prob",
                              id ="xgboost"),
                search_space = ps(
                  eta = p_dbl(lower = 0.1, upper = 1),
                  max_depth = p_int(lower = 1, upper = 10),
                  nrounds = p_int(lower = 1, upper = 16)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
 
  res_classif.xgboost<- 
    rr_classif.xgboost$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  ##5. RF
  rr_classif.ranger<- 
    tune_nested(learner = lrn("classif.ranger", 
                              predict_type = "prob",
                              id ="Random_Forest"),
                search_space = ps(
                  mtry = p_int(lower = 1, upper = task$ncol-1),
                  min.node.size = p_int(lower = 1, upper = 10),
                  num.trees = p_int(lower = 1, upper = 500)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
  
  
  res_classif.ranger<- 
    rr_classif.ranger$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  #6.glmnet
  rr_classif.glmnet<- 
    tune_nested(learner = lrn("classif.glmnet", 
                              predict_type = "prob",
                              id ="Glmnet"),
                search_space = ps(
                  lambda = p_dbl(lower = 0.001, upper = 1, logscale = TRUE),
                  alpha = p_dbl(lower = 0, upper = 1)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
  
 
  res_classif.glmnet<- 
    rr_classif.glmnet$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
 
  #6.LDA
  task_prep <- 
    po("filter", # 去除高度相关的列
       filter = mlr3filters::flt("find_correlation"), 
       filter.cutoff=0.3) %>>%
    po("scale", scale = F) %>>% # 中心化
    po("removeconstants") %>>% # 去掉零方差变量 
    po("encode")
  
  
  learner = task_prep%>>%lrn("classif.lda", predict_type = "prob")
  learner <- GraphLearner$new(learner)
  
  
  rr_classif.lda = resample(task, learner, outer_resample)
  
  
  res_classif.lda<- 
    rr_classif.lda$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  # 8.logreg
  learner = lrn("classif.log_reg", predict_type = "prob")
  
  rr_logreg = resample(task, learner, outer_resample)
  
  res_logreg<- 
    rr_logreg$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  # 9. naive_bayes
  rr_classif.naive_bayes<- 
    tune_nested(learner = lrn("classif.naive_bayes", predict_type = "prob"),
                search_space = ps(
                  laplace = p_dbl(lower = 0.1, upper = 1),
                  threshold = p_int(lower = 1, upper = 10)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
  
  res_classif.naive_bayes<- 
    rr_classif.naive_bayes$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  # 10.nnet
  rr_classif.nnet<- 
    tune_nested(learner = lrn("classif.nnet",predict_type = "prob"),
                search_space = ps(
                  size = p_int(lower = 1, upper = 10),
                  decay = p_dbl(lower = 0.1, upper = 0.9),
                  maxit = p_int(lower = 100, upper = 1000)),
                task = task,
                tuner = tuner,
                inner_resampling = inner_resample,
                outer_resampling = outer_resample,
                measure = measure)
  
  
  res_classif.nnet<- 
    rr_classif.nnet$score(final_score) |> 
    as.data.table() |> 
    as_tibble() |> 
    select(task_id,learner_id,iteration,contains("classif"))
  
  # final
  # Aggregate the results from all models into a single dataset for comparison  
  dat_aggregate.2 <- 
    rr_classif.rpart$aggregate(final_score) |> 
    bind_rows(rr_classif.svm$aggregate(final_score)) |> 
    bind_rows(rr_classif.xgboost$aggregate(final_score)) |> 
    bind_rows(rr_classif.kknn$aggregate(final_score)) |> 
    bind_rows(rr_classif.ranger$aggregate(final_score)) |> 
    bind_rows(rr_classif.glmnet$aggregate(final_score)) |> 
    bind_rows(rr_classif.lda$aggregate(final_score)) |> 
    bind_rows(rr_logreg$aggregate(final_score)) |> 
    bind_rows(rr_classif.nnet$aggregate(final_score)) |> 
    bind_rows(rr_classif.naive_bayes$aggregate(final_score)) |> 
    mutate(model = c("Classification Tree","SVM","Xgboost",
                     "KNN","Random Forest","Glmnet","LDA",
                     "logistic","NNET","Naive Bayes")) |> 
    column_to_rownames(var = "model") |> 
    as_tibble()
  
  # Rename the performance metric columns and round the values for visualization
  colnames(dat_aggregate.2) <-
    colnames(dat_aggregate.2) |> 
    str_replace_all("classif.","")
  
  
  # Create a long-format dataset for visualization purposes
  dat_aggregate_long.2 <- 
    rr_classif.rpart$aggregate(final_score) |> 
    bind_rows(rr_classif.svm$aggregate(final_score)) |> 
    bind_rows(rr_classif.xgboost$aggregate(final_score)) |> 
    bind_rows(rr_classif.kknn$aggregate(final_score)) |> 
    bind_rows(rr_classif.ranger$aggregate(final_score)) |> 
    bind_rows(rr_classif.glmnet$aggregate(final_score)) |> 
    bind_rows(rr_classif.lda$aggregate(final_score)) |> 
    bind_rows(rr_logreg$aggregate(final_score)) |> 
    bind_rows(rr_classif.nnet$aggregate(final_score)) |> 
    bind_rows(rr_classif.naive_bayes$aggregate(final_score)) |> 
    mutate(model = c("Classification Tree","SVM","Xgboost",
                     "KNN","Random Forest","Glmnet","LDA",
                     "logistic","NNET","Naive Bayes"))  |> 
    pivot_longer(!model,names_to = "measure", values_to = "value") 
  
  colnames(dat_aggregate_long.2) <-
    colnames(dat_aggregate_long.2) |> 
    str_replace_all("classif.","")
  
  # Create a heatmap visualization comparing the model performances using ggplot2
  dat_aggregate_long.2$value <- round(dat_aggregate_long.2$value,2)
  pdf("./heatmap.pdf",width = 9,height = 6)
  ggplot(dat_aggregate_long.2, aes(measure, model)) + 
    geom_tile(aes(fill = value), colour = "grey", size = 1)+
    scale_fill_gradient2(low = "#CEF0F1",mid = "#FFDFD3",high = "#BDBCEC",midpoint = 0.6) + 
    geom_text(aes(label=value),col ="black",size = 5) +
    theme_minimal() + 
    theme(axis.title.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"), 
          axis.text.y = element_text(size = 14, face = "bold")) + 
    labs(fill =paste0("Models Compare")) + 
    scale_x_discrete(position = "bottom")
  dev.off()
}

write.csv(as.data.frame(dat_aggregate_long.2),file = "dataggre.csv")


