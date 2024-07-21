##### Model Explanation #####
# Extract the AUC values from the aggregated data
AUC <- dat_aggregate_long.2[dat_aggregate_long.2$measure == "classif.auc", ]
# Order the AUC values in descending order and store the result
AUC <- AUC[order(AUC$value, decreasing = TRUE), ]
# Print the AUC values
AUC

# Train the decision tree model using the training data
rr_classif.rpart$learner$train(task, split$train)
# Predict on the test data
rr_classif.rpart$learner$predict(task, split$test)

# Set a seed for reproducibility
set.seed(123)
# Create a classification task from the expr2_imp dataset with the target variable "Group"
task <- as_task_classif(expr2_imp, target = "Group")
# Split the data into training and testing sets with a 70% training ratio
split <- partition(task, ratio = 0.7)
# Extract the feature data for the test set
dt_x = task$data(rows = split$test, cols = task$feature_names)
# Extract the target data for the test set
dt_y = task$data(rows = split$test, cols = task$target_names)

# Select the best neural network model based on the lowest cross-entropy
nn_final <- rr_classif.nnet$score(measures)[learner_id == 'classif.nnet.tuned', ][classif.ce == min(classif.ce)]$learner[[1]]
Train the selected neural network model
nn_final$train(task)

# The following lines are commented out and not executed
# Explain the neural network model using DALEXtra package (this requires the DALEXtra package and the model to be explained)
nn_exp = DALEXtra::explain_mlr3(nn_final,
                                 data = dt_x,
                                 y = as.numeric(dt_y$Group == 1),
                                 label = "nnet_model",
                                 colorize = FALSE)

# Load the iml library for model interpretation
library(iml)
# Create a new predictor object for the xgboost model
predictor = Predictor$new(nn_final, 
                          data = dt_x, 
                          y = dt_y)
# Create a feature importance object using the predictor and cross-entropy loss
importance <- FeatureImp$new(predictor, loss = "ce")
# Save the explanation as a PDF file
pdf(paste0("./explain/", name, ".pdf"))
plot(importance) + theme_bw()
dev.off()

# Save the feature importance results as a CSV file
write.csv(as.data.frame(importance$results), paste0("./explain/", name, ".csv"))
# Save the importance object and the glmnet model for future use
save(importance, rr_classif.glmnet, file = paste0("./explain/", name, ".RData"))

# Define a function to explain models and save the results
exp_mod <- function(mod, name) {
  # Similar to the above, select the best model based on the lowest cross-entropy
  nn_final <- mod$score()[learner_id == mod$score()$learner_id[1], ][classif.ce == min(classif.ce)]$learner[[1]]
  # Train the selected model
  nn_final$train(task)
  
  # Create a new predictor object for the model
  predictor = Predictor$new(nn_final, 
                            data = dt_x, 
                            y = dt_y)
  # Create a feature importance object
  importance <- FeatureImp$new(predictor, loss = "ce")
  # Save the explanation as a PDF file with specific dimensions
  p <- plot(importance) + theme_bw()
  ggsave(paste0("./explain/", name, ".pdf"), p, width = 8, height = 8)
  
  # Save the feature importance results as a CSV file
  write.csv(as.data.frame(importance$results), paste0("./explain/", name, ".csv"))
  # Save the importance object and the model for future use
  save(importance, mod, file = paste0("./explain/", name, ".RData"))
}

# Call the model explanation function with the glmnet model and a specified name for the output files
exp_mod(rr_classif.glmnet, name = "exp")