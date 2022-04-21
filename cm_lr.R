# library(glmnet)
# source("R/metrics/compute_metrics.R")

logistic_regression <- function(data.train, label.train, data.test, label.test,
                                classes, regularize = NA,
                                ...){
  # regularize <- "l2"
  # classes <- c("yes", "no")
  # model_name <- "logistic regression"
  # if (is.na(regularize)) {
  #   model_name <- paste("Simple", model_name)
  # } else if(regularize == 'l1') {
  #   model_name <- paste("L1 Regularized", model_name)
  # } else {
  #   model_name <- paste("L2 Regularized", model_name)
  # }
  
  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0) 
  
  try({
    label.train$Label <- ifelse(label.train$Label == classes[1], 0, 1)
    label.test$Label <- ifelse(label.test$Label == classes[1], 0, 1)
    
    
    if (!is.na(regularize)) {
      #alpha = 1 => l1 regularization (lasso)
      #alpha = 0 => l2 regularization (ridge)
      if(regularize == 'l1') {
        alpha <- 1
      }
      else {
        alpha <- 0
      }
      set.seed(1000)
      model <- glmnet::cv.glmnet(as.matrix(data.train), label.train$Label, 
                                 alpha = alpha, family = 'binomial', type.measure = 'mse')
      # plot(model)
      
      lambda_min <- model$lambda.min
      lambda_1se <- model$lambda.1se
      
      # best_acc <- -1
      best_cut_off <- 0.5
      cut_off <- 0.5
      # for(cut_off in seq(0.3, 0.7, 0.01)){
        # print(cut_off)
        pred_prob.train <- predict(model, newx = as.matrix(data.train), s = lambda_1se, type = 'response')
        pred.train <- ifelse(pred_prob.train > cut_off, 1, 0)
        acc <- mean(pred.train == label.train$Label)
        # print(acc)
        # if(acc > best_acc){
        #   best_acc <- acc
        #   best_cut_off <- cut_off
        # }
      # }
      
      
      pred_prob <- predict(model, newx = as.matrix(data.test), s = lambda_1se, type = 'response')
      pred <- ifelse(pred_prob > best_cut_off, 1, 0)
      # print(mean(pred == label.test$Label))
    }
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, true_label = label.test$Label, classes = c(0, 1))  
    # print(metrics)
    return(metrics)
    # result_df1 <- data.frame("TestDataClassName" = ifelse(label.test$Label == 0, classes[1], classes[2]),
    #                          "TestDataExpectedClassName" = ifelse(label.test$Label == 0, classes[1], classes[2]),
    #                          "TestDataClassId" = label.test$Label,
    #                          "Pred_prob" = pred_prob[,1],
    #                          "Prediction" = pred[,1])
    # 
    # result_df2 <- data.frame("TestDataClassName" = "PREREC",
    #                          "TestDataExpectedClassName" = "REC_TP", 
    #                          "TestDataClassId" = label.test2$Label,
    #                          "Pred_prob" = pred_prob2[,1],
    #                          "Prediction" = pred2[,1])
    # 
    # result_df <- rbind(result_df1, result_df2) 
    # colnames(result_df)[5] <- paste0("prediction_with_cutoff_", best_cut_off) 
    # 
    # # result_file_name <- "Data/prediction_result/transcriptomics.csv"
    # write.csv(result_df, result_file_name)  
  })
  
  
}


