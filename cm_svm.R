#library(e1071)
# source("R/metrics/compute_metrics.R")

svm_model <- function(data.train, label.train, data.test, label.test,
                      classes, kernel = "sigmoid", 
                      ...){
  
  kernel = "sigmoid"
  kernel_name <- paste(toupper(substring(kernel, 1, 1)), substring(kernel, 2), sep = "")
  model_name <- paste(kernel_name, "Kernel SVM")
  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0) 
  
  try({
    model <- e1071::svm(data.train, factor(label.train$Label, levels = classes), probability = TRUE, kernel = kernel)
    
    pred <- predict(model, data.test, probability = TRUE)
    pred_prob <- data.frame(attr(pred, 'probabilities'))[classes[2]]
    
    # result_df <- data.frame("TestDataClassName" = label.test$Label,
    #                          "TestDataExpectedClassName" = label.test$Label,
    #                          "Pred_prob" = pred_prob[,1],
    #                          "Prediction" = pred)
    # 
    # # result_file_name <- "Data/prediction_result/transcriptomics.csv"
    # write.csv(result_df, result_file_name)  
    
    
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, true_label = label.test$Label, classes = classes)    
    # print(metrics)
    return(metrics)
  })
  
  
}
