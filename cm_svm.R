#library(e1071)
# source("R/metrics/compute_metrics.R")

svm_model <- function(data.train, label.train, data.test = NA, label.test = NA,
                      classes, kernel = "sigmoid", 
                      ...){
  
  kernel = "sigmoid"
  kernel_name <- paste(toupper(substring(kernel, 1, 1)), substring(kernel, 2), sep = "")
  model_name <- paste(kernel_name, "Kernel SVM")
  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0, 0, 0) 
  
  try({
    set.seed(1000)
    model <- e1071::svm(data.train, factor(label.train$Label, levels = classes), probability = TRUE, kernel = kernel)
    
    if(!is.na(data.test) && !is.na(label.test)){
      pred <- predict(model, data.test, probability = TRUE)
      pred_prob <- data.frame(attr(pred, 'probabilities'))[classes[2]]   
      true_label = label.test$Label
    } else{
      pred <- predict(model, data.train, probability = TRUE)
      pred_prob <- data.frame(attr(pred, 'probabilities'))[classes[2]]  
      true_label = label.train$Label
    }  
    
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, true_label = true_label, classes = classes)    
    # print(metrics)
    return(metrics)
  })
  
  
}
