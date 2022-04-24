# library(randomForest)
# source("R/metrics/compute_metrics.R")

rf_model <- function(data.train, label.train, data.test = NA, label.test = NA, 
                     classes, random_seed = 1000, ...){
  model_name <- "Random Forest"
  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0, 0, 0)   
  
  try({
    set.seed(random_seed)
    if(sum(data.train - colMeans(data.train)) != 0){
      #ensure that atleast one column is not a constant vector
      #all columns constant causes the below line to run forever
      model <- randomForest::randomForest(x = data.train, y = factor(label.train$Label, levels = classes))
      
      if(!is.na(data.test) && !is.na(label.test)){
        pred_prob <- predict(model, data.test, type="prob")
        pred_prob <- data.frame(pred_prob)[classes[2]]
        pred <- ifelse(pred_prob > 0.5, classes[2], classes[1])
        true_label = label.test$Label
      } else{
        pred_prob <- predict(model, data.train, type="prob")
        pred_prob <- data.frame(pred_prob)[classes[2]]
        pred <- ifelse(pred_prob > 0.5, classes[2], classes[1])
        true_label = label.train$Label
      }  
      
      metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, 
                                 true_label = true_label, classes = classes)
      # print(metrics)    
    } else{
      print("data to RF : all fields constant!")
      print(dim(data.train))
    }
    
  })
  
  # return (list(model_name, metrics, feature_imp))
  return(metrics)
}
