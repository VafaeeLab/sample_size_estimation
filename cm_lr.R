logistic_regression <- function(data.train, label.train, data.test = NA, label.test = NA,
                                classes, regularize = NA,
                                ...){
  #if test data not provided then return metrics on train data
  
  #setting default value for metrics, to handle case where unable to train / execute classification model
  metrics <- c(0, 0, 0, 0) 
  
  try({
    label.train$Label <- ifelse(label.train$Label == classes[1], 0, 1)
    
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
      
      if(!is.na(data.test) && !is.na(label.test)){
        label.test$Label <- ifelse(label.test$Label == classes[1], 0, 1)    
        pred_prob <- predict(model, newx = as.matrix(data.test), s = lambda_1se, type = 'response')
        pred <- ifelse(pred_prob > best_cut_off, 1, 0)
        true_label <- label.test$Label
        # print(mean(pred == label.test$Label))
      } else{
        pred_prob <- pred_prob.train
        pred <- pred.train
        true_label <- label.train$Label
      }
      
    }
    metrics <- compute_metrics(pred = pred, pred_prob = pred_prob, 
                               true_label = true_label, classes = c(0, 1))  
    # print(metrics)
    return(metrics)
  })
  
  
}


