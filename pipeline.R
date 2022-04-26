library(tidyverse)
library(ranger)
library(mRMRe)
library(randomForest)
library(psdR)
library(umap)
source("cm_lr.R")
source("cm_svm.R")
source("cm_rf.R")
source("compute_metrics.R")

meta_data <- read.csv("data/meta_data.csv")

summary(factor(meta_data$subcohort))




pipeline <- function(data.train, label.train, 
                     data.test = NA, label.test = NA,
                     classes){
  #filter and transform
  filtered_features <- colSums(data.train) != 0
  data.train <- data.train[, filtered_features]
  # min(data.train[data.train != 0])
  #4e-05
  
  data.train <- data.train + 10^-5
  data.train <- log(data.train)
  
  #normalize
  normparam <- caret::preProcess(data.train) 
  data.train <- predict(normparam, data.train)
  # colSums(data.train)  
  
  #ranger
  random_seed <- 1000
  set.seed(random_seed)
  ranger_model <- ranger::ranger(x = data.train, y = factor(label.train$Label), 
                                 importance = "impurity_corrected")
  features <- which(ranger_model$variable.importance >= 0)
  data.train <- data.train[, features, drop = FALSE]
  print(assertthat::are_equal(rownames(data.train), label.train$sample))
  
  if(!is.na(data.test) && !is.na(label.test)){
    data.test <- data.test[, filtered_features]    
    data.test <- data.test + 10^-5
    data.test <- log(data.test)
    data.test <- predict(normparam, data.test) #normalizing test data using params from train data 
    data.test <- data.test[, features, drop = FALSE]
    print(assertthat::are_equal(rownames(data.test), label.test$sample))  
  }
  
  #classification model
  res_l2 <- logistic_regression(data.train, label.train, 
                      data.test, label.test,
                      classes, regularize = "l2")
  res_l2 <- res_l2[1:2]
  res_l1 <- logistic_regression(data.train, label.train, 
                                data.test, label.test,
                                classes, regularize = "l1")
  res_l1 <- res_l1[1:2]
  
  
  res_svmsig <- svm_model(data.train, label.train, data.test, label.test, 
                          classes, kernel = "sigmoid")
  res_svmsig <- res_svmsig[1:2]
  
  res_svmrad <- svm_model(data.train, label.train, data.test, label.test, 
                          classes, kernel = "radial")
  res_svmrad <- res_svmrad[1:2]
  
  res_rf <- rf_model(data.train, label.train, data.test, label.test,
                     classes)  
  res_rf <- res_rf[1:2]
  
  result_df <- data.frame(rbind(c("l2", res_l2),
                                c("l1", res_l1),
                                c("svmsig", res_svmsig),
                                c("svmrad", res_svmrad),
                                c("rf", res_rf)))  
  
  colnames(result_df) <- c("cm", "acc", "auc")
  
  return(result_df)
}

# classes = c("yes", "no")
# pipeline(data.train, label.train, data.test, label.test, classes)
# pipeline(data.train, label.train, NA, NA, classes)


#create smaller subsets of data with 0.98 size successively (total 30)
# and in each of those execute pipeline with whole data
#        and also run 50 different 80:20 train:test splits
run_pipeline_multiple_subsets <- function(data, output_labels, p = 0.98, ss_iter = 30){
  all_result_df <- data.frame(matrix(nrow = 0, ncol = 6, 
                                     dimnames = list(c(), 
                                                     c("size_iter", "samples",
                                                       "cm", 
                                                       "acc", "auc",
                                                       "iter"))
  )
  )
  for(ss_i in c(1:ss_iter)){
    # ss_i <- 2
    if(ss_i >= 2){
      train_index <- caret::createDataPartition(output_labels$Label, p = p, 
                                                list = FALSE, 
                                                times = 1)
      data <- data[train_index, ]
      output_labels <- output_labels[train_index, ]    
    }
    print(ss_i)
    print(dim(data))
    
    set.seed(1000)
    train_index <- caret::createMultiFolds(y = output_labels$Label, k = 5, times = 10)
    
    #result without train-test split
    result_df <- pipeline(data, output_labels, NA, NA, classes = c("yes", "no"))
    result_df <- result_df %>%
      mutate("iter" = NA)
    #result_df : cm, acc, auc, iter
    
    #results with train-test split
    for(i in c(1:50)){
      # i <- 1
      data.train <- data[train_index[[i]], ]
      label.train <- output_labels[train_index[[i]], ]
      
      data.test <- data[-train_index[[i]], ]
      label.test <- output_labels[-train_index[[i]], ]  
      
      # print(assertthat::are_equal(rownames(data.train), label.train$sample))
      # print(assertthat::are_equal(rownames(data.test), label.test$sample))
      
      result_df <- rbind(result_df,
                         pipeline(data.train, label.train, 
                                  data.test, label.test, 
                                  classes = c("yes", "no")) %>%
                           mutate("iter" = i))
    }
    
    result_df <- result_df %>%
      mutate(acc = as.double(acc), auc = as.double(auc))
    
    result_df <- result_df %>%
      mutate("size_iter" = ss_i, .before = "cm") %>%
      mutate("samples" = dim(data)[1], .after = "size_iter")
    
    all_result_df <- rbind(all_result_df,
                           result_df)
  } 
  return (all_result_df)
}


###################################################################

# with 164 patients

data <- read.csv("data/all_level_formatted_data.csv")
colnames(data)[1] <- "sample"
sum(is.na(data))
output_labels <- meta_data %>%
  # filter(!is.na(RECIST) & RECIST != "SD") %>%
  # filter(subcohort %in% c("PRIMM-UK")) %>%
  select(c(sample, ICIresponder)) 

output_labels <- output_labels %>%
  dplyr::rename(c("Label" = "ICIresponder"))

combined_data <- output_labels %>%
  inner_join(data)
missing1 <- output_labels %>%
  anti_join(data)
missing2 <- data %>%
  anti_join(output_labels)

data <- combined_data %>%
  select(-c(Label)) %>%
  column_to_rownames("sample")
assertthat::are_equal(rownames(data), output_labels$sample)

all_result_df <- run_pipeline_multiple_subsets(data, output_labels)
write.csv(all_result_df, "all_result_df_full_data_using_all_levels.csv", row.names = FALSE)




#with PRIMM-UK subcohort

# data <- read.csv("data/formatted_data.csv")

data <- read.csv("data/all_level_formatted_data.csv")
colnames(data)[1] <- "sample"
sum(is.na(data))
output_labels <- meta_data %>%
  # filter(!is.na(RECIST) & RECIST != "SD") %>%
  filter(subcohort %in% c("PRIMM-UK")) %>%
  select(c(sample, ICIresponder)) 

output_labels <- output_labels %>%
  dplyr::rename(c("Label" = "ICIresponder"))

combined_data <- output_labels %>%
  inner_join(data)
missing1 <- output_labels %>%
  anti_join(data)
missing2 <- data %>%
  anti_join(output_labels)

data <- combined_data %>%
  select(-c(Label)) %>%
  column_to_rownames("sample")
assertthat::are_equal(rownames(data), output_labels$sample)

# print(dim(data))
# for(ss_i in c(1:15)){
#   # ss_i <- 2
#   if(ss_i >= 2){
#     train_index <- caret::createDataPartition(output_labels$Label, p = 0.93, 
#                                               list = FALSE, 
#                                               times = 1)
#     data <- data[train_index, ]
#     output_labels <- output_labels[train_index, ]    
#   }
#   print(ss_i)
#   print(dim(data))
# }



all_result_df <- run_pipeline_multiple_subsets(data, output_labels, p = 0.93, ss_iter = 15)
write.csv(all_result_df, "all_result_df_PRIMMUK_data_using_all_levels.csv", row.names = FALSE)


full_data_result <- read.csv("all_result_df_full_data_using_all_levels.csv")
pu_result <- read.csv("all_result_df_PRIMMUK_data_using_all_levels.csv") 

# # ss = 1
# get_result_summary <- function(result_df, model, ss = 1){
#   # model <- "l2"
#   model_res <- result_df %>%
#     filter(cm == model)
#   if("size_iter" %in% colnames(model_res)){
#     model_res <- model_res %>%
#       filter(size_iter == ss)
#   }
#   print("AUC")
#   print(summary(as.double(model_res$auc)))
#   print("Acc")
#   print(summary(as.double(model_res$acc)))
# }
# 
# get_result_summary(all_result_df, "l1")
# get_result_summary(all_result_df, "l2")
# get_result_summary(all_result_df, "svmsig")
# get_result_summary(all_result_df, "svmrad")
# get_result_summary(all_result_df, "rf")
# 
# get_result_summary(all_result_df, "l1", 2)
# get_result_summary(all_result_df, "l2", 2)
# get_result_summary(all_result_df, "svmsig", 2)
# get_result_summary(all_result_df, "svmrad", 2)
# get_result_summary(all_result_df, "rf", 2)


result_df <- full_data_result

result <- "primm"
on_train_data = FALSE
metric_name = "auc"
model = "rf"

plot_variation <- function(result = "full", 
                           on_train_data = FALSE,
                           metric_name = "auc",
                           model = "rf"){
  
  if(result == "full"){
    result_df <- read.csv("all_result_df_full_data_using_all_levels.csv")
  } else{
    result_df <- read.csv("all_result_df_PRIMMUK_data_using_all_levels.csv") 
  }
  
  data_to_plot <- result_df %>%
    filter(cm == model)
  
  
  if(on_train_data){
    #create plot with fitted curve for training results
    data_to_plot <- data_to_plot %>%
      filter(is.na(iter)) 
    data_to_plot <- data_to_plot %>%
      select(samples, metric_name) %>%
      mutate("metric" = metric_name) %>%
      rename(c("val" = metric_name))    
  } else{
    #create plot with fitted curve for test results
    data_to_plot <- data_to_plot %>%
      filter(!is.na(iter)) %>%
      group_by(samples) %>%
      summarise(med_acc = median(acc), med_auc = median(auc)) %>%
      pivot_longer(!samples, names_to = "metric", values_to = "val")
    
    metric_name <- paste("med", metric_name, sep = "_")
    data_to_plot <- data_to_plot %>%
      filter(metric == metric_name)    
  }

  plot_title <- paste("Performance of", model, "on", 
                      if(on_train_data){
                        "train data"
                      } else{
                        "test data"
                      },
                      "with median", gsub("med_", "", metric_name))
  
  ggplot(data_to_plot) +
    geom_line(aes(x = samples, y = val, color = metric)) +
    geom_point(aes(x = samples, y = val, color = metric)) +
    xlab("Number of Samples") +
    ylab("Metric Value") +
    ggtitle(plot_title)
  
  file_name <- paste(result, model, on_train_data, metric_name, ".png", sep = "_")
  
  ggsave(file_name)
}

setwd("~/UNSW/LiverCancerGroup/sample_size_estimation/")

for(r in c("full", "PRIMMUK")){
  for(otc in c(FALSE, TRUE)){
    for(model in c("l1", "l2", "svmsig", "svmrad", "rf")){
      for(mn in c("acc", "auc")){
        plot_variation(result = r, on_train_data = otc, metric_name = mn, model = model)
      }
    }
  }
}
# plot_variation(result = "full", 
#                on_train_data = FALSE,
#                metric_name = "auc",
#                model = "rf")
# plot_variation(result = "full", 
#                on_train_data = TRUE,
#                metric_name = "auc",
#                model = "rf")


result = "PRIMMUK" 
on_train_data = FALSE
metric_name = "auc"
model = "rf"
plot_fit <- function(result = "PRIMMUK", 
                     on_train_data = FALSE,
                     metric_name = "auc",
                     model = "rf"){
  
  if(result == "full"){
    result_df <- read.csv("all_result_df_full_data_using_all_levels.csv")
  } else{
    result_df <- read.csv("all_result_df_PRIMMUK_data_using_all_levels.csv") 
  }
  
  data_to_plot <- result_df %>%
    filter(cm == model)
  
  
  if(on_train_data){
    #create plot with fitted curve for training results
    data_to_plot <- data_to_plot %>%
      filter(is.na(iter)) 
    data_to_plot <- data_to_plot %>%
      select(samples, metric_name) %>%
      mutate("metric" = metric_name) %>%
      rename(c("val" = metric_name))    
  } else{
    #create plot with fitted curve for test results
    data_to_plot <- data_to_plot %>%
      filter(!is.na(iter)) %>%
      group_by(samples) %>%
      summarise(med_acc = median(acc), med_auc = median(auc)) %>%
      pivot_longer(!samples, names_to = "metric", values_to = "val")
    
    metric_name <- paste("med", metric_name, sep = "_")
    data_to_plot <- data_to_plot %>%
      filter(metric == metric_name)    
  }
  
  
  #filter out values that come down
  data_to_plot <- data_to_plot[c(1:4, 9, 13),]
  
  
  
  file_name <- paste("fitted_curve",
    result, model, on_train_data, metric_name, ".png", sep = "_")

  x <- data_to_plot$samples
  y <- data_to_plot$val

  plot(x, y)

  m <- nls(y ~ a*log(x + b))
  m <- nls(y ~ a*log(x + b), start = list(a = 5,
                                          b = -10))

  m <- nls(y ~ (a * x) / (b + x))
  m <- nls(y ~ (a * x) / (b + x), start = list(a = 0.6,
                                               b = -5))
  # 
  m <- nls(y ~ a / exp(1/(x-3)))
  # 
  # m <- nls(y ~ a * x / (b + x), start = list(a = a_start,
  #                                            b = b_start))

  #get some estimation of goodness of fit
  
  m <- nls(y ~ a - exp(-x + b), start = list(a = 0.5, b = 30))
  # 0.9556264
  # Nonlinear regression model
  # model: y ~ a - exp(-x + b)
  # data: parent.frame()
  # a       b 
  # 0.7244 28.2211 
  # residual sum-of-squares: 0.002042
  
  m <- nls(y ~ a - exp(-x + b), start = list(a = 0.7, b = 30))
  
  m <- nls(y ~ - a*exp(-0.05*x + b), start = list(a = 0.9, b = 50))
  # m <- nls(y ~ a - exp(-x))
  # m <- nls(y ~ - exp(-x + b))
  
  cor(y,predict(m))
  #plot
  # plot(x,y)
  # lines(x,predict(m),lty=2,col="blue",lwd=3)
  
  ggplot() +
    geom_line(aes(x = c(30:100), y = predict(m, newdata = list(x = c(30:100)))), 
              linetype = "dashed", color = "blue") +
    geom_point(aes(x = x, y = y), color = "red") +
    xlab("Sample Size") +
    ylab("Median AUC")
  ggsave(file_name)
  
  
  #current best m
  # Nonlinear regression model
  # model: y ~ a * log(x + b)
  # data: parent.frame()
  # a      b 
  # 0.1819 0.2212 
  # residual sum-of-squares: 0.02912
  
  plot(x, log(10*y))
  
}

