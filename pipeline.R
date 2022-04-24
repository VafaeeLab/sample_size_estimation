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


#create smaller subsets of data with 0.9 size successively
# and in each of those execute pipeline with whole data
#        and also run 50 different 80:20 train:test splits
run_pipeline_multiple_subsets <- function(data, output_labels){
  all_result_df <- data.frame(matrix(nrow = 0, ncol = 6, 
                                     dimnames = list(c(), 
                                                     c("size_iter", "samples",
                                                       "cm", 
                                                       "acc", "auc",
                                                       "iter"))
  )
  )
  for(ss_i in c(1:10)){
    # ss_i <- 2
    if(ss_i >= 2){
      train_index <- caret::createDataPartition(output_labels$Label, p = .9, 
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

all_result_df <- run_pipeline_multiple_subsets(data, output_labels)
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
# 
# get_result_summary(all_result_df, "l1", 5)
# get_result_summary(all_result_df, "l2", 5)
# get_result_summary(all_result_df, "svmsig", 5)
# get_result_summary(all_result_df, "svmrad", 5)
# get_result_summary(all_result_df, "rf", 5)
# 
# write.csv(all_result_df, "all_result_df_PRIMMUK_data_using_all_level.csv",
#           row.names = FALSE)
# 
# 
# write.csv(result_df, "result_df_full_data_using_species.csv", row.names = FALSE)
# 
# write.csv(result_df, "result_df_PRIMMUK_data_using_species.csv", row.names = FALSE)
# 
# write.csv(result_df, "result_df_PRIMMUK_data_using_all_level.csv", row.names = FALSE)
# 
# 
# result_df <- read.csv("result_df_full_data_using_species.csv")
# get_result_summary(result_df, "l1")
# get_result_summary(result_df, "l2")
# get_result_summary(result_df, "svmsig")
# get_result_summary(result_df, "svmrad")
# get_result_summary(result_df, "rf")
# 
# 
# result_df <- read.csv("result_df_PRIMMUK_data_using_species.csv")
# get_result_summary(result_df, "l1")
# get_result_summary(result_df, "l2")
# get_result_summary(result_df, "svmsig")
# get_result_summary(result_df, "svmrad")
# get_result_summary(result_df, "rf")
# 
# 
# result_df <- read.csv("result_df_PRIMMUK_data_using_all_level.csv")
# get_result_summary(result_df, "l1")
# get_result_summary(result_df, "l2")
# get_result_summary(result_df, "svmsig")
# get_result_summary(result_df, "svmrad")
# get_result_summary(result_df, "rf")
# 
# 
# result_df <- read.csv("all_result_df_PRIMMUK_data_using_all_level.csv")
# get_result_summary(result_df, "l1")
# get_result_summary(result_df, "l2")
# get_result_summary(result_df, "svmsig")
# get_result_summary(result_df, "svmrad")
# get_result_summary(result_df, "rf")
# 
# get_result_summary(result_df, "l1", 2)
# get_result_summary(result_df, "l2", 2)
# get_result_summary(result_df, "svmsig", 2)
# get_result_summary(result_df, "svmrad", 2)
# get_result_summary(result_df, "rf", 2)
# 
# get_result_summary(result_df, "l1", 3)
# get_result_summary(result_df, "l2", 3)
# get_result_summary(result_df, "svmsig", 3)
# get_result_summary(result_df, "svmrad", 3)
# get_result_summary(result_df, "rf", 3)
# 
# get_result_summary(result_df, "l1", 5)
# get_result_summary(result_df, "l2", 5)
# get_result_summary(result_df, "svmsig", 5)
# get_result_summary(result_df, "svmrad", 5)
# get_result_summary(result_df, "rf", 5)
# 
# get_result_summary(result_df, "l1", 9)
# get_result_summary(result_df, "l2", 9)
# get_result_summary(result_df, "svmsig", 9)
# get_result_summary(result_df, "svmrad", 9)
# get_result_summary(result_df, "rf", 9)
# 
# get_result_summary(result_df, "l1", 10)
# get_result_summary(result_df, "l2", 10)
# get_result_summary(result_df, "svmsig", 10)
# get_result_summary(result_df, "svmrad", 10)
# get_result_summary(result_df, "rf", 10)

result_df <- full_data_result

model = "rf"
plot_title = "Performance of Random Forest model with median metrics"
plot_variation <- function(result_df, model = "rf",
                           plot_title = "Performance of Random Forest model with median metrics"){
  data_to_plot <- result_df %>%
    filter(cm == model)
  
  
  #create plot with fitted curve for training results
  data_to_plot <- data_to_plot %>%
    filter(is.na(iter))
  
  
  #create plot with fitted curve for test results
  data_to_plot <- data_to_plot %>%
    filter(!is.na(iter)) %>%
    group_by(samples) %>%
    summarise(med_acc = median(acc), med_auc = median(auc)) %>%
    pivot_longer(!samples, names_to = "metric", values_to = "val")
  
  metric_name <- "med_auc"
  data_to_plot <- data_to_plot %>%
    filter(metric == metric_name)
  
  data_to_plot <- data_to_plot %>%
    filter(samples > 32)
  
  
  ggplot(data_to_plot) +
    geom_line(aes(x = samples, y = val, color = metric)) +
    geom_point(aes(x = samples, y = val, color = metric)) +
    xlab("Number of Samples") +
    ylab("Metric Value") +
    ggtitle(plot_title)
  ggsave("variation_rf_median.png")
  
  
  
  x <- data_to_plot$samples
  y <- data_to_plot$val
  
  plot(x, y)
  
  m <- nls(y ~ log(x + a) + b)
  
  m <- nls(y ~ a * x / (b + x))
  
  m <- nls(y ~ a * x / (b + x), start = list(a = a_start,
                                             b = b_start))

  #get some estimation of goodness of fit
  cor(y,predict(m))
  #plot
  plot(x,y)
  lines(x,predict(m),lty=2,col="red",lwd=3)
  
}