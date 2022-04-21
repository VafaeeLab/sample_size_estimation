library(tidyverse)
library(ranger)
library(mRMRe)
library(randomForest)
library(psdR)
library(umap)

meta_data <- read.csv("data/meta_data.csv")

summary(factor(meta_data$subcohort))


data <- read.csv("data/formatted_data.csv")

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


for(ss_i in c(1:10)){
  if(ss_i >= 2){
    train_index <- caret::createDataPartition(output_labels$Label, p = .9, 
                                              list = FALSE, 
                                              times = 1)
    data <- data[train_index, ]
    output_labels <- output_labels[train_index, ]    
  }
  print(ss_i)
  print(dim(data))
}

set.seed(1000)
train_index <- caret::createMultiFolds(y = output_labels$Label, k = 5, times = 10)

result_df <- data.frame(matrix(nrow = 0, ncol = 4, dimnames = list(c(), 
                                                                   c("cm", "iter", "acc", "auc"))
                               )
                        )

for(i in c(1:50)){
  # i <- 1
  
  data.train <- data[train_index[[i]], ]
  label.train <- output_labels[train_index[[i]], ]
  
  data.test <- data[-train_index[[i]], ]
  label.test <- output_labels[-train_index[[i]], ]  
  
  # print(assertthat::are_equal(rownames(data.train), label.train$sample))
  # print(assertthat::are_equal(rownames(data.test), label.test$sample))
  
  filtered_features <- colSums(data.train) != 0
  # sum(filtered_features)
  data.train <- data.train[, filtered_features]
  data.test <- data.test[, filtered_features]
  
  # min(data.train)
  # max(data.train)
  # sum(data.train == 0)
  # dim(data.train)[1] * dim(data.train)[2]
  
  # min(data.train[data.train != 0])
  #4e-05
  
  data.train <- data.train + 10^-5
  data.test <- data.test + 10^-5
  
  # min(data.train)
  
  data.train <- log(data.train)
  data.test <- log(data.test)
  
  # data.train <- data.frame(psd(data.train))
  # data.test <- data.frame(psd(data.test))
  
  #filter and transform end
  
  normparam <- caret::preProcess(data.train) 
  data.train <- predict(normparam, data.train)
  data.test <- predict(normparam, data.test) #normalizing test data using params from train data 
  
  # colSums(data.train)
  
  
  
  
  #ranger
  random_seed <- 1000
  set.seed(random_seed)
  
  ranger_model <- ranger::ranger(x = data.train, y = factor(label.train$Label), 
                                 importance = "impurity_corrected")
  features <- which(ranger_model$variable.importance >= 0)
  data.train <- data.train[, features, drop = FALSE]
  data.test <- data.test[, features, drop = FALSE]
  
  print(assertthat::are_equal(rownames(data.train), label.train$sample))
  print(assertthat::are_equal(rownames(data.test), label.test$sample))
  
  
  #classification model
  classes = c("yes", "no")
  
  res_l2 <- logistic_regression(data.train, label.train, 
                      data.test, label.test,
                      classes, regularize = "l2")
  res_l2 <- c(i, res_l2[1:2])
  res_l1 <- logistic_regression(data.train, label.train, 
                      data.test, label.test,
                      classes, regularize = "l1")
  res_l1 <- c(i, res_l1[1:2])
  
  
  res_svmsig <- svm_model(data.train, label.train, data.test, label.test, 
            classes, kernel = "sigmoid")
  res_svmsig <- c(i, res_svmsig[1:2])
  
  res_svmrad <- svm_model(data.train, label.train, data.test, label.test, 
            classes, kernel = "radial")
  res_svmrad <- c(i, res_svmrad[1:2])
  
  res_rf <- rf_model(data.train, label.train, data.test, label.test, 
           classes)  
  res_rf <- c(i, res_rf[1:2])
  
  
  result_df <- rbind(result_df,
                     c("l2", res_l2),
                     c("l1", res_l1),
                     c("svmsig", res_svmsig),
                     c("svmrad", res_svmrad),
                     c("rf", res_rf))
}
colnames(result_df) <- c("cm", "iter", "acc", "auc")
result_df <- result_df %>%
  mutate(acc = as.double(acc), auc = as.double(auc))

get_result_summary <- function(result_df, model){
  # model <- "l2"
  model_res <- result_df %>%
    filter(cm == model)
  print("AUC")
  print(summary(as.double(model_res$auc)))
  print("Acc")
  print(summary(as.double(model_res$acc)))
}

get_result_summary(result_df, "l1")
get_result_summary(result_df, "l2")
get_result_summary(result_df, "svmsig")
get_result_summary(result_df, "svmrad")
get_result_summary(result_df, "rf")


write.csv(result_df, "result_df_full_data_using_species.csv", row.names = FALSE)

write.csv(result_df, "result_df_PRIMMUK_data_using_species.csv", row.names = FALSE)

write.csv(result_df, "result_df_PRIMMUK_data_using_all_level.csv", row.names = FALSE)


result_df <- read.csv("result_df_full_data_using_species.csv")
get_result_summary(result_df, "l1")
get_result_summary(result_df, "l2")
get_result_summary(result_df, "svmsig")
get_result_summary(result_df, "svmrad")
get_result_summary(result_df, "rf")


result_df <- read.csv("result_df_PRIMMUK_data_using_species.csv")
get_result_summary(result_df, "l1")
get_result_summary(result_df, "l2")
get_result_summary(result_df, "svmsig")
get_result_summary(result_df, "svmrad")
get_result_summary(result_df, "rf")


result_df <- read.csv("result_df_PRIMMUK_data_using_all_level.csv")
get_result_summary(result_df, "l1")
get_result_summary(result_df, "l2")
get_result_summary(result_df, "svmsig")
get_result_summary(result_df, "svmrad")
get_result_summary(result_df, "rf")
