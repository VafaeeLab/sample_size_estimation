library(tidyverse)
library(ranger)
library(mRMRe)

meta_data <- read.csv("data/meta_data.csv")

data <- read.csv("data/formatted_data.csv")
colnames(data)[1] <- "sample"
output_labels <- meta_data %>%
  select(c(sample, ICIresponder)) %>%
  rename(c("Label" = "ICIresponder"))

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

#############################

random_seed = 1000
set.seed(random_seed)
train_index <- caret::createDataPartition(output_labels$Label, p = .8, 
                                          list = FALSE, 
                                          times = 1)

data.train <- data[train_index, ]
label.train <- output_labels[train_index, ]

data.test <- data[-train_index, ]
label.test <- output_labels[-train_index, ]




#filter and transform

#filter features with overall abundance != 0

filtered_features <- colSums(data.train) != 0
sum(filtered_features)
data.train <- data.train[, filtered_features]
data.test <- data.test[, filtered_features]

min(data.train)
max(data.train)
sum(data.train == 0)
dim(data.train)[1] * dim(data.train)[2]

min(data.train[data.train != 0])
#4e-05

data.train <- data.train + 10^-5
data.test <- data.test + 10^-5

min(data.train)

data.train <- log(data.train)
data.test <- log(data.test)

#filter and transform end

normparam <- caret::preProcess(data.train) 
data.train <- predict(normparam, data.train)
data.test <- predict(normparam, data.test) #normalizing test data using params from train data 

colSums(data.train)


assertthat::are_equal(rownames(data.train), label.train$sample)
assertthat::are_equal(rownames(data.test), label.test$sample)



#ranger
set.seed(random_seed)

ranger_model <- ranger::ranger(x = data.train, y = factor(label.train$Label), 
                               importance = "impurity_corrected")

summary(ranger_model$variable.importance)

hist(ranger_model$variable.importance)

features <- which(ranger_model$variable.importance != 0)

features <- which(ranger_model$variable.importance >= quantile(ranger_model$variable.importance)[3])

data.train <- data.train[, features, drop = FALSE]
data.test <- data.test[, features, drop = FALSE]


assertthat::are_equal(rownames(data.train), label.train$sample)
assertthat::are_equal(rownames(data.test), label.test$sample)

#mrmr

classes = c("yes", "no")

mrmr.data.train <- mRMRe::mRMR.data(data = data.frame(
  target = factor(label.train$Label, levels = classes, ordered = TRUE),
  data.train))
filter <- mRMRe::mRMR.classic(data = mrmr.data.train, target_indices = c(1), feature_count = 400)

features <- mRMRe::solutions(filter)[[1]] - 1

data.train <- data.train[, features, drop = FALSE]
data.test <- data.test[, features, drop = FALSE]

#mrmr end



#classification

classes = c("yes", "no")

logistic_regression(data.train, label.train, 
                    data.test, label.test,
                    classes, regularize = "l2")
logistic_regression(data.train, label.train, 
                    data.test, label.test,
                    classes, regularize = "l1")
# [1] 0.5
# [1] 0.5757576
# [1] 0.5625
# [1] 0.5625 0.5000 0.0000 1.0000


# with ranger second option and l2 logreg
# [1] 0.5
# [1] 0.780303
# [1] 0.59375
# [1] 0.59375000 0.39285714 0.07142857 1.00000000


#with ranger second option and l1 logreg
# [1] 0.5
# [1] 0.5757576
# [1] 0.5625
# [1] 0.5625 0.5000 0.0000 1.0000


svm_model(data.train, label.train, data.test, label.test, 
          classes, kernel = "sigmoid")
svm_model(data.train, label.train, data.test, label.test, 
          classes, kernel = "radial")
