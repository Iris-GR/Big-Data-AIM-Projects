# Course MVE441, Gruop 4,
# Project 1, Question 2 - QDA assumptions data 

#### Load libraries ####
library(caret)
library(MASS)
library(rpart)
library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(corrplot)
library(e1071)
library(randomForest)
library(ranger)
library(pROC)
library(RColorBrewer)
#library(gplots)
#library(ggplot2)
library(lattice)

# Question 2: QDA and CART models trained on ONE QDA dataset (illustration) ####

## Generate QDA dataset having two variables x1 and x2
set.seed(345)
# Sample size (i.e number samples from a specific distribution)
n_qda <- 250

# Number of classes (number of distributions to sample from)
n_class_qda <- 4

# Mean value of distributions
mu_qda <- matrix(c(-1.5, 1.5, 1.5, -3, 1.5, 1.5, -1.5, -2.5),
                 nrow = n_class_qda)

# Covariance matrix for variables x1 and x2. Note difference in x1 and x2 
# direction
sigma_qda <- diag(c(1.7, 1))

# Generate variables x1 and x2
set.seed(12344)
x_qda <- do.call(rbind, lapply(1:4, function(i) MASS::mvrnorm(n_qda, 
                                                              mu_qda[i,], 
                                                              sigma_qda)))

# Assign the class/response variable labels
classes_qda <- c(rep("class_A", n_qda),
                 rep("class_B", n_qda), 
                 rep("class_C", n_qda),
                 rep("class_D", n_qda)) 

# Combine variables and classes (Note, ordered list.)
Data_qda <- tibble(x1 = x_qda[,1],
                   x2 = x_qda[,2],
                   Class = classes_qda)

# Randomize row order in Data_qda
set.seed(1245)
data_qda <- Data_qda[sample(1:nrow(Data_qda)), ]

# Plot the generated QDA dataset for overview 
ggplot(data_qda, aes(x = x1, y = x2)) +
  geom_point(aes(colour = as.factor(Class)))

xyplot(
  Data_qda$x2 ~ Data_qda$x1,
  group = Class,
  data = Data_qda,
  auto.key = list(space = "top", cex=1.7),
  xlab = list(label= "Feature 1", cex=2),
  ylab = list(label= "Feature 2", cex=2),
  jitter.x = TRUE,
  jitter.y = TRUE,
  reference.line = TRUE,
  par.settings = list(superpose.symbol = list(pch=21,
                                              cex = 1.1)))


# Split the datasets into test and training
inTraining_qda <- createDataPartition(data_qda$Class, p = .75, list = FALSE)
training_qda <- data_qda[ inTraining_qda,]
testing_qda  <- data_qda[-inTraining_qda,]

# Specify the type of resampling (10-fold CV, repeated 5 times)
fitControl_qda <- trainControl(method = "repeatedcv",
                               number = 10, 
                               repeats = 5)

# Fit the qda training data for QDA model using specified resampling method
qda_fit_qda_data <- train(Class ~ .,
                          data = training_qda,
                          method = "qda",
                          trControl = fitControl_qda)
qda_fit_qda_data

# Use QDA model fitted on QDA training data to predict QDA testing data class
qda_model_qda_test_data <- as.factor(predict(qda_fit_qda_data,
                                             newdata = testing_qda))

# (Instead of getting classes as output, get the probability for each class)
qda_predict_testing_qda_prob <- as.matrix(predict(qda_fit_qda_data,
                                                  newdata = testing_qda,
                                                  type = "prob"))

# Confusion matrix for QDA model fitted on QDA training data to predict QDA 
# testing data class (data = predicted by model, reference = observed/true)
CM_qda_model_qda_test_data <-
  confusionMatrix(data = qda_model_qda_test_data,
                  reference = as.factor(testing_qda$Class))
CM_qda_model_qda_test_data

# F1 scores per class
CM_qda_model_qda_test_data$byClass[,6]

## CART model!!
# Fit the qda training data for CART model using specified resampling method
# method = 'rpart' gives the complexity parameter as tuning parameter, and 
# method = 'rpart2' gives the tree max depth as tuning parameter
cart_fit_qda_data <- train(Class ~ .,
                           data = training_qda,
                           method = "rpart", # Note, rpart/rpart2?
                           trControl = fitControl_qda)
cart_fit_qda_data

# Use CART model fitted on QDA training data to predict QDA testing data class
cart_model_qda_test_data <- as.factor(predict(cart_fit_qda_data,
                                              newdata = testing_qda))

# Confusion matrix for CART model fitted on QDA training data to predict QDA 
# testing data class (data = predicted by model, reference = observed/true)
CM_cart_model_qda_test_data <-
  confusionMatrix(data = cart_model_qda_test_data,
                  reference = as.factor(testing_qda$Class))
CM_cart_model_qda_test_data # Observe how unstable CART is!!!

# F1 scores per class
CM_cart_model_qda_test_data$byClass[,6]

# Question 2: QDA and CART models trained on MANY QDA dataset ####

#Loop over multiple QAD datasets (for each i there is a k-fold-CV -> 1 confusion
# matrix)
N_datasets <- 10000

## Generate each QDA dataset having two variables x1 and x2

# Sample size (i.e number samples from a specific distribution)
n_qda <- 250

# Number of classes (number of distributions to sample from)
n_class_qda <- 4

# Mean value of distributions
mu_qda <- matrix(c(-1.5, 1.5, 1.5, -3, 1.5, -1.5, 1.5, -2.5),
                 nrow = n_class_qda)

# Covariance matrix for variables x1 and x2. Note difference in x1 and x2 
# direction
sigma_qda <- diag(c(1.7, 1))

# Initiate a list to collect the confusion matrix for each simulated dataset i
# where a QDA and CART models  has been used for prediction
CMs_qda_qda <- list()
CMs_cart_qda <- list()

# Initiate lists for F1 scores
F1s_class_qda_qda <- list()
F1s_mean_qda_qda <- list()
F1s_class_cart_qda <- list()
F1s_mean_cart_qda <- list()

for (i in 1:N_datasets){
  set.seed(i)
  # Generate variables x1 and x2
  x_qda <- do.call(rbind, lapply(1:4, function(j) MASS::mvrnorm(n_qda, 
                                                                mu_qda[j,], 
                                                                sigma_qda)))
  
  # Assign the class/response variable labels
  classes_qda <- c(rep("class_A", n_qda),
                   rep("class_B", n_qda), 
                   rep("class_C", n_qda),
                   rep("class_D", n_qda)) 
  
  # Combine variables and classes (Note, ordered list.)
  Data_qda <- tibble(x1 = x_qda[,1],
                     x2 = x_qda[,2],
                     Class = classes_qda)
  
  # Randomize row order in Data_qda
  set.seed(i+1)
  data_qda <- Data_qda[sample(1:nrow(Data_qda)), ]
  
  # Split the dataset into test and training
  inTraining_qda <- createDataPartition(data_qda$Class, p = .75, list = FALSE)
  training_qda <- data_qda[ inTraining_qda,]
  testing_qda  <- data_qda[-inTraining_qda,]
  
  # Specify the type of resampling (10-fold CV, repeated 5 times)
  fitControl_qda <- trainControl(method = "repeatedcv",
                                 number = 10, 
                                 repeats = 5)
  
  # Fit the qda training data for QDA model using specified resampling method
  qda_fit_qda_data <- train(Class ~ .,
                            data = training_qda,
                            method = "qda",
                            trControl = fitControl_qda)
  
  # Use QDA model fitted on QDA training data to predict QDA testing data class
  qda_model_qda_test_data <- as.factor(predict(qda_fit_qda_data,
                                               newdata = testing_qda))
  
  # Confusion matrix for QDA model fitted on QDA training data to predict QDA 
  # testing data class (data = predicted by model, reference = observed/true)
  CM_qda_model_qda_test_data <-
    confusionMatrix(data = qda_model_qda_test_data,
                    reference = as.factor(testing_qda$Class))
  
  # Store in list cinfusion matrix
  CMs_qda_qda[[i]] <- CM_qda_model_qda_test_data$table
  
  # Store in list F1 score per class, QDA model predict QDA test data
  F1_class_qda_qda <- as.matrix(CM_qda_model_qda_test_data$byClass[,6])
  colnames(F1_class_qda_qda ) <- c('F1')
  F1s_class_qda_qda[[i]]  <- F1_class_qda_qda 
  
  # Overall (arithmetic mean) F1 score for QDA model predict QDA test data
  F1s_mean_qda_qda[[i]] <- mean(CM_qda_model_qda_test_data$byClass[,6])
  
  ## CART model!!
  # Fit the qda training data for CART model using specified resampling method
  # method = 'rpart' gives the complexity parameter as tuning parameter, and 
  # method = 'rpart2' gives the tree max depth as tuning parameter
  cart_fit_qda_data <- train(Class ~ .,
                             data = training_qda,
                             method = "rpart", # Note, rpart/rpart2?
                             trControl = fitControl_qda)
  
  # Use CART model fitted on QDA training data to predict QDA testing data class
  cart_model_qda_test_data <- as.factor(predict(cart_fit_qda_data,
                                                newdata = testing_qda))
  
  # Confusion matrix for CART model fitted on QDA training data to predict QDA 
  # testing data class (data = predicted by model, reference = observed/true)
  CM_cart_model_qda_test_data <-
    confusionMatrix(data = cart_model_qda_test_data,
                    reference = as.factor(testing_qda$Class))
  # Store in list
  CMs_cart_qda[[i]] <- CM_cart_model_qda_test_data$table
  
  # Store in list F1 score per class, CART model predict QDA test data
  F1_class_cart_qda <- as.matrix(CM_cart_model_qda_test_data$byClass[,6])
  colnames(F1_class_qda_qda ) <- c('F1')
  F1s_class_cart_qda[[i]]  <- F1_class_cart_qda 
  
  # Overall (arithmetic mean) F1 score for CART model predict QDA test data
  F1s_mean_cart_qda[[i]] <- mean(CM_cart_model_qda_test_data$byClass[,6])
}
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# Compute mean and standard deviation of Confusion matrices 

# Mean confusion matrix of N_datasets simulated QDA datasets predicted by QDA 
# model and cart model respectively
Mean_CMs_qda_qda <- Reduce('+', CMs_qda_qda) / N_datasets
Mean_CMs_cart_qda <- Reduce('+', CMs_cart_qda) / N_datasets

# Variance of confusion matrices of N_datasets simulated QDA datasets predicted
# by QDA model and cart model respectively
ar_qda_qda <- array(unlist(CMs_qda_qda), 
                    c(dim(CMs_qda_qda[[1]]),
                      length(CMs_qda_qda)))
std_qda_qda <- round(apply(ar_qda_qda, c(1, 2), sd), 2)
dimnames(std_qda_qda) <- dimnames(Mean_CMs_qda_qda)

ar_cart_qda <- array(unlist(CMs_cart_qda), 
                     c(dim(CMs_cart_qda[[1]]),
                       length(CMs_cart_qda)))
std_cart_qda <- round(apply(ar_cart_qda, c(1, 2), sd), 2)
dimnames(std_cart_qda) <- dimnames(Mean_CMs_cart_qda)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# Mean F1 scores per class of N_datasets, simulated QDA datasets predicted by QDA 
# model and cart model respectively.
Mean_F1s_class_qda_qda<- Reduce('+', F1s_class_qda_qda) / N_datasets
Mean_F1s_class_cart_qda <- Reduce('+', F1s_class_cart_qda) / N_datasets

# Variance of F1 scores per class of N_datasets simulated QDA datasets predicted
# by QDA model and cart model respectively
ar_F1s_class_qda_qda <- array(unlist(F1s_class_qda_qda), 
                              c(dim(F1s_class_qda_qda[[1]]),
                                length(F1s_class_qda_qda)))
std_F1s_class_qda_qda <- round(apply(ar_F1s_class_qda_qda, c(1, 2), sd), 3)
dimnames(std_F1s_class_qda_qda) <- dimnames(Mean_F1s_class_qda_qda)

ar_F1s_class_cart_qda <- array(unlist(F1s_class_cart_qda), 
                               c(dim(F1s_class_cart_qda[[1]]),
                                 length(F1s_class_cart_qda)))
std_F1s_class_cart_qda <- round(apply(ar_F1s_class_cart_qda, c(1, 2), sd), 2)
dimnames(std_F1s_class_cart_qda) <- dimnames(Mean_F1s_class_cart_qda)


# Mean of mean F1 scores per class of N_datasets, simulated QDA datasets 
# predicted by QDA model and cart model respectively.
Mean_F1s_mean_qda_qda <- Reduce('+', F1s_mean_qda_qda) / N_datasets
Mean_F1s_mean_cart_qda <- Reduce('+', F1s_mean_cart_qda) / N_datasets

# Variance of mean F1 scores of N_datasets simulated QDA datasets predicted
# by QDA model and cart model respectively
std_F1s_mean_qda_qda <- sd(as.vector(unlist(F1s_mean_qda_qda)))
std_F1s_mean_cart_qda <- sd(as.vector(unlist(F1s_mean_cart_qda)))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# Compute mean and standard deviation of TEST missclassification rates for two 
# models also known as the error rate (error rate = 1 - accuracy)
error_rates_qda_qda <- matrix()
error_rates_cart_qda <- matrix()

for (i in 1:N_datasets){
  error_rates_qda_qda[i] <- rbind(1 - sum(diag(CMs_qda_qda[[i]]))/sum(CMs_qda_qda[[i]]))
  error_rates_cart_qda[i] <- rbind(1 - sum(diag(CMs_cart_qda[[i]]))/sum(CMs_cart_qda[[i]]))
}

mean_error_rates_qda_qda <- mean(error_rates_qda_qda)
sd_error_rates_qda_qda <- sd(error_rates_qda_qda)

mean_error_rates_cart_qda <- mean(error_rates_cart_qda)
sd_error_rates_cart_qda <- sd(error_rates_cart_qda)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Print out results
sprintf('The mean confusion matrix of %i simulated QDA datasets predicted by QDA model and its standard deviation', N_datasets)
round(Mean_CMs_qda_qda, 2)
std_qda_qda

sprintf('The F1 mean and per class of %i simulated QDA datasets predicted by QDA model and its standard deviation', N_datasets)
round(Mean_F1s_class_qda_qda, 2)
round(std_F1s_class_qda_qda, 2)

round(Mean_F1s_mean_qda_qda, 2)
std_F1s_mean_qda_qda 

sprintf('The test mean error rate of %i simulated QDA datasets predicted by QDA model and its standard deviation', N_datasets)
round(mean_error_rates_qda_qda, 2)
round(sd_error_rates_qda_qda, 2)

sprintf('The mean confusion matrix of %i simulated QDA datasets predicted by CART model and its standard deviation', N_datasets)
round(Mean_CMs_cart_qda, 2)
std_cart_qda

sprintf('The F1 mean and per class of%i simulated QDA datasets predicted by CART model and its standard deviation', N_datasets)
round(Mean_F1s_class_cart_qda, 2)
std_F1s_class_cart_qda

round(Mean_F1s_mean_cart_qda, 2)
round(std_F1s_mean_cart_qda, 2)

sprintf('The mean test error rate of %i simulated QDA datasets predicted by CART model and its standard deviation', N_datasets)
round(mean_error_rates_cart_qda, 2)
round(sd_error_rates_cart_qda, 2)
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
