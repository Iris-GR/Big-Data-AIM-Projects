# MVE441, Project 1, (Iris)

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
#library(tidyverse)

#### Question 1 ####
# What assumptions do the following methods make about the distribution of 
# features for a certain class? What does this mean abstractly but also 
# visually?
# 1. Linear discriminant analysis (distributions have same variance)
# 2. Quadratic discriminant analysis (distributions can have different variance)
# 3. kNN
# 4. CART
# Also think about the impact of k in kNN and the impact of different parameters 
# on the CART model.

#### Question 2 ####
# For the following questions focus on QDA and CART.
# 
# 1. Simulate training and test data that follows the assumptions for each of 
# the two models. (QDA, CART).
# 
# 2. What are your expectations on the test error if you were to apply QDA to 
# data simulated for CART and vice versa?
# 
# 3. Train the methods on both training data sets and compare their respective 
# test errors for each dataset. Were your expectations from Step 2 correct? Any 
# surprises? Discuss! (Confusion matrix (built in function))
#   
# 4. Taking what you have learned from the project so far, can you simulate 
# data for which the target method works best? (Clarifying example: Data 
# simulated for QDA should be classified best by QDA and less well CART, and the 
# other way around)

# Note! Make sure you repeat the simulation of the data many times. A result 
# based on a single training and test dataset could be pure chance, especially 
# due to the instability of CART. Write your code in such a fashion that you can
# easily create many training and test datasets repeat the comparison described 
# above. Always report average values but always with a measure of uncertainty 
# (e.g. standard deviation).
#

#### Question 2: QDA and CART models trained on ONE QDA dataset (illustration) ####

## Generate QDA dataset having two variables x1 and x2

# Sample size (i.e number samples from a specific distribution)
n_qda <- 250

# Number of classes (number of distributions to sample from)
n_class_qda <- 4

# Mean value of distributions
mu_qda <- matrix(c(-1.5, 1.5, 1.5, -1.5, 1.5, -1.5, 1.5, -1.5),
                 nrow = n_class_qda)

# Covariance matrix for variables x1 and x2. Note difference in x1 and x2 
# direction
sigma_qda <- diag(c(1.5, 1))

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
set.seed(124)
data_qda <- Data_qda[sample(1:nrow(Data_qda)), ]

# Plot the generated QDA dataset for overview 
ggplot(data_qda, aes(x = x1, y = x2)) +
  geom_point(aes(colour = as.factor(Class)))

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

### Question 2: QDA and CART models trained on MANY QDA dataset ####

#Loop over multiple QAD datasets (for each i there is a k-fold-CV -> 1 confusion
# matrix)
N_datasets <- 100

## Generate each QDA dataset having two variables x1 and x2

# Sample size (i.e number samples from a specific distribution)
n_qda <- 250

# Number of classes (number of distributions to sample from)
n_class_qda <- 4

# Mean value of distributions
mu_qda <- matrix(c(-1.5, 1.5, 1.5,-1.5, 1.5,-1.5, 1.5,-1.5),
                 nrow = n_class_qda)

# Covariance matrix for variables x1 and x2. Note difference in x1 and x2
# direction
sigma_qda <- diag(c(1.5, 1))

# Initiate a list to collect the confusion matrix for each simulated dataset i
# where a QDA model has been used for prediction
CMs_qda_qda <- list()

# ... and where a CART model has been used for prediction
CMs_cart_qda <- list()

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
  
  # Store in list
  CMs_qda_qda[[i]] <- CM_qda_model_qda_test_data$table
  
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
}

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

# Print out results
sprintf('The mean confusion matrix of %i simulated QDA datasets predicted by QDA model and its standard deviation', N_datasets)
Mean_CMs_qda_qda
std_qda_qda

sprintf('The mean confusion matrix of %i simulated QDA datasets predicted by CART model and its standard deviation', N_datasets)
Mean_CMs_cart_qda
std_cart_qda


#### Load and get overview of breast cancer dataset ####
# 1. Investigate if the class labels are balanced or imbalanced.
# 2. Investigate the features
#   a) Are they numerical or categorical? 
#   b) Do they have highly varying scales? Should data be scaled/normalized/
#      centred before being used in a classification method? (More on this in
#      lecture 5)
#   c) How does the correlation matrix between the features look like? Are there 
#      highly correlated features? Are there plausible reasons for the 
#      correlations?

# Factor names
names <- c("id_number", "diagnosis", "radius_mean",
           "texture_mean", "perimeter_mean", "area_mean",
           "smoothness_mean", "compactness_mean",
           "concavity_mean","concave_points_mean",
           "symmetry_mean", "fractal_dimension_mean",
           "radius_se", "texture_se", "perimeter_se",
           "area_se", "smoothness_se", "compactness_se",
           "concavity_se", "concave_points_se",
           "symmetry_se", "fractal_dimension_se",
           "radius_worst", "texture_worst",
           "perimeter_worst", "area_worst",
           "smoothness_worst", "compactness_worst",
           "concavity_worst", "concave_points_worst",
           "symmetry_worst", "fractal_dimension_worst")

# Load the breast cancer data set
uci_bc_data <- read_delim(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
  delim = ",",
  col_names = names,
  col_types = cols(
    .default = col_number(),
    id_number = col_integer(),
    diagnosis = col_factor()
  ))

# Define the response variable y 
y <- uci_bc_data %>% 
  mutate(diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1))) %>%
  select(diagnosis) %>%
  as_vector() %>%
  unname()

# Define the predictors (remove id and response variable y)
X <- uci_bc_data %>% 
  select(-id_number, -diagnosis) %>%
  as.matrix()

# Total number of observations
tot.nr.obs <- length(y)
tot.nr.obs

# Number observations for each class
table(y)
nr.class <- table(y)

# Proportion of observations in clas 1 (malignant tumor)
nr.class[2]/tot.nr.obs

# Are the features (predictors) numerical or categorical? They are numerical
lapply(X[1,], class)

# Spread
boxplot(X)
sprintf('The features area_mean and are_worst have much larger scale than the others')
sprintf('Probably be scaled/normalized/centered')

# Near zero variance predictors/freatures. Appears to be okay.
nzv <- nearZeroVar(X, saveMetrics= TRUE)

# Find linear combinations among features. Appears okay.
comboInfo <- findLinearCombos(X)
comboInfo

# Correlation matrix for features (X matrix). Yes, some appear highly positively 
# correlated (radii and area for ex.)
corr.x <- cor(X)
corrplot(corr.x)

# Other method of finding correlations
highCorr <- sum(abs(corr.x[upper.tri(corr.x)]) > .9)
#highCorr <- sum(abs(corr.x[upper.tri(corr.x)]) > .99)
#highCorr <- sum(abs(corr.x[upper.tri(corr.x)]) > .999)
summary(corr.x[upper.tri(corr.x)])

# Filter out features with high correlation
X_cor_high <- findCorrelation(corr.x, cutoff = .75)
X_high_cor_filtered <- X[,-X_cor_high]
cor_X_high_cor_filtered <- cor(X_high_cor_filtered)
summary(cor_X_high_cor_filtered[upper.tri(cor_X_high_cor_filtered)])

#### Question 3 #### 
# 1. Choose a classification methods (classifier) and test the effect of 
#    different classification metrics (such as accuracy, specificity, the F1 
#    score, etc.). What I mean by this is the following
#    - Split your data into a bunch of folds (do this without stratification for 
#      now)
#    - For each fold F:
#       - Train a model on the remaining training folds using the models 
#         standard mode of training (e.g. optimisation of a likelihood in the 
#         case of QDA, splitting on the Gini score in CART or Random Forest)
#       - Compute classification metrics on the test fold F.
#   - Report the average performance across folds
# 
# 2. Which classification metric(s?) is/are most suitable for our goal here? 
#    Explain why! (Could maybe be specificity, sensitivity, ROC, AUC)
#
# 3. Choose one additional classification method that is substantially different
#    (i.e. QDA and Random Forest, or logistic regression and CART. Not QDA and 
#    LDA, or logistic regression and LDA, they are too similar in their 
#    assumptions) and compare the two methods using the best classification 
#    metric(s) you determined in Step 2.
# 
# 4. Does using stratified cross-validation change/improve your results?

# Create dataset starting fresh from the uci_bc_data. Change diagnosis to 0-1
# factor and remove Id number
data <- mutate(uci_bc_data, diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1)))
data <- data[,-1]

# Dataset with 0-1 converted diagnosis, without id number, and features with 
# high correlation are removed (X_high_cor_filtered - calculated in previous 
# section). Very ugly, need to improve...
# data <- data %>% select(-(c(colnames(X_high_cor_filtered)))) 
# Tuning parameters seem a bit worse for data with correlated variables removed

### Split into training and testing data (random sampling WITHOUT stratification)
set.seed(122)
in_training <- sample(1:dim(data)[1],0.75*dim(data)[1], replace = FALSE)
data_train <- data[ in_training,]
data_test  <- data[-in_training,]


# Specify the type of sampling 
# fitControl <- trainControl(method = "cv",
#                            number = 10)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5)

#QDA----------------------------------------------------------------------------
# Fit the training data for QDA model for accuracy
qda_fit_data <- train(diagnosis ~ .,
                      data = data_train,
                      method = "qda",
                      metric = "Accuracy",
                      trControl = fitControl,
                      preProcess = c("center", "scale")) # could add corr correction already here

qda_fit_data
# Use QDA model fitted on training data to predict class of testing data
qda_model_test_data <- as.factor(predict(qda_fit_data,
                                              newdata = data_test))


# Confusion matrix for QDA model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_qda <-confusionMatrix(data = qda_model_test_data,
                  reference = as.factor(data_test$diagnosis))
CM_qda

# ROC curve---------------------------------------------------------------------

#.....

# Random Forest-----------------------------------------------------------------
# There are very many different types of random forests...
# Fit the training data for random forest model for accuracy
rf_fit_data <- train(diagnosis ~ .,
                      data = data_train,
                      method = 'ranger',
                      trControl = fitControl,
                      preProcess = c("center", "scale")) # could add corr correction already here

rf_fit_data
# Use QDA model fitted on training data to predict class of testing data
rf_model_test_data <- as.factor(predict(rf_fit_data,
                                         newdata = data_test))


# Confusion matrix for QDA model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_rf <-confusionMatrix(data = rf_model_test_data,
                         reference = as.factor(data_test$diagnosis))
CM_rf

# Add selcetion function... ?best

#-------------------------------------------------------------------------------


### Split into training and testing data (WITH stratification)
set.seed(123)
in_training_strat <- createDataPartition(data$diagnosis, p = .75, list = FALSE)
data_train_strat <- data[ in_training_strat,]
data_test_strat  <- data[-in_training_strat,]


# Specify the type of sampling 
# fitControl_strat <- trainControl(method = "cv",
#                                 number = 10)
fitControl_strat <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5)

#QDA----------------------------------------------------------------------------
# Fit the training data for QDA model for accuracy
qda_fit_data_strat <- train(diagnosis ~ .,
                      data = data_train_strat,
                      method = "qda",
                      metric = "Accuracy",
                      trControl = fitControl_strat,
                      preProcess = c("center", "scale")) # could add corr correction already here

qda_fit_data_strat
# Use QDA model fitted on training data to predict class of testing data
qda_model_test_data_strat <- as.factor(predict(qda_fit_data_strat,
                                         newdata = data_test_strat))


# Confusion matrix for QDA model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_qda_strat <-confusionMatrix(data = qda_model_test_data_strat,
                         reference = as.factor(data_test_strat$diagnosis))
CM_qda_strat
# Random Forest-----------------------------------------------------------------
# There are very many different types of random forests...
# Fit the training data for random forest model for accuracy
rf_fit_data_strat <- train(diagnosis ~ .,
                     data = data_train_strat,
                     method = 'ranger', # There is also 'rf'
                     trControl = fitControl_strat,
                     preProcess = c("center", "scale")) # could add corr correction already here

rf_fit_data_strat
# Use QDA model fitted on training data to predict class of testing data
rf_model_test_data_strat <- as.factor(predict(rf_fit_data_strat,
                                        newdata = data_test_strat))


# Confusion matrix for QDA model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_rf_strat <-confusionMatrix(data = rf_model_test_data_strat,
                        reference = as.factor(data_test_strat$diagnosis))
CM_rf_strat

#-------------------------------------------------------------------------------


#### Question 4 b) 

# Create dataset starting fresh from the uci_bc_data. Change diagnosis to 0-1
# factor and remove Id number
data <- mutate(uci_bc_data, diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1)))
data <- data[,-1]

# Sample two objects from the majority class at random
# Create a new X-feature by taking a weighted average of the two X-samples with 
# weights .95 and .05 for example, the point being that it should be a feature 
# close but not identical to an original sample.

# Generate imbalanced data out of the breast cancer dataset
set.seed(124)

# Filter out the data with class '0' as that is the majority class
data_class_0 <- data %>% filter(diagnosis == '0')

# Convert diagnosis from factor to numeric to be able to multiply with weight
data_class_0$diagnosis <- as.numeric(as.character(data_class_0$diagnosis))

# Sample two objects (two rows) from the data_class_0 and create new row
# Clear New_rows 
nr_new_rows <- 100
New_rows = NULL
for (i in 1:nr_new_rows) {
  set.seed(i)
  row.num <- sample(1:dim(data_class_0)[1], 2, replace = FALSE)
  new_2_rows <- data_class_0[row.num,] * c(0.95, 0.05, nrow(2))
  new_row <- new_2_rows[1,]+ new_2_rows[2,]
  New_rows <- rbind(new_row, New_rows)
}

# Convert back diagnosis from numeric to factor in the New rows
New_rows$diagnosis <- as.factor(as.numeric(New_rows$diagnosis))

# Construct a more imbalanced dataset by combining data and New_rows
imbal_data <- rbind(data, New_rows)


### Split imbalanced data into training and testing data (WITH stratification)
set.seed(125)
in_imbal_train <- createDataPartition(imbal_data$diagnosis, 
                                      p = .75, 
                                      list = FALSE)
imbal_train <- data[ in_imbal_train,]
imbal_test <- data[-in_imbal_train,]

# Down-sampling
down_train <- downSample(x = imbal_train[, -ncol(imbal_train)],
                         y = imbal_train$diagnosis)
table(down_train$diagnosis) 

# Up-sampling
up_train <- upSample(x = imbal_train[, -ncol(imbal_train)],
                     y = imbal_train$diagnosis)                         
table(up_train$diagnosis) 

# SMOTE
smote_train <- SMOTE(diagnosis ~ ., data  = imbal_train)                         
table(smote_train$diagnosis)

# ROSE
rose_train <- ROSE(diagnosis ~ ., data  = imbal_train)$data                         
table(rose_train$diagnosis) 


trl <- trainControl(method = "repeatedcv", repeats = 5,
                    classProbs = TRUE,
                    summaryFunction = twoClassSummary)

# Bagged classification and estimate the area under the ROC curve using five 
# repeats of 10-fold CV (straight from the documentation example so try out 
# other methods also...) 
set.seed(127)
orig_fit <- train(Class ~ ., data = imbal_train, 
                  method = "treebag",
                  nbagg = 50,
                  metric = "ROC",
                  trControl = ctrl)

down_outside <- train(Class ~ ., data = down_train, 
                      method = "treebag",
                      nbagg = 50,
                      metric = "ROC",
                      trControl = ctrl)

up_outside <- train(Class ~ ., data = up_train, 
                    method = "treebag",
                    nbagg = 50,
                    metric = "ROC",
                    trControl = ctrl)

rose_outside <- train(Class ~ ., data = rose_train, 
                      method = "treebag",
                      nbagg = 50,
                      metric = "ROC",
                      trControl = ctrl)


smote_outside <- train(Class ~ ., data = smote_train, 
                       method = "treebag",
                       nbagg = 50,
                       metric = "ROC",
                       trControl = ctrl)

outside_models <- list(original = orig_fit,
                       down = down_outside,
                       up = up_outside,
                       SMOTE = smote_outside,
                       ROSE = rose_outside)

# Compare methods and summarize
# outside_resampling <- resamples(outside_models)
# 
# test_roc <- function(model, data) {
#   roc_obj <- roc(data$Class, 
#                  predict(model, data, type = "prob")[, "Class1"],
#                  levels = c("Class2", "Class1"))
#   ci(roc_obj)
# }
# 
# outside_test <- lapply(outside_models, test_roc, data = imbal_test)
# outside_test <- lapply(outside_test, as.vector)
# outside_test <- do.call("rbind", outside_test)
# colnames(outside_test) <- c("lower", "ROC", "upper")
# outside_test <- as.data.frame(outside_test)
# 
# summary(outside_resampling, metric = "ROC")




















