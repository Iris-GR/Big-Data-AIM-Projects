# Question 3 - breast cancer data set analysis



# LIBRARIES:
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
library(MLmetrics)
library(caret)
library(doSNOW)

# OVERVIEW 
# 1. Investigate if the class labels are balanced or imbalanced.
# 2. Investigate the features
#   a) Are they numerical or categorical? 
#   b) Do they have highly varying scales? Should data be scaled/normalized/
#      centred before being used in a classification method? (More on this in
#      lecture 5)
#   c) How does the correlation matrix between the features look like? Are there 
#      highly correlated features? Are there plausible reasons for the 
#      correlations?

# Feature names
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

# Total number of observations: 569
tot.nr.obs <- length(y)
tot.nr.obs

# Number observations for each class
table(y)
nr.class <- table(y)

# Prop class 1 (malignant) = 38%
nr.class[2]/tot.nr.obs

# Numerical features
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

# Filter out features with high correlation (cutoff at 0.75)
X_cor_high <- findCorrelation(corr.x, cutoff = .75)
X_high_cor_filtered <- X[,-X_cor_high]
cor_X_high_cor_filtered <- cor(X_high_cor_filtered)
summary(cor_X_high_cor_filtered[upper.tri(cor_X_high_cor_filtered)])

X_new <- X_high_cor_filtered


# Create new dataset with filtered data
data_new <- as.data.frame(cbind(X_new, uci_bc_data$diagnosis))
rename(data_new, diagnosis = V13)
names(data_new)[13] <- 'diagnosis'



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

# Create dataset using filtered data. Change diagnosis to 0-1
data <- mutate(data_new, diagnosis = as.factor(case_when(diagnosis == "1" ~ 1, diagnosis == "2" ~ 0)))

# factor and remove Id number
#data <- mutate(uci_bc_data, diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1)))
#data <- data[,-1]
#rename(data, diagnosis = V13)
# Dataset with 0-1 converted diagnosis, without id number, and features with 
# high correlation are removed (X_high_cor_filtered - calculated in previous 
# section). Very ugly, need to improve...
# data <- data %>% select(-(c(colnames(X_high_cor_filtered)))) 
# Tuning parameters seem a bit worse for data with correlated variables removed


#------------------------WITHOUT STRATIFICATION---------------------------------
# Split data into training and testing 
set.seed(122)
in_training <- sample(1:dim(data)[1],0.75*dim(data)[1], replace = FALSE)
data_train <- data[ in_training,]
data_test  <- data[-in_training,]


# Function for calculating f1

f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- MLmetrics::F1_Score(y_pred = data$pred,
                                y_true = data$obs,
                                positive = lev[1])
  c(F1 = f1_val)
}

# Sampling type specification: 
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           summaryFunction = f1)

#--------------------------------QDA--------------------------------------------
## AVERAGE PER FOLD

# Fit the training data for QDA model
qda_fit_data <- train(diagnosis ~ .,
                      data = data_train,
                      method = "qda",
                      metric = "f1",
                      trControl = fitControl,
                      preProcess = c("center", "scale")) # could add corr correction already here

qda_fit_data

#-------------------------------------------------------------------------------
## TEST DATA

# Use QDA model fitted on training data to predict class of testing data
qda_model_test_data <- as.factor(predict(qda_fit_data,
                                         newdata = data_test))


# Confusion matrix for QDA model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_qda <-confusionMatrix(data = qda_model_test_data,
                         positive = "1",
                         reference = as.factor(data_test$diagnosis),
                         mode = "everything")
CM_qda

#--------------------------------RF---------------------------------------------
## AVERAGE PER FOLD

# Fit the training data for random forest model 
rf_fit_data <- train(diagnosis ~ .,
                     data = data_train,
                     method = 'rf',
                     metric = "f1",
                     trControl = fitControl,
                     preProcess = c("center", "scale")) # could add corr correction already here

rf_fit_data


#-------------------------------------------------------------------------------
## TEST DATA

# Use RF model fitted on training data to predict class of testing data
rf_model_test_data <- as.factor(predict(rf_fit_data,
                                        newdata = data_test))


# Confusion matrix for RF model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_rf <-confusionMatrix(data = rf_model_test_data,
                        positive = "1",
                        mode = "everything",
                        reference = as.factor(data_test$diagnosis))
CM_rf


#------------------------WITH STRATIFICATION------------------------------------
# Split into training and testing data 
set.seed(123)
in_training_strat <- createDataPartition(data$diagnosis, p = .75, list = FALSE)
data_train_strat <- data[ in_training_strat,]
data_test_strat  <- data[-in_training_strat,]


# Specify the type of sampling 
fitControl_strat <- trainControl(method = "repeatedcv",
                                 number = 10,
                                 repeats = 5)
                                 #summaryFunction = f1)

#--------------------------------QDA--------------------------------------------
## AVERAGE PER FOLD

# Fit the training data for QDA model 
qda_fit_data_strat <- train(diagnosis ~ .,
                            data = data_train_strat,
                            method = "qda",
                            metric = "f1",
                            trControl = fitControl_strat,
                            preProcess = c("center", "scale")) 

qda_fit_data_strat

#-------------------------------------------------------------------------------
## TEST DATA

# Use QDA model fitted on training data to predict class of testing data
qda_model_test_data_strat <- as.factor(predict(qda_fit_data_strat,
                                               newdata = data_test_strat))


# Confusion matrix for QDA model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_qda_strat <-confusionMatrix(data = qda_model_test_data_strat,
                               mode = "everything",
                               positive = "1",
                               reference = as.factor(data_test_strat$diagnosis))
CM_qda_strat

#--------------------------------RF---------------------------------------------
## AVERAGE PER FOLD

# Fit the training data for random forest model
rf_fit_data_strat <- train(diagnosis ~ .,
                           data = data_train_strat,
                           metric = "Specifity",
                           method = 'rf', # There is also 'rf'
                           trControl = fitControl_strat,
                           preProcess = c("center", "scale"))

rf_fit_data_strat

#-------------------------------------------------------------------------------
## TEST DATA

# Use RF model fitted on training data to predict class of testing data
rf_model_test_data_strat <- as.factor(predict(rf_fit_data_strat,
                                              newdata = data_test_strat))


# Confusion matrix for RF model fitted on training data to predict class of
# testing data (data = predicted by model, reference = observed/true)
CM_rf_strat <-confusionMatrix(data = rf_model_test_data_strat,
                              mode = "everything",
                              positive = "1",
                              reference = as.factor(data_test_strat$diagnosis))
CM_rf_strat

#-------------------------------------------------------------------------------

