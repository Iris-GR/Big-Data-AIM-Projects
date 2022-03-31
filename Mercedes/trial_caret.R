# copy of Irises code to gain understanding 

# Load libraries
library(caret)
library(MASS)
library(rpart)
library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(tidyverse)

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

#### Question 2: QDA and CART models trained on ONE QDA dataset ####

## Generate QDA dataset having two variables x1 and x2

# Sample size (i.e number samples from a specific distribution)
n_qda <- 250

# Number of classes (number of distributions to sample from)
n_class_qda <- 4

# Mean value of distributions
mu_qda <- matrix(c(-1.5, 1.5, 1.5, -1.5, 1.5, -1.5, 1.5, -1.5),
                 nrow = n_class_qda)

# Note: here we specify mean for both x1 and x2 for each of the distributions 
# corresponding to the n_class_qda classes we have. 

# Covariance matrix for variables x1 and x2
sigma_qda <- diag(c(1.5, 1))

# Note: for qda, covariance of x1 and x2 may vary

# Generate variables x1 and x2
set.seed(12344)
x_qda <- do.call(rbind, lapply(1:4, function(i) MASS::mvrnorm(n_qda, 
                                                              mu_qda[i,], 
                                                              sigma_qda)))
# Note: lapply here is used to build a list, where a function is used to generate
# the distributions with mean corresponding to the ith row in mu_qda and std dev
# corresponding to the values specified for x1 and x2 in sigma_qda.

# Assign the class/response variable labels
classes_qda <- c(rep("class_A", n_qda),
                 rep("class_B", n_qda), 
                 rep("class_C", n_qda),
                 rep("class_D", n_qda)) 

# Combine variables and classes (Note, ordered list.)
Data_qda <- tibble(x1 = x_qda[,1],
                   x2 = x_qda[,2],
                   Class = classes_qda)

# Tibble seems to be a convenient tool for creating a dataframe with specified
# column names and data input as seen in the case above. 

# Randomize row order in Data_qda
set.seed(124)
data_qda <- Data_qda[sample(1:nrow(Data_qda)), ] # The row order of the dataframe
# is randomised


# Plot the generated QDA dataset for overview 
ggplot(data_qda, aes(x = x1, y = x2)) +
  geom_point(aes(colour = as.factor(Class)))

# Split the datasets into test and training
inTraining_qda <- createDataPartition(data_qda$Class, p = .75, list = FALSE)
training_qda <- data_qda[ inTraining_qda,]
testing_qda  <- data_qda[-inTraining_qda,]


# Here the data seems to once again be randomised by the function 
# CreateDataPartition, as noted by:
#The function createDataPartition can be used to create balanced splits of 
# the data. If the y argument to this function is a factor, the random sampling 
# occurs within each class and should preserve the overall class distribution of 
# the data

# Specify the type of resampling (10-fold CV, repeated 5 times)
fitControl_qda <- trainControl(method = "repeatedcv",
                               number = 10, 
                               repeats = 5)
# Repeated cv means what it sounds like - you do repeated cvs and get several
# prediction errors (probably to increase certainty in result or something)
# Note: The function trainControl generates parameters that further control 
# how models are created

# Fit the qda training data for QDA model using specified resampling method
qda_fit_qda_data <- train(Class ~ .,
                          data = training_qda,
                          method = "qda",
                          trControl = fitControl_qda)
qda_fit_qda_data

# Note: To obtain predicted class probabilities within the resampling process, 
# the argument classProbs in trainControl must be set to TRUE. 
# This merges columns of probabilities into the predictions generated from 
# each resample (there is a column per class and the column names are the class names).

# Also note: it is through the argument trControl that we specify that we
# want to do repeated cv for qda 


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


#### Question 2: QDA and CART models trained on one CART dataset ####





#### Question 2: QDA and CART models trained on MANY QDA dataset ####
# Loop over multiple QAD datasets (for each i there is k-fold-CV)
N_datasets <- 1000
# for i in 1:N_datasets {}
set.seed(i)

#### Question 2: QDA and CART models trained on MANY CART dataset ####

#### Load and get overview of breast cancer dataset ####
# 1. Investigate if the class labels are balanced or imbalanced.
# 2. Investigate the features
#   a) Are they numerical or categorical? (lapply(dataset name, class))
#   b) Do they have highly varying scales? Should data be scaled/normalized/
#      centred before being used in a classification method? (More on this in
#      lecture 5)
#   c) How does the correlation matrix between the features look like? Are there 
#      highly correlated features? Are there plausible reasons for the 
#      correlations?

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

# Define the predictors
X <- uci_bc_data %>% 
  select(-id_number, -diagnosis) %>%
  as.matrix()

#### Question 4 b) ####



