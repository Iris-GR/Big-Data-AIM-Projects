# Question 3 - without stratification

# Breast cancer data 

library(readr)
library(corrplot)
library(tibble)
library(dplyr)
library(purrr)
library(caret)
library(tidyverse)
library(latex2exp) # Latex in ggplot2 labels
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour-blind friendly palette


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

uci_bc_data <- read_delim(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
  delim = ",",
  col_names = names,
  col_types = cols(
    .default = col_number(),
    id_number = col_integer(),
    diagnosis = col_factor()
  ))

# Define the response variable:
# Note - 1 = malignent and 0 = benign

y <- uci_bc_data %>% 
  mutate(diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1))) %>%
  select(diagnosis) %>%
  as_vector() %>%
  unname()

# Define the predictor
X <- uci_bc_data %>% 
  select(-id_number, -diagnosis) %>%
  as.matrix()

# Define a new data frame with data:

data <- mutate(uci_bc_data, diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1)))

# Randomize row order in Data_qda
set.seed(124)
data <- data[sample(1:nrow(data)), ] # The row order of the dataframe
# is randomised


# Plot the generated QDA dataset for overview 

ggplot(data, aes(x1 = radius_mean, x2 = diagnosis)) +
  geom_point(aes(colour = as.factor(diagnosis)))

# Split the datasets into test and training
set.seed(130)
inTraining <- sample(1:427, replace = FALSE)
training <- data[ inTraining,]
testing  <- data[-inTraining,]


# Here the data seems to once again be randomised by the function 
# CreateDataPartition, as noted by:
#The function createDataPartition can be used to create balanced splits of 
# the data. If the y argument to this function is a factor, the random sampling 
# occurs within each class and should preserve the overall class distribution of 
# the data

# Specify the type of resampling (10-fold CV, repeated 5 times)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10, 
                           repeats = 5)
# Repeated cv means what it sounds like - you do repeated cvs and get several
# prediction errors (probably to increase certainty in result or something)
# Note: The function trainControl generates parameters that further control 
# how models are created

# Fit the qda training data for QDA model using specified resampling method
fit_data <- train(diagnosis ~ .,
                  data = training,
                  method = "ranger",
                  trControl = fitControl)
fit_data

# Note: To obtain predicted class probabilities within the resampling process, 
# the argument classProbs in trainControl must be set to TRUE. 
# This merges columns of probabilities into the predictions generated from 
# each resample (there is a column per class and the column names are the class names).

# Also note: it is through the argument trControl that we specify that we
# want to do repeated cv for qda 


# Use QDA model fitted on QDA training data to predict QDA testing data class
model_test_data <- as.factor(predict(fit_data,
                                     newdata = testing))

# (Instead of getting classes as output, get the probability for each class)
predict_testing_prob <- as.matrix(predict(fit_data,
                                          newdata = testing))


# Confusion matrix for QDA model fitted on QDA training data to predict QDA 
# testing data class (data = predicted by model, reference = observed/true)
CM_model_test_data <-
  confusionMatrix(data = model_test_data,
                  reference = as.factor(testing$diagnosis))
CM_model_test_data

## CART model!!
# Fit the qda training data for CART model using specified resampling method
# method = 'rpart' gives the complexity parameter as tuning parameter, and 
# method = 'rpart2' gives the tree max depth as tuning parameter
qda_fit_data <- train(diagnosis ~ .,
                      data = training,
                      method = "qda", 
                      trControl = fitControl)
qda_fit_data

# Use CART model fitted on QDA training data to predict QDA testing data class
qda_model_test_data <- as.factor(predict(qda_fit_data,
                                         newdata = testing,
))

# Confusion matrix for CART model fitted on QDA training data to predict QDA 
# testing data class (data = predicted by model, reference = observed/true)
CM_qda_model_test_data <-
  confusionMatrix(data = qda_model_test_data,
                  reference = as.factor(testing$diagnosis))
CM_qda_model_test_data