# Breast cancer data 

library(readr)
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

# Question 2 
# Start by simulating normally distributed data 

set.seed(7334482)
# Simulate two nicely separated datasets,
# that does not seem relevant in the context of the classification problem
n <- 30
mu <- matrix(c(-1.5, 1, 1.5, 1.5, 1, -1.5), nrow = 3, ncol = 3)
Sigma <- diag(c(1.5, 1.5, 1.5))
# Generate variables
X_train <- MASS::mvrnorm(n,mu[1,],Sigma)
y_train <- rep(0,n)
X_train <- rbind(X_train,MASS::mvrnorm(n, mu[2,],Sigma))
y_train <- c(y_train, rep(1,n))
X_train <- rbind(X_train,MASS::mvrnorm(n, mu[3,],Sigma))
y_train <- c(y_train, rep(1,n))
### Note, you can change to different shapes Sigma for different classes and
### use different sample sizes as well of course.

# classic data frame use in R

data_train <- as.data.frame(cbind(X_train, y_train))
names(data_train) <- c("x1","x2","x3","class")

## ----2021-lda-qda-example, fig.width=5, fig.height=2.3, fig.align="center", echo=FALSE, warning=FALSE, results=FALSE, cache=TRUE, message=FALSE----
# Same data as for nearest centroid example
fit_qda <- MASS::qda(class ~ x1 + x2, data_train)

y_pred_qda <- predict(fit_qda, X_pred)$class

data_pred_da <- cbind(
  data_pred,
  tibble( class_qda = y_pred_qda)) %>%
  rename(class_centroid = class)

### Up to here it's pretty straight forward
## plotting principle same as above - try doing one at a time at home.








