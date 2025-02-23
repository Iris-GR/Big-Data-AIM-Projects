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
Sigma <- diag(c(1, 2, 3))
# Generate variables
X_train <- MASS::mvrnorm(n,mu[1,],Sigma)
y_train <- rep(1,n)
X_train <- rbind(X_train,MASS::mvrnorm(n, mu[2,],Sigma))
y_train <- c(y_train, rep(2,n))
X_train <- rbind(X_train,MASS::mvrnorm(n, mu[3,],Sigma))
y_train <- c(y_train, rep(3,n))


# Simulate two nicely separated datasets,
# that does not seem relevant in the context of the classification problem

#mu <- matrix(c(-1.5, 1.5, 1.5, -1.5), ncol = 2)
#Sigma <- diag(c(1.5, 1.5))
# Generate variables
#X_train <- MASS::mvrnorm(n,mu[1,],Sigma)
#y_train <- rep(0,n)
#X_train <- rbind(X_train,MASS::mvrnorm(n, mu[2,],Sigma))
#y_train <- c(y_train, rep(1,n))
### Note, you can change to different shapes Sigma for different classes and
### use different sample sizes as well of course.

# classic data frame use in R

data_train <- as.data.frame(cbind(X_train, y_train))
names(data_train) <- c("x1","x2","x3", "class")

## dplyr is great for this kind of housekeeping tasks, of course you can go the 
## long way around
#dd<-as.matrix(data_train[,-4]) # the features as numerical matrix
#data_centroid <- apply(dd[data_train[,1]=="setosa",],2,mean)
#data_centroid <- rbind(data_centroid, apply(dd[data_train[,1]=="versicolor",],2,mean))
#data_centroid <- rbind(data_centroid, apply(dd[data_train[,1]=="virginica",],2,mean))
#                # this shows that "setosa" is taken as a reference
 U <- data_train$class
 levels(U)[1] <- "low"
 levels(U)[2] <- "medium"
 levels(U)[3] <- "high"
 
 data_train <- as.data.frame(cbind(data_train,U))
 
 
 
 
fit <- nnet::multinom(U ~ x1, data_train)

# Linear predictor
n_pred <- 100
X_pred <- matrix(c(
  rep.int(1, n_pred),
  seq(4.2, 8.1, length.out = n_pred)),
  ncol = 2)
y_pred_lin <- X_pred %*% t(unname(coef(fit)))
# Transform to probabilities
y_pred <- cbind(
  rep.int(1, n_pred), exp(y_pred_lin)) /
  (1 + rowSums(exp(y_pred_lin)))
# Collect all predicted data for plotting
data_pred <- tibble(
  x = rep.int(X_pred[,2], 3),
  y = as.vector(y_pred),
  class = as.factor(
    rep(levels(data_train$U), each = n_pred)))


#test <- rep("A", 10)
#test <- rbind(test, rep("B", 10))
#test <- rbind(test, rep("C", 10))
#levels(test)


# Calculate centroids per class\
data_centroid <- data_train %>%
  group_by(class) %>%
  summarise(
    x1 = mean(x1),
    x2 = mean(x2),
    x3 = mean(x3))

# Classify with nearest centroid method
h <- 0.03
x1s <- seq(-3, 5, by = h)
x2s <- seq(-3, 5, by = h)
x3s <- seq(-3, 5, by = h)
X_pred <- expand.grid(x1s, x2s,x3s) # creates a grid from the two vectors - 
# can use rep and seq for this but the command is convenient.
colnames(X_pred) <- c("x1", "x2","x3")
n_pred <- dim(X_pred)[1]



## ----2021-lda-qda-example, fig.width=5, fig.height=2.3, fig.align="center", echo=FALSE, warning=FALSE, results=FALSE, cache=TRUE, message=FALSE----
# Same data as for nearest centroid example
fit_qda <- MASS::qda(class ~ x1 + x2 + x3, data_train)

y_pred_qda <- predict(fit_qda, X_pred)$U

data_pred_da <- cbind(
  data_pred,
  tibble(class_qda = y_pred_qda)) %>%
  rename(class_centroid = class)


# Plotting with ggplot


p_centroid <- ggplot() +
  geom_tile(
    aes(x = x1, y = x2, colour = class, fill = class_centroid),
    colour = "transparent",
    data = data_pred_da,
    width = h,
    height = h,
    alpha = 0.4) +
  ggtitle("Nearest Centroids")

p_lda <- ggplot() +
  geom_tile(
    aes(x = x1, y = x2, colour = class, fill = class_lda),
    colour = "transparent",
    data = data_pred_da,
    width = h,
    height = h,
    alpha = 0.4) +
  ggtitle("LDA")

p_qda <- ggplot() +
  geom_tile(
    aes(x = x1, y = x2, colour = class, fill = class_qda),
    colour = "transparent",
    data = data_pred_da,
    width = h,
    height = h,
    alpha = 0.4) +
  ggtitle("QDA")

plots <- lapply(list(p_centroid, p_lda, p_qda), function(p) {
  p + geom_jitter(
    aes(x = x1, y = x2, colour = class, shape = class),
    data = data_train, height = 0.1, width=0.1, size = 1) +
    scale_colour_manual(values = cbPalette[-1], guide = FALSE) +
    scale_fill_manual(values = cbPalette[-1], guide = FALSE) +
    scale_shape_discrete(guide = FALSE) +
    guides(  # Repeatedly using the class factor messes up the legend
      # Manual setup needed
      colour = guide_legend(
        title = "Species",
        override.aes = list(
          fill = "transparent",
          colour = cbPalette[2:4],
          shape = c(16, 17, 15),
          size = 2,
          linetype = 0))) +
    scale_x_continuous(name = "Sepal Length", expand = c(0, 0)) +
    scale_y_continuous(name = "Sepal Width", expand = c(0, 0)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(1, 1.5, 0, 0, "lines"),
      plot.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      legend.title = element_text(size = 9))
})

ggpubr::ggarrange(
  plotlist = plots, ncol = 3,
  common.legend = TRUE, legend = "bottom")








