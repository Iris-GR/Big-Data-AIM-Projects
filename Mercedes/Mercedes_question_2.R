
# General purpose packages (data handling, pretty plotting, ...)
library(tidyverse)
library(latex2exp) # Latex in ggplot2 labels
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour-blind friendly palette

## ----2021-mulit-class-logistic-regression-example, fig.width=5, fig.height=2.8, fig.align="center", echo=FALSE, warning=FALSE, results=FALSE, cache=TRUE, message=FALSE----
# Keep all classes this time
data_train <- as_tibble(datasets::iris)

# Fit multinomial logistic model == multi-class logistic model
# nnet::multinom takes a factor on the left hand side and takes the
# first factor level as a reference class
# levels(data_train$Species) # returns the factor levels in order
#                # this shows that "setosa" is taken as a reference
fit <- nnet::multinom(Species ~ Sepal.Length, data_train)

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
  Species = as.factor(
    rep(levels(data_train$Species), each = n_pred)))

## Of course you can use the more simple plotting commands here too. par(mfrow=c(3,1))
## is a way to organize subpanels for a plot (here 3 rows and 1 column)
## However, ggplot does produce nicer plots - just take your time to ensure you
## plot what you want in the right place if you reuse this kind of code. Careful
## with what you put into a data frame and reuse of variable names etc.


## ----2021-nearest-centroids-example, fig.width=3.5, fig.height=3.3, fig.align="center", echo=FALSE, warning=FALSE, results=FALSE, cache=TRUE, message=FALSE----
# Keep all classes this time
data_train <- as_tibble(datasets::iris) %>%
  dplyr::select(Species, Sepal.Length, Sepal.Width) %>%
  rename(class = Species, x1 = Sepal.Length, x2 = Sepal.Width)


# Calculate centroids per class
data_centroid <- data_train %>%
  group_by(class) %>%
  summarise(
    x1 = mean(x1),
    x2 = mean(x2))

## dplyr is great for this kind of housekeeping tasks, of course you can go the 
## long way around
dd<-as.matrix(data_train[,-1]) # the features as numerical matrix
data_centroid <- apply(dd[data_train[,1]=="setosa",],2,mean)
data_centroid <- rbind(data_centroid, apply(dd[data_train[,1]=="versicolor",],2,mean))
data_centroid <- rbind(data_centroid, apply(dd[data_train[,1]=="virginica",],2,mean))

# Classify with nearest centroid method
h <- 0.03
x1s <- seq(4, 8, by = h)
x2s <- seq(1.8, 4.8, by = h)
X_pred <- expand.grid(x1s, x2s) # creates a grid from the two vectors - 
# can use rep and seq for this but the command is convenient.
colnames(X_pred) <- c("x1", "x2")
n_pred <- dim(X_pred)[1]


data_centroid <- as.data.frame(cbind(c(1,2,3),data_centroid))
names(data_centroid)<-c("class","x1","x2")
centroids <- as.matrix(data_centroid[,2:3]) 


y_pred2 <- apply(
  apply(centroids, 1,
        function(c) {
          # Calculate squared Euclidean norm for each vector in the grid
          # for the current centroid c
          rowSums(
            (X_pred - matrix(
              rep.int(c, n_pred), nrow = n_pred, byrow =TRUE)) ^ 2)
        }),
  1, which.min) # Determine if the vector is closest to centroid 1, 2, or 3


## ----2021-lda-qda-example, fig.width=5, fig.height=2.3, fig.align="center", echo=FALSE, warning=FALSE, results=FALSE, cache=TRUE, message=FALSE----
# Same data as for nearest centroid example
fit_qda <- MASS::qda(class ~ x1 + x2, data_train)

y_pred_qda <- predict(fit_qda, X_pred)$class

data_pred_da <- cbind(
  data_pred,
  tibble(class_qda = y_pred_qda)) %>%
  rename(class_centroid = class)

### Up to here it's pretty straight forward
## plotting principle same as above - try doing one at a time at home.



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

