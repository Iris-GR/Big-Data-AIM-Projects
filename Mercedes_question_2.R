
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

p1 <- ggplot(
  data_pred,
  aes(x = x, y = y, colour = Species)) +
  geom_line() +
  scale_x_continuous(
    name = NULL,
    labels = NULL,
    lim = c(4.2, 8.1),
    expand = c(0, 0)) +
  scale_y_continuous(name = "Probability", breaks = c(0, 0.5, 1)) +
  scale_colour_manual(values = cbPalette[-1], guide = FALSE) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(
      margin = margin(0, 3.5, 0, 0, "lines")),
    plot.margin = margin(.25, .5, 0, 0.1, "lines"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9))

p2 <- ggplot(
  tibble(
    x = rep.int(X_pred[,2], 3),
    y = c(rep.int(0, n_pred), as.vector(y_pred_lin)),
    Species = as.factor(
      rep(levels(data_train$Species), each = n_pred))),
  aes(x = x, y = y, colour = Species)) +
  geom_line() +
  scale_x_continuous(
    name = NULL,
    labels = NULL,
    lim = c(4.2, 8.1),
    expand = c(0, 0)) +
  scale_y_continuous(name = "Linear predictor") +
  scale_colour_manual(values = cbPalette[-1], guide = FALSE) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(
      margin = margin(0, 3.5, 0, 0, "lines")),
    plot.margin = margin(.25, .5, 0, 0.1, "lines"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9))

p3 <- ggplot(
  data_train,
  aes(x = Sepal.Length,
      y = as.numeric(Species),
      colour = Species)) +
  geom_jitter(height = 0, width=0.1) +
  scale_x_continuous(
    name = "Sepal Length",
    lim = c(4.2, 8.1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Species",
    breaks = c(1, 2, 3),
    labels = levels(data_train$Species)) +
  scale_colour_manual(values = cbPalette[-1], guide = FALSE) +
  theme_minimal() +
  theme(
    plot.margin = margin(0, .5, 0, 0.1, "lines"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9))

ggpubr::ggarrange(
  p1, p2, p3, nrow = 3, heights = c(1, 1, 1))

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

## Long and ugly way around...

dd<-function(xin,cs) {
  sum((xin-cs)^2)
}
DD<-apply(X_pred,1,dd,cs=centroids[1,])
DD<-cbind(DD,apply(X_pred,1,dd,cs=centroids[2,]))
DD<-cbind(DD,apply(X_pred,1,dd,cs=centroids[3,]))
y_pred<-apply(DD,1,which.min)

y_pred<-apply(X_pred,1,closestcentroid,cs=centroids)  

colvec <- cbPalette[1:3]
ycol <- y_pred
ycol[y_pred==1]<-colvec[1]
ycol[y_pred==2]<-colvec[2]
ycol[y_pred==3]<-colvec[3]
#
plot(X_pred,col=ycol,pch="S")
colvec2 <- cbPalette[4:6]
ycol2 <- rep(0,length(data_train$class))
ycol2[data_train$class=="setosa"]<-colvec2[1]
ycol2[data_train$class=="versicolor"]<-colvec2[2]
ycol2[data_train$class=="virginica"]<-colvec2[3]
points(data_train[,2:3],col=ycol2,pch=16)
points(centroids,col=colvec2,pch=4,lwd=2,cex=2.5)




## ----2021-lda-qda-example, fig.width=5, fig.height=2.3, fig.align="center", echo=FALSE, warning=FALSE, results=FALSE, cache=TRUE, message=FALSE----
# Same data as for nearest centroid example
fit_lda <- MASS::lda(class ~ x1 + x2, data_train)
fit_qda <- MASS::qda(class ~ x1 + x2, data_train)

y_pred_lda <- predict(fit_lda, X_pred)$class
y_pred_qda <- predict(fit_qda, X_pred)$class

data_pred_da <- cbind(
  data_pred,
  tibble(class_lda = y_pred_lda, class_qda = y_pred_qda)) %>%
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

