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

#Toy data


y <- as.factor(rep(c(1,2,3,4), each = 25))
X <- rbind(
  matrix(rnorm(100) - 1, ncol = 2),
  matrix(rnorm(100) + 1, ncol = 2)
)
dev.off()
options(repr.plot.width=6, repr.plot.height=6)
ggplot(tibble(y = y, x1 = X[,1], x2 = X[,2])) +
  geom_point(aes(x1, x2, colour = y)) +
  coord_equal() +
  scale_colour_discrete("Label") +
  theme_minimal() +
  theme(text = element_text(size=16))



fit <- MASS::qda(y ~ X)
labels_pred <- predict(fit, as.data.frame(X))$class
levels(labels_pred)

caret::confusionMatrix(y, labels_pred, positive = "1")

