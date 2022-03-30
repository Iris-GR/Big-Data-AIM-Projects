library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(caret)
library(tidyverse)
library(latex2exp)
library(rpart)
library(rpart.plot)
library(MASS)

#Test error
R_te <- NULL

#Loop for calc. test error on many different samples from a single distribution
#(The speed might improve by removing unecessary plots and calculations)
for(l in 1:100){

  #Generating data set for CART:
  n = 1000
  X<-NULL
  
  for (i in 1:n) {
    x = runif(1,0,4)
    y = runif(1,0,4)
    class = -1
    
    if(x<2){
      if(y<2){
        class = 0
      } else{
        class = 1
      }
    } else if(x>3){
      if(y>3){
        class = 0
      } else{
        class = 1
      }
    } else {
      class = 1
    }
    
    X = rbind(X, c(x,y,class))
  }
  
  #Plot the generated data
  dfX <- as.data.frame(X)
  attach(dfX); plot(V1, V2, col=c("red","blue")[V3+1]); detach(dfX)
  
  #Train CART and plot the decision tree (on all data)
  dfX$V3 <- as.factor(dfX$V3)
  z = rpart(V3 ~ ., data = dfX, method = "class", cp =0, minsplit = 2, minbucket = 0)
  prp(z, type = 1, extra = 1, under = TRUE)
  
  #Splitting data into train and test
  ix = createDataPartition(dfX$V3, p = 0.75, list = FALSE)
  train_set = dfX[ix,]
  test_set = dfX[-ix,] #Don't really understand this syntax "-ix", but it works
  
  #Try training QDA-classifier and see how well it works
  model_QDA = qda(V3 ~ ., data = train_set)
  
  #Calculate test error
  predictions_QDA = data.frame(predict(model_QDA, test_set))
  
  #Concatenate dataframes
  predictions_QDA = cbind(test_set, predictions_QDA)
  
  #Score = #correct predictions / #total predictions
  #predictions_QDA %>%
  #  summarize(score = mean(class == V3))
  
  #Percentage of correct predictions
  score = mean(predictions_QDA$class == predictions_QDA$V3)
  
  #Test error using 0-1 loss
  R_te = append(R_te,1 - score)

}

hist(R_te)
mean(R_te)
sd(R_te)


#Another (more flexible) implementation of training the model
#Can add folding to this piece of code
model <- train(V3 ~ ., train_set,
               method = "qda")

predictions_QDA = data.frame(predict(model, test_set))
predictions_QDA = cbind(test_set, predictions_QDA)

#Percentage of correct predictions
score = mean(predictions_QDA$predict.model..test_set. == predictions_QDA$V3)

#Test error using 0-1 loss
R_te = 1 - score



# -------- PLOTTING DECISION BOUNDARY -----------------

#plot decision boundary
# decisionplot <- function(model, data, class = NULL, predict_type = "class",
#                          resolution = 100, showgrid = TRUE, ...) {
# 
#   if(!is.null(class)) cl <- data[,class] else cl <- 1
#   data <- data[,1:2]
#   k <- length(unique(cl))
# 
#   plot(data, col = as.integer(cl)+1L, pch = as.integer(cl)+1L, ...)
# 
#   # make grid
#   r <- sapply(data, range, na.rm = TRUE)
#   xs <- seq(r[1,1], r[2,1], length.out = resolution)
#   ys <- seq(r[1,2], r[2,2], length.out = resolution)
#   g <- cbind(rep(xs, each=resolution), rep(ys, time = resolution))
#   colnames(g) <- colnames(r)
#   g <- as.data.frame(g)
# 
#   ### guess how to get class labels from predict
#   ### (unfortunately not very consistent between models)
#   p <- predict(model, g, type = predict_type)
#   if(is.list(p)) p <- p$class
#   p <- as.factor(p)
# 
#   if(showgrid) points(g, col = as.integer(p)+1L, pch = ".")
# 
#   z <- matrix(as.integer(p), nrow = resolution, byrow = TRUE)
#   contour(xs, ys, z, add = TRUE, drawlabels = FALSE,
#           lwd = 2, levels = (1:(k-1))+.5)
# 
#   invisible(z)
# }
# 
# decisionplot(model_QDA, test_set, class = "V3")
