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

predictions_QDA = cbind(test_set, predictions_QDA)

predictions_QDA %>%
  count(class, V3)

predictions_QDA %>%
  summarize(score = mean(class == V3))

