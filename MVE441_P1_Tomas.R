library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(caret)
library(tidyverse)
library(latex2exp)
library(rpart)
library(rpart.plot)

#Generating data set for CART:
X<-NULL

for (i in 1:100) {
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

#Train CART and plot the decision tree
dfX$V3 <- as.factor(dfX$V3)
z = rpart(V3 ~ ., data = dfX, method = "class", cp =0, minsplit = 2, minbucket = 0)
prp(z, type = 1, extra = 1, under = TRUE)

#Try splitting data into train and test and see how well it works

#Try training QDA-classifier and see how well it works

#Generalize to be able to try on many different datasets

