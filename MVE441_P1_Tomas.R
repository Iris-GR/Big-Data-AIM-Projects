library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(caret)
library(tidyverse)
library(latex2exp)

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


dfX <- as.data.frame(X)
attach(dfX); plot(V1, V2, col=c("red","blue")[V3+1]); detach(dfX)


