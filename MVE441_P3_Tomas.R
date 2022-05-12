library(glmnet)
library(ggplot2)
library(MASS)
library(base)

MSE_min = NULL
MSE_1se = NULL

for(i in 1:10){
  #Skapa lite data
  x = cbind(runif(100),runif(100))
  y = x[,1] + rnorm(100,0,0.1)
  
  x_train = x[1:50,]
  x_test = x[51:100,]
  
  y_train = y[1:50]
  y_test = y[51:100]
  
  #Default settings:
  #alpha = 1 (Gives Lasso penalty)
  #nfolds = 10 (10 fold CV)
  fit_lasso = cv.glmnet(x_train, y_train)
  
  predicted_min = predict(fit_lasso, x_test, s = "lambda.min")
  predicted_1se = predict(fit_lasso, x_test, s = "lambda.1se")
  
  #Calculate MSE
  MSE_min_new = sum((y_test - predicted_min)^2)/length(y_test)
  MSE_1se_new = sum((y_test - predicted_1se)^2)/length(y_test)
  
  #Append to global list
  MSE_min = append(MSE_min, MSE_min_new)
  MSE_1se = append(MSE_1se, MSE_1se_new)
}

#NOTE! QUESTION: <------------------------------------------ OBSERVE!!
#How does MSE behave for different n, and sparsity levels?
#Make multiplot (with boxplots) where vertically we have n
#                                     horizontally we have sparsity
boxplot(cbind(MSE_min,MSE_1se))


#Not necessary for prediction, maybe for feature selection
#Get coefficients
beta_min = as.matrix(coef(fit_lasso, s = "lambda.min"))
beta_1se = as.matrix(coef(fit_lasso, s = "lambda.1se"))

#Create logic vector: TRUE if != 0 and FALSE if == 0
beta_min_logic = beta_min != 0
beta_1se_logic = beta_1se != 0

beta_logic = c(TRUE, FALSE, TRUE)

#Fixing so that there are no errors
#Remember to divide by length(beta_logic)-2 when computing
#sensitivity and specificity because of this
beta_min_logic = append(beta_min_logic, c(TRUE, FALSE))
beta_1se_logic = append(beta_1se_logic, c(TRUE, FALSE))
beta_logic = append(beta_logic, c(FALSE, TRUE))

tb_min = table(beta_min_logic, beta_logic)
tb_1se = table(beta_1se_logic, beta_logic)

p = length(beta_logic)-2
sens_min = tb_min[2,2]/p
spec_min = tb_min[1,1]/p

sens_1se = tb_1se[2,2]/p
spec_1se = tb_1se[1,1]/p

#Plots for different values of "n" and "sparsity":
#Make one plot for average values
#Make one plot for standard deviations
