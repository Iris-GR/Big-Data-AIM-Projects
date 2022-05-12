library(glmnet)
library(ggplot2)
library(MASS)
library(base)
library(viridis)

################################################################################
############ Function for simulating data PART 1 ###############################
################################################################################
#' Simulate data for Project 3, Part 1
#' 
#' @param n Number of samples
#' @param p Number of features
#' @param sparsity Percentage of zero elements in simulated regression coefficients
#' @param SNR Signal-to-noise ratio (see explanation above)
#' @param beta_scale Scaling for the coefficient to make sure they are large
#' 
#' @return A list containing
#'     {X}{n x p matrix of features}
#'     {y}{n vector of responses}
#'     {beta}{p vector of regression coefficients}
simulate_data <- function(n, p, sparsity = 0.8, SNR = 2, beta_scale = 5) {
  X <- matrix(rnorm(n * p), ncol = p)
  
  q <- ceiling((1 - sparsity) * p)
  beta <- rep.int(0, p)
  beta[1:q] <- beta_scale * rnorm(q)
  
  sigma <- sqrt(sum((X %*% beta) ^ 2) / (n - 1)) / SNR
  
  y <- X %*% beta + sigma * rnorm(n)
  
  # Shuffle columns so that non-zero features appear
  # not simply in the first (1 - sparsity) * p columns
  idx_col <- sample(p, p, replace = FALSE)
  
  list(
    X = X[, idx_col],
    y = y,
    beta = beta[idx_col]
  )
}
###############################################################################
###############################################################################

# Settings for simulated data:
# n = [200, 700, 750]
# p = 1000
# sparsity = [0.75, 0.9, 0.95, 0.99]
# SNR = 3
# beta_scale = 5


#----------------- set up -----------------------------
n = c(200, 700, 750)
p = 1000
sparsity = c(0.75, 0.9, 0.95, 0.99)
SNR = 3
beta_scale = 5

n_test = 1000

#Arrays for saving computed MSE's
#length(n) = 3, length(sparsity) = 4, and 5 runs each
MSE_min = array(0, c(3,4,5))
MSE_1se = array(0, c(3,4,5))

#Arrays for saving sensitivity and specificity for lambda_min and lambda_1se
#for each pair of n and sparsity over 5 runs. In 4th dimension 1=sens, 2=spec.
sens_spec_min = array(0, c(3,4,5,2))
sens_spec_1se = array(0, c(3,4,5,2))

#----------------- simulation part --------------------
i = 1
j = 1
k = 1


for(i in 1:3){
  for(j in 1:4){
    for(k in 1:5){
      
      #Simulate data and split into training and testing sets
      df = simulate_data(n[i] + n_test,p,sparsity[j],SNR,beta_scale)
      df_train = list(
        X = df$X[1:n[i],],
        y = df$y[1:n[i],]
      )
      df_test = list(
        X = df$X[-(1:n[i]),],
        y = df$y[-(1:n[i])]
      )
      
      #Fit lasso model to training data using cross validation
      fit_lasso = cv.glmnet(df_train$X, df_train$y)
      
      #--------------------- Computation of MSE ---------------------
      
      #Predict outcomes from model on tets data
      predicted_min = predict(fit_lasso, df_test$X, s = "lambda.min")
      predicted_1se = predict(fit_lasso, df_test$X, s = "lambda.1se")
      
      #Calculate MSE
      MSE_min_new = sum((df_test$y - predicted_min)^2)/n_test
      MSE_1se_new = sum((df_test$y - predicted_1se)^2)/n_test
      
      #Save MSE in global array for later use
      MSE_min[i,j,k] = MSE_min_new
      MSE_1se[i,j,k] = MSE_1se_new
      
      #----------- Computation of sensitivity and specificity --------------------
      
      #Extract coefficients from fitted model
      beta_min = as.matrix(coef(fit_lasso, s = "lambda.min"))
      beta_1se = as.matrix(coef(fit_lasso, s = "lambda.1se"))
      #Remove intercept <--------------------------------------- Questionable!! OBS I don't know if this is 100% correct
      #But there seems to be no intercept in the simulated data (i.e. intercept == 0)
      #Motivation for removing: We're only interested in the coefficients that
      #are related to the features
      beta_min = beta_min[-1]
      beta_1se = beta_1se[-1]
      
      #Create logic vector: TRUE if != 0 and FALSE if == 0
      beta_min_logic = beta_min != 0
      beta_1se_logic = beta_1se != 0
      
      beta_logic = df$beta != 0
      
      #Fixing so that there are no errors
      #Remember to divide by length(beta_logic)-2 when computing
      #sensitivity and specificity because of this
      beta_min_logic = append(beta_min_logic, c(TRUE, FALSE))
      beta_1se_logic = append(beta_1se_logic, c(TRUE, FALSE))
      beta_logic = append(beta_logic, c(FALSE, TRUE))
      
      #Confusion matrix
      tb_min = table(beta_min_logic, beta_logic)
      tb_1se = table(beta_1se_logic, beta_logic)
      
      #Calculations of sensitivity and specificity
      #(Note p = #features, defined above)
      sens_min = tb_min[2,2]/p
      spec_min = tb_min[1,1]/p
      
      sens_1se = tb_1se[2,2]/p
      spec_1se = tb_1se[1,1]/p
      
      #Save metrics into global array for later use
      sens_spec_min[i,j,k,] = c(sens_min,spec_min)
      sens_spec_1se[i,j,k,] = c(sens_1se,spec_1se)
    }  
  }
}

MSE_min_plotdata = NULL
MSE_1se_plotdata = NULL
MSE_difference_plotdata = NULL

for(i in 1:3){
  for(j in 1:4){
    MSE_min_plotdata = rbind(MSE_min_plotdata,c(n[i],sparsity[j],log(mean(MSE_min[i,j,]))))
    MSE_1se_plotdata = rbind(MSE_1se_plotdata,c(n[i],sparsity[j],log(mean(MSE_1se[i,j,]))))
    MSE_difference_plotdata = rbind(MSE_difference_plotdata,c(n[i],sparsity[j],mean(MSE_1se[i,j,])-mean(MSE_min[i,j,])))
  }
}

MSE_min_plotdata = as.data.frame(MSE_min_plotdata)
MSE_1se_plotdata = as.data.frame(MSE_1se_plotdata)
MSE_difference_plotdata = as.data.frame(MSE_difference_plotdata)


#Plotting the logarithm of the difference between MSE_1se and MSE_min
#NOTE: No error bars included.
#Conclusions:
# - Difference is positive in all cases
# - Difference is greater the smaller the sparsity
# - Difference is greater the higher n/p
names_sparsity = c('0.75','0.90','0.95','0.99')
windows()
barplot(V3 ~., MSE_difference_plotdata,
        #legend = TRUE,
        xlab = "Sparsity",
        ylab = "Difference in MSE",
        names.arg = names_sparsity,
        legend.text = c("n = 200", "n = 700", "n = 750"))



# Gradient color
ggplot(MSE_difference_plotdata, aes(V1, V2))+
  geom_point(aes(color = V3)) +
  scale_color_viridis(option = "D")+
  theme_minimal() +
  theme(legend.position = "bottom")

#Presenting results:

#NOTE! QUESTION: <------------------------------------------ OBSERVE!!
#How does MSE behave for different n, and sparsity levels?
#Make multiplot (with boxplots) where vertically we have n
#                                     horizontally we have sparsity

#Plots for different values of "n" and "sparsity":
#Make one plot for average values
#Make one plot for standard deviations

