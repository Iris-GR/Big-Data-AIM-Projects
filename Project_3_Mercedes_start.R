
library(glmnet)

# Feature selection and regularised regression


#Goal: give us an impression of what these methods can be used for, but
# also what limitations they have.


#WARNING: Both parts of assignment require some amount of repeated simulation/
# estimation of models which can take some time to run!

## ---- PART 1: Prediction quality vs feature selection ------------------------

# Lasso encodes variable selection through penalisation. 
#install.packages("tidyverse")
#install.packages("glmnet", repos = "https://cran.us.r-project.org")
library(tidyverse)
library(latex2exp)

# Simulating data: 

#' 
#' @param n Number of samples 
#' @param p Number of features
#' @param sparsity Percentage of zero elements in simulated regression coefficients
#' @param SNR Signal-to-noise ratio (see explanation above)
#' @param scaling.beta Scaling for the coefficient to make sure they are large
#' 
#' @return A list containing
#     {X}{n x p matrix of features}
#     {y}{n vector of responses}
#     {beta}{p vector of regression coefficients}


sim.data <- function(n, p, sparsity, SNR = 2, scaling.beta = 5){
  
  X <- matrix(rnorm(n*p), ncol = p) # nXp design matrix 
  
  non.zero.coeff <- ceiling((1-sparsity)*p) # Number of non-zero coeff
  
  # ceiling will take a numeric number (such as 1.1) and return the integer
  # closest and higher in value to this number (so 2 in my example)
  
  beta <- rep.int(0,p) #initialise 1Xp vector with 0
  
  beta[1:non.zero.coeff] <- scaling.beta * rnorm(non.zero.coeff) # Insert q scaled coeff from sd normal dist
  
  # (sd is beta_scale for the betas)
  
  
  # I think what we are doing here is basically creating some noise/
  # a sample
  # beta*X is the noise-less response
  # beta*X + sigma*epsilon is the actual response
  # sigma = sd_noise = sd_signal/SNR - we generate a reasonable noise sd
  
  # (So basically like we did with beta_scale above, we want to rescale
  # our random error so that we get reasonable results)
  sd.noise <- sqrt(sum((X %*% beta)^2)/(n-1))/SNR 
  
  
  y <- X%*%beta + sd.noise * rnorm(n) # Response 
  
  # We do not want the non-zero features in the first (1-sparsity)p 
  # number of columns only, so we reorder the columns:
  
  reord <- sample(p, p, replace = F)
  
  list(
    
    X = X[, reord],
    y = y,
    beta = beta[reord]
  )
  
}

# Test n - the data we test with (see Rebeckas notes)
n.test <- 500

p = 1000
# Generate n observations from Gaussian, isometric dist (they are theoretically 
# uncorrelated)
# Rebeckas way of doing this is by creating a function that does this:

# Different ns to train with
number <- c(200,500,750)

# Different sparsity
sparse <- c(0.75,0.9, 0.95, 0.99)

return <- list()

# Number of repeated calculations for each n and sparsity
loops <- 10

mean.squared.error <- matrix(0, nrow =12, ncol = 2*loops)

beta_est <- matrix(0,nrow = p+1, ncol = 2*loops*12)

# Something to note is that I don't think true beta contains intercept
# whilst beta_est (the estimated beta) does



beta_true <- matrix(0, nrow = p, ncol = 120)

for (l in 1:loops) {
  
  for (i in 1:3) {
    
    #for (i in 1:3) {}
    #for (j in 1:4){}
    
    # CASE 1 - train your data on n observations
    n.train = number[i] 
    n = n.train + n.test
    
    for (j in 1:4) {
      sparsity = sparse[j]
      
      return <- sim.data(n,p,sparsity)
      
      # 2. Determine hyperparameters
      
      attach(return)
    
      
      
      
      
      # Define MSE function (from: https://stackoverflow.com/questions/39482436/why-calculating-mse-in-lasso-regression-gives-different-outputs)
      
      MSE <- function(x,y) { mean((x-y)^2)}
      
      # sample size
      
      cv.sample.size <- floor(n.train)
      
      # set seed to make reproductible - probably not what we want to do
      
      #set.seed(907)
      
      train.ind <- sample(seq_len(n), size = cv.sample.size)
      
      # Training set
      
      train.X <- X[train.ind, ]
      
      train.y <- y[train.ind]
      
      # Test set
      
      test.X <- X[-train.ind, ]
      
      test.y <- y[-train.ind]
      
      # Fit linear model with lasso on train data
      
      probs.lasso <- cv.glmnet(train.X, train.y,
                               type.measure = "mse",
                               keep = T,
                               alpha =1
      )
      
      
      
      
      
      lasso.lambda <- probs.lasso$lambda.min
      
      id.lambda.min <- which(probs.lasso$lambda == probs.lasso$lambda.min)
      
      # get pred values on training folds with lambda.min 
      
      MSE.1 <- MSE(probs.lasso$fit[,id.lambda.min], train.y)
      
      cat("MSE (method 1):", MSE.1, "\n")
      
      
      
      
      id.lambda.1se <- which(probs.lasso$lambda == probs.lasso$lambda.1se)
      
      # get pred values on training folds with lambda.min 
      
      MSE.2 <- MSE(probs.lasso$fit[,id.lambda.1se], train.y)
      
      cat("MSE (method 2):", MSE.2, "\n")
      
      if (i == 1 || i == 3) {
        
        mean.squared.error[i*i+j-1,l] <-  MSE.1
        mean.squared.error[i*i+j-1,l+10] <- MSE.2
        
        
        beta_true[,i*i+(j-1)+12*(l-1)] <- as.vector(return$beta)
        
        beta_est[,i*i+(j-1)+12*(l-1)] <- as.vector(coef(probs.lasso, s = "lambda.min"))
        
        beta_est[,i*i+(j-1)+12*(l-1)+120] <- as.vector(coef(probs.lasso, s = "lambda.1se"))
        
      } else if (i == 2) {
        
        mean.squared.error[i*i+j,l] <-  MSE.1
        mean.squared.error[i*i+j,l+10] <- MSE.2
        
        beta_true[,i*i+(j)+12*(l-1)] <- as.vector(return$beta)
        
        
        beta_est[,i*i+(j)+12*(l-1)] <- as.vector(coef(probs.lasso, s = "lambda.min"))
        
        beta_est[,i*i+(j)+12*(l-1)+120] <- as.vector(coef(probs.lasso, s = "lambda.1se"))
        
      }
        
      # } else if (i == 3) {
      #   
      #   mean.squared.error[i*i+j-1,l] <-  MSE.1
      #   mean.squared.error[i*i+j-1,l+10] <- MSE.2
      #   
      #   beta_true[,i*i+(j-1)+12*(l-1)] <- return$beta
      #   
      #   beta_est[,i*i+(j-1)+12*(l-1)] <- as.vector(coef(probs.lasso, s = "lambda.min"))
      #   
      #   beta_est[,i*i+(j-1)+12*(l-1)+120] <- as.vector(coef(probs.lasso, s = "lambda.1se"))
      #   
     #  }
      
    }
    
  }
  
}



# Compute mean for each row (the first 10 col are lambda_min)
# the last 10 col are lambda_1se
mean.lambda.min <- apply(mean.squared.error[,c(1:10)], 1, mean)

mean.lambda.1se <- apply(mean.squared.error[,c(11:20)], 1, mean)


# Compute respective sd
sd.lambda.min <- apply(mean.squared.error[,c(1:10)], 1, sd)

sd.lambda.1se <- apply(mean.squared.error[,c(11:20)], 1, sd)

# Plot mean and sd for lambda_min and lambda_1se

par(mfrow = c(2,2))

plot(mean.lambda.min,
     ylab = "lambda_min",
     main = "MSE mean lambda_min")
plot(mean.lambda.1se,
     ylab = "lambda_1se",
     main = "MSE mean lambda_1se")

plot(sd.lambda.min,
     ylab = "lambda_min",
     main = "MSE sd lambda_min")
plot(sd.lambda.min,
     ylab = "lambda_1se",
     main = "MSE sd lambda_1se")




# Plot with mean and sd as bars instead:

# Load ggplot2
library(ggplot2)

df_min<-data.frame(Mean=mean.lambda.min,
                   sd=sd.lambda.min)

ggplot(df_min, aes(x = c(1:12), y=Mean)) + geom_point() +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2)


df_1se<-data.frame(Mean=mean.lambda.1se,
                   sd=sd.lambda.1se)

ggplot(df_1se, aes(x = c(1:12), y=Mean)) + geom_point() +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2)


## ---- Feature selection ------------------------------------------------------


# Make your coeff matrices binary


 make.01 <- function(v, param, loop){
   
   x_01 <- matrix(2, nrow = param+1, ncol = 2*loop*12)
   
   for (i in 1:(param+1)) {
     for (j in 1:(2*loop*12)) {
       
       if (round(v[i,j],2) != 0) # < 0 || beta_est[i,j] > 0){
       {
         x_01[i,j]  <- 1
         
       } else if (round(v[i,j],2) == 0){
         x_01[i,j]  <- 0
       }
     }
   }
 return <- x_01
}

beta_01_est <- make.01(beta_est, p, loops )

beta_01_true <- make.01(beta_true, p-1, loops/2)


prop.1.est <- apply(beta_01_est[-1,], 1, mean)

prop.1.true <- apply(beta_01_true, 1, mean)


par(mfrow = c(1,2))

plot(prop.1.est)

plot(prop.1.true)

# seems like there are barely any features being selected 1/5 of times...
# probably something is not correct


## ---- OLD STUFF --------------------------------------------------------------


#xx <- matrix(c(1.5,0,-2,0),c(3,0,2,-2), nrow = 10+1, ncol = 2*4*12)

# beta_01_est <- matrix(2, nrow = p+1, ncol = 2*loops*12)
# 
# for (i in 1:p+1) {
#   for (j in 1:2*loops*12) {
#     
#     if (round(beta_est[i,j],2) != 0.00) # < 0 || beta_est[i,j] > 0){
# {
#       beta_01_est[i,j]  <- 1
# 
#     } else if (round(beta_est[i,j],2) == 0.00){
#       beta_01_est[i,j]  <- 0
# }
#   }
# }





#par(mfrow = c(1,2))

#plot(mean.squared.error[,1],
# ylab = "lambda_min",
# main = "MSE lambda_min")

#plot(mean.squared.error[,2],
# ylab = "lambda_1se",
# main = "MSE lambda_1se")







# 
# X.train <- X[1:n.train,]
# 
# y.train <- y[1:n.train]
# 
# X.test <- X[(n.train+1):n,]
# 
# y.test <- y[(n.train+1):n]
# 
# sd.y <- sqrt(var(y.train)*((n)-1)/n)
# 
# # Just using the methods from 
# # https://glmnet.stanford.edu/articles/glmnet.html
# 
# # Fit lasso regression using glmnet
# fit.lasso.train <- glmnet(X.train,y.train)
# 
# plot(fit.lasso.train)
# 
# print(fit.lasso.train)
# 
# # Obtain model coeff at different lambda ranges within sequence
# 
# coef(fit.lasso.train, s = 0.1)
# 
# # Built-in cv tool:
# 
# # note that X here has been generated with the n not n+n_test
# 
# 
# 
# 
# # 
#  TEST <- matrix(0, nrow = p+1, ncol = 2)
# # 
#  cv.fit.lasso.train <- cv.glmnet(X.train,y.train, 
#                            nfolds = 10,
#                            alpha = 1
#                            )
# # 
# # plot(cv.fit.lasso.train)
# # 
# # # Built-in tool to get lambda_min + coeff 
# # 
#  TEST[,1] <- as.vector(coef(cv.fit.lasso.train, s = "lambda.min"))
#  TEST[,2] <- as.vector(coef(cv.fit.lasso.train, s = "lambda.1se"))

   
   
   
   
   #cv.fit.lasso.train$lambda.min
# 
# coef(cv.fit.lasso.train, s = "lambda.min")
# 
# #predict(cv.fit.lasso.train, newx = X[1:5,],s="lambda.min" )
# 
# #matrix.lambda.min <- as.matrix(coef(cv.fit.lasso, s = "lambda.min"))
# 
# #non.zero.coeff <- matrix.lambda.min>0
# 
# #matrix.lambda.min <- matrix.lambda.min[non.zero.coeff]
# 
# # Built-in tool to get lambda_1se
# 
# cv.fit.lasso.train$lambda.1se
# 
# coef(cv.fit.lasso.train, s = "lambda.1se")
# 
# # Use n_test to calculate MSE 






















# Generate data with test observations instead

#r <- sim.data(n_test,p,sparsity)
# mse.lasso.min <- sum((y.test -
#                         predict(
#                           fit.lasso.train,
#                           X.test,
#                           s = sd.y * cv.fit.lasso.train$lambda.min/n,
#                           exact = T,
#                           x = X.train,
#                           y = y.train))^2)
# 
# mse.lasso.1se <- sum((y.test -
#                         predict(
#                           fit.lasso.train,
#                           X.test,
#                           s = sd.y * cv.fit.lasso.train$lambda.1se/n,
#                           exact = T,
#                           x = X.train,
#                           y = y.train))^2)




















