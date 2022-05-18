library(glmnet)
library(ggplot2)
library(MASS)
library(base)
library(viridis)
library(cowplot)

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
      #Note - 1 because of error removal thing in definiing beta_min_logic etc.
      sens_min = tb_min[2,2]/(tb_min[2,2] + tb_min[1,2] - 1)
      spec_min = tb_min[1,1]/(tb_min[1,1] + tb_min[2,1] - 1)
      
      sens_1se = tb_1se[2,2]/(tb_min[2,2] + tb_min[1,2] - 1)
      spec_1se = tb_1se[1,1]/(tb_min[1,1] + tb_min[2,1] - 1)
      
      #Save metrics into global array for later use
      sens_spec_min[i,j,k,] = c(sens_min,spec_min)
      sens_spec_1se[i,j,k,] = c(sens_1se,spec_1se)
    }  
  }
}

#---------------------- Plotting results of MSE-calculations -----------------

MSE_min_plotdata = NULL
MSE_1se_plotdata = NULL
MSE_difference_plotdata = NULL

MSE_min_se = NULL #Standard error
MSE_1se_se = NULL

for(i in 1:3){
  for(j in 1:4){
    MSE_min_plotdata = rbind(MSE_min_plotdata,c(n[i],sparsity[j],mean(MSE_min[i,j,])))
    MSE_1se_plotdata = rbind(MSE_1se_plotdata,c(n[i],sparsity[j],mean(MSE_1se[i,j,])))
    MSE_difference_plotdata = rbind(MSE_difference_plotdata,c(n[i],sparsity[j],mean(MSE_1se[i,j,])-mean(MSE_min[i,j,])))
    
    #Calculate standard errors
    MSE_min_se = rbind(MSE_min_se,sqrt(var(MSE_min[i,j,])/5))
    MSE_1se_se = rbind(MSE_1se_se,sqrt(var(MSE_1se[i,j,])/5))
  }
}


MSE_min_plotdata = as.data.frame(MSE_min_plotdata)
MSE_1se_plotdata = as.data.frame(MSE_1se_plotdata)
MSE_difference_plotdata = as.data.frame(MSE_difference_plotdata)

#Plotting each MSE for lambda ___min___:
limits <- aes(ymax = V3 + MSE_min_se,
              ymin = V3 - MSE_min_se)

p_min <- ggplot(data = MSE_min_plotdata, aes(x = factor(V2), y = V3,
                               fill = factor(V1)))
p_min = p_min + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "Sparsity", y = "MSE") +
  ggtitle("MSE using lambda.min") +
  scale_fill_discrete(name = "No. samples")

#Plotting each MSE for lambda ___1se___:
limits <- aes(ymax = V3 + MSE_1se_se,
              ymin = V3 - MSE_1se_se)

p_1se <- ggplot(data = MSE_1se_plotdata, aes(x = factor(V2), y = V3,
                                         fill = factor(V1)))
p_1se = p_1se + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "Sparsity", y = "MSE") +
  ggtitle("MSE using lambda.1se") +
  scale_fill_discrete(name = "No. samples")

plot_grid(p_min,p_1se)

#Plotting the difference between MSE_1se and MSE_min
#NOTE: No error bars included.
#Conclusions:
# - Difference is positive in all cases
# - Difference is greater the smaller the sparsity
# - Difference is greater the higher n/p
names_sparsity = c('0.75','0.90','0.95','0.99')
windows()
barplot(V3 ~., MSE_difference_plotdata,
        beside = TRUE,
        xlab = "Sparsity",
        ylab = "Difference in MSE",
        names.arg = names_sparsity,
        legend.text = c("n = 200", "n = 700", "n = 750"))

# --------------- Plotting results from feature selection calculations --------

#Fattar inte vad som ska plottas i en scatterplot...
#Gör samma sak som ovan med en barplot? En för min och en för 1se?

# SETUP FOR PLOTTING:
# windows()
# par(mfcol=c(2,2))

# ----- For lambda.min

sens_min_plotdata = NULL
spec_min_plotdata = NULL

sens_min_se = NULL
spec_min_se = NULL

for(i in 1:3){
  for(j in 1:4){
    sens_min_plotdata = rbind(sens_min_plotdata,
                                   c(n[i],
                                     sparsity[j],
                                     mean(sens_spec_min[i,j,,1])))
    spec_min_plotdata = rbind(spec_min_plotdata,
                              c(n[i],
                                sparsity[j],
                                mean(sens_spec_min[i,j,,2])))
    
    sens_min_se = rbind(sens_min_se,sqrt(var(sens_spec_min[i,j,,1])/5))
    spec_min_se = rbind(spec_min_se,sqrt(var(sens_spec_min[i,j,,2])/5))
  }
}

sens_min_plotdata = as.data.frame(sens_min_plotdata)
spec_min_plotdata = as.data.frame(spec_min_plotdata)

#Plotting sens for lambda min:
limits <- aes(ymax = V3 + sens_min_se,
              ymin = V3 - sens_min_se)

p1 <- ggplot(data = sens_min_plotdata, aes(x = factor(V2), y = V3,
                                         fill = factor(V1)))
p1 = p1 + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "Sparsity", y = "Sensitivity") +
  ggtitle("Sensitivity using lambda.min") +
  scale_fill_discrete(name = "No. samples")

#Plotting spec for lambda min:
limits <- aes(ymax = V3 + spec_min_se,
              ymin = V3 - spec_min_se)

p2 <- ggplot(data = spec_min_plotdata, aes(x = factor(V2), y = V3,
                                          fill = factor(V1)))
p2 = p2 + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "Sparsity", y = "Specificity") +
  ggtitle("Specificity using lambda.min") +
  scale_fill_discrete(name = "No. samples")

# barplot(V3 ~., sens_min_plotdata,beside = TRUE,
#         xlab = "Sparsity",
#         ylab = "Sensitivity (lambda_min)",
#         names.arg = names_sparsity,
#         legend.text = c("n = 200", "n = 700", "n = 750"),
#         args.legend = list(x = "top",
#                            inset = c(0, -0.25)))
# barplot(V3 ~., spec_min_plotdata,beside = TRUE,
#         xlab = "Sparsity",
#         ylab = "Specificity (lambda_min)",
#         names.arg = names_sparsity,
#         legend.text = c("n = 200", "n = 700", "n = 750"),
#         args.legend = list(x = "top",
#                            inset = c(0, -0.25)))

# ------ For lambda.1se

sens_1se_plotdata = NULL
spec_1se_plotdata = NULL

sens_1se_se = NULL
spec_1se_se = NULL

for(i in 1:3){
  for(j in 1:4){
    sens_1se_plotdata = rbind(sens_1se_plotdata,
                              c(n[i],
                                sparsity[j],
                                mean(sens_spec_1se[i,j,,1])))
    spec_1se_plotdata = rbind(spec_1se_plotdata,
                              c(n[i],
                                sparsity[j],
                                mean(sens_spec_1se[i,j,,2])))
    
    sens_1se_se = rbind(sens_1se_se,sqrt(var(sens_spec_1se[i,j,,1])/5))
    spec_1se_se = rbind(spec_1se_se,sqrt(var(sens_spec_1se[i,j,,2])/5))
  }
}

sens_1se_plotdata = as.data.frame(sens_1se_plotdata)
spec_1se_plotdata = as.data.frame(spec_1se_plotdata)

#Plotting sens for lambda 1se:
limits <- aes(ymax = V3 + sens_1se_se,
              ymin = V3 - sens_1se_se)

p3 <- ggplot(data = sens_1se_plotdata, aes(x = factor(V2), y = V3,
                                          fill = factor(V1)))
p3 = p3 + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "Sparsity", y = "Sensitivity") +
  ggtitle("Sensitivity using lambda.1se") +
  scale_fill_discrete(name = "No. samples")

#Plotting spec for lambda 1se:
limits <- aes(ymax = V3 + spec_1se_se,
              ymin = V3 - spec_1se_se)

p4 <- ggplot(data = spec_1se_plotdata, aes(x = factor(V2), y = V3,
                                          fill = factor(V1)))
p4 = p4 + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "Sparsity", y = "Specificity") +
  ggtitle("Specificity using lambda.1se") +
  scale_fill_discrete(name = "No. samples")

# barplot(V3 ~., sens_1se_plotdata,beside = TRUE,
#         xlab = "Sparsity",
#         ylab = "Sensitivity (lambda_1se)",
#         names.arg = names_sparsity,
#         legend.text = c("n = 200", "n = 700", "n = 750"),
#         args.legend = list(x = "top",
#                            inset = c(0, -0.25)))
# barplot(V3 ~., spec_1se_plotdata,beside = TRUE,
#         xlab = "Sparsity",
#         ylab = "Specificity (lambda_1se)",
#         names.arg = names_sparsity,
#         legend.text = c("n = 200", "n = 700", "n = 750"),
#         args.legend = list(x = "top",
#                            inset = c(0, -0.25)))

#plotting all in same pic
plot_grid(p1,p3,p2,p4)


# ----- For difference (lambda.1se - lambda.min)

sens_diff_plotdata = NULL
spec_diff_plotdata = NULL

for(i in 1:3){
  for(j in 1:4){
    sens_diff_plotdata = rbind(sens_diff_plotdata,
                              c(n[i],
                                sparsity[j],
                                mean(sens_spec_1se[i,j,,1]) - mean(sens_spec_min[i,j,,1])))
    spec_diff_plotdata = rbind(spec_diff_plotdata,
                              c(n[i],
                                sparsity[j],
                                mean(sens_spec_1se[i,j,,2]) - mean(sens_spec_min[i,j,,2])))
  }
}

sens_diff_plotdata = as.data.frame(sens_diff_plotdata)
spec_diff_plotdata = as.data.frame(spec_diff_plotdata)

windows()

par(mfrow = c(2,1))

barplot(V3 ~., sens_diff_plotdata,beside = TRUE,
        xlab = "Sparsity",
        ylab = "Sensitivity (1se minus min)",
        names.arg = names_sparsity,
        legend.text = c("n = 200", "n = 700", "n = 750"),
        args.legend = list(x = "top",
                           inset = c(0, -0.55)))

barplot(V3 ~., spec_diff_plotdata,beside = TRUE,
        xlab = "Sparsity",
        ylab = "Specificity (1se minus min)",
        names.arg = names_sparsity,
        legend.text = c("n = 200", "n = 700", "n = 750"),
        args.legend = list(x = "top",
                           inset = c(0, -0.55)))

#Conclusions:
#lambda.1se model is worse at including features that were actually in the model
#But it is better att not including feature that were not in the model
#However, the difference in sensitivity is lower than the difference in specificity
#And the difference in sensitivity becomes smaller with a less sparse data set






################################################################################
###                               THE END                                    ###
################################################################################