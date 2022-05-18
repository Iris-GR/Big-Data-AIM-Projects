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
n = c(200, 500, 700, 900)
p = 1000
sparsity = c(0.75, 0.9, 0.95, 0.99)
SNR = c(2,3,4,5)
beta_scale = 5

n_test = 1000


#----------------- simulation part --------------------
i = 1
j = 1
k = 1

#Keep 2 features constant, and vary the third
#1. 1-4 Vary n
#2. 5-8 vary sparsity
#3. 9-12 vary SNR

n = c(n,rep(500,8))
sparsity = c(rep(0.9,4),sparsity,rep(0.9,4))
SNR = c(rep(3,8),SNR)

beta_counter = matrix(0,1000,12)
beta_true = matrix(0,1000,12)

for(i in 1:12){
  
  #Simulate data and split into training and testing sets
  df = simulate_data(n[i],p,sparsity[i],SNR[i],beta_scale)
  
  #1000 features bootstrap samples
  beta_sum = matrix(0,1000,1)
  
  #Repeat M = 50 times:
  for(m in 1:50){
    
    #Create bootstrap sample of filtered data (X_filtered) with 801 samples in each
    ix = sample(1:n[i], n[i], replace = TRUE)
    X_m = df$X[ix,]
    y_m = df$y[ix]
    
    #Perform feature selection, i.e. determine which coeffs in fit that are non-zero
    #10 folds is default
    fit_m <- cv.glmnet(X_m, y_m)
    
    #Record in an array of length(beta) the accumulated counts thus far
    beta_m = coef(fit_m, s="lambda.1se")
    beta_m = as.matrix(beta_m)
    beta_m = beta_m[-1,] #Remove intercept
    beta_m = beta_m != 0 #Make logical
    
    #Behöver reda ut vilka features som ska få en count
    #Kolla i beta_m på namnen på raderna, blir det samma ordning varje gång?
    #Isf, bara gör en count vector och summera succesivt för varje run
    beta_sum = beta_sum + beta_m
  }
  
  beta_counter[,i] = beta_sum
  beta_true[,i] = df$beta
}


for(i in 1:12){
  
  sf = which(beta_counter[,i]/50 > 0.8)
  beta_e = rep(0,1000)
  beta_e[sf] = 1
  
  beta_logic = beta_true[,i] != 0
  
  #Confusion matrix
  tb = table(beta_e, beta_logic)
  
  #Calculations of sensitivity and specificity
  #Note - 1 because of error removal thing in definiing beta_min_logic etc.
  sens = tb[2,2]/(tb[2,2] + tb[1,2])
  spec = tb[1,1]/(tb[1,1] + tb[2,1])
  
  cat("sens = ", sens, "n = ", n[i], "spar = ", sparsity[i], "SNR = ",SNR[i],"\n")
  cat("spec =", spec, "n = ", n[i], "spar = ", sparsity[i], "SNR = ", SNR[i],"\n")
  cat("\n")
}


#---------------Plotting histograms for feature selection:--------------------

#------ Varying n
for(i in 1:4){
  windows()
  
  counter_plot = as.data.frame(beta_counter[,i])
  counter_plot = counter_plot/50
  counter_plot$feature = row.names(counter_plot)
  colnames(counter_plot) <- c('V1','feature')
  
  print(ggplot(counter_plot) +
        geom_bar( aes(x=reorder(feature,-V1), y=V1), stat="identity", fill="skyblue", alpha=0.7)+
        coord_cartesian(xlim = c(1, 100)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )

}

#------- Varying sparsity
for(i in 5:8){
  windows()
  
  counter_plot = as.data.frame(beta_counter[,i])
  counter_plot = counter_plot/50
  counter_plot$feature = row.names(counter_plot)
  colnames(counter_plot) <- c('V1','feature')
  
  print(ggplot(counter_plot) +
          geom_bar( aes(x=reorder(feature,-V1), y=V1), stat="identity", fill="skyblue", alpha=0.7)+
          coord_cartesian(xlim = c(1, 100)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )
  
}

#------- Varying signal strength
for(i in 9:12){
  windows()
  
  counter_plot = as.data.frame(beta_counter[,i])
  counter_plot = counter_plot/50
  counter_plot$feature = row.names(counter_plot)
  colnames(counter_plot) <- c('V1','feature')
  
  print(ggplot(counter_plot) +
          geom_bar( aes(x=reorder(feature,-V1), y=V1), stat="identity", fill="skyblue", alpha=0.7)+
          coord_cartesian(xlim = c(1, 100)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )
  
}
