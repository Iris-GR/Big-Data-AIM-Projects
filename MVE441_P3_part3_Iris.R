# MVE441 Project 3, Question 3, (Iris)

# This question is about selecting features with confidence. The dataset used 
# is the gene expression dataset.

# Load libraries 
library(tibble)
library(glmnet)
library(tidyverse)
library(ggplot2)

# Simulate data ----------------------------------------------------------------
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
simulate_data <- function(n, p, sparsity, SNR = 3, beta_scale = 5) {
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

##############

# Number of datasets simulated 
nr.ds <- 100

lambda.names <- c("L.small", "L.1se", "L.big")
stab.cut.off <- c(0.5, 0.7, 0.9)

# Create a data frame to gather all the results 
# # (nrow = nr lambda * nr cutoffs * nr data simulations)
# sel.acc <- as.data.frame(matrix(nrow = length(lambda.names)*length(stab.cut.off)*nr.ds, 
#                                 ncol = 4))
# (nrow = nr lambda * nr cutoffs * nr data simulations)
sel.acc <- as.data.frame(matrix(ncol = 4))
colnames(sel.acc) <- c("sensitivity", "speceficity", "lambda", "stab.cut.off")

###############

for (jj in 1:nr.ds) {
  # Generate the data X, y, and beta for the true features -----------------------
  # Number of observations (rows) 
  n <- 750
  
  # Number of features (columns)
  p <- 500
  
  # The sparsity of true features
  sparsity <- 0.8
  
  # Call data simulation function
  data <- simulate_data(n, p, sparsity)
  
  # X-data
  X.data <- data.matrix(data$X)
  colnames(X.data) <- paste(rep("P", p), c(1:p), sep = "_")
  rownames(X.data) <- paste(rep("n", n), c(1:n), sep = "_")
  
  # y-data
  y.data <- data.matrix(data$y)
  rownames(y.data) <- paste(rep("n", n), c(1:n), sep = "_")
  
  # True betas
  b.true <- as.data.frame(data$beta)
  rownames(b.true) <- paste(rep("b_true_P", p), c(1:p), sep = "_")
  
  # Collect the indexes of the true betas
  true.feat.ind <- which(b.true != 0)
  
  # Subset the true beta values for the true features
  true.feat.bval <- subset(b.true, data$beta != 0)
  
  # Binary notation of true features
  b.true.bin <- b.true
  b.true.bin[b.true.bin != 0] <- 1
  
  # Prepare data frame for collection of feature selection frequencies -----------
  sel.freq <- as.data.frame(matrix(0, nrow = p + 1, ncol = 3))
  rownames(sel.freq) <- c("(Intercept)", colnames(X.data))
  colnames(sel.freq) <- c("L.small", "L.1se", "L.big")
  sel.freq
  
  #* Boot strap ----------------------------------------------------------------
  
  # Number of repeated bootstraps (M=100 is rather time consuming)
  M = 100
  
  # Factor that lambda.1se is changed by [Look into an appropriate number!!!]
  L.fac <- 5
  
  # For each bootstrap
  for (j in 1:M) {
    # Take a sample with replacement from the X data matrix and the same
    # from the y vector
    boot.rows <- sample(1:nrow(X.data), nrow(X.data), replace = TRUE)
    X.boot <- X.data[boot.rows, ]
    y.boot <- y.data[boot.rows, ]
    
    # 10-fold CV of lasso
    cv.boot <- cv.glmnet(X.boot,
                         y.boot,
                         type.measure = "mse",
                         alpha = 1)
    
    # Extract the feature coefficients for different lambda
    pred.coef <- predict(
      cv.boot,
      newx = X.boot,
      s = c(
        cv.boot$lambda.1se / L.fac,
        cv.boot$lambda.1se,
        cv.boot$lambda.1se * L.fac
      ),
      type = "coef"
    )
    colnames(pred.coef) <- c("L.small", "L.1se", "L.big")
    
    # Count frequency of features that were selected
    for (i in 1:length(pred.coef[, 1])) {
      if (pred.coef[i, 1] != 0) {
        sel.freq$L.small[i] = sel.freq$L.small[i] + 1
      }
      
      if (pred.coef[i, 2] != 0) {
        sel.freq$L.1se[i] = sel.freq$L.1se[i] + 1
      }
      
      if (pred.coef[i, 3] != 0) {
        sel.freq$L.big[i] = sel.freq$L.big[i] + 1
      }
    }
  }
  
  # Total counts of feature selection
  sel.freq
  
  # Feature selection stability across bootstrap samples
  sel.stab <- sel.freq / M
  
  # How many and which features are left after cutoff at different selection
  # stabilites 70, 80, 90%?
  
  # Remove intercept from sel.stab and reset the row indexes
  sel.stab.feat <- sel.stab[-1, ]
  
  # Calculate sensitivity and specificity --------------------------------------
  
  # Condition positive (CP), the number of real positive cases in the data
  CP <- length(true.feat.ind)
  
  # Condition negative (CN), the number of real negative cases in the data
  CN <- p-length(true.feat.ind)
  
  for (ii in 1:length(lambda.names)) {
    
    for (kk in 1:length(stab.cut.off)) {
      
      # The indexes of the selected features that are left after stability cutoff
      feat.ind <- which(sel.stab.feat[ii] > stab.cut.off[kk])
      
      # True positives. Number of by selected features (estimates) that are also 
      # true features
      TP <- length(Reduce(intersect, list(true.feat.ind, feat.ind)))
      
      # False positive (FP). A test result which wrongly indicates that a 
      # particular condition or attribute is present
      if (length(feat.ind) > CP) {
        FP <- length(feat.ind) - CP
      } else {
        FP <- 0
      }
      
      # Calculate sensitivity
      sens <- TP/CP
      
      # Calculate specificity
      spec <- 1 - FP/CN
      
      # Insert into sel.acc matrix
      newrow <- c(round(sens, digits = 3), 
                  round(spec, digits = 3), 
                  lambda.names[ii],  
                  stab.cut.off[kk])
      sel.acc <- rbind(sel.acc, newrow) # Not good with rbind
      
    }
  }
}

# Remove starting row (not nice)
sel.acc <- sel.acc[-1,]
rownames(sel.acc) <- NULL
sel.acc

result <- sel.acc


# Make the variables into factors
result$lambda <- factor(result$lambda)
result$speceficity <- as.numeric(result$speceficity)
result$sensitivity <- as.numeric(result$sensitivity)
result$stab.cut.off<- as.numeric(result$stab.cut.off)

result
# Look at pairs plot
plot(result$stab.cut.off, result$sensitivity)
plot(result$stab.cut.off, result$speceficity)
plot(result$lambda, result$sensitivity)
plot(result$lambda, result$speceficity)
plot(result$speceficity, result$sensitivity)
plot(result$sensitivity, result$speceficity)

# Fit model to specificity (fail)
fit.spec.gussian <- glm(speceficity ~ lambda+stab.cut.off,
                        data = result,
                        family = 'gaussian')

summary(fit.spec.gussian)
plot(fit.spec.gussian)

# Fit model to sensitivity (fail)
fit.sens.gussian <- glm(sensitivity ~ lambda+stab.cut.off,
                        data = result,
                        family = 'gaussian')
summary(fit.sens.gussian)
plot(fit.sens.gussian)















