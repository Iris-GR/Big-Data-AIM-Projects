# This is my first attempt at looking at the dataset for p2

## ---- New Packages -----------------------------------------------------------

install.packages("data.table")
install.packages("irlba")
install.packages("plot3D")
install.packages("rgl")
install.packages("purr")
install.packages("tidyverse")
install.packages("latex2exp")
install.packages("kableExtre")
install.packages("NbClust")



## ---- Loading Packages Into Memory -------------------------------------------
library("readr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("stringr")
library("purr")
library("forcats")
library("grDevices")
library("data.table")
library("tidyverse")
library("latex2exp")
library("irlba") 
library("plot3D")
library("rgl")
library("kableExtra")
library("NbClust")


## ---- Colour Palette ---------------------------------------------------------
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour-blind friendly palette

## ---- Data Set ---------------------------------------------------------------

# Gene expression cancer RNA-Seq Data Set
# https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq

# loading dataset 
data <- as_tibble(read.csv("TCGA-PANCAN-HiSeq-801x20531//data.csv",
                           header = TRUE,
                           row.names = 1))


## ---- Q.1.1 - Overview -------------------------------------------------------

## ---- TASK 1 -----------------------------------------------------------------
## ---- Overview of Data Set ---------------------------------------------------


head(data)                  # Overview

summary(data)               # Summary 

not_missing <- is.na(data)  # Note, this does not work as missing values
# are not coded as NA.

table(data == 0)            # Potentially, missing values coded as 0, 
# if that is the case, there are 2338124 missing

means <- apply(data, 2, mean) # Calculate means of features, note that some 
# means are 0, indicating constant variables

vars <- apply(data, 2, var)   # Calculate variance features, same note as above

hist(means, 
     main = paste("Histogram of feature means")) # Histogram of means
hist(vars,
     main = paste("Histogram of feature variances"),
     xlab = "Variance")                          # Histogram of variances

# 0 variance = constant values
constants <- which(vars == 0)  # In total 267 constants                   

## ---- TASK 2 -----------------------------------------------------------------
## ---- High Variance Filtering ------------------------------------------------

HVF <- which(vars>= 2.237) # Filtered out features with low variance
# leaving exactly 5000 features

#new_data <- rbind(vars,data)
#print("Original Dataframe:")
#print(new_data)

#rownames <- rownames(new_data)
#print("Original Rownames:")
#print(rownames)

#rownames(new_data) <- c("Variance", 1:801)


filtered_hv_data <- data %>%
  select(where(~ var(.) >= 2.237))           # Create new tibble with hvf 

means_hv <- apply(filtered_hv_data, 2, mean) # Calculate means of features, note that some 
# means are 0, indicating constant variables

vars_hv <- apply(filtered_hv_data, 2, var)   # Calculate variance features, same note as above

hist(means_hv, 
     main = paste("Histogram of high-variance feature means"),
     xlab = "means")                         # Histogram of means
hist(vars_hv,
     main = paste("Histogram of high-variance feature variances"),
     xlab = "Variance")                      # Histogram of variances

## ---- TASK 3 -----------------------------------------------------------------
## ---- PCA --------------------------------------------------------------------

# Calculate PCA for raw data and standardised data (using code from canvas L5)
# Data is centred so that SVD can be utilised

# Compute singular-value decomposition for the raw data (not scaled)
X_svd_raw <- svd(scale(as.matrix(filtered_hv_data), scale = FALSE))

# Compute singular-value decomposition for the standardised data (scaled)
X_svd_std <- svd(scale(as.matrix(filtered_hv_data)))

# Compute the principal components for the raw data
X_proj_raw <- scale(as.matrix(filtered_hv_data), scale = FALSE) %*%
  X_svd_raw$v


# Compute the principal components for the standardised data
X_proj_std <- scale(as.matrix(filtered_hv_data)) %*% X_svd_std$v

# Create tibble with first three pc for raw data
filtered_hv_proj_raw <- tibble(
  PC1 = X_proj_raw[,1],
  PC2 = X_proj_raw[,2],
  PC3 = X_proj_raw[,3])

# Create tibble with first three pc for standardised data
filtered_hv_proj_std <- tibble(
  PC1 = X_proj_std[,1],
  PC2 = X_proj_std[,2],
  PC3 = X_proj_std[,3])

# Plot pc1 vs pc2 raw
p3 <- ggplot(filtered_hv_proj_raw) +
  geom_point(aes(x = PC1, y = PC2), size = 0.6) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)) +
  ggtitle("Raw data pair plot of PC1 vs PC2")

# Plot pc1 vs pc2 standardised
p4 <- ggplot(filtered_hv_proj_std) +
  geom_point(aes(x = PC1, y = PC2), size = 0.6) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)) +
  ggtitle("Standardised data pair plot of PC1 vs PC2")

# Plot pc1 vs pc3 raw
p5 <- ggplot(filtered_hv_proj_raw) +
  geom_point(aes(x = PC1, y = PC3), size = 0.6) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)) +
  ggtitle("Raw data pair plot of PC1 vs PC3")

# Plot pc1 vs pc3 standardised
p6 <- ggplot(filtered_hv_proj_std) +
  geom_point(aes(x = PC1, y = PC3), size = 0.6) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)) +
  ggtitle("Standardised data pair plot of PC1 vs PC3") 


# Arrange multiple ggplots on the same page:
ggpubr::ggarrange(
  p3, p4, p5, p6, ncol = 2, nrow = 2, heights = c(1.1, 1))

# Plot pc3 vs pc2 raw
p7 <- ggplot(filtered_hv_proj_raw) +
  geom_point(aes(x = PC2, y = PC3), size = 0.6) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)) +
  ggtitle("Raw data pair plot of PC2 vs PC3")

# Plot pc3 vs pc2 raw
p8 <- ggplot(filtered_hv_proj_std) +
  geom_point(aes(x = PC2, y = PC3), size = 0.6) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)) +
  ggtitle("Standardised data pair plot of PC2 vs PC3") 


# Arrange multiple ggplots on the same page:
ggpubr::ggarrange(
  p7, p8, ncol = 2, nrow = 2, heights = c(1.1, 1))


## ---- PCA From the Internet --------------------------------------------------

# Perform PCA
results <- prcomp(filtered_hv_data, scale = TRUE)

# Calculate total variance explained by each pc
var_explained = results$sdev^2 / sum(results$sdev^2)

# Note that sdev are the standard deviations of the principal components
# i.e. the square roots of the eigenvalues of the covariance/correlation matrix


# Create scree plot
qplot(c(1:10), var_explained[1:10]) +
  geom_line() +
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0,1)

# Note, qplot stands for quick plot. It is a convenient wrapper for creating a 
# number of different types of plots using a convenient calling scheme.
# It is a quick fix for the more complex ggplot.


## ---- Q.1.2 - Clustering -----------------------------------------------------

# We want to determine potential cluster count.
# Can use internal clustering indices to determine quality of clustering results







