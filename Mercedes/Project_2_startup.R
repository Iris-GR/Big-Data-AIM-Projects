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
install.packages("MASS")
install.packages("class")
install.packages("FNN")
install.packages("cluster")
install.packages("mclust")


## ---- Loading Packages Into Memory -------------------------------------------
library(readr)
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
library(factoextra)


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


labels <- as_tibble(read.csv("TCGA-PANCAN-HiSeq-801x20531//labels.csv",
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
constants <- vars == 0  # In total 267 constants    

# Data set with constants only:
test <- data[constants]

# Data set without constants:
data <- data[!constants]


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

## ---- PCA using FactoMineR ---------------------------------------------------

install.packages("FactoMineR")

library(FactoMineR)

X <- filtered_hv_data


res.pca <- PCA(X, scale.unit = TRUE, ncp = 5, graph = FALSE)

barplot(res.pca$eig[,1], 
        main = "Eigenvalues")
#names.arg = 1:nrow(res.pca$eig))
plot(res.pca, 
     choix = "ind",
     habillage = 801,
     label = "none",
     xlab = "PC1 (15.46%)",
     ylab = "PC2 (11.34%)")

# Get eigenvalues using factoextra package
# The eigenvalues measure the amount of variation
# retained by each principal component. 
eig.val <- get_eigenvalue(res.pca)
eig.val

# Eigenvalues are large for the first pc and small for the latter (because
# the first components are the directions containing the most of the datasets 
# variation)


# Scree plot using factoextra package
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0,50))

# From this plot we see that we might want to stop at the 10th principal
# component, as more than 50% of variation is explained by these components. 

# Graph of variables
vari <- get_pca_var(res.pca)

# To plot variables:
fviz_pca_var(res.pca, col.var = "black", label = "none")

# Colour by groups
# Create a grouping variable using kmeans
# Trying 5 centers
set.seed(123)
res.km <- kmeans(vari$coord, centers = 5, nstart = 25)
grp <- as.factor(res.km$cluster)

#colour variables by groups (from labels)

# Add labels to data
data_newest <- cbind(labels, data)

# Factorise categorical data (to be able to use col.ind)
C <- as.factor(data_newest$Class)

# Get individuals
ind <- get_pca_ind(res.pca)

# Plot individuals grouped by colour
fviz_pca_ind(res.pca, 
             geom.ind = "point",
             col.ind = C,
             palette = c(cbPalette[6],
                         cbPalette[4],
                         cbPalette[5],
                         cbPalette[7],
                         cbPalette[8]),
             addEllipses = TRUE,
             legend.title = "Groups",
             xlab = "PC1 (15.5%)",
             ylab = "PC2 (11.3%)")

# Without ellipses
fviz_pca_ind(res.pca, 
             geom.ind = "point",
             col.ind = C,
             palette = c(cbPalette[6],
                         cbPalette[4],
                         cbPalette[5],
                         cbPalette[7],
                         cbPalette[8]),
             addEllipses = FALSE,
             legend.title = "Groups",
             xlab = "PC1 (15.5%)",
             ylab = "PC2 (11.3%)")


## ---- Q.1.2 - Clustering -----------------------------------------------------

# We want to determine potential cluster count.
# Can use internal clustering indices to determine quality of clustering results

## ---- Libraries for Clustering -----------------------------------------------
library("NbClust")
library("class")
library("tidyverse")
library("latex2exp")
library("MASS")
library("FNN")
library("cluster")
library("mclust") # Good for GMM clustering/ BIC
library(factoextra)

# NBclust provides 30 indices for determining the number of clusters and 
# proposes to users the best scheme from the different results obtained by
# varying all combinations of number of clusters, distance measures, and 
# clustering methods. 


# mclust is a package for model-based clusterin, classification, and density 
# estimation based on finite normal mixture modelling. 
# It provides functions for parameter estimation via the EM algorithm for normal 
# mixture models with a variety of covariance structures, and functions for 
# simulation of these models. 
# Also included are functions that combine 


# Note: There are internal and external indices 
# Aim: Achieve high between cluster scatter and low within cluster scatter

X <- filtered_hv_data

## ---- TASK 4 ----------------------------------------------------------------- 

## ---- Within Cluster Scatter -----------------------------------------------

# Note that I start by assuming that we should expect 5 clusters (for each 
# cancer type) - however, this is not necessarily true ig

# Simple solution to get k from 
# https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/
# Idea is to compute k-means clustering using different values of clusters k
# Next, wss (within sum of squares) is drawn according to the number of clusters
# The location of a bend (knee) in the plot is generally considered as an 
# indicator of the appropriate number of clusters

# From https://rpkgs.datanovia.com/factoextra/reference/fviz_nbclust.html
# fviz_nbclust() determines and visualise the optimale number of clusters using
# different methods:
# - Within cluster sum of squares
# - Average silhouette
# - gap statistics

# fviz_gap_stat() visulise the gap statistic generated by the function
# clusGap() (in cluster package)
# The optimal number of clusters is specified using the "firstmax" method


# Optimal number of clusters in the data
library(cluster)
## ---- Looking at the knees/ wss ----------------------------------------------
# The lower wss, the better
if (FALSE){
  # kmeans:
  fviz_nbclust(X, kmeans, method = "wss") + 
    geom_vline(xintercept = 3, linetype = 2)
  
  # pam:
  fviz_nbclust(X, cluster::pam, method = "wss") + 
    geom_vline(xintercept = 3, linetype = 2)
  
  # hierarchical clustering
  fviz_nbclust(X, hcut, method = "wss") + 
    geom_vline(xintercept = 3, linetype = 2)
  
}

## ---- Average Clusters Silhouette --------------------------------------------
# The higher the better.
# Average silhouette

if (FALSE){
  # kmeans:
  fviz_nbclust(X, kmeans, method = "silhouette") + 
    geom_vline(xintercept = 3, linetype = 2)
  
  # pam:
  fviz_nbclust(X, pam, method = "silhouette") + 
    geom_vline(xintercept = 3, linetype = 2)
  
  # hierarchical clustering
  fviz_nbclust(X, hcut, method = "silhouette") + 
    geom_vline(xintercept = 3, linetype = 2)
  
}


## ---- Gap statistic ----------------------------------------------------------
# The higher, the better

set.seed(123)
# Compute gap statistic for kmeans
# Recommended B = 500
gap_stat <- clusGap(X, FUN = kmeans, nstart = 25, K.max = 10, B = 10)
print(gap_stat, method = "firstmax")

fviz_gap_stat(gap_stat)

# Gap statistic for hierarchical clustering
gap_stat_hclust <- clusGap(X, FUN = hcut, K.max = 10, B = 10)
fviz_gap_stat(gap_stat_hclust)



## ---- Using NbClust(): 30 indices for choosing the bet number of clusters ----
# data: matrix
# diss: dissimilarity matrix to be used. Default, diss = NULL
# distance: the distance measure to be used to compute the dissimilarity matrix
# e.g. "euclidean", "manhattan", "NULL"
# min.nc, max.nc: mininal and maximal number of clusters, respectively.
# method: the cluster analysis method to be used,
# including "ward.D", "ward.D2", "single", "complete", "average", "kmeans", etc

# NbClust() directly returns the optimal number of clustering based on the 
# frequency distribution histogram.

res.nbclust <- NbClust(X, distance = "euclidean",
                       min.nc = 2, 
                       max.nc = 10,
                       method = "kmeans",
                       index = "all")

# The Hubert index is a graphical method of determining the number of clusters.
# In the plot of the Hubert index, we seek a significant knee that corresponds 
# to a significant increase of the value of the measure, i.e. the 
# significant peak in Hubert index second differences plot.

# The D index is a graphical method of determining the number of clusters. 
# In the plot of the D index, we seek a  significant knee (the significant
# peak in D-index second differences plot) that corresponds to a 
# significant increase of the value of the measure. 

# Visualise the optimal number of clusters distribution
idx <- seq(from = 1, to = 52, by = 2)

cluster_number <- res.nbclust[["Best.nc"]][idx]

#ftable <- data.frame(matrix(, nrow = ))



## ---- K-means ----------------------------------------------------------------

# From http://varianceexplained.org/r/kmeans-free-lunch/
# K-means is a widely used method in cluster analysis.
# This method does NOT require ANY assumptions
# i.e. given a data set and a pre-specified number of clusters,k,
# you can just apply this algorithm which minimise the SSE, 
# SSE = within cluster square error

# K-means is therefore essentially an optimisation problem. 

# Drawbacks: 
# k-means assume:
# 1. The variance of the distribution of each attribute
# (variable) is spherical;

# 2. All variables has the same variancec;

# 3. The prior probability of all k clusters are the same, 
# i.e. each cluster has roughly equal number of observations;

# If any of these 3 assumptions is violated, then it means k-means
# will fail. 

## INTERNET 
# - https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/

# Compute k-means with k = 5
set.seed(123)
km.res <- kmeans(X, 5, nstart = 25)

# nstart = 25 means that R will try 25 different random starting assignments
# and then select the best results corresponding to the one with the lowest 
# within cluster variation.
# Default nstart value is 1, but strongly recommended to compute
# k-means clustering with a large value of nstart such as 25 or 50 
# in order to have more stable results.

# print the results:
print(km.res)

# Displays: the cluster means or centers: 
# a matrix, which rows are cluster numbers and columns are variables
# the clustering vector: A vector of integers (from 1:k) indicating the
# the cluster to which point is allocated. 

# It is possible to compute the mean of each variables by clusters using the 
# original data:
cluster_means <- aggregate(X, by = list(cluster = km.res$cluster),mean)


# If you want to add point classifications to the original data:
dd <- cbind(X, cluster = km.res$cluster)
head(dd)

# Cluster size:
km.res$size

# Cluster means:
km.res$centers


# Visualisation
# Reduce dimensions by PCA and plot data according to the first 
# two principal component coordinates.

# The function fviz_cluster() can be used to easily visualise k-means clusters.
# It takes k-means results and the original data as arguments.
# In the resulting plot, observations are represented by points, using 
# principal components if the number of variables are greater than 2. 

# Also possible to draw concentration ellipse around each cluster. 

# All the fviz_code below is taken from the side:
# https://rpkgs.datanovia.com/factoextra/reference/fviz_cluster.html

# With only points:
fviz_cluster(km.res, X, 
             ellipse.type = "convex",
             geom = "point",
             main = "k-means cluster plot (k = 5)",
             xlab = "PC1 (15.5%)",
             ylab = "PC2 (11.3%)")

# With points and text (default):
fviz_cluster(km.res, X, 
             ellipse.type = "norm",
             geom = c("point","text"),
             main = "k-means cluster plot (k = 5)",
             xlab = "PC1 (15.5%)",
             ylab = "PC2 (11.3%)")

# Use repel = TRUE to avoid overplotting (whatever that is)

# Try on scaled data:
X.scaled <- scale(X)

km.res.scaled <- kmeans(X.scaled,5,nstart=25)

fviz_cluster(km.res.scaled,X.scaled,ellipse.type = "norm")

# No difference - maybe the data was already scaled heh...

#Change colour palette and theme:
fviz_cluster(km.res, X, palette = "Set2", ggtheme = theme_minimal())

if (FALSE) {
  # Show points only
  fviz_cluster(km.res, X, geom = "point")
  # Show text only
  fviz_cluster(km.res, X, geom = "text")
  
  # PAM  clustering
  require(cluster)
  pam.res <- pam(X, 5)
  # Visualise pam clustering
  fviz_cluster(pam.res,geom="point", ellipse.type = "norm")
  
  # Hierarchical clustering
  # Use hcut() which compute hclust and cut the tree
  hc.cut <- hcut(X, k = 5, hc_method = "complete")
  # Visualise dendrogram
  fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)
  
  # Visualise cluster 
  fviz_cluster(hc.cut, 
               ellipse.type = "convex",
               geom = "point",
               main = "Hierarchical clustering (k = 5)", 
               xlab = "PC1 (15.5%)",
               ylab = "PC2 (11.3%)")
}



## ---- Silhouette Width -------------------------------------------------------





## ---- the Davies-Bouldin index -----------------------------------------------

## ---- the Gap  Statistic -----------------------------------------------------

## ---- GMM Clustering ---------------------------------------------------------

# Bayesion Information Criterion using mclust:

X <- filtered_hv_data

## ---- mclust BIC -------------------------------------------------------------
BIC <- mclustBIC(X)

plot(BIC)

summary(BIC)


# Gaussian Mixture model fitted to EM algorithm:
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)

#windows()
#plot(mod1, what = "classification") # Impossible to do with 5000 genes ofc...

## ---- mclust ICL -------------------------------------------------------------
ICL <- mclustICL(X)
summary(ICL)


## ---- Some other





