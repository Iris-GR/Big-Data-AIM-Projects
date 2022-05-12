library(MASS)
library(car)
library(e1071); library(ggplot2)
library(tidyverse)
library(base)
library(matlib)
library(ggfortify)
library(cowplot)

library(fpc)
library(NbClust)
library(mclust)
library("FactoMineR")
library(factoextra)
library(cluster)
library(ClusterR)
library(ggpubr)
library(cclust)
library(clusterSim)
library(Polychrome)

#Read data
df = readr::read_csv("TCGA-PANCAN-HiSeq-801x20531//data.csv")
df = df[,-1] #Remove first column containing row names
df_labels = readr::read_csv("TCGA-PANCAN-HiSeq-801x20531//labels.csv")
df_labels = df_labels[,-1]


#Variance filtering for 5000 features
s = apply(df,2,sd) #Calculate standard deviation of each column (feature)
ix_filtered = order(s, decreasing =TRUE)
ix_filtered = ix_filtered[1:5000]

df_filtered = df[,ix_filtered]

#PCA on variance filtered data, and save first 10 PC's
df_filtered.pca = prcomp(df_filtered[,-ncol(df_filtered)], scale. = TRUE)
df_filtered_PCA = as.data.frame(df_filtered.pca$x[,1:10])


#------------------Performing hierarchial clustering--------------------------
#Can change some values like "euclidean", "average" and K.

#Distance matrix
distance_mat = dist(df_filtered_PCA, method = "euclidean")
#NOTE! Kanske ska vara squared?? <-------------------------------------???

#Fit hierarchial clustering model to data
Hierar_cl <- hclust(distance_mat, method = "single")

#Plot dendogram
plot(Hierar_cl)


#Number of clusters
K = 5

# Cutting tree by K number of clusters
fit <- cutree(Hierar_cl, k = K)

#Plot clusters:
ind.coord <- as.data.frame(get_pca_ind(df_filtered.pca)$coord)
ind.coord$cluster <- factor(fit)
ind.coord$Cancer <- df_labels$Class

#Plot for many clusters, becomes gray
# P50 = createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
# 
# ggscatter(
#   ind.coord, x = "Dim.1", y = "Dim.2",
#   color = "cluster", palette = P50, ellipse = TRUE, ellipse.type = "convex",
#   shape = "Cancer", size = 1.5
# ) +
#   xlab("PC1")+
#   ylab("PC2")

#Plot for few number of clusters
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2",
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Cancer", size = 1.5,  legend = "right", ggtheme = theme_bw()
) +
  stat_mean(aes(color = cluster), size = K)+
  xlab("PC1")+
  ylab("PC2")

#------------Calculating indices for hierarchial clustering--------------------
#This gives hierarchial clustering with "euclidean distance" and "single linkage"
#Calculates indices for clusters from 2 to max_clusters.
#Can change "euclidean" and "single" in NbClust for different versions of
#hierarchial clustering. Default is "eculidean" and "single" with
#max_clusters = 100


max_clusters = 100
index_silhouette = NbClust(data = df_filtered_PCA,
                           diss= NULL,
                           distance = "euclidean",
                           min.nc=2,
                           max.nc=max_clusters,
                           method = "single",
                           index = c("silhouette"),
                           alphaBeale = 0.1)

index_db = NbClust(data = df_filtered_PCA,
                   diss= NULL,
                   distance = "euclidean",
                   min.nc=2,
                   max.nc=max_clusters,
                   method = "single",
                   index = c("db"),
                   alphaBeale = 0.1)

index_WSS = NbClust(data = df_filtered_PCA,
                    diss= NULL,
                    distance = "euclidean",
                    min.nc=2,
                    max.nc=max_clusters,
                    method = "single",
                    index = c("ball"),
                    alphaBeale = 0.1)

index_silhouette = index_silhouette$All.index
index_db = index_db$All.index
index_WSS = index_WSS$All.index

n_clusters = 2:max_clusters
d = cbind(n_clusters, index_silhouette, index_db, index_WSS)

#Plotting the indices as functions of number of clusters
ggplot(d, aes(x = n_clusters)) +
  geom_point(aes(y = index_silhouette), color = "#426cf5") +
  geom_line(aes(y = index_silhouette), color = "#426cf5") +
  xlab("Number of clusters") +
  ylab("Average silhouette width")

ggplot(d, aes(x = n_clusters)) +
  geom_point(aes(y = index_db), color = "#426cf5") +
  geom_line(aes(y = index_db), color = "#426cf5") +
  xlab("Number of clusters") +
  ylab("Davis-Bouldin index")

ggplot(d, aes(x = n_clusters)) +
  geom_point(aes(y = index_WSS), color = "#426cf5") +
  geom_line(aes(y = index_WSS), color = "#426cf5") +
  xlab("Number of clusters") +
  ylab("Within cluster sum of squares")


