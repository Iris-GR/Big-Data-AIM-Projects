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

#-------------------------------------------------------------------------
##########################################################################
#---------------------------PROJECT 2-------------------------------------

df = readr::read_csv("TCGA-PANCAN-HiSeq-801x20531//data.csv")
df = df[,-1] #Remove first column containing row names


#How are missing values encoded?
#-Not a clue

#Calculate mean and sd of the features

#------Can try plotting density curves for different genes------
# Kernel density estimation, using gaussian kernel
d <- density(df$gene_1)

# Kernel density plot
plot(d, lwd = 2, main = "Default kernel density plot")


#-----TASK1-----
#Beräkna mean av varje feature?
m = NULL
s = NULL
m = apply(df,2,mean)
s = apply(df,2,sd)
hist(m,50)
d <- density(m)
plot(d, lwd = 2, main = "Default kernel density plot") #"does not vary so much"
hist(s,50)
d <- density(s)
plot(d, lwd = 2, main = "Default kernel density plot")
#Some constant features? Since variance \approx 0.

#-----TASK2-----
#Variance filtering to obtain 5000 features with most variance.
ix_filtered = order(s, decreasing =TRUE)
ix_filtered = ix_filtered[1:5000]

m_filtered = m[ix_filtered]
s_filtered = s[ix_filtered]
hist(m_filtered,50)
d <- density(m_filtered)
plot(d, lwd = 2, main = "Default kernel density plot")
hist(s_filtered,50)
d <- density(s_filtered)
plot(d, lwd = 2, main = "Default kernel density plot")

#Feature means are now quite evenly spread around 5

#Most feature variances are between 1 and 3.

#Need for center and standardize data?
#-Probably less need than before. But it can depend on what method
# one is going to apply? I think it would be best to try both methods
# and compare the results.

#-----TASK3-----
#Perform PCA on the variance filtered dataset (using centering and sclaing)
df_filtered = df[,ix_filtered]
df_filtered.pca = prcomp(df_filtered[,-ncol(df_filtered)], scale. = TRUE)

df_filtered.pca.plot <- autoplot(df_filtered.pca,
                          data = df_filtered,
                          x = 1,
                          y = 2)

df_filtered.pca.plot

#Plotta flera i ett grid eller skita i det?

#calculate total variance explained by each principal component
var_explained = df_filtered.pca$sdev^2 / sum(df_filtered.pca$sdev^2)

n_PC = 75

#plot first 10 points in Scree-plot
qplot(c(1:n_PC), var_explained[1:n_PC]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.25)

cum_var_explained <- NULL
cum_var_explained[1] = var_explained[1]
for(i in 2:n_PC){
  cum_var_explained[i] = cum_var_explained[i-1]+var_explained[i]
}

qplot(c(1:n_PC), cum_var_explained[1:n_PC]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Cumulative Variance Explained") +
  ggtitle("Cumulative Variance Explained Plot") +
  ylim(0, 1)

x = 1:n_PC
y = var_explained[1:n_PC]
z = cum_var_explained[1:n_PC]
d = cbind(x,y,z)
ggplot(d, aes(x = x)) +
  geom_point(aes(y = y), color = "#338BFF") +
  geom_line(aes(y = y), color = "#338BFF") +
  geom_point(aes(y = z), color = "#FF5233") +
  geom_line(aes(y = z), color = "#FF5233") +
  ggtitle("Scree plot/Cumulative variance explained") +
  xlab("Principal component") +
  ylab("Variance explained") +
  labs(colour = "y")


#Looking at plot, seems like there is an "elbow" at #pc = 5
#suggesting K = 5. With K = 5 we explain 46.5% of the variance.
#We take K = 10, and thus explain 56% of the variance.
#P.S. the elbow thing maybe doesn't mean anything here. I might be confusing
#it with the index plots.


#>>>NOTE<<<: If it is not working well, we can come back here and use more PC's.
#If taking 75 PC's we can explain 75% of the variance

#-----TASK4-----
#First select the PC's decided in the previous task
df_filtered_PCA = as.data.frame(df_filtered.pca$x[,1:10])

#Cluster indices to use:
# >Within cluster scatter (mostly for k-means)
# >silhouette width
# >Davis-Bouldin index
# >Gap statistic

#Cluster algorithm: k-means
kk <- kmeans(df_filtered_PCA,5)

plot(df_filtered_PCA[c("PC1", "PC2")])
plot(df_filtered_PCA[c("PC1", "PC2")], 
     col = kk$cluster)

#Use within cluster scatter to determine how "good" the cluster is
NbClust(data = df_filtered_PCA, diss= NULL, distance = "euclidean",
        min.nc=2, max.nc=15, method = "kmeans",
        index = "all", alphaBeale = 0.1)

#Looking at the summary of the various indices, a clear majority proposed
#5 as the best number of clusters. (Although, maybe some indices are more
#suitable depending on what the clusters actually look like. But if we do not
#know anything about the data, maybe it is reasonable to just look at what
#the majority says.)

#Furthermore, looking at the index vs #clusters plot, we see the elbow at 5
#clusters, indicating that introducing more clusters does not improve the
#results significantly, suggesting that we are just splitting already found
#true clusters in half.

#Best number of clusters when using k-means: 5

#Comparing with true labels:
#Load label data
df_labels = readr::read_csv("TCGA-PANCAN-HiSeq-801x20531//labels.csv")
df_labels = df_labels[,-1]
barplot(table(df_labels$Class))

#We see that indeed there were only 5 types of cancer in the data, in accordance
#with our findings!

cm = table(df_labels$Class, kk$cluster)
cm


#Plot av clusters (sämre än den nedan)
# fviz_cluster(kk, data = df_filtered,
#              palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FF00FF", "#00FFFF"), 
#              geom = "point",
#              ellipse.type = "convex", 
#              ggtheme = theme_bw()
# )


ind.coord <- as.data.frame(get_pca_ind(df_filtered.pca)$coord)
ind.coord$cluster <- factor(kk$cluster)
ind.coord$Cancer <- df_labels$Class

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Cancer", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  # xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  # ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 5)

#OBS! Ser ut som att clustrerna överlappar, trots att k-means gör att datapunkterna
#assignas till närmaste centroid, i.e. borde inte kunna finnas en punkt av en
#viss färg mitt bland punkter av andra färger. Men kom då ihåg att vi har använt
#10 PC's vid beräkningen i k-means, och projicerar det nu bara på två av PC'na.
#Så det ser ut som att det har "blivit fel" om man tänker att k-means bara har
#kört i 2 dim, medan det (förmodligen) har blivit rätt bara att klustrerna är
#hypersfärer som inte går att visualisera.

#Ser att k-means har stora problem! Varför?


### Test av hierarchial clustering (hclust)

distance_mat = dist(df_filtered_PCA, method = "euclidean")

#Fit hierarchial clustering model to data
Hierar_cl <- hclust(distance_mat, method = "average")

#Plot dendogram
plot(Hierar_cl)


#Number of clusters
K = 6

# Cutting tree by k number of clusters
fit <- cutree(Hierar_cl, k = K)

#table(fit)
#rect.hclust(Hierar_cl, k = K, border = "green")

dis = dist(df_filtered_PCA)^2
#sil = silhouette(fit,dis)
# windows()
# plot(sil)


ind.coord <- as.data.frame(get_pca_ind(df_filtered.pca)$coord)
ind.coord$cluster <- factor(fit)
ind.coord$Cancer <- df_labels$Class

#Plot för många kluster
# P50 = createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
# 
# ggscatter(
#   ind.coord, x = "Dim.1", y = "Dim.2",
#   color = "cluster", palette = P50, ellipse = TRUE, ellipse.type = "convex",
#   shape = "Cancer", size = 1.5
# ) +
#   xlab("PC1")+
#   ylab("PC2")

# Plot för få antal kluster
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
#Calculate in turn the silhouette, db, and WSS- index and append to indices
max_clusters = 100
index_silhouette = NbClust(data = df_filtered_PCA,
                diss= NULL,
                distance = "euclidean",
                min.nc=2,
                max.nc=max_clusters,
                method = "average",
                index = c("silhouette"),
                alphaBeale = 0.1)

index_db = NbClust(data = df_filtered_PCA,
                    diss= NULL,
                    distance = "euclidean",
                    min.nc=2,
                    max.nc=max_clusters,
                    method = "average",
                    index = c("db"),
                    alphaBeale = 0.1)

index_WSS = NbClust(data = df_filtered_PCA,
                    diss= NULL,
                    distance = "euclidean",
                    min.nc=2,
                    max.nc=max_clusters,
                    method = "average",
                    index = c("ball"),
                    alphaBeale = 0.1)

index_silhouette = index_silhouette$All.index
index_db = index_db$All.index
index_WSS = index_WSS$All.index

n_clusters = 2:max_clusters
d = cbind(n_clusters, index_silhouette, index_db, index_WSS)


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

#------------------------------------------------

NbClust(data = df_filtered_PCA,
                           diss= NULL,
                           distance = "euclidean",
                           min.nc=2,
                           max.nc=15,
                           method = "average",
                           index = "all",
                           alphaBeale = 0.1)

Hierar_cl <- hclust(distance_mat, method = "single")

K = 37
  
fit = cutree(Hierar_cl, K)

cm = table(df_labels$Class, fit)
cm

##########################################################################
#-------------------------TASK 5------------------------------------------

library(ConsensusClusterPlus)



#-------------------------------------------------

# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

#-------HOW TO PLOT MULTIPLE LINES-----------------

a = 1:5
b = 1:5
c = (1:5)^2
d = cbind(a,b,c)

ggplot(d, aes(x = a)) +
  geom_point(aes(y = b), color = "red") +
  geom_line(aes(y = b), color = "red") +
  geom_point(aes(y = c), color = "blue") +
  geom_line(aes(y = c), color = "blue")


#-------LOGISTIC REGRESSION TEST-------------------

n = 100;

x = rnorm(n, 0, 0.5)
y = rnorm(n, 0, 0.5)
p0 = cbind(x,y,0)

x = rnorm(n, 1, 0.5)
y = rnorm(n, 1, 0.5)
p1 = cbind(x,y,1)

p = rbind(p0,p1)
df = as.data.frame(p)
colnames(df)[3] <- "class"

#qplot(x,y, colour = class, data = df)

model <- glm(class ~.,family=binomial(link='logit'),data=df)

beta = model$coefficients

x1 = -1
y1 = -(beta[1] + beta[2]*x1)/beta[3]

x2 = 2
y2 = -(beta[1] + beta[2]*x2)/beta[3]

sp = ggplot(data = df, aes(x = x, y = y)) +
  geom_point(aes(colour = as.factor(class)))
sp + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2))

#----------------------------------------------------------------------

#######################################################################

#-----------------PCA TEST---------------------------------------------


iris.pca = prcomp(iris[,-ncol(iris)], center = TRUE, scale. = TRUE)
#NOTE: iris.pca$x returnerar den transformerade datan efter basbytet!

#biplot(prcomps, scale = 0)

iris.pca.plot <- autoplot(iris.pca,
                          data = iris,
                          colour = 'Species')

iris.pca.plot

plot.iris.pca <- plot(iris.pca, type="l")
plot.iris.pca

#-----------------------------------------------------------------------

Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2, byrow = TRUE)
df = mvrnorm(n = 100, c(0, 1), Sigma)
df = as.dataframe(df)
colnames(df) <- c("x","y")

#Plotta se om det ser bra ut
ggplot(data = df, aes(x = x, y = y)) + geom_point()

#Prova använda PCA och se om det ser bra ut
df.pca = prcomp(df, scale. = TRUE)
df.pca.plot <- autoplot(df.pca,
                        data = df)
df.pca.plot
biplot(df.pca, scale = 0)