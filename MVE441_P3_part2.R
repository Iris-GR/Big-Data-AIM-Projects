# MVE441 Project 3, Question 2

# This question is about selecting features with confidence. The dataset used 
# is the gene expression dataset.

# Load libraries 
library(tibble)
library(glmnet)
library(tidyverse)
library(ggplot2)

# Task 1, Feature reduction----------------------------------------------------

# Load data
data <- tibble::as_tibble(read.csv("TCGA-PANCAN-HiSeq-801x20531//data.csv", 
                                   header = TRUE, 
                                   row.names = 1))

labels <- tibble::as_tibble(read.csv("TCGA-PANCAN-HiSeq-801x20531//labels.csv", 
                                     header = TRUE))

# Convert the data from tibble to a data frame and the labels into factors
X <- as.data.frame(data)
y.df <- as.data.frame(labels)
y <- as.factor(y.df$Class)

## Feature reduction by selecting top 200 features using the F statistic

# Remove constant features
X.no.const <- X[, apply(X, 2, var) > 0]

# Fitt a linear model to every non-constant column, do ANOVA and extract the 
# F-value. Store the F-values in 
F.value <- apply(X.no.const, 2, function(x) anova(lm(x~y))$`F value`[1])

# Extract the top 200 features/columns with the highest correlation to y 
# (the labels)
X.filtered <- X.no.const[, order(F.value, decreasing = TRUE)[1:200]]

# Task 2, (a mess) Performing feature selection--------------------------------

# alpha = 0, Ridge Regression, ungrouped
res.1.1 <- glmnet(
  X.filtered,
  y,
  family = "multinomial",
  alpha = 0,
  nlambda = 100,
  type.multinomial = "ungrouped"
)

# alpha = 1, Lasso Regression, ungrouped
res.1.2 <- glmnet(
  X.filtered,
  y,
  family = "multinomial",
  alpha = 1,
  nlambda = 100,
  type.multinomial = "ungrouped"
)

# alpha = 0.5, a 50/50 mixture of Ridge and Lasso Regression, ungrouped
res.1.3 <- glmnet(
  X.filtered,
  y,
  family = "multinomial",
  alpha = 0.5,
  nlambda = 100,
  type.multinomial = "ungrouped"
)

plot(res.1.1, label = TRUE)
plot(res.1.2, label = TRUE)
plot(res.1.3, label = TRUE)

# alpha = 0, Ridge Regression, grouped
res.2.1 <- glmnet(
  X.filtered,
  y,
  family = "multinomial",
  alpha = 0,
  nlambda = 100,
  type.multinomial = "grouped"
)

# alpha = 1, Lasso Regression, grouped
res.2.2 <- glmnet(
  X.filtered,
  y,
  family = "multinomial",
  alpha = 1,
  nlambda = 100,
  type.multinomial = "grouped"
)

# alpha = 0.5, a 50/50 mixture of Ridge and Lasso Regression, grouped
res.2.3 <- glmnet(
  X.filtered,
  y,
  family = "multinomial",
  alpha = 0.5,
  nlambda = 100,
  type.multinomial = "grouped"
)
plot(res.2.1, label = TRUE)
plot(res.2.2, label = TRUE)
plot(res.2.3, label = TRUE)


# 10-fold cross validation, alpha = 1, Lasso Regression
res.3 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "mse",
                   alpha = 1, 
                   type.multinomial = "ungrouped")
res.4 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "deviance",
                   alpha = 1, 
                   type.multinomial = "ungrouped")
res.5 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "class",
                   alpha = 1, 
                   type.multinomial = "ungrouped")
res.6 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "auc",
                   alpha = 1, 
                   type.multinomial = "ungrouped")

par(mfrow = c(2,2))
plot(res.3)
plot(res.4)
plot(res.5)
plot(res.6)

# 10-fold cross validation, alpha = 1, Ridge Regression
res.7 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "mse",
                   alpha = 0, 
                   type.multinomial = "ungrouped")
res.8 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "deviance",
                   alpha = 0, 
                   type.multinomial = "ungrouped")
res.9 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "class",
                   alpha = 0, 
                   type.multinomial = "ungrouped")
res.10 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "auc",
                   alpha = 0, 
                   type.multinomial = "ungrouped")

par(mfrow = c(2,2))
plot(res.7)
plot(res.8)
plot(res.9)
plot(res.10)

# 10-fold cross validation, alpha = 0.5, a 50/50 mixture of Ridge and Lasso 
# Regression, ungrouped
res.11 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "mse",
                   alpha = 0.5, 
                   type.multinomial = "ungrouped")
res.12 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "deviance",
                   alpha = 0.5, 
                   type.multinomial = "ungrouped")
res.13 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "class",
                   alpha = 0.5, 
                   type.multinomial = "ungrouped")
res.14 <- cv.glmnet(data.matrix(X.filtered),
                   y,
                   family = "multinomial", 
                   type.measure = "auc",
                   alpha = 0.5, 
                   type.multinomial = "ungrouped")

par(mfrow = c(2,2))
plot(res.11)
plot(res.12)
plot(res.13)
plot(res.14)

predict(res.11, newx = data.matrix(X.filtered),
        s = "lambda.min",
        type = "coef")

pred.a <- predict(res.11, newx = data.matrix(X.filtered),
        s = "lambda.min")

predict(res.11, newx = data.matrix(X.filtered)[1:10,],
        s = "lambda.1se",
        type = "class")

# predict(res.2,newx = x, s='lambda.min')
# predict(res.3, newx = x[1:10,], s = "lambda.min", type = "class")

res.4.pred.min <- predict(res.4, newx = data.matrix(X.filtered),
        s = "lambda.min",
        type = "coef")
res.4.pred.min$PRAD

res.4.pred.se <- predict(res.4, newx = data.matrix(X.filtered),
                          s = "lambda.1se",
                          type = "coef")
res.4.pred.se$PRAD

# predict(res.4, newx = data.matrix(X.filtered),
#         s = "lambda.min",
#         type = "class")
# 
# predict(res.4, newx = data.matrix(X.filtered),
#         s = "lambda.min",
#         type = "coef")




# Task 3 Deviance measure and lambda.1se --------------------------------------

# Make an empty list with five classes where each class has 200 + 1 long vectors
sel.freq.dev.1se <- as.data.frame(matrix(0, 
                                         nrow = 201,
                                         ncol = 5))

# Add gene names as row names and the classes as column names
rownames(sel.freq.dev.1se) <- c("(Intercept)",colnames(X.filtered))
colnames(sel.freq.dev.1se) <- c("BRCA.sel.freq",
                                "COAD.sel.freq",
                                "KIRC.sel.freq",
                                "LUAD.sel.freq",
                                "PRAD.sel.freq")

print(sel.freq.dev.1se)

# Convert the filtered data matrix to data.matrix format
dat.mat.filt <- data.matrix(X.filtered)

#* Boot strap ----------------------------------------------------------------

# Number of repeated bootstraps 
M = 200

# For each bootstrap
for (j in 1:M) {
  set.seed(j)
  
  # Take a sample with replacement from the filtered data matrix and the same
  # from the labels vector
  boot.rows <- sample(1:nrow(dat.mat.filt), nrow(dat.mat.filt), replace = TRUE)
  X.filtered.boot <- dat.mat.filt[boot.rows,]
  y.boot <- y[boot.rows]
  
  # 10-fold CV of multi class multinomial logistic, for deviance measure
  cv.boot.dev <- cv.glmnet(X.filtered.boot,
                           y.boot,
                           family = "multinomial", 
                           type.measure = "deviance",
                           alpha = 1, 
                           type.multinomial = "ungrouped")
  
  # Extract the feature coefficients for lambda within a std of the min
  pred.coef.se <- predict(cv.boot.dev, newx = X.filtered.boot,
                          s = "lambda.1se",
                          type = "coef")
  
  # Count frequency of features that were selected
  for (i in 1:length(pred.coef.se$BRCA)) {
    if(pred.coef.se$BRCA[i] != 0) {
      sel.freq.dev.1se$BRCA.sel.freq[i] = sel.freq.dev.1se$BRCA.sel.freq[i] + 1}
    
    if(pred.coef.se$COAD[i] != 0) {
      sel.freq.dev.1se$COAD.sel.freq[i] = sel.freq.dev.1se$COAD.sel.freq[i] + 1}
    
    if(pred.coef.se$KIRC[i] != 0) {
      sel.freq.dev.1se$KIRC.sel.freq[i] = sel.freq.dev.1se$KIRC.sel.freq[i] + 1}
    
    if(pred.coef.se$LUAD[i] != 0) {
      sel.freq.dev.1se$LUAD.sel.freq[i] = sel.freq.dev.1se$LUAD.sel.freq[i] + 1}
    
    if(pred.coef.se$PRAD[i] != 0) {
      sel.freq.dev.1se$PRAD.sel.freq[i] = sel.freq.dev.1se$PRAD.sel.freq[i] + 1}
  }
  
}

sel.freq.dev.1se

#* Prep of class data frames -------------------------------------------------

# Calculate the standard deviation of the frequency of selection for each class
sel.sd.dev.1se <- as.data.frame(matrix(NA, 
                                       nrow = length(sel.freq.dev.1se[,1]),
                                       ncol = length(sel.freq.dev.1se[1,]),
                                       dimnames = list(rownames(sel.freq.dev.1se),
                                                       c("BRCA.sel.sd", 
                                                         "COAD.sel.sd", 
                                                         "KIRC.sel.sd", 
                                                         "LUAD.sel.sd", 
                                                         "PRAD.sel.sd"))))

for (i in 1:length(sel.sd.dev.1se[,1])) {
  for (j in 1:length(sel.sd.dev.1se[1,])){
    sel.sd.dev.1se[i,j] <- sd(c(rep(0, M-sel.freq.dev.1se[i,j]), rep(1,sel.freq.dev.1se[i,j])))
  }
}

# Gather the coefficient name, selection frequency, stability and standard 
# deviation in one data frame
BRCA.sel.dev.1se <- data.frame(name = rownames(sel.freq.dev.1se),
                               sel.freq = sel.freq.dev.1se$BRCA.sel.freq,
                               sel.stab = sel.freq.dev.1se$BRCA.sel.freq/M,
                               sel.sd = sel.sd.dev.1se$BRCA.sel.sd)
COAD.sel.dev.1se <- data.frame(name = rownames(sel.freq.dev.1se),
                               sel.freq = sel.freq.dev.1se$COAD.sel.freq,
                               sel.stab = sel.freq.dev.1se$COAD.sel.freq/M,
                               sel.sd = sel.sd.dev.1se$COAD.sel.sd)
KIRC.sel.dev.1se <- data.frame(name = rownames(sel.freq.dev.1se),
                               sel.freq = sel.freq.dev.1se$KIRC.sel.freq,
                               sel.stab = sel.freq.dev.1se$KIRC.sel.freq/M,
                               sel.sd = sel.sd.dev.1se$KIRC.sel.sd)
LUAD.sel.dev.1se <- data.frame(name = rownames(sel.freq.dev.1se),
                               sel.freq = sel.freq.dev.1se$LUAD.sel.freq,
                               sel.stab = sel.freq.dev.1se$LUAD.sel.freq/M,
                               sel.sd = sel.sd.dev.1se$LUAD.sel.sd)
PRAD.sel.dev.1se <- data.frame(name = rownames(sel.freq.dev.1se),
                               sel.freq = sel.freq.dev.1se$PRAD.sel.freq,
                               sel.stab = sel.freq.dev.1se$PRAD.sel.freq/M,
                               sel.sd = sel.sd.dev.1se$PRAD.sel.sd)

# Order after decreasing coefficient selection frequency
BRCA.sel.dev.1se <- BRCA.sel.dev.1se[order(BRCA.sel.dev.1se$sel.freq,
                                           decreasing = TRUE), ]
COAD.sel.dev.1se <- COAD.sel.dev.1se[order(COAD.sel.dev.1se$sel.freq,
                                           decreasing = TRUE), ]
KIRC.sel.dev.1se <- KIRC.sel.dev.1se[order(KIRC.sel.dev.1se$sel.freq,
                                           decreasing = TRUE), ]
LUAD.sel.dev.1se <- LUAD.sel.dev.1se[order(LUAD.sel.dev.1se$sel.freq,
                                           decreasing = TRUE), ]
PRAD.sel.dev.1se <- PRAD.sel.dev.1se[order(PRAD.sel.dev.1se$sel.freq,
                                           decreasing = TRUE), ]

# Select rows with non-zero selection frequency
BRCA.sel.nz.dev.1se <- subset(BRCA.sel.dev.1se, BRCA.sel.dev.1se$sel.freq > 0)
COAD.sel.nz.dev.1se <- subset(COAD.sel.dev.1se, COAD.sel.dev.1se$sel.freq > 0)
KIRC.sel.nz.dev.1se <- subset(KIRC.sel.dev.1se, KIRC.sel.dev.1se$sel.freq > 0)
LUAD.sel.nz.dev.1se <- subset(LUAD.sel.dev.1se, LUAD.sel.dev.1se$sel.freq > 0)
PRAD.sel.nz.dev.1se <- subset(PRAD.sel.dev.1se, PRAD.sel.dev.1se$sel.freq > 0)

# For ordering of bars in plot, set the gene names as factors
BRCA.sel.nz.dev.1se$name <- factor(BRCA.sel.nz.dev.1se$name, 
                                   levels = BRCA.sel.nz.dev.1se$name)
COAD.sel.nz.dev.1se$name <- factor(COAD.sel.nz.dev.1se$name, 
                                   levels = COAD.sel.nz.dev.1se$name)
KIRC.sel.nz.dev.1se$name <- factor(KIRC.sel.nz.dev.1se$name, 
                                   levels = KIRC.sel.nz.dev.1se$name)
LUAD.sel.nz.dev.1se$name <- factor(LUAD.sel.nz.dev.1se$name, 
                                   levels = LUAD.sel.nz.dev.1se$name)
PRAD.sel.nz.dev.1se$name <- factor(PRAD.sel.nz.dev.1se$name, 
                                   levels = PRAD.sel.nz.dev.1se$name)

#* Plot stability of non-zero coefficients -----------------------------------
# BRCA
ggplot(BRCA.sel.nz.dev.1se) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: BRCA",
       subtitle = paste0("(CV measure: deviance, Lambda type: 1se, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# COAD
ggplot(COAD.sel.nz.dev.1se) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: COAD",
       subtitle = paste0("(CV measure: deviance, Lambda type: 1se, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# KIRC
ggplot(KIRC.sel.nz.dev.1se) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: KIRC",
       subtitle = paste0("(CV measure: deviance, Lambda type: 1se, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# LUAD
ggplot(LUAD.sel.nz.dev.1se) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: LUAD",
       subtitle = paste0("(CV measure: deviance, Lambda type: 1se, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# PRAD
ggplot(PRAD.sel.nz.dev.1se) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: PRAD",
       subtitle = paste0("(CV measure: deviance, Lambda type: 1se, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))


beep(5)

#* Same feature selected in multiple classes? ----------------------------------

# Class simultaneous selection of features for all non-zero features
mult.feat.sel.dev.1se <- c(BRCA.sel.nz.dev.1se$name,
                           COAD.sel.nz.dev.1se$name,
                           KIRC.sel.nz.dev.1se$name,
                           LUAD.sel.nz.dev.1se$name,
                           PRAD.sel.nz.dev.1se$name)
mult.feat.sel.dev.1se[duplicated(mult.feat.sel.dev.1se)]

# Class simultaneous selection of features for top 5 stable features per class
# (Don't include the intercept)
top.mult.feat.sel.dev.1se <- c(BRCA.sel.nz.dev.1se$name[2:6],
                               COAD.sel.nz.dev.1se$name[2:6],
                               KIRC.sel.nz.dev.1se$name[2:6],
                               LUAD.sel.nz.dev.1se$name[2:6],
                               PRAD.sel.nz.dev.1se$name[2:6])
top.mult.feat.sel.dev.1se[duplicated(top.mult.feat.sel.dev.1se)] 

# Class simultaneous selection of features for stability cutoffs per class
# (Don't include the intercept)
cut.offs <- seq(from = 0, to = 1, by = 0.01)
nr.dup.names.1se <- vector()

for (i in 1:length(cut.offs)) {
  dup.names.1se <- c(subset(BRCA.sel.nz.dev.1se$name, 
                          (BRCA.sel.nz.dev.1se$sel.stab > cut.offs[i]) 
                           & (BRCA.sel.nz.dev.1se$name != "(Intercept)")),
                   subset(COAD.sel.nz.dev.1se$name, 
                          (COAD.sel.nz.dev.1se$sel.stab > cut.offs[i]) 
                          & (COAD.sel.nz.dev.1se$name != "(Intercept)")),
                   subset(KIRC.sel.nz.dev.1se$name, 
                          (KIRC.sel.nz.dev.1se$sel.stab > cut.offs[i]) 
                          & (KIRC.sel.nz.dev.1se$name != "(Intercept)")),
                   subset(LUAD.sel.nz.dev.1se$name, 
                          (LUAD.sel.nz.dev.1se$sel.stab > cut.offs[i]) 
                          & (LUAD.sel.nz.dev.1se$name != "(Intercept)")),
                   subset(PRAD.sel.nz.dev.1se$name, 
                          (PRAD.sel.nz.dev.1se$sel.stab > cut.offs[i]) 
                          & (PRAD.sel.nz.dev.1se$name != "(Intercept)")))

  nr.dup.names.1se[i] <- length(dup.names.1se[duplicated(dup.names.1se)])
}

nr.dup.names.1se
plot(nr.dup.names.1se)


# Task 3 Deviance measure and lambda.min --------------------------------------
# Make an empty list with five classes where each class has 200 + 1 long vectors
sel.freq.dev.min <- as.data.frame(matrix(0, 
                                         nrow = 201,
                                         ncol = 5))

# Add gene names as row names and the classes as column names
rownames(sel.freq.dev.min) <- c("(Intercept)",colnames(X.filtered))
colnames(sel.freq.dev.min) <- c("BRCA.sel.freq",
                                "COAD.sel.freq",
                                "KIRC.sel.freq",
                                "LUAD.sel.freq",
                                "PRAD.sel.freq")

print(sel.freq.dev.min)

# Convert the filtered datamatrix to data.matrix format
dat.mat.filt <- data.matrix(X.filtered)

#* Boot strap ------------------------------------------------------------------

# Number of repeated bootstraps 
M = 200

# For each bootstrap
for (j in 1:M) {
  set.seed(j)
  
  # Take a sample with replacement from the filtered data matrix and the same
  # from the labels vector
  boot.rows <- sample(1:nrow(dat.mat.filt), nrow(dat.mat.filt), replace = TRUE)
  X.filtered.boot <- dat.mat.filt[boot.rows,]
  y.boot <- y[boot.rows]
  
  # 10-fold CV of multiclass multinomial logistic, for deviance measure
  cv.boot.dev <- cv.glmnet(X.filtered.boot,
                           y.boot,
                           family = "multinomial", 
                           type.measure = "deviance",
                           alpha = 1, 
                           type.multinomial = "ungrouped")
  
  # Extract the freature coefficients for minimum lambda
  pred.coef.se <- predict(cv.boot.dev, newx = X.filtered.boot,
                          s = "lambda.min",
                          type = "coef")
  
  # Count frequency of features that were selected
  for (i in 1:length(pred.coef.se$BRCA)) {
    if(pred.coef.se$BRCA[i] != 0) {
      sel.freq.dev.min$BRCA.sel.freq[i] = sel.freq.dev.min$BRCA.sel.freq[i] + 1}
    
    if(pred.coef.se$COAD[i] != 0) {
      sel.freq.dev.min$COAD.sel.freq[i] = sel.freq.dev.min$COAD.sel.freq[i] + 1}
    
    if(pred.coef.se$KIRC[i] != 0) {
      sel.freq.dev.min$KIRC.sel.freq[i] = sel.freq.dev.min$KIRC.sel.freq[i] + 1}
    
    if(pred.coef.se$LUAD[i] != 0) {
      sel.freq.dev.min$LUAD.sel.freq[i] = sel.freq.dev.min$LUAD.sel.freq[i] + 1}
    
    if(pred.coef.se$PRAD[i] != 0) {
      sel.freq.dev.min$PRAD.sel.freq[i] = sel.freq.dev.min$PRAD.sel.freq[i] + 1}
  }
  
}

sel.freq.dev.min

#* Prep of class data frames -------------------------------------------------

# Calculate the standard deviation of the frequency of selection for each class
sel.sd.dev.min <- as.data.frame(matrix(NA, 
                                       nrow = length(sel.freq.dev.min[,1]),
                                       ncol = length(sel.freq.dev.min[1,]),
                                       dimnames = list(rownames(sel.freq.dev.min),
                                                       c("BRCA.sel.sd", 
                                                         "COAD.sel.sd", 
                                                         "KIRC.sel.sd", 
                                                         "LUAD.sel.sd", 
                                                         "PRAD.sel.sd"))))

for (i in 1:length(sel.sd.dev.min[,1])) {
  for (j in 1:length(sel.sd.dev.min[1,])){
    sel.sd.dev.min[i,j] <- sd(c(rep(0, M-sel.freq.dev.min[i,j]), rep(1,sel.freq.dev.min[i,j])))
  }
}

# Gather the coefficient name, selection frequency, stability and standard 
# deviation in one data frame
BRCA.sel.dev.min <- data.frame(name = rownames(sel.freq.dev.min),
                               sel.freq = sel.freq.dev.min$BRCA.sel.freq,
                               sel.stab = sel.freq.dev.min$BRCA.sel.freq/M,
                               sel.sd = sel.sd.dev.min$BRCA.sel.sd)
COAD.sel.dev.min <- data.frame(name = rownames(sel.freq.dev.min),
                               sel.freq = sel.freq.dev.min$COAD.sel.freq,
                               sel.stab = sel.freq.dev.min$COAD.sel.freq/M,
                               sel.sd = sel.sd.dev.min$COAD.sel.sd)
KIRC.sel.dev.min <- data.frame(name = rownames(sel.freq.dev.min),
                               sel.freq = sel.freq.dev.min$KIRC.sel.freq,
                               sel.stab = sel.freq.dev.min$KIRC.sel.freq/M,
                               sel.sd = sel.sd.dev.min$KIRC.sel.sd)
LUAD.sel.dev.min <- data.frame(name = rownames(sel.freq.dev.min),
                               sel.freq = sel.freq.dev.min$LUAD.sel.freq,
                               sel.stab = sel.freq.dev.min$LUAD.sel.freq/M,
                               sel.sd = sel.sd.dev.min$LUAD.sel.sd)
PRAD.sel.dev.min <- data.frame(name = rownames(sel.freq.dev.min),
                               sel.freq = sel.freq.dev.min$PRAD.sel.freq,
                               sel.stab = sel.freq.dev.min$PRAD.sel.freq/M,
                               sel.sd = sel.sd.dev.min$PRAD.sel.sd)

# Order after decreasing coefficient selection frequency
BRCA.sel.dev.min <- BRCA.sel.dev.min[order(BRCA.sel.dev.min$sel.freq,
                                           decreasing = TRUE), ]
COAD.sel.dev.min <- COAD.sel.dev.min[order(COAD.sel.dev.min$sel.freq,
                                           decreasing = TRUE), ]
KIRC.sel.dev.min <- KIRC.sel.dev.min[order(KIRC.sel.dev.min$sel.freq,
                                           decreasing = TRUE), ]
LUAD.sel.dev.min <- LUAD.sel.dev.min[order(LUAD.sel.dev.min$sel.freq,
                                           decreasing = TRUE), ]
PRAD.sel.dev.min <- PRAD.sel.dev.min[order(PRAD.sel.dev.min$sel.freq,
                                           decreasing = TRUE), ]

# Select rows with non-zero selection frequency
BRCA.sel.nz.dev.min <- subset(BRCA.sel.dev.min, BRCA.sel.dev.min$sel.freq > 0)
COAD.sel.nz.dev.min <- subset(COAD.sel.dev.min, COAD.sel.dev.min$sel.freq > 0)
KIRC.sel.nz.dev.min <- subset(KIRC.sel.dev.min, KIRC.sel.dev.min$sel.freq > 0)
LUAD.sel.nz.dev.min <- subset(LUAD.sel.dev.min, LUAD.sel.dev.min$sel.freq > 0)
PRAD.sel.nz.dev.min <- subset(PRAD.sel.dev.min, PRAD.sel.dev.min$sel.freq > 0)

# For ordering of bars in plot, set the gene names as factors
BRCA.sel.nz.dev.min$name <- factor(BRCA.sel.nz.dev.min$name, 
                                   levels = BRCA.sel.nz.dev.min$name)
COAD.sel.nz.dev.min$name <- factor(COAD.sel.nz.dev.min$name, 
                                   levels = COAD.sel.nz.dev.min$name)
KIRC.sel.nz.dev.min$name <- factor(KIRC.sel.nz.dev.min$name, 
                                   levels = KIRC.sel.nz.dev.min$name)
LUAD.sel.nz.dev.min$name <- factor(LUAD.sel.nz.dev.min$name, 
                                   levels = LUAD.sel.nz.dev.min$name)
PRAD.sel.nz.dev.min$name <- factor(PRAD.sel.nz.dev.min$name, 
                                   levels = PRAD.sel.nz.dev.min$name)

#* Plot stability of non-zero coefficients -----------------------------------
# BRCA
ggplot(BRCA.sel.nz.dev.min) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: BRCA",
       subtitle = paste0("(CV measure: deviance, Lambda type: min, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# COAD
ggplot(COAD.sel.nz.dev.min) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: COAD",
       subtitle = paste0("(CV measure: deviance, Lambda type: min, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# KIRC
ggplot(KIRC.sel.nz.dev.min) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: KIRC",
       subtitle = paste0("(CV measure: deviance, Lambda type: min, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# LUAD
ggplot(LUAD.sel.nz.dev.min) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: LUAD",
       subtitle = paste0("(CV measure: deviance, Lambda type: min, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

# PRAD
ggplot(PRAD.sel.nz.dev.min) +
  geom_bar(aes(x = name, 
               y = sel.stab),
           stat = "identity", 
           fill = "skyblue", 
           alpha = 0.7) +
  geom_errorbar(aes(x = name, 
                    ymin = sel.stab - sel.sd/sqrt(M), 
                    ymax = sel.stab + sel.sd/sqrt(M)), 
                width = 0.4, 
                colour = "orange", 
                alpha = 0.9, 
                size = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  labs(title = "Feature Selection Stability, Class: PRAD",
       subtitle = paste0("(CV measure: deviance, Lambda type: min, Nr. bootstrap: ",
                         M, ")"),
       x = NULL, 
       y = "Stability") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))

#* Same feature selected in multiple classes? ----------------------------------

# Class simultaneous selection of features for all non-zero features
mult.feat.sel.dev.min <- c(BRCA.sel.nz.dev.min$name,
                           COAD.sel.nz.dev.min$name,
                           KIRC.sel.nz.dev.min$name,
                           LUAD.sel.nz.dev.min$name,
                           PRAD.sel.nz.dev.min$name)
mult.feat.sel.dev.min[duplicated(mult.feat.sel.dev.min)]

# Class simultaneous selection of features for top 5 stable features per class
# (Don't include the intercept)
top.mult.feat.sel.dev.min <- c(BRCA.sel.nz.dev.min$name[2:6],
                               COAD.sel.nz.dev.min$name[2:6],
                               KIRC.sel.nz.dev.min$name[2:6],
                               LUAD.sel.nz.dev.min$name[2:6],
                               PRAD.sel.nz.dev.min$name[2:6])
top.mult.feat.sel.dev.min[duplicated(top.mult.feat.sel.dev.min)] 

