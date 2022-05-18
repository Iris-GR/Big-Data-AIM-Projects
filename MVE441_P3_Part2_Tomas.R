library(glmnet)
library(ggplot2)
library(MASS)
library(base)
library(viridis)

#Load data
X = readr::read_csv("TCGA-PANCAN-HiSeq-801x20531//data.csv")
X = X[,-1] #Remove first column containing row names

#Load labels
y = readr::read_csv("TCGA-PANCAN-HiSeq-801x20531//labels.csv")
y = y[,-1] #Remove first column containing row names
y = unlist(y)

#Filtering using F-statistic (using features with top 200 largest F-value)
X_no_const <- X[, apply(X, 2, var) > 0]
fval <- apply(X_no_const, 2, function(x) anova(lm(x ~ y))$`F value`[1])
X_filtered <- X_no_const[, order(fval, decreasing = TRUE)[1:200]]
X_filtered = as.matrix(X_filtered)

#----------------- Create single selection of features from one run -----------
#Fitting multi-class logistic regression model regularized with the L_1-norm
#of the regression coefficients.
#NOTE!!! We use lambda.min for our model
fit <- cv.glmnet(X_filtered,y,family = "multinomial")
beta = coef(fit, s="lambda.min")

#----------------- Check how good the selection was ---------------------------

beta_counter= matrix(0,200,5)

#Repeat M = 50 times:
for(m in 1:50){
  
  #Create bootstrap sample of filtered data (X_filtered) with 801 samples in each
  ix = sample(1:801, 801, replace = TRUE)
  X_m = X_filtered[ix,]
  y_m = y[ix]
  
  #Perform feature selection, i.e. determine which coeffs in fit that are non-zero
  #10 folds is default
  #CV loss function: Deviance
  fit_m <- cv.glmnet(X_m, y_m, family = "multinomial", type.measure = "deviance")
  
  #Record in an array of length(beta) the accumulated counts thus far
  beta_m = coef(fit_m, s="lambda.min")
  beta_m = cbind(as.matrix(beta_m$'BRCA'),
                 as.matrix(beta_m$'COAD'),
                 as.matrix(beta_m$'KIRC'),
                 as.matrix(beta_m$'LUAD'),
                 as.matrix(beta_m$'PRAD'))
  beta_m = beta_m[-1,] #Remove intercept
  beta_m = beta_m != 0 #Make logical
  
  #Behöver reda ut vilka features som ska få en count
  #Kolla i beta_m på namnen på raderna, blir det samma ordning varje gång?
  #Isf, bara gör en count vector och summera succesivt för varje run
  beta_counter = beta_counter + beta_m
}

beta_counter = as.data.frame(beta_counter)
colnames(beta_counter) <- c('BRCA','COAD','KIRC','LUAD','PRAD')
beta_counter= beta_counter/50

beta_counter = data.frame(
  gene = row.names(beta_counter),
  count = beta_counter
)

#Can try to plot the other ones too.
#This result is essentially the same as iris's
ggplot(beta_counter) +
  geom_bar( aes(x=reorder(gene,-count.BRCA), y=count.BRCA), stat="identity", fill="skyblue", alpha=0.7)+
  coord_cartesian(xlim = c(1, 20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#QUESTIONS:
#Q: 5 most important features for each class?
#A: BRCA: See plot
#Q: Are there features that occur appear often in several classes simultaneously?
#A: Write code to count how many times the 20 most important features appear from each class
#   Make histogram

