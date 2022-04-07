library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(corrplot)
library(caret)
library(CRAN)
library(smotefamily)

names <- c("id_number", "diagnosis", "radius_mean",
           "texture_mean", "perimeter_mean", "area_mean",
           "smoothness_mean", "compactness_mean",
           "concavity_mean","concave_points_mean",
           "symmetry_mean", "fractal_dimension_mean",
           "radius_se", "texture_se", "perimeter_se",
           "area_se", "smoothness_se", "compactness_se",
           "concavity_se", "concave_points_se",
           "symmetry_se", "fractal_dimension_se",
           "radius_worst", "texture_worst",
           "perimeter_worst", "area_worst",
           "smoothness_worst", "compactness_worst",
           "concavity_worst", "concave_points_worst",
           "symmetry_worst", "fractal_dimension_worst")

uci_bc_data <- read_delim(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
  delim = ",",
  col_names = names,
  col_types = cols(
    .default = col_number(),
    id_number = col_integer(),
    diagnosis = col_factor()
  ))

y <- uci_bc_data %>% 
  mutate(diagnosis = as.factor(case_when(diagnosis == "B" ~ 0, diagnosis == "M" ~ 1))) %>%
  select(diagnosis) %>%
  as_vector() %>%
  unname()
X <- uci_bc_data %>% 
  select(-id_number, -diagnosis) %>%
  as.matrix()

# -------------------- TOMAS CODE BELOW -----------------------

#Remove id_number
uci_bc_data = subset(uci_bc_data, select = -c(id_number))

#Summary of task
#M = 212, B = 357 i.e. 37% M, 63% B
#Adjusted imbalances:
#M = 212, B = 212 i.e. 50/50 (-145 of B)
#M = 212, B = 318 i.e. 40/60 (-39)
#M = 212, B = 495 i.e. 30/70 (+138)
#M = 212, B = 848 i.e. 20/80 (+491)
#M = 212, B = 1908 i.e. 10/90 (+1551)

#Help functions:
error_metric=function(CM,beta)
{
  
  TN =CM[1,1]
  TP =CM[2,2]
  FP =CM[1,2]
  FN =CM[2,1]
  precision =(TP)/(TP+FP)
  accuracy_model  =(TP+TN)/(TP+TN+FP+FN)
  recall_model = (TP)/(TP+FN)
  F_beta = (1+beta^2)*(precision*recall_model)/(beta^2 * precision + recall_model)
  # print(paste("Recall value of the model: ",round(recall_model,2)))
  # print(paste("F_beta of the model: ",round(F_beta,2)))
  # print(paste("Accuracy of the model: ",round(accuracy_model,2)))
  return(c(recall_model,F_beta,accuracy_model))
}

mean_metric_values_QDA=function(DF, nfolds)
{
  ix_folds = createFolds(DF$diagnosis, k = nfolds)
  
  metric_values <- NULL
  
  for(i in 1:nfolds){
    test_set = DF[unlist(ix_folds[i]),]
    train_set = DF[-unlist(ix_folds[i]),]
    
    model_QDA <- train(diagnosis ~ ., train_set,
                       method = "qda")
    
    predicted_QDA = data.frame(predict(model_QDA, test_set))
    colnames(predicted_QDA)<-c('class')
    
    CM = confusionMatrix(predicted_QDA$class, test_set$diagnosis)
    metric_value = error_metric(CM$table, 1.5)
    metric_values = rbind(metric_values, metric_value);
  }
  
  #Calculate the average values and the respective standard deviations:
  avg_metric_values = colMeans(metric_values)
  sd_metric_values = sapply(as.data.frame(metric_values), sd)
  return(c(avg_metric_values, sd_metric_values))
}

#Initialize list to keep track of means for the different imbalanced sets
metrics_noadjust <- NULL
metrics_oversample <- NULL

#50/50 --------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
ix_remove = sample(ix_B,145)
df = uci_bc_data[-ix_remove,]
df %>% count(diagnosis)

#No adjustment
metrics_noadjust = rbind(metrics,mean_metric_values_QDA(df,10))

#Oversample
df = SMOTE() #Fattar inte syntaxen

#Class weight

  #40/60 ---------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
ix_remove = sample(ix_B,39)
df = uci_bc_data[-ix_remove,]
df %>% count(diagnosis)

#30/70 ---------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
df = uci_bc_data

for(i in 1:138){
  ix = sample(ix_B,2)
  ncols = ncol(uci_bc_data)
  
  #Convex combination of the two points
  new_data = 0.95*uci_bc_data[ix[1],2:ncols] + 0.05*uci_bc_data[ix[2],2:ncols]
  
  #Appending the class "B"
  new_data = cbind(diagnosis = factor("B"), new_data)
  
  #Appending to dataframe
  df = rbind(df, new_data)
}

#Shuffle dataframe (so that not all "artificial points" are at the end)
df = df[sample(1:nrow(df)),]

df %>% count(diagnosis)

#20/80 ------------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
df = uci_bc_data

for(i in 1:491){
  ix = sample(ix_B,2)
  ncols = ncol(uci_bc_data)
  
  #Convex combination of the two points
  new_data = 0.95*uci_bc_data[ix[1],2:ncols] + 0.05*uci_bc_data[ix[2],2:ncols]
  
  #Appending the class "B"
  new_data = cbind(diagnosis = factor("B"), new_data)
  
  #Appending to dataframe
  df = rbind(df, new_data)
}

#Shuffle dataframe (so that not all "artificial points" are at the end)
df = df[sample(1:nrow(df)),]

df %>% count(diagnosis)

#10/90 -----------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
df = uci_bc_data

for(i in 1:1551){
  ix = sample(ix_B,2)
  ncols = ncol(uci_bc_data)
  
  #Convex combination of the two points
  new_data = 0.95*uci_bc_data[ix[1],2:ncols] + 0.05*uci_bc_data[ix[2],2:ncols]
  
  #Appending the class "B"
  new_data = cbind(diagnosis = factor("B"), new_data)
  
  #Appending to dataframe
  df = rbind(df, new_data)
}

#Shuffle dataframe (so that not all "artificial points" are at the end)
df = df[sample(1:nrow(df)),]

df %>% count(diagnosis)



