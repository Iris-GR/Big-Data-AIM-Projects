library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(corrplot)
library(caret)
library(CRAN)
library(smotefamily)
library(groupdata2)

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

Calc_metrics = function(fit_model, df_te)
{
  predicted_QDA = data.frame(predict(fit_model,df_te))
  colnames(predicted_QDA)<-c('class')
  
  CM = confusionMatrix(predicted_QDA$class, df_te$diagnosis)
  metric_value = error_metric(CM$table, 1.5)
  return(metric_value)
}


#Initialize list to keep track of means for the different imbalanced sets
metrics_noadjust <- NULL
metrics_upsample <- NULL
metrics_downsample <- NULL
metrics_Weighted <- NULL

#50/50 --------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
ix_remove = sample(ix_B,145)
df = uci_bc_data[-ix_remove,]
df %>% count(diagnosis)

#Partition data into training and test set (with stratification):
ix_tr = createDataPartition(df$diagnosis, p = 0.75, list = FALSE)
df_tr = df[ix_tr,]
df_te = df[-ix_tr,]

#Noadjust -----------------------|

fit_noadjust <- train(diagnosis ~.,
                      data = df_tr,
                      methods = 'qda')

metric_value = Calc_metrics(fit_noadjust, df_te)
metrics_noadjust = rbind(metrics_noadjust, metric_value)

#Upsample -----------------------|
df_upsample = upsample(df_tr, cat_col = "diagnosis")
df_upsample %>% count(diagnosis)

fit_upsample <- train(diagnosis ~.,
                      data = df_upsample,
                      method = 'qda')

metric_value = Calc_metrics(fit_upsample, df_te)
metrics_upsample = rbind(metrics_upsample, metric_value)


#Downsample -----------------------|
df_downsample = downsample(df_tr, cat_col = "diagnosis")
df_downsample %>% count(diagnosis)

fit_downsample <- train(diagnosis ~.,
                        data = df_downsample,
                        method = 'qda')


metric_value = Calc_metrics(fit_downsample, df_te)
metrics_downsample = rbind(metrics_downsample, metric_value)

#Class wieghts -------------------|
#Class weight (#ifelse(test, value if true, value if false))
model_weights <- ifelse(df_tr$diagnosis == "M",
                        (1/table(df_tr$diagnosis)[1]) * 0.5,
                        (1/table(df_tr$diagnosis)[2]) * 0.5)

fit_wieghted <- train(diagnosis ~.,
                      data = df_tr,
                      method = 'qda',
                      wieghts = model_weights)

metric_value = Calc_metrics(fit_wieghted, df_te)
metrics_Weighted = rbind(metrics_Weighted, metric_value)


  #40/60 ---------------------------------------------------------------------
ix_B = which(uci_bc_data$diagnosis == "B")
ix_remove = sample(ix_B,39)
df = uci_bc_data[-ix_remove,]
df %>% count(diagnosis)

#Partition data into training and test set (with stratification):
ix_tr = createDataPartition(df$diagnosis, p = 0.75, list = FALSE)
df_tr = df[ix_tr,]
df_te = df[-ix_tr,]

#Noadjust -----------------------|

fit_noadjust <- train(diagnosis ~.,
                      data = df_tr,
                      methods = 'qda')

metric_value = Calc_metrics(fit_noadjust, df_te)
metrics_noadjust = rbind(metrics_noadjust, metric_value)

#Upsample -----------------------|
df_upsample = upsample(df_tr, cat_col = "diagnosis")
df_upsample %>% count(diagnosis)

fit_upsample <- train(diagnosis ~.,
                      data = df_upsample,
                      method = 'qda')

metric_value = Calc_metrics(fit_upsample, df_te)
metrics_upsample = rbind(metrics_upsample, metric_value)


#Downsample -----------------------|
df_downsample = downsample(df_tr, cat_col = "diagnosis")
df_downsample %>% count(diagnosis)

fit_downsample <- train(diagnosis ~.,
                      data = df_downsample,
                      method = 'qda')


metric_value = Calc_metrics(fit_downsample, df_te)
metrics_downsample = rbind(metrics_downsample, metric_value)

#Class wieghts -------------------|
#Class weight (#ifelse(test, value if true, value if false))
model_weights <- ifelse(df_tr$diagnosis == "M",
                        (1/table(df_tr$diagnosis)[1]) * 0.5,
                        (1/table(df_tr$diagnosis)[2]) * 0.5)

fit_wieghted <- train(diagnosis ~.,
                      data = df_tr,
                      method = 'qda',
                      wieghts = model_weights)

metric_value = Calc_metrics(fit_wieghted, df_te)
metrics_Weighted = rbind(metrics_Weighted, metric_value) #Varför ger denna samma som metrics upsample ??????


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

#Partition data into training and test set (with stratification):
ix_tr = createDataPartition(df$diagnosis, p = 0.75, list = FALSE)
df_tr = df[ix_tr,]
df_te = df[-ix_tr,]

#Noadjust -----------------------|

fit_noadjust <- train(diagnosis ~.,
                      data = df_tr,
                      methods = 'qda')

metric_value = Calc_metrics(fit_noadjust, df_te)
metrics_noadjust = rbind(metrics_noadjust, metric_value)

#Upsample -----------------------|
df_upsample = upsample(df_tr, cat_col = "diagnosis")
df_upsample %>% count(diagnosis)

fit_upsample <- train(diagnosis ~.,
                      data = df_upsample,
                      method = 'qda')

metric_value = Calc_metrics(fit_upsample, df_te)
metrics_upsample = rbind(metrics_upsample, metric_value)


#Downsample -----------------------|
df_downsample = downsample(df_tr, cat_col = "diagnosis")
df_downsample %>% count(diagnosis)

fit_downsample <- train(diagnosis ~.,
                        data = df_downsample,
                        method = 'qda')


metric_value = Calc_metrics(fit_downsample, df_te)
metrics_downsample = rbind(metrics_downsample, metric_value)

#Class wieghts -------------------|
#Class weight (#ifelse(test, value if true, value if false))
model_weights <- ifelse(df_tr$diagnosis == "M",
                        (1/table(df_tr$diagnosis)[1]) * 0.5,
                        (1/table(df_tr$diagnosis)[2]) * 0.5)

fit_wieghted <- train(diagnosis ~.,
                      data = df_tr,
                      method = 'qda',
                      wieghts = model_weights)

metric_value = Calc_metrics(fit_wieghted, df_te)
metrics_Weighted = rbind(metrics_Weighted, metric_value)

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

#Partition data into training and test set (with stratification):
ix_tr = createDataPartition(df$diagnosis, p = 0.75, list = FALSE)
df_tr = df[ix_tr,]
df_te = df[-ix_tr,]

#Noadjust -----------------------|

fit_noadjust <- train(diagnosis ~.,
                      data = df_tr,
                      methods = 'qda')

metric_value = Calc_metrics(fit_noadjust, df_te)
metrics_noadjust = rbind(metrics_noadjust, metric_value)

#Upsample -----------------------|
df_upsample = upsample(df_tr, cat_col = "diagnosis")
df_upsample %>% count(diagnosis)

fit_upsample <- train(diagnosis ~.,
                      data = df_upsample,
                      method = 'qda')

metric_value = Calc_metrics(fit_upsample, df_te)
metrics_upsample = rbind(metrics_upsample, metric_value)


#Downsample -----------------------|
df_downsample = downsample(df_tr, cat_col = "diagnosis")
df_downsample %>% count(diagnosis)

fit_downsample <- train(diagnosis ~.,
                        data = df_downsample,
                        method = 'qda')


metric_value = Calc_metrics(fit_downsample, df_te)
metrics_downsample = rbind(metrics_downsample, metric_value)

#Class wieghts -------------------|
#Class weight (#ifelse(test, value if true, value if false))
model_weights <- ifelse(df_tr$diagnosis == "M",
                        (1/table(df_tr$diagnosis)[1]) * 0.5,
                        (1/table(df_tr$diagnosis)[2]) * 0.5)

fit_wieghted <- train(diagnosis ~.,
                      data = df_tr,
                      method = 'qda',
                      wieghts = model_weights)

metric_value = Calc_metrics(fit_wieghted, df_te)
metrics_Weighted = rbind(metrics_Weighted, metric_value)

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

#Partition data into training and test set (with stratification):
ix_tr = createDataPartition(df$diagnosis, p = 0.75, list = FALSE)
df_tr = df[ix_tr,]
df_te = df[-ix_tr,]

#Noadjust -----------------------|

fit_noadjust <- train(diagnosis ~.,
                      data = df_tr,
                      methods = 'qda')

metric_value = Calc_metrics(fit_noadjust, df_te)
metrics_noadjust = rbind(metrics_noadjust, metric_value)

#Upsample -----------------------|
df_upsample = upsample(df_tr, cat_col = "diagnosis")
df_upsample %>% count(diagnosis)

fit_upsample <- train(diagnosis ~.,
                      data = df_upsample,
                      method = 'qda')

metric_value = Calc_metrics(fit_upsample, df_te)
metrics_upsample = rbind(metrics_upsample, metric_value)


#Downsample -----------------------|
df_downsample = downsample(df_tr, cat_col = "diagnosis")
df_downsample %>% count(diagnosis)

fit_downsample <- train(diagnosis ~.,
                        data = df_downsample,
                        method = 'qda')


metric_value = Calc_metrics(fit_downsample, df_te)
metrics_downsample = rbind(metrics_downsample, metric_value)

#Class wieghts -------------------|
#Class weight (#ifelse(test, value if true, value if false))
model_weights <- ifelse(df_tr$diagnosis == "M",
                        (1/table(df_tr$diagnosis)[1]) * 0.5,
                        (1/table(df_tr$diagnosis)[2]) * 0.5)

fit_wieghted <- train(diagnosis ~.,
                      data = df_tr,
                      method = 'qda',
                      wieghts = model_weights)

metric_value = Calc_metrics(fit_wieghted, df_te)
metrics_Weighted = rbind(metrics_Weighted, metric_value)


# ---------------------------------------------
# summarize data:
# Recall:

i = 1;
rec_beta_acc = cbind(metrics_noadjust[,i],metrics_upsample[,i],metrics_downsample[,i],metrics_Weighted[,i])
matplot(rec_beta_acc, type = c("b"),pch=1,col = 1:4)

i = 2;
rec_beta_acc = cbind(metrics_noadjust[,i],metrics_upsample[,i],metrics_downsample[,i],metrics_Weighted[,i])
matplot(rec_beta_acc, type = c("b"),pch=1,col = 1:4)

i = 3;
rec_beta_acc = cbind(metrics_noadjust[,i],metrics_upsample[,i],metrics_downsample[,i],metrics_Weighted[,i])
matplot(rec_beta_acc, type = c("b"),pch=1,col = 1:4)

# Potential discard ----------------------------------------------------------

# mean_metric_values_QDA=function(DF, nfolds)
# {
#   ix_folds = createFolds(DF$diagnosis, k = nfolds)
#   
#   metric_values <- NULL
#   
#   for(i in 1:nfolds){
#     test_set = DF[unlist(ix_folds[i]),]
#     train_set = DF[-unlist(ix_folds[i]),]
#     
#     model_QDA <- train(diagnosis ~ ., train_set,
#                        method = "qda")
#     
#     predicted_QDA = data.frame(predict(model_QDA, test_set))
#     colnames(predicted_QDA)<-c('class')
#     
#     CM = confusionMatrix(predicted_QDA$class, test_set$diagnosis)
#     metric_value = error_metric(CM$table, 1.5)
#     metric_values = rbind(metric_values, metric_value);
#   }
#   
#   #Calculate the average values and the respective standard deviations:
#   avg_metric_values = colMeans(metric_values)
#   sd_metric_values = sapply(as.data.frame(metric_values), sd)
#   return(c(avg_metric_values, sd_metric_values))
# }