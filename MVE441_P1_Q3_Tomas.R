## The following R code loads the data in the same fashion as the Python code above

library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(corrplot)
library(caret)
library(CRAN)

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

#Pre-questions:

#1.
#Check if the class labels (B = Benign, M = Malignant) are imbalanced
#-->> M = 212, B = 357 i.e. 37% M, 63% B (slightly unbalanced)
uci_bc_data %>% count(diagnosis)

#2.A
#All parameters are numerical
str(uci_bc_data) #is there a way to summarize #parameters with certain data type?

#2.B
#Yes the data have different scales
boxplot(uci_bc_data[,3:NCOL(uci_bc_data)])

#2.C
#We see some strong correlations
#Plausible reasons for this? <---- have not checked
corr_matrix = cor(dplyr::select(uci_bc_data, c(-id_number,-diagnosis)))
corrplot(corr_matrix)

#Q3:

#Remove the identification number of patient, since it is not necessary here
uci_bc_data = subset(uci_bc_data, select = -c(id_number))

#Q3.1 <-----------------------------------
#Metrics:
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


#Skapa 10 folds:
ix_folds = createFolds(uci_bc_data$diagnosis, k = 10)

#Lista för att spara metrik-värdena
metric_values <- NULL

for(i in 1:10){
  test_set = uci_bc_data[unlist(ix_folds[i]),]
  train_set = uci_bc_data[-unlist(ix_folds[i]),]
  
  model_rf <- train(diagnosis ~ ., train_set,
                 method = "rf")
  
  predicted_rf = data.frame(predict(model_rf, test_set))
  colnames(predicted_rf)<-c('class')
  
  CM = confusionMatrix(predicted_rf$class, test_set$diagnosis)
  metric_value = error_metric(CM$table, 1.5)
  metric_values = rbind(metric_values, metric_value);
}

#Calculate the average values and the respective standard deviations:
avg_metric_values = colMeans(metric_values)
sd_metric_values = sapply(as.data.frame(metric_values), sd)
#averages = 0.9618497 0.9673302 0.9629386
#sd = 0.02633959 0.02137206 0.02440664


#Q3.2 <----------------------------------------------
#We want to ensure that someone with a malignant tumour is not classified
#as someone with a benign tumour, in other words we want to minimize the
#number of false negatives (and maximize the number of true positives)
#That is the same thing as maximizing the recall.

#At the same time, however, we also want people with benign tumours not to
#be classified as someone with a malignant tumour, since they then would
#undergo unecessary dangerous treatment. In other words we would like to
#minimize the number of false positives (and maximizing the number of true
#positives). That is the same thing as maximizing the precision.

#To compromise between the two, we are thus interested in the F_beta-score,
#which is somewhere inbetween, depending on the value of beta.

#Q3.3 <-----------------------------------------------

ix_folds = createFolds(uci_bc_data$diagnosis, k = 10)

metric_values <- NULL

for(i in 1:10){
  test_set = uci_bc_data[unlist(ix_folds[i]),]
  train_set = uci_bc_data[-unlist(ix_folds[i]),]
  
  model_QDA <- train(diagnosis ~ ., train_set,
                    method = "qda")
  
  predicted_rf = data.frame(predict(model_rf, test_set))
  colnames(predicted_rf)<-c('class')
  
  CM = confusionMatrix(predicted_rf$class, test_set$diagnosis)
  metric_value = error_metric(CM$table, 1.5)
  metric_values = rbind(metric_values, metric_value);
}

#Calculate the average values and the respective standard deviations:
avg_metric_values = colMeans(metric_values)
sd_metric_values = sapply(as.data.frame(metric_values), sd)
#averages = 0.9972222 0.9980603 0.9982143
#sd = 0.008784105 0.006133728 0.005646924

#Conclusion: QDA is working alot better! Don't know exactly what the settings
#of 'rf', can't find it in the documentation. But maybe it could produce better
#results by tweaking it, or using a different type of random forest.

#Q3.4 <------------------------------------------------------

#Stratified sampling?
