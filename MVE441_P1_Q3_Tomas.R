## The following R code loads the data in the same fashion as the Python code above

library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(corrplot)

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

#Prova random forest
