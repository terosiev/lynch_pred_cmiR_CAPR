# c-miR risk sum cancer risk prediction model

# Script is written by Matti Hyvärinen & Tero Sievänen (2023)

# Load libraries
library(glmnet)
library(caret)
library(car)
library(tidyverse)
library(survival)
library(survminer)
library(magrittr)

# Import meta data
meta.data <- read.csv("phenoData.csv",sep="\t", header=TRUE)

#Rename columns
meta.data <- meta.data %>% 
  dplyr::rename("status" = "Survival_status")

meta.data <- meta.data %>% 
  dplyr::rename("time" = "Age_surveillance_end")

# Create model in full sample ----

# Create surv-object
y <- Surv(meta.data$time, meta.data$status)

# Extract predictors
x <- as.matrix(meta.data[, c(14:50)]) %>% # c-miRs
  scale(center = T) %>%
  na.omit()

# Set seed for reproducibility
set.seed(1234)

# Inspect LASSO fit lambdas
fit <- glmnet(x, y, family ="cox")
print(fit)

# Plot
plot(fit, label = T, xvar="lambda")

# Create folds
nfolds <- 10
foldid <- sample(rep(seq(nfolds), length.out = nrow(x)))

# Run Lasso
cvfit <- cv.glmnet(x,y,family = "cox", alpha = 1, nfolds = nfolds, foldid = foldid, maxit = 1000)
plot(cvfit)

# Get coefficients (miRs)
coefs <- coef(cvfit, s = "lambda.min")
coefs

# Predictors: hsa.miR.10b.5p, hsa.miR.125b.5p, hsa.miR.200a.3p ja hsa.miR.3613.5p and hsa-miR.3615

# Fit model in full sample
full.cox <- coxph(Surv(time, status) ~ hsa.miR.10b.5p + hsa.miR.125b.5p + 
                    hsa.miR.200a.3p + hsa.miR.3613.5p + hsa.miR.3615, data = meta.data)
summary(full.cox)
confint(full.cox)

# The full model
#                 coef    exp(coef) se(coef)  z    p
#hsa.miR.10b.5p  1.8840    6.5797   0.7776 2.423 0.0154
#hsa.miR.125b.5p 0.9344    2.5456   0.7579 1.233 0.2176
#hsa.miR.200a.3p 0.8394    2.3149   0.6431 1.305 0.1918
#hsa.miR.3613.5p 1.0443    2.8415   0.5753 1.815 0.0695
#hsa.miR.3615    1.5757    4.8343   1.0037 1.570 0.1164

# Compute c-miR risk sum score (expression * Cox regression coefficient ..)
meta.data$full_score <- 
  meta.data$hsa.miR.10b.5p*1.8840 + 
  meta.data$hsa.miR.125b.5p*0.9344 +
  meta.data$hsa.miR.200a.3p*0.8394 + 
  meta.data$hsa.miR.3613.5p*1.0443 + 
  meta.data$hsa.miR.3615*1.5757

# Fit c-miR risk sum score model (unadjusted)
full.cox.score <- coxph(Surv(time, status) ~ full_score, data = meta.data)
summary(full.cox.score)

# Score HR 2.718, p=0.000113
# Model C-index: 0.722 (se=0.076)

# Inspect Schoenfeld residuals
scho.cox <- cox.zph(full.cox.score)

# Plot Schoenfeld residuals
scho.residuals <- ggcoxzph(scho.cox)

# Plot dfbetas
ggcoxdiagnostics(full.cox.score, type = "dfbeta",
                 linear.predictions = F, 
                 ggtheme = theme_classic())

# C-miR risk prediction model cross-validation ----

# Set seed
set.seed(1234)

# Number of folds
folds <- 5
n <- length(meta.data$Filename)

# Split data into N folds
fold <- sample(rep(seq(folds), length.out = nrow(meta.data)))
meta.data$fold <- fold

# For loop for 5x cross-validation of the full model

for(i in 1:folds) { 
  
  message("Create training and validation data")
  int <- which(meta.data$fold == i)
  train.data <- meta.data[-int, ]
  test.data <- meta.data[int, ]
  
  # New df for coxph in the training data set
  x.train <- train.data[, c("hsa.miR.10b.5p", "hsa.miR.125b.5p","hsa.miR.200a.3p", "hsa.miR.3613.5p", "hsa.miR.3615" )]
  x.train$time <- train.data$time
  x.train$status <- train.data$status
  
  # New df for coxph in the validation data set
  x.test <- test.data[, c("hsa.miR.10b.5p", "hsa.miR.125b.5p","hsa.miR.200a.3p", "hsa.miR.3613.5p", "hsa.miR.3615" )]
  x.test$time <- test.data$time
  x.test$status <- test.data$status
  
  # Fit the multivariate Cox regression model
  message ("Fit Cox regression model")
  res.cox <- coxph(Surv(time, status) ~  .
                   , data = x.train)
  
  # Loop to compute risk sum score
  message("Compute c-miR risk sum score for training and validation data:")
  score <- numeric(length(x.train[,1]))
  score.test <- numeric(length(x.test[,1]))
  
  for(i in 1:length(res.cox$coefficients)) { 
    score <- score + (x.train[,i]*res.cox$coefficients[[i]]) # Score for training set
    score.test <- score.test + (x.test[,i]*res.cox$coefficients[[i]]) # Score for validation set
  }
  
  # Training data
  message("Training data risk score:")
  x.train$score <- score
  print(x.train$score)
  
  # Validation data
  message("Validation data risk score:")
  x.test$score <- score.test
  print(x.test$score)
  
  # Create a model in the train set
  message("Fit Cox regression model in the training dataset using risk sum score")
  message("Model summary:")
  res.cox.score <- coxph(Surv(time, status) ~  score, data = x.train)
  print(summary(res.cox.score))
  
  # Make predictions in the validation set
  message ("Making predictions with validation dataset")
  message ("C-index:")
  y2 <- Surv(x.test$time, x.test$status)
  length <- length(x.test)
  x2 <- x.test[, c((length-2):length)]
  cox_preds <- predict(res.cox.score, newdata = x2)
  
  # Calculate C-index
  Cindex <- Cindex(cox_preds, y2)
  print(Cindex)
}
