# Multiple imputation

# Data setup ----

# Load libraries
library(mice)
library(miceafter)
library(foreign)
library(tidyverse)
library(survival)
library(survminer)
library(car)
library(caret)
library(miceadds)

# Import meta data
Meta_data <- read.csv("phenoDataMir11.csv",sep="\t", header=TRUE)

# Split data to training and test cohors
set.seed(12)
training.samples <- Meta_data$Survival_status %>% 
  createDataPartition(p = 0.5, list = FALSE)
train.data  <- Meta_data[training.samples, ]
test.data <- Meta_data[-training.samples, ]

# Inspect missing REG IDs
joo <- which(is.na(Meta_data$REG_ID))
Meta_data$REG_ID[48] <- 9999 #Re-assing

# Save REG IDs as separate column
Ids <- test.data$REG_ID %>% data.frame()
rownames(Ids) <- Ids$.

# Create a new variable "Group" with info which
Meta_data$Group <- numeric(110)
testi <- Meta_data$REG_ID %in% rownames(Ids)
Meta_data$Group[which(testi==TRUE)] <- 1

# Create new dfs for multiple imputation
meta.data <- Meta_data[, c(8,9,12,16,18,19,20,21,23, 24,47,49,53:90)]

# Multiple imputation----

# Inspect NAs
sapply(meta.data, function(x) sum(is.na(x)))

# Assing variables
meta.data <- meta.data %>%
  #mutate(Survival_status = as.factor(Survival_status)) %>%
  #mutate(NSAID = as.factor(NSAID)) %>%
  mutate(Sex = as.factor(Sex)) %>%
  mutate(Variant = as.factor(Variant))

# Inspect the data structure
str(meta.data)

# Imputation
set.seed(1224)
init = mice(meta.data, maxit=0) 
meth = init$method
predM = init$predictorMatrix
predM # inspect
predM[,49] <- 0
imputed = mice(meta.data, method=meth, predictorMatrix=predM, m=50, maxit = 50)

# Save dataset
save(imputed, file="meta.data.imputed.RData")

# Check
imputed$imp

# Plot & inspect
plot(imputed)

# Cox model fitting and pooling ----


## Fit the model with the training group

# Imputed data sets n=50
nimp <- 50

# List to save all the models
models1 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  
  # Select only the cases in the training group
  sub <- subset(complete, Group == 0)
  
  
  # Fit coxph with lifestyle variables to get regression coefficients for c-miRs
  model1 <- coxph(Surv(Age_surveillance_end, Survival_status) ~ hsa.miR.10b.5p + hsa.miR.125b.5p + hsa.miR.200a.3p + hsa.miR.3613.5p + hsa.miR.3615 
                  + Sex, 
                  data = sub)
  
  # Save the 50 coxph models
  models1[i] <- list(model1)
  
  
}

# Pool and inspect
mirafit1 <- as.mira(models1)
pooled1 <- pool(mirafit1)
S1 <- summary(pooled1, conf.int = TRUE)
pooled1$glanced # Inspect

# miR coefficients adjusted with sex for computing risk sum score
# hsa.miR.10b.5p    2.1050
# hsa.miR.125b.5p   0.9368
# hsa.miR.200a.2p   1.1032
# hsa.miR.3613.5p   0.9866
# hsa.miR.3615      0.7875

# miR coefficients for computing risk sum score
# hsa.miR.10b.5p    1.8840
# hsa.miR.125b.5p   0.9344
# hsa.miR.200a.2p   0.8394
# hsa.miR.3613.5p   1.0443
# hsa.miR.3615      1.5757

### Fit a model to training group with risk sum score

# List to save all the models
models2 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.1050 + complete$hsa.miR.125b.5p*0.9368 + complete$hsa.miR.200a.3p*1.1032 
  + complete$hsa.miR.3613.5p*0.9866 + complete$hsa.miR.3615*0.7875
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select only the cases in the validation group
  sub <- subset(complete, Group == 0)
  
  # Fit coxph with only the risk sum score
  model2 <- coxph(Surv(Age_surveillance_end, Survival_status) ~ score,
                  data = sub)
  
  # Tallenenetaan mallit
  models2[i] <- list(model2)
  
}

# Pool and inspect
mirafit2 <- as.mira(models2)
pooled2 <- pool(mirafit2)
S2 <- summary(pooled2, conf.int = TRUE)
pooled2$glanced

# Make it a dataframe for plotting later on
model_table2 <- data.frame(matrix(ncol = 7 , nrow =1))
colnames(model_table2) <- c(' ', 'Estimate', "HR", "HR 95% lb", "HR 95% ub", "p_value", "C-index")
model_table2[,1] <- as.character(S2$term)
model_table2$Estimate <- S2$estimate
model_table2$HR <- exp(S2$estimate)
model_table2$p_value <- S2$p.value
model_table2$`HR 95% lb` <- exp(S2$`2.5 %`)
model_table2$`HR 95% ub` <- exp(S2$`97.5 %`)
model_table2$`C-index` <- 0.81

# Inspect
model_table2
write.table(x = model_table2, "training_group_cox_2.tsv", sep="\t")

### Test the model in training cohort with lifestyle variables

# List to save all the models
models3 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.813105 + complete$hsa.miR.3613.5p*1.077253  + complete$hsa.miR.200a.3p*1.976044  + complete$hsa.miR.3615*2.254473
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select only the cases in the validation group
  sub <- subset(complete, Group == 0)
  
  # Fit coxph with only the risk sum score
  model3 <- coxph(Surv(Age_surveillance_end, Survival_status) ~ score + Sex,
                  data = sub)
  
  # Save the models
  models3[i] <- list(model3)
  
}

# Pool and inspect
mirafit3 <- as.mira(models3)
pooled3 <- pool(mirafit3)
S3 <- summary(pooled3, conf.int = TRUE)
pooled3$glanced

c.indeces <- pooled3$glanced[13]
avg.c <- mean(c.indeces$concordance)
avg.c

# Make it a dataframe for plotting later on
model_table3 <- data.frame(matrix(ncol = 6 , nrow =2))
colnames(model_table3) <- c(' ', 'Estimate', "HR", "HR 95% lb", "HR 95% ub", "p_value")
model_table3[,1] <- as.character(S3$term)
model_table3$Estimate <- S3$estimate
model_table3$HR <- exp(S3$estimate)
model_table3$p_value <- S3$p.value
model_table3$`HR 95% lb` <- exp(S3$`2.5 %`)
model_table3$`HR 95% ub` <- exp(S3$`97.5 %`)
#model_table3$`C-index` <- 0.82

# Inspect
model_table3
write.table(x = model_table3, "training_group_cox_2_lifestyle.tsv", sep="\t")

### Test the fitted model with the validation group

# List to save all the models
models4 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.813105 + complete$hsa.miR.3613.5p*1.077253  + complete$hsa.miR.200a.3p*1.976044  + complete$hsa.miR.3615*2.254473
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select only the cases in the validation group
  sub <- subset(complete, Group == 1)
  
  # Fit coxph with only the risk sum score
  model4 <- coxph(Surv(Age_surveillance_end, Survival_status) ~ score,
                  data = sub)
  
  # Save the models
  models4[i] <- list(model4)
  
}

# Pool and inspect
mirafit4 <- as.mira(models4)
pooled4 <- pool(mirafit4)
S4 <- summary(pooled4, conf.int = TRUE)
pooled4$glanced

# Make it a dataframe for plotting later on
model_table4 <- data.frame(matrix(ncol = 7 , nrow =1))
colnames(model_table4) <- c(' ', 'Estimate', "HR", "HR 95% lb", "HR 95% ub", "p_value", "C-index")
model_table4[,1] <- as.character(S4$term)
model_table4$Estimate <- S4$estimate
model_table4$HR <- exp(S4$estimate)
model_table4$p_value <- S4$p.value
model_table4$`HR 95% lb` <- exp(S4$`2.5 %`)
model_table4$`HR 95% ub` <- exp(S4$`97.5 %`)
model_table4$`C-index` <- 0.60

# Inspect
model_table4
write.table(x = model_table4, "validation_group_cox_2.tsv", sep="\t")


### Test the models in the validation cohort with lifestyle variables

# List to save all the models
models5 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.813105 + complete$hsa.miR.3613.5p*1.077253  + complete$hsa.miR.200a.3p*1.976044  + complete$hsa.miR.3615*2.254473
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select only the cases in the validation group
  sub <- subset(complete, Group == 1)
  
  # Fit coxph with only the risk sum score
  model5 <- coxph(Surv(Age_surveillance_end, Survival_status) ~ score + Sex,
                  data = sub)
  
  # Save the models
  models5[i] <- list(model5)
  
}

# Pool and inspect
mirafit5 <- as.mira(models5)
pooled5 <- pool(mirafit5)
S5 <- summary(pooled5, conf.int = TRUE)
pooled5$glanced

c.indeces <- pooled5$glanced[13]
avg.c <- mean(c.indeces$concordance)
avg.c

# Make it a dataframe for plotting later on
model_table5 <- data.frame(matrix(ncol = 6 , nrow =2))
colnames(model_table5) <- c(' ', 'Estimate', "HR", "HR 95% lb", "HR 95% ub", "p_value")
model_table5[,1] <- as.character(S5$term)
model_table5$Estimate <- S5$estimate
model_table5$HR <- exp(S5$estimate)
model_table5$p_value <- S5$p.value
model_table5$`HR 95% lb` <- exp(S5$`2.5 %`)
model_table5$`HR 95% ub` <- exp(S5$`97.5 %`)
#model_table5$`C-index` <- 0.60

# Inspect
model_table5
write.table(x = model_table5, "validation_group_cox_2_lifestyle.tsv", sep="\t")


# Correlations of risk sum score and lifestyle variables ----

# PA

# Imputed data sets n=50
nimp <- 50

# List to save all the models
models6 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.1050 + complete$hsa.miR.125b.5p*0.9368 + complete$hsa.miR.200a.3p*1.1032 
  + complete$hsa.miR.3613.5p*0.9866 + complete$hsa.miR.3615*0.7875
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))

  
  # Select all cases
  sub <- subset(complete, Group <= 0)
  
  # Fit linear regression to see correlations
  model6 <- with(sub, lm(scale(score) ~ scale(PA_met_d)))

  # Save the models
  models6[i] <- list(model6)
  
}

# Pool and inspect
mirafit6 <- as.mira(models6)
pooled6 <- pool(mirafit6)
S6 <- summary(pooled6, conf.int = TRUE)
S6
#res_PA <- pool.r.squared(pooled6, adjusted = F)
#res_PA


# BMI

# List to save all the models
models7 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.1050 + complete$hsa.miR.125b.5p*0.9368 + complete$hsa.miR.200a.3p*1.1032 
  + complete$hsa.miR.3613.5p*0.9866 + complete$hsa.miR.3615*0.7875
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select all cases
  sub <- subset(complete, Group <= 0)
  
  # Fit linear regression to see correlations
  model7 <- with(sub, lm(scale(score) ~ scale(BMI_nayte_tai_kysely )))
  #model7 <- with(sub, lm(score ~ BMI_nayte_tai_kysely ))

  
  # Save the models
  models7[i] <- list(model7)
  
}

# Pool and inspect
mirafit7 <- as.mira(models7)
pooled7 <- pool(mirafit7)
S7 <- summary(pooled7, conf.int = TRUE)
S7
#res_BMI <- pool.r.squared(pooled7, adjusted = F)
#res_BMI


# Dietary fiber

# List to save all the models
models8 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.1050 + complete$hsa.miR.125b.5p*0.9368 + complete$hsa.miR.200a.3p*1.1032 
  + complete$hsa.miR.3613.5p*0.9866 + complete$hsa.miR.3615*0.7875
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select all cases
  sub <- subset(complete, Group <= 0)
  
  # Fit linear regression to see correlations
  model8 <- with(sub, lm(scale(score) ~ scale(fibc )))
  
  # Save the models
  models8[i] <- list(model8)
  
}

# Pool and inspect
mirafit8 <- as.mira(models8)
pooled8 <- pool(mirafit8)
S8 <- summary(pooled8, conf.int = TRUE)
S8
#res_fiber <- pool.r.squared(pooled8, adjusted = F)
#res_fiber


# NSAID usage

# List to save all the models
models9 <- list()

# For loop to fit coxph model for all the imputed datasets
for (i in 1:nimp) {
  
  # Assing complete dataset (no missing values)
  complete <- complete(imputed, i)
  
  # Compute risk sum score and score classificator for all the cases
  complete$score <- complete$hsa.miR.10b.5p*2.1050 + complete$hsa.miR.125b.5p*0.9368 + complete$hsa.miR.200a.3p*1.1032 
  + complete$hsa.miR.3613.5p*0.9866 + complete$hsa.miR.3615*0.7875
  tert <- quantile(complete$score, c(0:3/3))
  complete$risk_total <- with(complete, 
                              cut(score, tert, 
                                  include.lowest = T,
                                  labels=c("low", "medium", "high")))
  
  # Select all cases
  sub <- subset(complete, Group <= 0)
  
  # Fit linear regression to see correlations
  #model9 <- with(sub, lm(scale(score) ~ scale(NSAID )))
  model9 <- with(sub, lm(scale(score) ~ scale(Age )))
  
  # Save the models
  models9[i] <- list(model9)
  
}

# Pool and inspect
mirafit9 <- as.mira(models9)
pooled9 <- pool(mirafit9)
S9 <- summary(pooled9, conf.int = TRUE)
S9
#res_nsaid <- pool.r.squared(pooled9, adjusted = F)
#res_nsaid

