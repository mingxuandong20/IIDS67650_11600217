#Data import and initial exploration
##Read data

library(tidyverse)

df <- read_csv("C:/Users/mingx/Desktop/dissertation/Framingham Dataset.csv")

###Check out the previous lines
head(df)

###View data structures
str(df)

###Check the data dimensions
dim(df)

###Basic description of the data
summary(df)

###Look at the situation of the variables in each PERIOD to determine that it is longitudinal data
table(df$PERIOD)

##Extract the start and stop to construct the survival interval
library(dplyr)
data_sorted <- arrange(df, RANDID, TIME)
data_grouped <- group_by(data_sorted, RANDID)
data_mutated <- mutate(
  data_grouped,
  start = TIME, 
  next_time = lead(TIME),
  stop_raw = pmin(next_time, TIMECVD, na.rm = TRUE),
  stop = ifelse(is.na(stop_raw), TIMECVD, stop_raw),
  event = ifelse(stop >= TIMECVD & first(CVD) == 1, 1, 0)
)
data_ungrouped <- ungroup(data_mutated)
data_filtered <- filter(data_ungrouped, stop > start)
data <- select(data_filtered, -next_time, -stop_raw)

head(data, 10)

###Select the appropriate variables and reassemble the new dataset
###Pick out the variables you need and generate a new dataset
useable_vars <- c("RANDID", "AGE", "SEX", "SYSBP", "BMI", "TOTCHOL", 
                  "CURSMOKE", "DIABETES", "TIMECVD", "PERIOD", "start", "stop", "event")

data_selected <- data[, useable_vars]

colnames(data_selected)[colnames(data_selected) == "event"] <- "CVD"

head(data_selected)

###Check the new dataset structure
str(data_selected)
summary(data_selected)

###Query for missing status
###Missing scale table
colSums(is.na(data_selected)) / nrow(data_selected)

###The proportion of missing values is very small, so delete them directly
data_final <- data_selected[complete.cases(data_selected), ]

###Basic description of the data
summary(data_final)

###See how variables are doing in each PERIOD
table(data_final$PERIOD)

##Handle outliers
###List of continuous variables
cont_vars <- c("AGE", "SYSBP", "BMI", "TOTCHOL")

###Draw a boxplot
par(mfrow=c(2,2))
for (v in cont_vars) {
  boxplot(data_final[[v]], main=paste("Boxplot of", v))
}
par(mfrow=c(1,1))

summary(data_final[, c("AGE", "SYSBP", "BMI", "TOTCHOL")])

###Define the outlier value function
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR
  upper <- Q3 + 1.5 * IQR
  x >= lower & x <= upper
}

###Remove any variable that is marked as an extreme value
outlier_flag <- rep(TRUE, nrow(data_final))
for (v in cont_vars) {
  outlier_flag <- outlier_flag & remove_outliers(data_final[[v]])
}
data_final <- data_final[outlier_flag, ]

###Review the results
summary(data_final[, cont_vars])
par(mfrow = c(2, 2))
for (v in cont_vars) {
  boxplot(data_final[[v]], main = paste("Boxplot of", v))
}
par(mfrow = c(1, 1))

###Check the new dataset structure
str(data_final)
summary(data_final)
################################################################################
##Divide the training set and the test set
set.seed(123)
###Get all unique patient IDs
ids <- unique(data_final$RANDID)

###The training set ID is randomly selected in a 7:3 ratio
train_id <- sample(ids, size = round(length(ids) * 0.7))

###Divide the dataset and keep only the rows of RANDID in the data in the train_id
train_data <- data_final[data_final$RANDID %in% train_id, ]
test_data_mix <- data_final[!data_final$RANDID %in% train_id, ]

###The test set randomly selects one data for each patient
library(dplyr)
set.seed(2025)
test_data <- test_data_mix %>% group_by(RANDID) %>% slice_sample(n = 1) %>% ungroup()

###check
length(unique(train_data$RANDID))  ###Number of patients in the training set
length(unique(test_data$RANDID))   ###Number of patients in the test set
intersect(train_data$RANDID, test_data$RANDID) ###No overlapping data

################################################################################

##################################
#Simulate missing mechanisms and proportions
##Simulate the absence of the entire MCAR follow-up record
##0.2
set.seed(100)
train_data_mcar1 <- train_data

###Randomly select 20% of the rows as "missing whole rows"
drop_n <- round(0.2 * nrow(train_data_mcar1))
idx_drop <- sample(1:nrow(train_data_mcar1), size = drop_n)

###Delete these lines directly
train_data_mcar1 <- train_data_mcar1[-idx_drop, ]

###See the remaining sample size
nrow(train_data_mcar1)


##0.4
set.seed(200)
train_data_mcar2 <- train_data

###Randomly select 40% of the rows as "whole row missing"
drop_n <- round(0.4 * nrow(train_data_mcar2))
idx_drop <- sample(1:nrow(train_data_mcar2), size = drop_n)

###Delete these lines directly
train_data_mcar2 <- train_data_mcar2[-idx_drop, ]

###See the remaining sample size
nrow(train_data_mcar2)


##0.6
set.seed(300)
train_data_mcar3 <- train_data

###Randomly select 60% of the rows as "whole row missing"
drop_n <- round(0.6 * nrow(train_data_mcar3))
idx_drop <- sample(1:nrow(train_data_mcar3), size = drop_n)

###Delete these lines directly
train_data_mcar3 <- train_data_mcar3[-idx_drop, ]

###See the remaining sample size
nrow(train_data_mcar3)


##0.8
set.seed(400)
train_data_mcar4 <- train_data

###Randomly select 80% of the rows as "whole row missing"
drop_n <- round(0.8 * nrow(train_data_mcar4))
idx_drop <- sample(1:nrow(train_data_mcar4), size = drop_n)

###Delete these lines directly
train_data_mcar4 <- train_data_mcar4[-idx_drop, ]

###See the remaining sample size
nrow(train_data_mcar4)


###############################
##Simulate the absence of MAR
###The probability of missing univariate condition is 60% for smokers
set.seed(500)
train_data_mar1 <- train_data

for (i in 1:nrow(train_data_mar1)) {
  if (train_data_mar1$CURSMOKE[i] == 1) {
    if (runif(1) < 0.6) train_data_mar1$BMI[i] <- NA  ###60% of smokers are missing
  }
  ###Non-smokers do not do the missing work
}

###Check for missing items
cat("Proportion of overall BMI missing:", round(mean(is.na(train_data_mar1$BMI)), 3), "\n")
cat("Smoker BMI missing proportion:", round(mean(is.na(train_data_mar1$BMI[train_data_mar1$CURSMOKE == 1])), 3), "\n")
cat("Proportion of BMI missing in non-smokers:", round(mean(is.na(train_data_mar1$BMI[train_data_mar1$CURSMOKE == 0])), 3), "\n")


##Simulate multivariate MAR missingness
##The probability of BMI loss in elderly diabetic patients is 60%
set.seed(600)
train_data_mar2 <- train_data

for (i in 1:nrow(train_data_mar2)) {
  if (train_data_mar2$AGE[i] > 60 && train_data_mar2$DIABETES[i] == 1) {
    if (runif(1) < 0.6) train_data_mar2$BMI[i] <- NA  ###60% of elderly and diabetic patients are missing
  }
}

###Check for missing items
cat("Proportion of overall BMI missing:", round(mean(is.na(train_data_mar2$BMI)), 3), "\n")
cat("Proportion of older age and diabetes BMI missing:", 
    round(mean(is.na(train_data_mar2$BMI[train_data_mar2$AGE > 60 & train_data_mar2$DIABETES == 1])), 3), "\n")
cat("Proportion of BMI missing in other groups:", 
    round(mean(is.na(train_data_mar2$BMI[!(train_data_mar2$AGE > 60 & train_data_mar2$DIABETES == 1)])), 3), "\n")

###Remove the missing lines of BMI before building the model
train_data_mar2_cc <- train_data_mar2[!is.na(train_data_mar2$BMI), ]

################################################################################
##Random forest model construction
library(randomForest)

###Construct a random forest model
set.seed(123)

#Ensure that the CVD of the training and test sets is factor
train_data$CVD <- as.factor(train_data$CVD)
test_data$CVD <- as.factor(test_data$CVD)

rf_model <- randomForest(
  as.factor(CVD) ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data,
  ntree = 500,
  importance = TRUE
)

###Review model importance
print(rf_model)
importance(rf_model)

###Predict probabilities in the test set
rf_pred <- predict(rf_model, newdata = test_data, type = "prob")[, 2]


#First, two candidate points (days) are given according to "clinically commonly used": 5 years and 15 years
times_days <- round(c(5, 15) * 365.25)

###calibration
library(riskRegression)

# Make sure the CVD is a 0/1 value
if (is.factor(test_data$CVD)) {
  test_data$CVD <- as.numeric(as.character(test_data$CVD))
}

#Generate riskRegression assessment results (AUC/Brier + Calibration)
sc_rf_td <- Score(
  object   = list("RF" = rf_pred),  
  formula  = Surv(stop, CVD) ~ 1,              
  data     = test_data,
  times    = times_days,  
  entry    = "start",
  metrics  = c("auc","brier"),
  plots    = c("ROC","calibration"),
  conf.int = TRUE,
  B        = 200,
  null.model = FALSE
)

##View AUC / Brier
sc_rf_td$AUC
sc_rf_td$Brier


cols <- c("#2E86C1", "#E74C3C") 
plotROC(sc_rf_td, times = times_days[1], col = cols[1], lwd = 2,
        main = "Time-dependent ROC of Logistic model")
plotROC(sc_rf_td, times = times_days[2], add = TRUE, col = cols[2], lwd = 2)

legend("bottomright",
       legend = c("5-year ROC", "15-year ROC"),
       col = cols, lwd = 2)


plotCalibration(sc_rf_td, times = times_days[1], legend = FALSE,
                xlab = "Predicted risk @5y", ylab = "Observed risk @5y", lwd = 2)

plotCalibration(sc_rf_td, times = times_days[2], legend = FALSE,
                xlab = "Predicted risk @15y", ylab = "Observed risk @15y", lwd = 2)

################################
##MSimulate MCAR missingness
###0.2
train_data_mcar1$CVD <- as.factor(train_data_mcar1$CVD)
rf_model1 <- randomForest(
  CVD ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data_mcar1,
  ntree = 500
)
rf_pred1 <- predict(rf_model1, newdata = test_data, type = "prob")[,2]

###0.4
train_data_mcar2$CVD <- as.factor(train_data_mcar2$CVD)
rf_model2 <- randomForest(
  CVD ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data_mcar2,
  ntree = 500
)
rf_pred2 <- predict(rf_model2, newdata = test_data, type = "prob")[,2]

###0.6
train_data_mcar3$CVD <- as.factor(train_data_mcar3$CVD)
rf_model3 <- randomForest(
  CVD ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data_mcar3,
  ntree = 500
)
rf_pred3 <- predict(rf_model3, newdata = test_data, type = "prob")[,2]

###0.8
train_data_mcar4$CVD <- as.factor(train_data_mcar4$CVD)
rf_model4 <- randomForest(
  CVD ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data_mcar4,
  ntree = 500
)
rf_pred4 <- predict(rf_model4, newdata = test_data, type = "prob")[,2]


par(mfrow = c(1, 3))
###calibration

preds <- list(
  "None"    = rf_pred,
  "MCAR-20" = rf_pred1,
  "MCAR-40" = rf_pred2,
  "MCAR-60" = rf_pred3,
  "MCAR-80" = rf_pred4
)

#Scenario-by-scenario evaluation and plotting (per scenario: one overlay ROC + two calibration drawings)
scores <- lapply(names(preds), function(lbl){
  sc <- Score(
    object   = list(lbl = preds[[lbl]]),
    formula  = Surv(stop, CVD) ~ 1,
    data     = test_data,
    times    = times_days,
    entry    = "start",
    metrics  = c("auc","brier"),
    plots    = c("ROC","calibration"),
    conf.int = TRUE, B = 200
  )

  cols <- c("#2E86C1", "#E74C3C")
  plotROC(sc, times = times_days[1], col = cols[1], lwd = 2,
          main = paste0("Time-dependent ROC - ", lbl))
  plotROC(sc, times = times_days[2], add = TRUE, col = cols[2], lwd = 2)
  legend("bottomright", legend = c("5-year ROC", "15-year ROC"),
         col = cols, lwd = 2)
  

  plotCalibration(sc, times = times_days[1], legend = FALSE, lwd = 2,
                  xlab = paste0("Predicted risk @5y - ", lbl),
                  ylab = "Observed risk @5y")
  plotCalibration(sc, times = times_days[2], legend = FALSE, lwd = 2,
                  xlab = paste0("Predicted risk @15y - ", lbl),
                  ylab = "Observed risk @15y")
  
  sc
})
names(scores) <- names(preds)
##View the AUC / Brier for each scene
lapply(scores, function(sc) sc$AUC)
lapply(scores, function(sc) sc$Brier)

#####################################
##MAR simulation is missing
##Simulate univariate MAR: 60% BMI missing for smokers
train_data_mar1_cc <- train_data_mar1[!is.na(train_data_mar1$BMI), ]
train_data_mar1_cc$CVD <- as.factor(train_data_mar1_cc$CVD)
rf_model_mar1 <- randomForest(
  CVD ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data_mar1_cc,
  ntree = 500
)
rf_pred_mar1 <- predict(rf_model_mar1, newdata = test_data, type = "prob")[,2]

##Absence of multivariate conditions (advanced age and 60% absence of diabetes BMI)
train_data_mar2_cc <- train_data_mar2[!is.na(train_data_mar2$BMI), ]
train_data_mar2_cc$CVD <- as.factor(train_data_mar2_cc$CVD)
rf_model_mar2 <- randomForest(
  CVD ~ AGE + SEX + SYSBP + BMI + TOTCHOL + CURSMOKE + DIABETES,
  data = train_data_mar2_cc,
  ntree = 500
)
rf_pred_mar2 <- predict(rf_model_mar2, newdata = test_data, type = "prob")[,2]

###calibration
preds_mar <- list(
  "MAR-Smoke60"     = rf_pred_mar1,
  "MAR-Older&DM60"  = rf_pred_mar2
)

#Scenario-by-scenario evaluation and plotting (per scenario: one overlay ROC (5y+15y) + two calibration drawings)
scores_mar <- lapply(names(preds_mar), function(lbl){
  sc <- Score(
    object   = list(lbl = preds_mar[[lbl]]),
    formula  = Surv(stop, CVD) ~ 1,    
    data     = test_data,
    times    = times_days,
    entry    = "start",
    metrics  = c("auc","brier"),
    plots    = c("ROC","calibration"),
    conf.int = TRUE, B = 200
  )
  
  cols <- c("#2E86C1", "#E74C3C")
  plotROC(sc, times = times_days[1], col = cols[1], lwd = 2,
          main = paste0("Time-dependent ROC - ", lbl))
  plotROC(sc, times = times_days[2], add = TRUE, col = cols[2], lwd = 2)
  legend("bottomright", legend = c("5-year ROC", "15-year ROC"),
         col = cols, lwd = 2)
  
  plotCalibration(sc, times = times_days[1], legend = FALSE, lwd = 2,
                  xlab = paste0("Predicted risk @5y - ", lbl),
                  ylab = "Observed risk @5y")
  plotCalibration(sc, times = times_days[2], legend = FALSE, lwd = 2,
                  xlab = paste0("Predicted risk @15y - ", lbl),
                  ylab = "Observed risk @15y")
  
  sc
})
names(scores_mar) <- names(preds_mar)

##View AUC / Brier
lapply(scores_mar, function(sc) sc$AUC)
lapply(scores_mar, function(sc) sc$Brier)