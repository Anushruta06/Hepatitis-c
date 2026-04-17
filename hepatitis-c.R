#----library needed----
library(ggplot2)
library(caTools)
library(VGAM)
library(dplyr)
library(reshape2)
#----setting dataset----

setwd("C:/Users/aparu/OneDrive/Desktop/Hepatitis C")
raw_data <- read.csv("hcvdat0.csv")

data <- na.omit(raw_data)
data <- data[ ,-c(1,4)]
n_runs <- 1000
acc <- numeric(n_runs)
b_acc <- numeric(n_runs)
f1_all <- matrix(0, nrow = n_runs, ncol = 4)
sens_all <- matrix(0, nrow = n_runs, ncol = 4)
spec_all <- matrix(0, nrow = n_runs, ncol = 4)


#deal with categories

table(data$Category)
prop.table(table(data$Category))

data$Category[data$Category == "0s=suspect Blood Donor"] <- "0=Blood Donor"
data$Category_num <- as.numeric(factor(data$Category,
                                        levels = c("0=Blood Donor","1=Hepatitis","2=Fibrosis","3=Cirrhosis"))) - 1

data<- data[ ,-1]
data[,-12]<- scale(data[,-12])

for(i in 1:n_runs){
set.seed(273)

split <- sample.split(data$Category_num, SplitRatio = 0.75)
X_tr <- subset(data, split == "TRUE")
X_ts <- subset(data, split == "FALSE")

X_balanced <- X_tr %>%
  group_by(Category_num) %>%
  sample_n(min(table(X_tr$Category_num)))
Y <- data[ c(554,556,557,558,559,560,561,562,563,564,565),]

#----fitting model----
class_counts <- table(X_tr$Category_num)
weights <- as.numeric(1 / class_counts[as.character(X_tr$Category_num)])
weights <- weights / mean(weights)


model <- vglm(Category_num ~ Age + ALB + ALT + AST + BIL + CHE + CHOL + CREA + GGT + ALP + PROT,
              family = cumulative(parallel = TRUE),
              data = X_tr,
              weights = weights)

summary(model)


model_fixed <- vglm(Category_num ~Age+ ALB + ALT + AST + BIL +CHOL + CREA +PROT,
                    family = cumulative(parallel = TRUE),
                    data = X_tr,
                    weights = weights)

summary(model_fixed)

prob <- predict(model, newdata = X_tr, type = "response")
pred <- apply(prob, 1, which.max) - 1
table(True = X_tr$Category_num, Predicted = pred)
mean(pred == X_tr$Category_num)
head(prob)

prob_new <- predict(model_fixed, newdata = X_tr, type = "response")
pred_new <- apply(prob_new, 1, which.max) - 1
table(True = X_tr$Category_num, Predicted = pred_new)
mean(pred_new == X_tr$Category_num)
head(prob_new)


#----test the model----

# prob_test <- predict(model, newdata = X_ts, type = "response")
# pred_test <- apply(prob_test, 1, which.max) - 1
# table(True = X_ts$Category_num, Predicted = pred_test)
# mean(pred_test == X_ts$Category_num)
# head(prob_test)

prob_test <- predict(model_fixed, newdata = X_ts, type = "response")
pred_test <- apply(prob_test, 1, which.max) - 1
table(True = X_ts$Category_num, Predicted = pred_test)
mean(pred_test == X_ts$Category_num)
head(prob_test)

# prob_test<- predict(model_fixed, newdata = X_balanced, type = "response")
# pred_test <- apply(prob_test, 1, which.max) - 1
# table(True = X_balanced$Category_num, Predicted = pred_test)
# mean(pred_test == X_balanced$Category_num)
# head(prob_test)
# 
# 
# prob_test <- predict(model_fixed, newdata = Y, type = "response")
# pred_test <- apply(prob_test, 1, which.max) - 1
# table(True = Y$Category_num, Predicted = pred_test)
# mean(pred_test == Y$Category_num)
# head(prob_test)
# # 
# X_ts_last10 <- tail(X_ts, 10)
# prob_test <- predict(model_fixed, newdata = X_ts_last10, type = "response")
# pred_test <- apply(prob_test, 1, which.max) - 1
# table(True = X_ts_last10$Category_num, Predicted = pred_test)
# mean(pred_test == X_ts_last10$Category_num)
# head(prob_test)
# 
# data_last10 <- tail(data, 56)
# prob_test <- predict(model_fixed, newdata = data_last10, type = "response")
# pred_test <- apply(prob_test, 1, which.max) - 1
# table(True = data_last10$Category_num, Predicted = pred_test)
# mean(pred_test == data_last10$Category_num)
# head(prob_test)

cm <- table(
  factor(X_ts$Category_num, levels = 0:3),
  factor(pred_test, levels = 0:3)
)
# ensure all classes exist
cm <- as.matrix(cm)
classes <- 0:3

for(k in classes){
  
  TP <- cm[as.character(k), as.character(k)]
  
  FN <- sum(cm[as.character(k), ]) - TP
  FP <- sum(cm[, as.character(k)]) - TP
  TN <- sum(cm) - TP - FN - FP
  
  sens_all[i, k+1] <- TP / (TP + FN)
  spec_all[i, k+1] <- TN / (TN + FP)
  
}



# Accuracy
acc[i] <- mean(pred_test == X_ts$Category_num)

# ---- F1 SCORE CALCULATION ----

f1_per_class <- numeric(length(classes))

for (k in classes) {
  
  TP <- cm[as.character(k), as.character(k)]
  FN <- sum(cm[as.character(k), ]) - TP
  FP <- sum(cm[, as.character(k)]) - TP
  
  precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  recall    <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  
  f1_per_class[k + 1] <- ifelse(
    (precision + recall) == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
}

# store per run (optional but recommended)
# you should create this before loop:
# f1_all <- matrix(0, nrow = n_runs, ncol = 4)

f1_all[i, ] <- f1_per_class
}

# Results
mean(acc)   # average accuracy
sd(acc)     # variability
min(acc)    # worst case
max(acc)    # best case

ci_lower <- mean(acc) - 1.96 * sd(acc) / sqrt(n_runs)
ci_upper <- mean(acc) + 1.96 * sd(acc) / sqrt(n_runs)

c(ci_lower, ci_upper) #CI


colMeans(sens_all)
colMeans(spec_all)


tab <- table(True = X_ts$Category_num, Predicted = pred_test)
sensitivity <- diag(tab) / rowSums(tab)
sensitivity

specificity <- sapply(1:nrow(tab), function(i) {
  TP <- tab[i, i]
  FN <- sum(tab[i, ]) - TP
  FP <- sum(tab[, i]) - TP
  TN <- sum(tab) - TP - FN - FP
  TN / (TN + FP)
})
specificity
f1_all[i, ]
colMeans(f1_all)

balanced_acc <- rowMeans(sens_all)
mean_balanced_acc <- mean(balanced_acc)
cm_prop <- prop.table(cm, margin = 1)
cm_df <- melt(cm_prop)
colnames(cm_df) <- c("True", "Predicted", "Value")
ggplot(cm_df, aes(x = Predicted, y = True, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), size = 5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Confusion Matrix Heatmap",
       x = "Predicted Class",
       y = "True Class") +
  theme_minimal()
