library(readxl)
library(randomForest)
library(pROC)
library(caret)

set.seed(123)

data_matrix <- read_excel("data_matrix.xlsx", col_names = FALSE)
head(data_matrix)

colnames(data_matrix) <- c("t0", "t1", "t3", "t4", "mechProb", "outbreak_label")
head(data_matrix)

n=nrow(data_matrix)
n

# Random forest model - subset

# data_matrix <- data_matrix[1:100000,]
index=sample(1:n, size = 0.1*n, replace = FALSE)
data_matrix=data_matrix[index,]

n=nrow(data_matrix)
n

trainIndex=sample(1:n, size=0.7*n, replace=FALSE)
train=data_matrix[trainIndex,]
test=data_matrix[-trainIndex,]

#trainIndex <- createDataPartition(data_matrix$outbreak_label, p = 0.7, list = FALSE)
#train <- data_matrix[trainIndex,]
#test <- data_matrix[-trainIndex,]

dim(train)
dim(test)
head(train)

rf_model <- randomForest(as.factor(outbreak_label) ~ ., data = train[,-5], 
                         ntree = 700, mtry = 3, importance = TRUE, probability = TRUE)

# Check feature importance
importance(rf_model)

# Predict probabilities
pred_probs <- predict(rf_model, newdata = test[,1:4], type = "prob")[, 2]

# Binarise predictions using a threshold (e.g., 0.5)
pred_labels <- ifelse(pred_probs > 0.5, 1, 0)

# Evaluate performance
confusion_matrix <- table(Predicted = pred_labels, Actual = test$outbreak_label)
print(confusion_matrix)

# Calculate ROC-AUC
roc_curve <- roc(test$outbreak_label, pred_probs)
auc <- auc(roc_curve)
print(auc)
plot(roc_curve)

test[,5]

mechProb <- test$mechProb

mechPred_labels <- ifelse(mechProb > 0.5, 1, 0)

# Evaluate performance
mech_confusion_matrix <- table(Predicted = mechPred_labels, Actual = test$outbreak_label)
print(mech_confusion_matrix)

mech_roc_curve <- roc(test$outbreak_label, mechProb)
mech_auc <- auc(mech_roc_curve)
print(mech_auc)
plot(mech_roc_curve)


library(ggplot2)
comparison_df <- data.frame(
  Classifier_Prob = pred_probs,
  External_Prob = test$mechProb,
  Label = test$outbreak_label
)

ggplot(comparison_df, aes(x = Classifier_Prob, y = External_Prob, colour = as.factor(Label))) +
  geom_point(alpha = 0.5) +
  labs(title = "Probability Comparison", x = "Classifier Probability", y = "External Probability") +
  theme_minimal()

plot(roc_curve, col = "blue", main = "ROC Curve Comparison", lwd = 2)
plot(mech_roc_curve, col = "red", add = TRUE, lwd = 2)
legend("bottomright", legend = c("Classifier", "External (mechProb)"), col = c("blue", "red"), lwd = 2)
