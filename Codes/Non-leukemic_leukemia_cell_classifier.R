library(Seurat)
library(dplyr)
library(pROC)
library(gelnet)
library(caret)
rm(list=ls())
# set the path
setwd('./')

# read data and labels
use_data <- readRDS( './NT_classifier_data.rds')
use_labels <- readRDS('./NT_classifier_labels.rds')
table(use_labels)
idx <- createDataPartition(use_labels, p = 0.8, list = FALSE)

# generate training and testing datasets
use_data_train <- as.matrix(use_data[,as.vector(idx)])
use_label_train <- use_labels[as.vector(idx)]
use_data_test <- as.matrix(use_data[,as.vector(-idx)])
use_label_test <- use_labels[as.vector(-idx)]

# train the model by classical binary logistic regression
mm <- gelnet( t(use_data_train), use_label_train, 0, 1 )
# model in paper:
# mm <- readRDS('./data/oclr_mm.rds')

# Apply Non-leukemic/leukemia cell classifier on training dataset
print("train")
df_oclr <- apply( use_data_train, 2, function(z) {sum(mm$w*z)+mm$b} )
label_pred <- 1:length(df_oclr)
label_pred[df_oclr>0] <- 'N'
label_pred[df_oclr<0] <- 'T'
label_pred <- as.factor(label_pred)

print("Accuracy:")
sum(label_pred==use_label_train)/length(label_pred)
print("Confusion matrix:")
table(data.frame(truth = use_label_train,pred = label_pred))
F1 <- 2*true_pred[1, 1]/(2*true_pred[1, 1]+true_pred[2, 1]+true_pred[1, 2])
print("F1:")
F1
  
# Apply Non-leukemic/leukemia cell classifier on testing dataset
print("test")
df_oclr <- apply( use_data_test, 2, function(z) {sum(mm$w*z)+mm$b} )
label_pred <- 1:length(df_oclr)
label_pred[df_oclr>0] <- 'N'
label_pred[df_oclr<0] <- 'T'
label_pred <- as.factor(label_pred)

print("Accuracy:")
sum(label_pred==use_label_test)/length(label_pred)
print("Confusion matrix:")
table(data.frame(truth = use_label_test,pred = label_pred))
true_pred <- table(data.frame(truth =use_label_test,pred = label_pred))
F1 <- 2*true_pred[1, 1]/(2*true_pred[1, 1]+true_pred[2, 1]+true_pred[1, 2])
print("F1:")
F1
