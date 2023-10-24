rm(list = ls())
library(Biostrings)
library(DNAshapeR)
library(devtools)
library(caret)
library(glmnet)
setwd("/Users/anikmitra/Documents/Documents/Documents/Grad School/USC/Coursework & Teaching/Fall 2023/QBIO-481/optional_HW/QBIO481-master/gcPBM/")

#######################

MadBinding <- read.table("Mad.txt", header = FALSE)

fv_1mer_Mad <- encodeSeqShape(fastaFileName = "Mad.txt.fa",featureNames = "1-mer")
colnames(fv_1mer_Mad) <- c(sprintf("n%02d", seq(1,ncol(fv_1mer_Mad))))

pred_shape_Mad <- getShape("Mad.txt.fa")

fv_1mer_shape_Mad <- encodeSeqShape(fastaFileName = "Mad.txt.fa",shapeMatrix = pred_shape_Mad, featureNames = "1-shape")
colnames(fv_1mer_shape_Mad) <- c(sprintf("n%02d", seq(1,ncol(fv_1mer_shape_Mad))))

set.seed(1234)
trainControlParameter <- trainControl(method = "cv", number = 10)
Model_1mer_Mad <- train(x = fv_1mer_Mad,
                        y = MadBinding[,2],
                        method = "glmnet",
                        trControl = trainControlParameter,
                        tuneGrid = expand.grid(alpha = 0, lambda = seq(0, 0.01, length = 10)))
print(Model_1mer_Mad)
rsquare_1mer_Mad <- max(Model_1mer_Mad$result)
print(paste("Rsquared_1mer_Mad:", rsquare_1mer_Mad))

set.seed(1234)
trainControlParameter <- trainControl(method = "cv", number = 10)
Model_1mer_shape_Mad <- train(x = fv_1mer_shape_Mad,
                        y = MadBinding[,2],
                        method = "glmnet",
                        trControl = trainControlParameter,
                        tuneGrid = expand.grid(alpha = 0, lambda = seq(0, 0.01, length = 10)))
print(Model_1mer_shape_Mad)
rsquare_1mer_shape_Mad <- max(Model_1mer_shape_Mad$result)
print(paste("Rsquared_1mer_shape_Mad:", rsquare_1mer_shape_Mad))

#########################

MaxBinding <- read.table("Max.txt", header = FALSE)

fv_1mer_Max <- encodeSeqShape(fastaFileName = "Max.txt.fa",featureNames = "1-mer")
colnames(fv_1mer_Max) <- c(sprintf("n%02d", seq(1,ncol(fv_1mer_Max))))

pred_shape_Max <- getShape("Max.txt.fa")

fv_1mer_shape_Max <- encodeSeqShape(fastaFileName = "Max.txt.fa",shapeMatrix = pred_shape_Max, featureNames = "1-shape")
colnames(fv_1mer_shape_Max) <- c(sprintf("n%02d", seq(1,ncol(fv_1mer_shape_Max))))

set.seed(1234)
trainControlParameter <- trainControl(method = "cv", number = 10)
Model_1mer_Max <- train(x = fv_1mer_Max,
                        y = MaxBinding[,2],
                        method = "glmnet",
                        trControl = trainControlParameter,
                        tuneGrid = expand.grid(alpha = 0, lambda = seq(0, 0.01, length = 10)))
print(Model_1mer_Max)
rsquare_1mer_Max <- max(Model_1mer_Max$result)
print(paste("Rsquared_1mer_Max:", rsquare_1mer_Max))

set.seed(1234)
trainControlParameter <- trainControl(method = "cv", number = 10)
Model_1mer_shape_Max <- train(x = fv_1mer_shape_Max,
                              y = MaxBinding[,2],
                              method = "glmnet",
                              trControl = trainControlParameter,
                              tuneGrid = expand.grid(alpha = 0, lambda = seq(0, 0.01, length = 10)))
print(Model_1mer_shape_Max)
rsquare_1mer_shape_Max <- max(Model_1mer_shape_Max$result)
print(paste("Rsquared_1mer_shape_Max:", rsquare_1mer_shape_Max))

#######################

MycBinding <- read.table("Myc.txt", header = FALSE)

fv_1mer_Myc <- encodeSeqShape(fastaFileName = "Myc.txt.fa",featureNames = "1-mer")
colnames(fv_1mer_Myc) <- c(sprintf("n%02d", seq(1,ncol(fv_1mer_Myc))))

pred_shape_Myc <- getShape("Myc.txt.fa")

fv_1mer_shape_Myc <- encodeSeqShape(fastaFileName = "Myc.txt.fa",shapeMatrix = pred_shape_Myc, featureNames = "1-shape")
colnames(fv_1mer_shape_Myc) <- c(sprintf("n%02d", seq(1,ncol(fv_1mer_shape_Myc))))

set.seed(1234)
trainControlParameter <- trainControl(method = "cv", number = 10)
Model_1mer_Myc <- train(x = fv_1mer_Myc,
                        y = MycBinding[,2],
                        method = "glmnet",
                        trControl = trainControlParameter,
                        tuneGrid = expand.grid(alpha = 0, lambda = seq(0, 0.01, length = 10)))
print(Model_1mer_Myc)
rsquare_1mer_Myc <- max(Model_1mer_Myc$result)
print(paste("Rsquared_1mer_Myc:", rsquare_1mer_Myc))

set.seed(1234)
trainControlParameter <- trainControl(method = "cv", number = 10)
Model_1mer_shape_Myc <- train(x = fv_1mer_shape_Myc,
                              y = MycBinding[,2],
                              method = "glmnet",
                              trControl = trainControlParameter,
                              tuneGrid = expand.grid(alpha = 0, lambda = seq(0, 0.01, length = 10)))
print(Model_1mer_shape_Myc)
rsquare_1mer_shape_Myc <- max(Model_1mer_shape_Myc$result)
print(paste("Rsquared_1mer_shape_Myc:", rsquare_1mer_shape_Myc))

Rsquared_values <- data.frame(models=c("1mer_Mad","1mer_Max","1mer_Myc","1mer+shape_Mad","1mer_shape_Max","1mer_shape_Myc"),
                              rsquared=c(rsquare_1mer_Mad,rsquare_1mer_Max,rsquare_1mer_Myc,rsquare_1mer_shape_Mad,rsquare_1mer_shape_Max,rsquare_1mer_shape_Myc))

write.table(Rsquared_values, file="Q4_R_squared.txt",row.names = FALSE, col.names = TRUE,fileEncoding = "UTF-8")


#Q7

library(fields)
setwd("/Users/anikmitra/Documents/Documents/Documents/Grad School/USC/Coursework & Teaching/Fall 2023/QBIO-481/optional_HW/QBIO481-master/CTCF/")
pred_shape_bound <- getShape("bound.fa")
pred_shape_unbound <- getShape("unbound.fa")

plotShape(pred_shape_bound$MGW)
plotShape(pred_shape_bound$ProT)
plotShape(pred_shape_bound$Roll)
plotShape(pred_shape_bound$HelT)

plotShape(pred_shape_unbound$MGW)
plotShape(pred_shape_unbound$ProT)
plotShape(pred_shape_unbound$Roll)
plotShape(pred_shape_unbound$HelT)

heatShape(pred_shape_bound$MGW, 30)
heatShape(pred_shape_bound$ProT, 30)
heatShape(pred_shape_bound$Roll, 30)
heatShape(pred_shape_bound$HelT, 30)

heatShape(pred_shape_unbound$MGW, 30)
heatShape(pred_shape_unbound$ProT, 30)
heatShape(pred_shape_unbound$Roll, 30)
heatShape(pred_shape_unbound$HelT, 30)
