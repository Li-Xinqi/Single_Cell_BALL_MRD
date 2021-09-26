library(Seurat)
library(caret)
library(gelnet)
rm(list=ls())

# set the path
dir <- './'
setwd(dir)

# read Seurat object of H1
H1.obj <- readRDS('./H1.rds' )

# 479 feature genes
markers <- readRDS('./healthy_bone_markers.rds')
top75 <- markers %>% group_by(cell_type) %>% top_n(n = 75, wt = avg_logFC)
use_features_unique <- unique(top75[top75$cell_type!='cycling',]$gene)
use_features <- setdiff(use_features_unique,c(cc.genes[[1]], cc.genes[[2]]))
length(use_features)

# data and labels 
labels <- H1.obj$label
labels <- factor(labels, levels = c('HSC/LMPP','CLP','proB','preBI','preBII','immatureB','matureB','activatedB','cycling'))
use_data <- GetAssayData(object = H1.obj, slot = 'data')[use_features,]
use_data <- use_data[,labels!='cycling']
use_labels <- labels[labels!='cycling']

# generate training and testing datasets
idx <- createDataPartition(use_labels, p = 0.8, list = FALSE)
use_data_train <- use_data[,as.vector(idx)]
use_label_train <- use_labels[as.vector(idx)]
use_data_test <- use_data[,as.vector(-idx)]
use_label_test <- use_labels[as.vector(-idx)]

# training
for (title in unique(use_label_train)){
  use_train_stage <- as.matrix(use_data_train[,use_label_train==title])
  # train the model by one-class logistic regression classifier
  mm <- gelnet( t(use_train_stage), NULL, 0, 1 )
  df_oclr <- data.frame(w = mm$w)
  if (title=='HSC/LMPP'){
    title <- 'HSC-LMPP'
  }
  rownames(df_oclr) <- use_features
  write.table(df_oclr,paste0(dir, '/data/',title,'_marker_train.rds'),row.names = T, col.names =F)
  
}

# testing
df_oclr <- as.data.frame(matrix(NA,nrow = ncol(use_data_test), ncol = 8))
titles <- c('HSC/LMPP','CLP','proB','preBI','preBII','immatureB','matureB','activatedB')
colnames(df_oclr) <- titles

for (title in titles){
  title_file <- title
  if (title=='HSC/LMPP'){
    title_file <- 'HSC-LMPP'
  }
  # model in paper:
  # w_df <- read.table(paste0(dir, '/data/',title_file,'_marker.rds'),row.names = 1, header = F)
  w_df <- read.table(paste0(dir, '/data/',title_file,'_marker_train.rds'),row.names = 1, header = F)
  use_data <- use_data_test[use_features,]
  w <- w_df[use_features,1]
  df_oclr[[title]] <- apply( use_data, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
  
}

colname <- colnames(df_oclr)
pred_type <- apply(df_oclr,1,function(z){colname[which.max(z)]})
pred_type <- factor(pred_type, levels = c('HSC/LMPP','CLP','proB','preBI','preBII','immatureB','matureB','activatedB'))

# testing data of confusion matrix
table(data.frame(truth = use_label_test,pred = pred_type))


##########################################################################
# test on H2
H2.obj <- readRDS('./H2.rds')
data <- GetAssayData(object = H2.obj,slot = 'data')
df_oclr <- as.data.frame(matrix(NA,nrow = dim(data)[2], ncol = 8))

titles <- c('HSC/LMPP','CLP','proB','preBI','preBII','immatureB','matureB','activatedB')
colnames(df_oclr) <- titles

for (title in titles){
  title_file <- title
  if (title=='HSC/LMPP'){
    title_file <- 'HSC-LMPP'
  }
  genes <-rownames(H2.obj)
  
  # model in paper:
  # w_df <- read.table(paste0(dir, '/data/',title_file,'_marker.rds'),row.names = 1, header = F)
  w_df <- read.table(paste0(dir, '/data/',title_file,'_marker_train.rds'),row.names = 1, header = F)
  
  use_genes <- intersect(genes,rownames(w_df))
  use_data <- data[use_genes,]
  w <- w_df[use_genes,1]
  
  df_oclr[[title]] <- apply( use_data, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )})
  
}

colname <- colnames(df_oclr)
pred_type <- apply(df_oclr,1,function(z){colname[which.max(z)]})
pred_type <- factor(pred_type, levels = c('HSC/LMPP','CLP','proB','preBI','preBII','immatureB','matureB','activatedB'))
H2.obj <- AddMetaData(H2.obj,pred_type,'pred')

# UMAP showing the predicting results
use_stage_colors <- c("#EB5B55","#EF833B","#FEEB4F","#3A8345","#4CA06E","#7BCDDF","#2771A1","#992E5D","#764D21")
p1 <- DimPlot(H2.obj,label = T,group.by = 'pred')+scale_color_manual(values = use_stage_colors)
print(p1)

