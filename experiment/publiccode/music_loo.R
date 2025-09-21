
rm(list=ls())
gene_expressions <- read.csv("data/train_data.csv", row.names=1)
metadata <- read.csv("data/metadata.csv", row.names=1)
library(stringr)
library(dplyr)
library(MuSiC)
library(deconvolution)
library(SingleCellExperiment)

# extract ages
ages <- sapply(colnames(gene_expressions),
               function(cname) str_extract_all(cname, "\\d+")[[1]]) |>
  as.integer()

gene_expressions <- t(gene_expressions) |> as.data.frame()
age_filter <- ages >= 2 & ages <= 23
gene_expressions <- gene_expressions[age_filter, ]
ages <- ages[age_filter]


# load marker genes
markergene_info <- read.csv("markergene/markergene_anova.csv")
# markergene_info_subset <- markergene_info %>% group_by(Cluster) %>% arrange(pval) %>%
#     slice_head(prop=0.3)
markergene <- markergene_info$Gene
gene_expressions <- gene_expressions[, markergene]


unique_ages <- unique(ages)



weights_mat_train <- matrix(0, nrow=nrow(gene_expressions), ncol=length(unique_ages))
predicted_train <- rep(0, nrow(gene_expressions))
rownames(weights_mat_train) <- rownames(gene_expressions)
colnames(weights_mat_train) <- sprintf("Age_%d", unique_ages)


for (j in 1:nrow(gene_expressions)){
  print(j)
  X_expressions <- gene_expressions[-j, ]
  X_expressions$Age <- ages[-j]
  X_expressions_mean <- X_expressions %>% group_by(Age) %>% summarise_all(mean)
  X_expressions$Age <- NULL
  X_expressions_mean$Age <- NULL
  
  Y_expressions <- gene_expressions[j, ] |> as.matrix() |> as.vector()
  
  # normalized_data <- preprocess(X=X_expressions_mean, Y=Y_expressions)
  # normalized_X <- normalized_data$normalized_X
  # normalized_Y <- normalized_data$normalized_Y
  
  signature_mat <- as.matrix(log(X_expressions_mean)) |> t()
  colnames(signature_mat) <- sprintf("Age_%d", unique_ages)
  rownames(signature_mat) <- colnames(X_expressions_mean)
  
  log_Y <- log(Y_expressions)
  names(log_Y) <- colnames(X_expressions_mean)
  
  size_vec <- rep(1, length(unique_ages))
  names(size_vec) <- sprintf("Age_%d", unique_ages)
  Sigma_mat <- matrix(0, nrow=nrow(signature_mat), ncol=ncol(signature_mat))
  rownames(Sigma_mat) <- rownames(signature_mat)
  colnames(Sigma_mat) <- colnames(signature_mat)
  
  
  solmusic <- music.iter(Y=log_Y, D=signature_mat,
                         S=size_vec, Sigma=Sigma_mat)
  
  
  # deconv_result <- l2_solve(X=normalized_X, Y=normalized_Y,
  #                           labels=unique_ages)
  
  weights_mat_train[j, ] <- solmusic$p.weight
  predicted_train[j] <- sum(solmusic$p.weight * unique_ages)
  
}

train_prediction_df <- data.frame(Age=ages,
                                  Prediction=predicted_train)
weights_train_df <- as.data.frame(weights_mat_train)

write.csv(train_prediction_df, "planarian/train_prediction_MuSiC.csv",
          row.names=T)
write.csv(weights_train_df, "planarian/train_coefficients_MuSiC.csv",
          row.names=T)



# col_metadata <- data.frame(Sample="Sample1",
#                            Age=sprintf("Age_%d", ages))
# rownames(col_metadata) <- rownames(gene_expressions)
# cellsize_df <- data.frame(Age=sprintf("Age_%d", unique(ages)),
#                           Size=1)
# geneexp_experiment <- SingleCellExperiment(
#   assays=list(counts=t(gene_expressions)),
#   colData=col_metadata
# )
# 
# 
# data_summary <- music_basis(x=geneexp_experiment,
#                             clusters="Age", samples="Sample",
#                             cell_size=cellsize_df)
# clusters <- cellsize_df$Age
# samples <- rep("Sample1", nrow(cellsize_df))
# S <- sapply(unique(clusters), function(ct){
#   my.rowMeans(sapply(unique(samples), function(sid){
#     y = counts(geneexp_experiment)[, clusters %in% ct & samples %in% sid]
#     if(is.null(dim(y))){
#       return(sum(y))
#     }else{
#       return(sum(y)/ncol(y))
#     }
#   }), na.rm = TRUE)
# })
# 
# 
# Sigma <- sapply(unique(clusters), function(ct){
#   apply(sapply(unique(samples), function(sid){
#     y = counts(geneexp_experiment)[,clusters %in% ct & samples %in% sid]
#     if(is.null(dim(y))){
#       return(y/sum(y))
#     }else{
#       return(rowSums(y)/sum(y))
#     }
#   }), 1, mean, na.rm = TRUE)
# })
# 
# 
# values <- sapply(unique(samples), function(sid){
#   y = counts(geneexp_experiment)[,clusters %in% c("Age_2") & samples %in% sid]
#   if(is.null(dim(y))){
#     return(y/sum(y))
#   }else{
#     return(rowSums(y)/sum(y))
#   }
# })




