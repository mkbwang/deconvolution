
rm(list=ls())
gene_expressions <- read.csv("data/train_data.csv", row.names=1)
metadata <- read.csv("data/metadata.csv", row.names=1)
library(stringr)
library(dplyr)
library(DWLS)

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
  signature_mat <- as.matrix(log(X_expressions_mean)) |> t()
  
  Y_expressions <- gene_expressions[j, ] |> as.matrix() |> as.vector()
  log_Y <- log(Y_expressions)
  
  solDWLS <- solveDampenedWLS(signature_mat, log_Y)
  # normalized_data <- preprocess(X=X_expressions_median, Y=Y_expressions)
  # normalized_X <- normalized_data$normalized_X
  # normalized_Y <- normalized_data$normalized_Y
  
  # deconv_result <- l2_solve(X=normalized_X, Y=normalized_Y,
  #                           labels=unique_ages)
  
  weights_mat_train[j, ] <- solDWLS
  predicted_train[j] <- sum(solDWLS * unique_ages)
  
}

train_prediction_df <- data.frame(Age=ages,
                                  Prediction=predicted_train)
weights_train_df <- as.data.frame(weights_mat_train)

write.csv(train_prediction_df, "planarian/train_prediction_DWLS.csv",
          row.names=T)
write.csv(weights_train_df, "planarian/train_coefficients_DWLS.csv",
          row.names=T)




