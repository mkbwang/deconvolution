


rm(list=ls())
gene_expressions <- read.csv("experiment/train_data.csv", row.names=1)
metadata <- read.csv("experiment/metadata.csv", row.names=1)
library(stringr)
library(dplyr)

# extract ages
ages <- sapply(colnames(gene_expressions),
               function(cname) str_extract_all(cname, "\\d+")[[1]]) |>
    as.integer()

gene_expressions <- t(gene_expressions) |> as.data.frame()
age_filter <- ages >= 2 & ages <= 23
gene_expressions <- gene_expressions[age_filter, ]
ages <- ages[age_filter]

# load marker genes
markergene_info <- read.csv("experiment/markergene_anova.csv")
markergene_info_subset <- markergene_info %>% group_by(Cluster) %>% arrange(pval) %>%
    slice_head(prop=0.3)
markergene <- markergene_info_subset$Gene
gene_expressions <- gene_expressions[, markergene]


unique_ages <- unique(ages)

weights_mat <- matrix(0, nrow=nrow(gene_expressions), ncol=length(unique_ages))
rownames(weights_mat) <- rownames(gene_expressions)
colnames(weights_mat) <- sprintf("Age_%d", unique_ages)


for (j in 1:nrow(gene_expressions)){
    X_expressions <- gene_expressions[-j, ]
    X_expressions$Age <- ages[-j]
    X_expressions_median <- X_expressions %>% group_by(Age) %>% summarise_all(median)
    X_expressions$Age <- NULL
    X_expressions_median$Age <- NULL


    Y_expressions <- gene_expressions[j, ] |> as.matrix() |> as.vector()

    normalized_data <- preprocess(X=X_expressions_median, Y=Y_expressions)
    normalized_X <- normalized_data$normalized_X
    normalized_Y <- normalized_data$normalized_Y

    deconv_result <- l2_solve(X=normalized_X, Y=normalized_Y,
                              labels=unique_ages)

    weights_mat[j, ] <- deconv_result$weights
}



# Bayesian fitting

X_expressions <- gene_expressions[-5, ]
X_expressions$Age <- ages[-5]
X_expressions_median <- X_expressions %>% group_by(Age) %>% summarise_all(median)
X_expressions$Age <- NULL
X_expressions_median$Age <- NULL


Y_expressions <- gene_expressions[5, ] |> as.matrix() |> as.vector()

normalized_data <- preprocess(X=X_expressions_median, Y=Y_expressions)
normalized_X <- normalized_data$normalized_X
normalized_Y <- normalized_data$normalized_Y

bayes_deconv_result <- deconv_bayes_2(X=normalized_X, Y=normalized_Y,
                                      n.chains=4, n.iter=3e4, n.burnin=1.5e4, n.thin=4)

bayes_deconv_result$gelman
bayes_deconv_result$tau1
bayes_deconv_result$tau2
bayes_deconv_result$gamma

hist(bayes_deconv_result$Z)

# load test data



# visualization
library(ggplot2)
library(patchwork)




combined_weight_plots <- list()

for (j in 1:nrow(gene_expressions_subset)){

    y_max <- max(c(outcome_nopenalty$weight_mat[j, ],
                   outcome_penalty_1$weight_mat[j, ],
                   outcome_penalty_2$weight_mat[j, ],
                   outcome_penalty_3$weight_mat[j, ]))

    plot_0 <- outcome_nopenalty$weight_plots[[j]] + ylim(0, y_max+0.05)
    plot_1 <- outcome_penalty_1$weight_plots[[j]] + ylim(0, y_max+0.05)
    plot_2 <- outcome_penalty_2$weight_plots[[j]] + ylim(0, y_max+0.05)
    plot_3 <- outcome_penalty_3$weight_plots[[j]] + ylim(0, y_max+0.05)

    combined_weight_plots[[j]] <- (plot_0 | plot_1) /
        (plot_2 | plot_3)

    output_file <- sprintf("experiment/plots_new/weights_%d.pdf", j)

    ggsave(output_file,
           plot = combined_weight_plots[[j]],
           width = 10,
           height = 6,
           units = "in")

}


