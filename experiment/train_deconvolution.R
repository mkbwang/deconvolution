

gene_expressions <- read.csv("experiment/train_data.csv", row.names=1)
library(stringr)

# extract ages
ages <- sapply(colnames(gene_expressions),
               function(cname) str_extract_all(cname, "\\d+")[[1]]) |>
    as.integer()

gene_expressions <- t(gene_expressions) |> as.data.frame()

# retain prevalent genes
prevalences <- colMeans(gene_expressions > 0)
subset_genes <- names(prevalences[prevalences == 1])
gene_expressions <- gene_expressions[, subset_genes]


outcome <- f_test_filter(geneexp_df=gene_expressions,
                         labels=ages,
                         min_label=2,
                         max_label=23)


subset_genes <- outcome$subset_genes
variance_df <- outcome$variance_subset
dendrogram <- hcluster_genes(geneexp_df=gene_expressions,
                             labels=ages,
                             subset_genes=subset_genes,
                             min_label=2,
                             max_label=23)


plot(dendrogram, labels = FALSE, xlab="Genes")

clusters <- cutree(dendrogram, k = 15)
variance_df$Cluster <- clusters

library(dplyr)
variance_df <- variance_df %>% group_by(Cluster) %>%
    arrange(pval)

# select top 2.5% of features
variance_subset_df <- variance_df %>% group_by(Cluster) %>%
    slice_head(prop = 1/40)
write.csv(variance_subset_df, "experiment/selected_genes.csv", row.names=FALSE,
          quote=FALSE)

marker_genes <- variance_subset_df$Gene


ages_subset <- ages[ages >=2 & ages <= 23]
gene_expressions_subset <- gene_expressions[ages >=2 & ages <= 23, marker_genes] |>
    as.matrix()



# first fit deconvolution among the samples with no penalties

fit_deconv <- function(predictors, true_labels, p0=0, p1=0){


    estimated_age <-rep(0, length(true_labels))
    estimated_weights <- matrix(0, nrow=length(true_labels),
                                ncol=length(true_labels))
    for (j in 1:length(true_labels)){

        predictors_train <- predictors[-j, ]
        predictors_test <- predictors[j, ]

        ages_train <- true_labels[-j]
        ages_test <- true_labels[j]

        l2_result <- deconvolution(X=predictors_train, Y=predictors_test,
                                   labels=ages_train, type="l2",
                                   p0=p0, p1=p1)

        estimated_weights[j, -j] <- l2_result$weights
        estimated_age[j] <- l2_result$estim

    }

    mae <- mean(abs(estimated_age - true_labels))
    plot_title_1 <- sprintf("P0 = %d; P1= %d, Mean Absolute Error=%.3f",
                            p0, p1, mae)

    prediction_plot <- viz_predict(truth=true_labels, predicted=estimated_age,
                                   diagonal=T, title=plot_title_1)

    weight_plot_list <- list()
    for (j in 1:nrow(gene_expressions_subset)){

        true_age <- ages_subset[j]
        estimate <- round(estimated_age[j], digits=3)
        plot_title_2 <- sprintf("Truth=%d; Estimated=%.3f; P0 = %d; P1= %d",
                                true_age, estimate, p0, p1)
        weight_plot <- viz_weights(labels=ages_subset, weights=estimated_weights[j, ],
                                   truth=ages_subset[j], exclude=j, title=plot_title_2)
        weight_plot_list[[j]] <- weight_plot

    }


    return(list(result = data.frame(Truth=true_labels, Predicted=estimated_age),
                weight_mat = estimated_weights,
                prediction_plot=prediction_plot,
                weight_plots=weight_plot_list))


}


outcome_nopenalty <- fit_deconv(predictors=gene_expressions_subset,
                                  true_labels=ages_subset,
                                p0=0, p1=0)


outcome_penalty_1 <- fit_deconv(predictors=gene_expressions_subset,
                                  true_labels=ages_subset,
                                p0=50, p1=0)


outcome_penalty_2 <- fit_deconv(predictors=gene_expressions_subset,
                                true_labels=ages_subset,
                                p0=0, p1=50)


outcome_penalty_3 <- fit_deconv(predictors=gene_expressions_subset,
                                true_labels=ages_subset,
                                p0=50, p1=50)




library(ggplot2)
library(patchwork)

# combined_prediction_plot <- (outcome_penalty_0$prediction_plot | outcome_penalty_5$prediction_plot) /
#     (outcome_penalty_10$prediction_plot | outcome_penalty_20$prediction_plot)
#
# ggsave("experiment/plots/predictions_plot.pdf",
#        plot = combined_prediction_plot,
#        width = 8,
#        height = 6,
#        units = "in")


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


