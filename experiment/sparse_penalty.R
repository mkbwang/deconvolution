

library(dplyr)
rm(list=ls())

data <- read.csv("experiment/train_data_subset.csv")
ncols <- ncol(data)




fit_deconv <- function(data, p0=0, p1=0, p2=0){

    nsamples <- nrow(data)
    nfeatures <- ncol(data)

    unique_ages <- unique(data$Age)
    true_ages <- data$Age
    weights_mat <- matrix(0, nrow=nsamples, ncol=length(unique_ages))
    predicted_ages <- rep(0, nsamples)

    for(j in 1:nsamples){
        target_sample <- data[j, -nfeatures] |> as.matrix() |> as.vector()
        log_target <- log(target_sample)
        reference_samples <- data[-j, ]
        reference_samples[, -nfeatures] <- log(reference_samples[, -nfeatures])
        reference_samples <- reference_samples %>% group_by(Age) %>% summarise_all(mean)
        ages <- reference_samples$Age
        log_reference <- reference_samples[, -1] |> as.matrix()

        scaled_log_reference <- scale(log_reference)

        mean_log_reference <- colMeans(log_reference)
        sd_log_reference <- apply(log_reference, 2, sd)

        scaled_log_target <- (log_target - mean_log_reference) / sd_log_reference

        result <- deconvolution_BFGS(X=scaled_log_reference, y=scaled_log_target, labels=ages,
                                     p0=p0, p1=p1, p2=p2)
        weights_mat[j, ] <- result$weights
        predicted_ages[j] <- result$estimate
    }

    # plot prediction
    mae <- mean(abs(predicted_ages - true_ages))
    plot_title_1 <- sprintf("P1 = %d; P2= %d, Mean Absolute Error=%.3f",
                            p1, p2, mae)
    prediction_plot <- viz_predict(truth=true_ages, predicted=predicted_ages,
                                   diagonal=T, title=plot_title_1)

    # plot weights
    weight_plot_list <- list()
    for (j in 1:nsamples){

        true_age <- true_ages[j]
        estimate <- round(predicted_ages[j], digits=3)
        plot_title_2 <- sprintf("Truth=%d; Estimated=%.3f; P1 = %d; P2= %d",
                                true_age, estimate, p1, p2)
        weight_plot <- viz_weights(labels=ages, weights=weights_mat[j, ], title=plot_title_2)
        weight_plot_list[[j]] <- weight_plot

    }

    return(list(result = data.frame(Truth=true_ages, Predicted=predicted_ages),
                weight_mat = weights_mat,
                prediction_plot=prediction_plot,
                weight_plots=weight_plot_list))

}

results_0 <- fit_deconv(data=data, p0=0, p1=0, p2=0)

results_10 <- fit_deconv(data=data, p0=0, p1=0, p2=10)

results_20 <- fit_deconv(data=data, p0=0, p1=0, p2=20)

results_30 <- fit_deconv(data=data, p0=0, p1=0, p2=30)


library(ggplot2)
library(patchwork)
combined_weight_plots <- list()

for (j in 1:nrow(data)){

    y_max <- max(c(results_0$weight_mat[j, ],
                   results_10$weight_mat[j, ],
                   results_20$weight_mat[j, ],
                   results_30$weight_mat[j, ]))

    plot_0 <- results_0$weight_plots[[j]] + ylim(0, y_max+0.05)
    plot_1 <- results_10$weight_plots[[j]] + ylim(0, y_max+0.05)
    plot_2 <- results_20$weight_plots[[j]] + ylim(0, y_max+0.05)
    plot_3 <- results_30$weight_plots[[j]] + ylim(0, y_max+0.05)

    combined_weight_plots[[j]] <- (plot_0 | plot_1) /
        (plot_2 | plot_3)

    output_file <- sprintf("experiment/plots_penalty/weights_%d.pdf", j)

    ggsave(output_file,
           plot = combined_weight_plots[[j]],
           width = 10,
           height = 6,
           units = "in")

}



