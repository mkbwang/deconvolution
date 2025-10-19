



fit_deconv <- function(train_data, test_data, train_labels, test_label,
                       min_scale=0.01, samplename){

  lambda_ranges <- c(0, 2000, 4000)

  # tuning different smoothness penalty and sparseness penalty
  ## first sparseness
  # predicted_age_sparseness <- data.frame(Shrinkage=lambda_ranges,
  #                                        Predicted=0)

  cat("Tune sparseness penalty...")
  fitted_weights_sparseness <- list()
  for (j in 1:length(lambda_ranges)){
    deconv_result <- iter_deconv(X=train_data,
                                 Y=test_data,
                                 labels=train_labels,
                                 lambda1=0, lambda2=lambda_ranges[j], log=FALSE,
                                 min_scale=min_scale,
                                 max_iter=20)
    weights_df <- data.frame(Age=train_labels, Weight=deconv_result$weights,
                             normalized_weights=deconv_result$normalized_weights,
                             Shrinkage=lambda_ranges[j])
    weights_df$Predicted <- deconv_result$estim
    fitted_weights_sparseness[[j]] <- weights_df
  }
  fitted_weights_sparseness_df <- do.call(rbind, fitted_weights_sparseness)
  predicted_age_sparseness <- unique(fitted_weights_sparseness_df[, c("Shrinkage", "Predicted")])
  fitted_weights_sparseness_df$Shrinkage <- factor(fitted_weights_sparseness_df$Shrinkage,
                                                   levels=lambda_ranges)

  prediction_plot_sparseness <- ggplot(predicted_age_sparseness, aes(x=Shrinkage, y=Predicted)) +
    geom_point(size=1.3) + geom_line() + xlab("Shrinkage") + ylab("Prediction") +
    scale_x_continuous(breaks=predicted_age_sparseness$Shrinkage) +
    scale_y_continuous(limits=c(min(train_labels), max(train_labels)),
                       breaks=seq(min(train_labels), max(train_labels), 10))+
    geom_hline(yintercept=test_label, linetype="dashed", linewidth=0.7)+
    ggtitle("Sparseness") + theme_bw(base_size = 10)



  weight_plot_sparseness <- ggplot(fitted_weights_sparseness_df,
                                   aes(x=Age, y=Weight, color=Shrinkage)) +
    geom_line(alpha=0.6) + geom_point(size=1.3, alpha=0.7) +
    scale_x_continuous(limits=c(min(train_labels), max(train_labels)),
                       breaks=seq(min(train_labels), max(train_labels), 5)) +
    geom_vline(xintercept=test_label, linetype="dashed", linewidth=0.7)+
    scale_color_manual(values=c("#E69F00", "#0072B2", "#222222"))+
    xlab("Age") + ylab("Weights") + ggtitle("Sparseness") +
    theme_bw(base_size = 10) + theme(plot.title = element_blank())


  normalized_weight_plot_sparseness <- ggplot(fitted_weights_sparseness_df,
                                              aes(x=Age, y=normalized_weights, color=Shrinkage)) +
    geom_line(alpha=0.6) + geom_point(size=1.3, alpha=0.7) +
    scale_x_continuous(limits=c(min(train_labels), max(train_labels)),
                       breaks=seq(min(train_labels), max(train_labels), 5)) +
    geom_vline(xintercept=test_label, linetype="dashed", linewidth=0.7)+
    scale_color_manual(values=c("#E69F00", "#0072B2", "#222222"))+
    xlab("Age") + ylab("Normalized Weights") + ggtitle("Sparseness") +
    theme_bw(base_size = 10) + theme(plot.title = element_blank())


  ## then smoothness
  cat("Tune smoothness penalty...")
  fitted_weights_smoothness <- list()
  for (j in 1:length(lambda_ranges)){

    deconv_result <- iter_deconv(X=train_data,
                                 Y=test_data,
                                 labels=train_labels,
                                 lambda1=lambda_ranges[j], lambda2=0, log=FALSE,
                                 min_scale=min_scale,
                                 max_iter=20)
    weights_df <- data.frame(Age=train_labels, Weight=deconv_result$weights,
                             normalized_weights=deconv_result$normalized_weights,
                             Shrinkage=lambda_ranges[j])
    weights_df$Predicted <- deconv_result$estim
    fitted_weights_smoothness[[j]] <- weights_df

  }

  fitted_weights_smoothness_df <- do.call(rbind, fitted_weights_smoothness)

  predicted_age_smoothness <- fitted_weights_smoothness_df[, c("Shrinkage", "Predicted")] |> unique()
  fitted_weights_smoothness_df$Shrinkage <- factor(fitted_weights_smoothness_df$Shrinkage,
                                                   levels=lambda_ranges)

  prediction_plot_smoothness <- ggplot(predicted_age_smoothness, aes(x=Shrinkage, y=Predicted)) +
    geom_point(size=1.3) + geom_line() + xlab("Shrinkage") + ylab("Prediction") +
    scale_x_continuous(breaks=predicted_age_smoothness$Shrinkage) +
    scale_y_continuous(limits=c(min(train_labels), max(train_labels)),
                       breaks=seq(min(train_labels), max(train_labels), 10))+
    geom_hline(yintercept=test_label, linetype="dashed", linewidth=0.7)+
    scale_color_manual(values=c("#E69F00", "#0072B2", "#222222"))+
    ggtitle("Smoothness") + theme_bw(base_size = 10)


  weight_plot_smoothness <- ggplot(fitted_weights_smoothness_df,
                                   aes(x=Age, y=Weight, color=Shrinkage)) +
    geom_line(alpha=0.6) + geom_point(size=1.3, alpha=0.7) +
    scale_x_continuous(limits=c(min(train_labels), max(train_labels)),
                       breaks=seq(min(train_labels), max(train_labels), 5)) +
    geom_vline(xintercept=test_label, linetype="dashed", linewidth=0.7)+
    scale_color_manual(values=c("#E69F00", "#0072B2", "#222222"))+
    xlab("Age") + ylab("Weights") + ggtitle("Smoothness") +
    theme_bw(base_size = 10) + theme(plot.title = element_blank())


  normalized_weight_plot_smoothness <- ggplot(fitted_weights_smoothness_df,
                                              aes(x=Age, y=normalized_weights, color=Shrinkage)) +
    geom_line(alpha=0.6) + geom_point(size=1.3, alpha=0.7) +
    scale_x_continuous(limits=c(min(train_labels), max(train_labels)),
                       breaks=seq(min(train_labels), max(train_labels), 5)) +
    geom_vline(xintercept=test_label, linetype="dashed", linewidth=0.7)+
    scale_color_manual(values=c("#E69F00", "#0072B2", "#222222"))+
    xlab("Age") + ylab("Normalized Weights") + ggtitle("Smoothness") +
    theme_bw(base_size = 10)+ theme(plot.title = element_blank())



  combined_plots <- (prediction_plot_sparseness + prediction_plot_smoothness) /
    (weight_plot_sparseness + weight_plot_smoothness) /
    (normalized_weight_plot_sparseness + normalized_weight_plot_smoothness) +
    plot_annotation(title=samplename)+
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  sparseness_optimal_id <- which.min(abs(predicted_age_sparseness$Predicted - test_label))
  sparseness_penalty <- predicted_age_sparseness$Shrinkage[sparseness_optimal_id]
  sparseness_prediction <- predicted_age_sparseness$Predicted[sparseness_optimal_id]

  smoothness_optimal_id <- which.min(abs(predicted_age_smoothness$Predicted - test_label))
  smoothness_penalty <- predicted_age_smoothness$Shrinkage[smoothness_optimal_id]
  smoothness_prediction <- predicted_age_smoothness$Predicted[smoothness_optimal_id]

  output <- list(visualization=combined_plots,
                 sparseness_penalty=sparseness_penalty,
                 sparseness_prediction=sparseness_prediction,
                 smoothness_penalty=smoothness_penalty,
                 smoothness_prediction=smoothness_prediction)

  return(output)

}

