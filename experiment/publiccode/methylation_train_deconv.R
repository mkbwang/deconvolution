
library(arrow)
library(dplyr)
library(ggplot2)
library(patchwork)
rm(list=ls())
data_folder <- "experiment/publicdata"
code_folder <- "experiment/publiccode"


# list all the file names
files <- list.files(path=data_folder, pattern="*.feather")
meta_files <- list.files(path=data_folder, pattern="*.csv")


fname <- files[1] #TODO: change here
metaname <- meta_files[1]
study_name <- gsub("[metadata_]|[.csv]", "", metaname)

# load all the metadata and the biomarker data
metadata <- read.csv(file.path(data_folder, metaname), row.names=1)
biomarker <- read_feather(file.path(data_folder, fname)) |> as.data.frame()
rownames(biomarker) <- rownames(metadata)


# I have fit elastic net regression on it, get the train/test split
en_prediction <- read.csv(file.path(data_folder, "en_model", sprintf("%s_pred.csv", study_name)),
                          row.names=1)
## retrieve the training data and training ages
en_prediction_train <- en_prediction %>% filter(Batch == "Train")
metadata_train <- metadata[rownames(en_prediction_train), ]
train_ages <- metadata_train$age
biomarker_train <- biomarker[rownames(en_prediction_train), ]
# for each feature, take the median at every age
biomarker_train$age <- train_ages
biomarker_train_avg <- biomarker_train %>% group_by(age) %>% summarise_all(median)

unique_train_ages <- biomarker_train_avg$age ## unique ages as training labels
biomarker_train_avg$age <- NULL


## retrieve the test data
en_prediction_test <- en_prediction %>% filter(Batch=="Test")

test_samplenames <- rownames(en_prediction_test)
metadata_test <- metadata[test_samplenames, ]
deconv_prediction_test <- en_prediction_test
deconv_prediction_test$sparseness <- 0
deconv_prediction_test$sparseness_prediction <- 0
deconv_prediction_test$smoothness <- 0
deconv_prediction_test$smoothness_prediction <- 0
deconv_prediction_test$Prediction <- NULL

visualization_list <- list()

biomarker_test <- biomarker[test_samplenames, ]
## robust normalization of test data based on the training data (age averaged)
normalized_data <- preprocess(X=biomarker_train_avg, Y=biomarker_test,
                              takelog=FALSE)
biomarker_train_normalized <- normalized_data$normalized_X
biomarker_test_normalized <- normalized_data$normalized_Y


# function for fitting and evaluating deconvolution
fit_deconv <- function(train_data, test_data, train_labels, test_label){

    lambda_ranges <- c(0, 2000, 4000)

    # tuning different smoothness penalty and sparseness penalty
    ## first sparseness
    predicted_age_sparseness <- data.frame(Shrinkage=lambda_ranges,
                                           Predicted=0)
    fitted_weights_sparseness <- list()

    cat("Tune sparseness penalty...")
    for (j in 1:length(lambda_ranges)){
        deconv_result <- l2_solve(X=train_data,
                                  Y=test_data,
                                  labels=train_labels,
                                  lambda1=0, lambda2=lambda_ranges[j])

        predicted_age_sparseness$Predicted[j] <- deconv_result$estim
        weights_df <- data.frame(Age=train_labels, Weight=deconv_result$weights,
                                 normalized_weights=deconv_result$normalized_weights,
                                 Shrinkage=lambda_ranges[j])
        fitted_weights_sparseness[[j]] <- weights_df

    }
    fitted_weights_sparseness_df <- do.call(rbind, fitted_weights_sparseness)
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
    predicted_age_smoothness <- data.frame(Shrinkage=lambda_ranges,
                                           Predicted=0)
    fitted_weights_smoothness <- list()
    cat("Tune smoothness penalty...")
    for (j in 1:length(lambda_ranges)){
        deconv_result <- l2_solve(X=train_data,
                                  Y=test_data,
                                  labels=train_labels,
                                  lambda1=lambda_ranges[j], lambda2=0)

        predicted_age_smoothness$Predicted[j] <- deconv_result$estim
        weights_df <- data.frame(Age=train_labels, Weight=deconv_result$weights,
                                 normalized_weights=deconv_result$normalized_weights,
                                 Shrinkage=lambda_ranges[j])
        fitted_weights_smoothness[[j]] <- weights_df

    }

    fitted_weights_smoothness_df <- do.call(rbind, fitted_weights_smoothness)
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

# fit deconvolution one by one
for (j in 1:nrow(metadata_test)){

    print(rownames(biomarker_test_normalized)[j])
    biomarker_target <- biomarker_test_normalized[j, ]
    metadata_target <- metadata_test[j, ]
    age_target <- metadata_target$age
    result <- fit_deconv(train_data=biomarker_train_normalized,
                         train_labels=unique_train_ages,
                         test_data=biomarker_target,
                         test_label=age_target)

    deconv_prediction_test$sparseness[j] <- result$sparseness_penalty
    deconv_prediction_test$sparseness_prediction[j] <- result$sparseness_prediction
    deconv_prediction_test$smoothness[j] <- result$smoothness_penalty
    deconv_prediction_test$smoothness_prediction[j] <- result$smoothness_prediction

    visualization_list[[j]] <- result$visualization

}

plots_fname <- sprintf("experiment/publicdata/deconv_model/%s_deconv.pdf",
                       study_name)
pdf(plots_fname, width = 9, height = 9)
for (plot in visualization_list) {
    print(plot)
}
dev.off()


prediction_fname <- sprintf("experiment/publicdata/deconv_model/%s_deconv.csv",
                           study_name)
write.csv(deconv_prediction_test, prediction_fname, quote=FALSE)



