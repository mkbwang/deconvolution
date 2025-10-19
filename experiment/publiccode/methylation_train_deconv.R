
library(arrow)
library(dplyr)
library(ggplot2)
library(patchwork)
library(parallelly)
library(foreach)
library(deconvolution)
library(doParallel)
rm(list=ls())

experiment_folder <- "experiment/publiccode"
methylation_folder <- "experiment/publicdata"
# studies <- dir(file.path(methylation_folder, "performance", "deconvolution"))


library(optparse)
option_list <- list(make_option(c("-i", "--index"), type="integer", default=1,
                                help="index [default=1]"),
                    make_option(c("-n", "--name"), type="character", default="",
                                help="study name [default=]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

study_name <- opt$name
id <- opt$index
study_name <- "TGCA_KIRC"
id <- 10

fname <- sprintf("methylation_%s.feather", study_name)
metaname <- sprintf("metadata_%s.csv", study_name)

# load all the metadata and the biomarker data
metadata <- read.csv(file.path(methylation_folder, metaname), row.names=1)
biomarker <- read_feather(file.path(methylation_folder, fname)) |> as.data.frame()
rownames(biomarker) <- rownames(metadata)


# I have fit elastic net regression on it, get the train/test split
en_prediction <- read.csv(file.path(methylation_folder, "en_model",
                                    sprintf("%s_pred.csv", study_name)), row.names=1)
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
# deconv_prediction_test <- en_prediction_test
# deconv_prediction_test$sparseness <- 0
# deconv_prediction_test$sparseness_prediction <- 0
# deconv_prediction_test$smoothness <- 0
# deconv_prediction_test$smoothness_prediction <- 0
# deconv_prediction_test$Prediction <- NULL

# visualization_list <- list()

biomarker_test <- biomarker[test_samplenames, ]


source(file.path(experiment_folder, "deconv_utils.R"))

# fit deconvolution one by one
# for (j in 1:nrow(metadata_test)){

sname <- rownames(biomarker_test)[id]
biomarker_target <- biomarker_test[id, ]
metadata_target <- metadata_test[id, ]
age_target <- metadata_target$age

scaling_params <- robust_scale(as.matrix(biomarker_train_avg))

preprocessed_data <- preprocess(X=as.matrix(biomarker_train_avg),
                                Y=as.vector(as.matrix(biomarker_target)),
                                takelog=FALSE,
                                min_scale=0.01)
standardized_X <- preprocessed_data$normalized_X
standardized_Y <- preprocessed_data$normalized_Y

begin <- proc.time()
result <- l2_solve(X=standardized_X, Y=standardized_Y, labels=unique_train_ages,
                   lambda1=0, lambda2=0)
end <- proc.time()


result <- fit_deconv(train_data=as.matrix(biomarker_train_avg),
                     train_labels=unique_train_ages,
                     test_data=as.vector(as.matrix(biomarker_target)),
                     test_label=age_target,
                     samplename=sname)



deconv_prediction_test$sparseness[j] <- result$sparseness_penalty
deconv_prediction_test$sparseness_prediction[j] <- result$sparseness_prediction
deconv_prediction_test$smoothness[j] <- result$smoothness_penalty
deconv_prediction_test$smoothness_prediction[j] <- result$smoothness_prediction

    # visualization_list[[j]] <- result$visualization

# }




plots_fname <- file.path(methylation_folder, "performance", "deconvolution", study_name,
                         sprintf("%s_deconv.pdf", study_name))
pdf(plots_fname, width = 9, height = 9)
for (plot in visualization_list) {
    print(plot)
}
dev.off()



prediction_fname <- file.path(methylation_folder, "performance", "deconvolution", study_name,
                              sprintf("%s_pred.csv", study_name))
write.csv(deconv_prediction_test, prediction_fname, quote=FALSE)

