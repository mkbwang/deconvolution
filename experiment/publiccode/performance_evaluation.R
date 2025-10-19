
library(dplyr)
library(patchwork)
library(ggplot2)
folder <- "publicdata"

en_files <- list.files("experiment/publicdata/en_model", pattern="*_pred.csv")
deconv_files <- list.files("experiment/publicdata/deconv_model", pattern="*_pred.csv")

study_name <- strsplit(en_files[2], split="_")[[1]][1]
en_prediction_df <- read.csv(file.path("experiment/publicdata/en_model", en_files[2]))
en_prediction_test <- en_prediction_df %>% filter(Batch == "Test")

mae_en <- mean(abs(en_prediction_test$Truth - en_prediction_test$Prediction))
cor_en <- cor(en_prediction_test$Truth, en_prediction_test$Prediction)
title_en <- sprintf("Elastic Net (MAE %.2f, Correlation %.2f)",
                 mae_en, cor_en)

en_performance_plot <- viz_predict(truth=en_prediction_test$Truth,
            predicted=en_prediction_test$Prediction,
            title=title_en)


deconv_prediction_test <- read.csv(file.path("experiment/publicdata/deconv_model", deconv_files[2]))
mae_deconv <- mean(abs(deconv_prediction_test$Truth - deconv_prediction_test$sparseness_prediction))
cor_deconv <- cor(deconv_prediction_test$Truth, deconv_prediction_test$sparseness_prediction)
title_deconv <- sprintf("Deconvolution (MAE %.2f, Correlation %.2f)",
                    mae_deconv, cor_deconv)


deconv_performance_plot <- viz_predict(truth=deconv_prediction_test$Truth,
                                       predicted=deconv_prediction_test$sparseness_prediction,
                                       title=title_deconv)


combined_performance_plot <- en_performance_plot + deconv_performance_plot
plot_name <- sprintf("experiment/publicdata/%s_comparison.pdf", study_name)

ggsave(filename=plot_name, plot=combined_performance_plot,
       width=8, height=4)
