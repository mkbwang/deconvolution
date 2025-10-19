


library(arrow)
library(dplyr)
rm(list=ls())
data_folder <- "experiment/publicdata"
code_folder <- "experiment/publiccode"

# list all the file names
files <- list.files(path=data_folder, pattern="*.rds")
extract_studyname <- function(filename){
  filename_nosuffix <- strsplit(filename, split="[.]")[[1]][1]
  study <- strsplit(filename_nosuffix, split="[_]")[[1]][2]
  return(study)
}
organs <- sapply(files, extract_studyname)
meta_files <- sprintf("metadata_%s.csv", organs)

source(file.path(code_folder, "en_utils.R"))


fname <- files[3] #TODO: change here
metaname <- meta_files[3]
organ <- organs[3]

metadata <- read.csv(file.path(data_folder, metaname), row.names=1)
biomarker <- readRDS(file.path(data_folder, fname))  |> t()
rownames(metadata) <- metadata$SRR.ID



set.seed(2025)
test_id <- sample(nrow(metadata), nrow(metadata)*0.25)
train_id <- setdiff(seq(1, nrow(metadata)), test_id)


train_meta <- metadata[train_id, ]
test_meta <- metadata[test_id, ]

train_biomarker <- biomarker[train_id, ]
test_biomarker <- biomarker[test_id, ]

train_ages <- train_meta$Age
test_ages <- test_meta$Age


scaling_params <- robust_scale(train_biomarker)


if(sum(scaling_params$scale_vals < 1e-3) > 0){

  cols_remove <- which(scaling_params$scale_vals < 1e-3)
  scaling_params$median_vals <- scaling_params$median_vals[-cols_remove]
  scaling_params$scale_vals <- scaling_params$scale_vals[-cols_remove]

  train_biomarker <- train_biomarker[, -cols_remove]
  test_biomarker <- test_biomarker[, -cols_remove]

}



train_biomarker_scaled <- scale_transform(value_mat=train_biomarker,
                                median_vals=scaling_params$median_vals,
                                scale_vals=scaling_params$scale_vals)

test_biomarker_scaled <- scale_transform(value_mat=test_biomarker,
                               median_vals=scaling_params$median_vals,
                               scale_vals=scaling_params$scale_vals)



en_result <- en_fit_func(train_data=train_biomarker,
                         train_ages=train_ages,
                         test_data=test_biomarker,
                         test_ages=test_ages)
en_coef_df <-  en_result$coefs
en_pred_df <- en_result$prediction
en_pred_df_test <- en_pred_df %>% filter(Batch=="Test")
# write.csv(en_coef_df, sprintf("/nfs/turbo/sph-ligen/wangmk/GeneExp_Aging/experiment/prediction/coef_raw_%s.csv", selected_study_name),
#           row.names=FALSE, quote=FALSE)
# write.csv(en_pred_df, sprintf("/nfs/turbo/sph-ligen/wangmk/GeneExp_Aging/experiment/prediction/prediction_raw_%s.csv", selected_study_name),
#           quote=FALSE)

en_result_scaled <- en_fit_func(train_data=train_biomarker_scaled,
                                train_ages=train_ages,
                                test_data=test_biomarker_scaled,
                                test_ages=test_ages)

en_coef_scaled_df <- en_result_scaled$coefs
en_pred_scaled_df <- en_result_scaled$prediction
en_pred_scaled_df_test <- en_pred_scaled_df %>% filter(Batch == "Test")
write.csv(en_coef_scaled_df, file.path(data_folder, "en_model", sprintf("%s_coef.csv", organ)),
          row.names=FALSE, quote=FALSE)
write.csv(en_pred_scaled_df, file.path(data_folder, "en_model", sprintf("%s_pred.csv", organ)),
          quote=FALSE)






