

library(stringr)
rm(list=ls())
# load all the data
folder <- "experiment/planarian_data"
marker_file <- "AgingClock_hisat_tpm.89.redSum.combat.csv"
metadata_file <- "AgingClock_meta.89.csv"



marker_data <- read.csv(file.path(folder, marker_file), row.names = 1)
metadata <- read.csv(file.path(folder, metadata_file), row.names=1)
sample_names <- colnames(marker_data)

train_samples <- sample_names[8:40] # skip samples at age 2 for now
test_samples <- sample_names[grepl("MK", sample_names)]


train_data <- marker_data[, train_samples]
test_data <- marker_data[, test_samples]
prevalence <- rowMeans(train_data > 0)
train_data <- train_data[prevalence == 1, ]
test_data <- test_data[prevalence == 1, ]


# fit cubic spline to all the features
train_ages <- sapply(train_samples, function(sname){
    str_extract_all(sname, "\\d+")[[1]] |> as.integer()
})
log_train_data <- log(train_data) |> as.matrix()

x <- unique(train_ages) # knots
new_t <- seq(6, 23, 0.2)

predicted_values <- matrix(0, nrow=nrow(train_data), ncol=length(new_t))
se_values <- matrix(0, nrow=nrow(train_data), ncol=length(new_t))
rownames(predicted_values) <- rownames(se_values) <- rownames(train_data)
colnames(predicted_values) <- colnames(se_values) <- sprintf("Age_%.1f", new_t)


feature_spline_summary <- data.frame(Feature=rownames(train_data),
                                EDF=0,
                                Fstat=0,
                                Pval=0)

for (j in 1:nrow(train_data)){

    if (j %% 100 == 0){
        print(j)
    }

    input_y <- log_train_data[j, ]
    spline_fit <- fit_cspline(t=train_ages,
                              y=input_y,x=x)
    feature_spline_summary$EDF[j] <- spline_fit$EDF
    feature_spline_summary$Fstat[j] <- spline_fit$Fstat
    feature_spline_summary$Pval[j] <- spline_fit$pval

    prediction <- predict_cspline(t=new_t, fitted_result=spline_fit)
    predicted_values[j, ] <- prediction$y_hat
    se_values[j, ] <- prediction$se_pred

}


feature_spline_summary$padj <- p.adjust(feature_spline_summary$Pval, method="BH")

viz_prediction <- function(observed_t, observed_y, new_t, y_hat, se_hat, new_y=NULL){

    plot(observed_t, observed_y)
    lines(new_t, y_hat)
    lines(new_t, y_hat + se_hat, lty=2)
    lines(new_t, y_hat - se_hat, lty=2)
    if (!is.null(new_y)){
        lines(new_t, rep(new_y, length(new_t)), lty=2, col="red")
    }

}

id <- 1145
viz_prediction(observed_t = train_ages, observed_y = log_train_data[id, ],
               new_t=new_t, y_hat=predicted_values[id, ], se_hat=se_values[id, ])

## select features whose adjusted p values are smaller than 0.02

marker_features <- feature_spline_summary$Feature[feature_spline_summary$padj < 0.02]

marker_predict_values <- predicted_values[marker_features, ]
marker_se_values <- se_values[marker_features, ]
marker_original_values <- log_train_data[marker_features, ]

# given a new test sample


## filter nonzero entries and take log
target <- test_data[marker_features, "MK3"]
nonzero_indices <- which(target > 0)
target_nonzero <- target[nonzero_indices]
log_target_mat <- matrix(rep(log(target_nonzero), length(new_t)), nrow=length(target_nonzero),
                      ncol=length(new_t))
marker_subset_mean <- marker_predict_values[nonzero_indices, ]
marker_subset_se <- marker_se_values[nonzero_indices, ]
marker_subset_original <- marker_original_values[nonzero_indices, ]
rownames(log_target_mat) <- rownames(marker_subset_mean)

## remove feature values that are much larger or smaller than the fitted training spline range
upper_bound <- apply(marker_subset_mean + marker_subset_se*1, 1, max)
lower_bound <- apply(marker_subset_mean - marker_subset_se*1, 1, min)

abnormal_filter <- (log_target_mat[, 1] > lower_bound) & (log_target_mat[, 1] < upper_bound)
marker_subset_mean <- marker_subset_mean[abnormal_filter, ]
marker_subset_se <- marker_subset_se[abnormal_filter, ]
marker_subset_original <- marker_subset_original[abnormal_filter, ]
log_target_mat <- log_target_mat[abnormal_filter, ]


## calculate posterior

posterior_mat <- exp(-(log_target_mat - marker_subset_mean)^2 / (2*marker_subset_se^2)) / marker_subset_se
posterior_mat <- sweep(posterior_mat, 1, rowSums(posterior_mat), FUN="/")




## Simpson's index, plot posterior probabilities
simpson <- rowSums(posterior_mat^2)
top_examples <- simpson[simpson > 0.04 & simpson < 0.05]

View(feature_spline_summary %>% filter(Feature %in% names(top_examples)))

par(mfrow = c(1, 2))
for (j in 1:6){
    example <- names(top_examples)[j]
    viz_prediction(observed_t=train_ages, observed_y=marker_subset_original[example, ],
                   new_t = new_t, y_hat=marker_subset_mean[example, ], se_hat=marker_subset_se[example, ],
                   new_y = log_target_mat[example, 1])
    plot(new_t, posterior_mat[example, ])
}


peak <- 12

bases <- monotone_bases(labels=new_t, peak=peak)
l_bases <- bases$l_bases[1:length(new_t), ]
g_bases <- bases$g_bases[1:length(new_t), ]
n_l_bases <- bases$n_l_bases
n_g_bases <- bases$n_g_bases


#TODO: nlopt for optimization of rayleigh quotient


