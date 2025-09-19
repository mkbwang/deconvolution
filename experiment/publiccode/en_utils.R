

# calculating scaling parameter
robust_scale <- function(value_mat){
  # take median and report scale based on IQR
  quantile_vals <- apply(value_mat, MARGIN=2, FUN=function(vec){
    quantile(vec, c(0.25, 0.5, 0.75))
  })
  
  med_vals <- quantile_vals[2, ]
  scale_vals <- quantile_vals[3, ] - quantile_vals[1, ]
  
  return(list(median_vals=med_vals,
              scale_vals=scale_vals))
}

# center and scale the data
scale_transform <- function(value_mat, median_vals, scale_vals){
  if (!is.null(dim(value_mat))){
    t_scaled_mat <- (t(value_mat) - median_vals) / scale_vals
    return(t(t_scaled_mat))
  } else{
    scaled_vec <- (value_mat - median_vals) / scale_vals
    return(scaled_vec)
  }
}



# elastic net
library(glmnet)

en_fit_func <- function(train_data, train_ages, test_data, test_ages){
  
  en_fit <- cv.glmnet(x=as.matrix(train_data), y=train_ages, alpha=0.5,
                      type.measure="mae")
  en_coefs <- coef(en_fit, s="lambda.min")
  en_coefs <- en_coefs[-1]
  nonzero_indices <- which(en_coefs != 0)

  coef_df <- data.frame(Feature=colnames(train_data)[nonzero_indices],
                        Value=en_coefs[nonzero_indices])
  
  en_train_prediction <- predict(en_fit, newx=as.matrix(train_data),
                                 s="lambda.min")
  en_test_prediction <- predict(en_fit, newx=as.matrix(test_data),
                                s="lambda.min")
  
  prediction_df <- data.frame(Truth=c(train_ages, test_ages),
                              Prediction=c(en_train_prediction, en_test_prediction))
  
  prediction_df$Batch <- c(rep("Train", length(train_ages)), rep("Test", length(test_ages)))
  rownames(prediction_df) <- c(rownames(train_data), rownames(test_data))
  output <- list(coefs=coef_df,
                 prediction=prediction_df)
  
  return(output)
  
}


