


#' Elastic Net Regression Prediction
#'
#' @param train_data feature values of training samples (nsample*nfeature)
#' @param train_labels response values of training samples
#' @param test_data feature values of test samples (nsample*nfeature)
#' @param test_labels response values of test samples
#'
#' @returns Estimated coefficients and prediction
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef predict
#' @export
elasticnet <- function(train_data, train_labels, test_data, test_labels){

    # cross validation
    en_fit <- cv.glmnet(x=as.matrix(train_data), y=train_labels, alpha=0.5,
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

    prediction_df <- data.frame(Truth=c(train_labels, test_labels),
                                Prediction=c(en_train_prediction, en_test_prediction))

    prediction_df$Batch <- c(rep("Train", length(train_labels)), rep("Test", length(test_labels)))
    rownames(prediction_df) <- c(rownames(train_data), rownames(test_data))
    output <- list(coefs=coef_df,
                   prediction=prediction_df)

    return(output)

}



