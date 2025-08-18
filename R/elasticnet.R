


#' Elastic Net Regression Prediction
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param log whether to take log of the values
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated coefficients and prediction
#'
#' @importFrom caret createFolds
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef predict
#' @export
elasticnet <- function(X, Y, labels, log=T, verbose=FALSE){

    # replace zeros in X with pseudocount
    impute_zero <- function(xvec){
        if(sum(xvec == 0) > 0){
            xvec[xvec == 0] <- min(xvec[xvec > 0])/2
        }
        return(xvec)
    }

    if (log){ # take log
        X <- apply(X, 2, impute_zero) # impute zeros in X
        minimum_X <- apply(X, 2, min)
        if (sum(Y == 0) > 0){
            Y[Y == 0] <- minimum_X[Y == 0]
        }
        X <- log(X)
        Y <- log(Y)
    }

    # standardize X and Y
    mean_X <- colMeans(X)
    repeat_means <- matrix(mean_X, nrow=length(labels),
                           ncol=length(mean_X),
                           byrow=TRUE)
    sd_X <- apply(X, 2, sd)
    repeat_sd <- matrix(sd_X, nrow=length(labels),
                        ncol=length(sd_X),
                        byrow=TRUE)
    standardized_X <- (X - repeat_means) / repeat_sd
    standardized_Y <- (Y - mean_X) / sd_X


    # create cross validation folds
    folds <- createFolds(factor(labels), k = 4)
    # Convert list of folds into a numeric vector where each element indicates the fold ID
    foldid <- rep(NA, length(labels))
    for(i in seq_along(folds)) {
        foldid[folds[[i]]] <- i
    }

    # cross validation
    cv_fit <- cv.glmnet(x=standardized_X,
                        y=labels,
                        alpha=0.5, type.measure = "mse",
                        foldid=foldid)

    en_coefs <- coef(cv_fit, s = "lambda.min") |> as.vector()
    predicted_Y <- predict(cv_fit, s="lambda.min", newx=standardized_Y)
    result <- list(coefs=en_coefs, prediction=predicted_Y)
    return(result)

}


