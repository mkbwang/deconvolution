


if (getRversion() >= "2.15.1") {
    utils::globalVariables(".data")
}


#' take log and normalize both training data and test data
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param takelog whether to take log of X and Y before normalization
#' @returns normalized training data and test data
#'
#' @export
preprocess <- function(X, Y, takelog=T){

    # replace zeros in X with pseudocount
    impute_zero <- function(xvec){
        if(sum(xvec == 0) > 0){
            xvec[xvec == 0] <- min(xvec[xvec > 0])/2
        }
        return(xvec)
    }

    if (takelog){ # take log
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
    repeat_means <- matrix(mean_X, nrow=nrow(X),
                           ncol=ncol(X),
                           byrow=TRUE)
    sd_X <- apply(X, 2, sd)
    repeat_sd <- matrix(sd_X, nrow=nrow(X),
                        ncol=ncol(X),
                        byrow=TRUE)
    standardized_X <- (X - repeat_means) / repeat_sd
    standardized_Y <- (Y - mean_X) / sd_X

    return(list(normalized_X=standardized_X,
                normalized_Y=standardized_Y))

}

