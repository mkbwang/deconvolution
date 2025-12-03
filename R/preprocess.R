


if (getRversion() >= "2.15.1") {
    utils::globalVariables(".data")
}


#' Calculate median value and IQR of each feature
#' @param value_mat data matrix (nsample * nfeature)
#' @param margin calculate scale for each row (1) or each column (2)
#' @returns median value vector and scale vector
#'
#' @importFrom stats quantile
#' @export
robust_scale <- function(value_mat, margin=1){
    # take median and report scale based on IQR
    quantile_vals <- apply(value_mat, MARGIN=margin, FUN=function(vec){
        quantile(vec, c(0.25, 0.5, 0.75))
    })

    med_vals <- quantile_vals[2, ]
    scale_vals <- quantile_vals[3, ] - quantile_vals[1, ]

    return(list(median_vals=med_vals,
                scale_vals=scale_vals))
}

#'  Robust scaling of a matrix or a vector
#' @param value_mat data matrix (nsample*nfeature) or vector (nfeature)
#' @param median_vals median value vector (nfeature)
#' @param scale_vals scale value vector (nfeature)
#' @param margin calculate scale for each row (1) or each column (2)
#' @export
scale_transform <- function(value_mat, median_vals, scale_vals, margin=1){

    if (!is.null(dim(value_mat))){
        if (margin == 1){
            scaled_mat <- (value_mat - median_vals) / scale_vals
        } else{
            t_scaled_mat <- (t(value_mat) - median_vals) / scale_vals
            scaled_mat <- t(t_scaled_mat)
        }
        return(scaled_mat)
    } else{
        scaled_vec <- (value_mat - median_vals) / scale_vals
        return(scaled_vec)
    }

}

#' take log and normalize both training data and test data
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param takelog whether to take log of X and Y before normalization
#' @param margin calculate scale for each row (1) or each column (2)
#' @param min_scale features whose IQR value smaller than min_scale will be removed
#' @returns normalized training data and test data
#'
#' @export
preprocess <- function(X, Y, takelog=T, margin=1, min_scale=0.01){

    # replace zeros in X with pseudocount
    impute_zero <- function(xvec){
        if(sum(xvec == 0) > 0){
            xvec[xvec == 0] <- min(xvec[xvec > 0])/2
        }
        return(xvec)
    }

    if (takelog){ # take log
        X <- apply(X, margin, impute_zero) # impute zeros in X
        minimum_X <- apply(X, margin, min)
        if (sum(Y == 0) > 0){
            Y[Y == 0] <- minimum_X[Y == 0]
        }
        X <- log(X)
        Y <- log(Y)
    }

    # standardize X and Y
    scaling_params <- robust_scale(X, margin=margin)
    mask <- scaling_params$scale_vals > min_scale
    if (margin == 1){
        X <- X[mask, ]
    } else{# margin == 2
        X <- X[, mask]
    }

    if (is.null(dim(Y))){
        Y <- Y[mask]
    } else{
        if (margin == 1){
            Y <- Y[mask, ]
        } else{ # margin == 2
            Y <- Y[, mask]
        }
    }

    X_scaled <- scale_transform(X, scaling_params$median_vals[mask],
                                scaling_params$scale_vals[mask],
                                margin=margin)
    Y_scaled <- scale_transform(Y, scaling_params$median_vals[mask],
                                scaling_params$scale_vals[mask],
                                margin=margin)


    return(list(normalized_X=X_scaled,
                normalized_Y=Y_scaled))

}

