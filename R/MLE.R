




#' Generate posterior density for each feature given the cubic spline fit result (mean and se)
#'
#' @param mean_mat posterior mean matrix, nfeature*ntimestamps
#' @param se_mat posterior mean matrix, nfeature*ntimestamps
#' @param y feature values of a new sample, a vector of length nfeature
#' @param exclude_factor if a feature's value is larger than max_mean+factor*max_se or smaller than
#'  min_mean-factor*max_se, then the feature is considered abnormal and not considered for prediction
#'  @return posterior density for each feature and a validity mask
#' @export
likelihood <- function(mean_mat, se_mat, y, exclude_factor=1){

    max_range <- apply(mean_mat + exclude_factor*se_mat, 1, max)
    min_range <- apply(mean_mat - exclude_factor*se_mat, 1, min)

    validity <- (y < max_range) & (y > min_range) # features with abnormal values are not considered

    lik <- exp(-(y-mean_mat)^2 / (2*se_mat^2))

    return(list(lik=lik, mask=validity))

}


#' Given posterior densities fo all the valid features, generate a consensus posterior density and use it for final prediction
#' @param lik likelihood values across all labels for each feature
#' @param feature_significance a vector with length nfeature representing the usefulness of each feature, based on adjustedR2 of cubic spline fit
#' @param mask a binary vector with length equal to number of features, exclude features whose values are beyond the training data range
#' @param labels label values
#' @return joint likelihood curve and mle estimate
#' @export
mle <- function(lik, feature_significance=rep(1, nrow(lik)), mask, labels){



    weight_vec <- mask * feature_significance
    names(weight_vec) <- rownames(lik)

    lik_square <- t(lik) %*% (weight_vec * lik)

    eigen_result <- eigen(lik_square, symmetric=TRUE)

    evector <- eigen_result$vectors[, 1] # guarantee to be of the same sign given Perron Frobenius theorem
    stopifnot(all(evector >= 0) | all(evector <= 0))
    evector <- abs(evector)

    joint_likelihood <- data.frame(Label=labels,
                                   lik=evector)

    mle_estimate <- labels[which.max(evector)]

    output <- list(lik=joint_likelihood,
                   estimate=mle_estimate)

    return(output)

}



