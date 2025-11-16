

#' Generate posterior density for each feature given the cubic spline fit result (mean and se)
#'
#' @param postmean posterior mean matrix, nfeature*ntimestamps
#' @param postse posterior mean matrix, nfeature*ntimestamps
#' @param y feature values of a new sample, a vector of length nfeature
#' @param exclude_factor if a feature's value is larger than max_mean+factor*max_se or smaller than
#'  min_mean-factor*max_se, then the feature is considered abnormal and not considered for prediction
#'  @return posterior density for each feature and a validity mask
#' @export
posterior <- function(postmean, postse, y, exclude_factor=1.5){

    max_range <- apply(postmean + exclude_factor*postse, 1, max)
    min_range <- apply(postmean - exclude_factor*postse, 1, min)

    validity <- (y < max_range) & (y > min_range) # features with abnormal values are not considered

    post_density <- 1 / postse * exp(-(y-postmean)^2 / (2*postse^2)) # not normalized yet
    normalized_density <- post_density / rowSums(post_density)

    return(list(post_density=normalized_density, mask=validity))

}

#' Given posterior densities fo all the valid features, generate a consensus posterior density and use it for final prediction
#' @param post_density posterior density matrix, nfeature * nlabels, each row sum to one
#' @param validity a binary vector with length equal to number of features
#' @param labels label values
#' @return posterior probability, MAP, posterior mean and the importance of each feature
#' @export
eigen_posterior <- function(post_density, validity, labels){

    # weigh each feature by simpson index - 1/length(labels) (based on Jensen inequality)
    # take first eigen vector


}



