

#' Generate posterior density for each feature given the cubic spline fit result (mean and se)
#'
#' @param postmean posterior mean matrix, nfeature*ntimestamps
#' @param postse posterior mean matrix, nfeature*ntimestamps
#' @param y feature values of a new sample, a vector of length nfeature
#' @param prior prior distribution
#' @param exclude_factor if a feature's value is larger than max_mean+factor*max_se or smaller than
#'  min_mean-factor*max_se, then the feature is considered abnormal and not considered for prediction
#'  @return posterior density for each feature and a validity mask
#' @export
posterior <- function(postmean, postse, y, prior=rep(1, ncol(postse)), exclude_factor=1){

    max_range <- apply(postmean + exclude_factor*postse, 1, max)
    min_range <- apply(postmean - exclude_factor*postse, 1, min)

    prior_mat <- matrix(prior, nrow=nrow(postmean), ncol=length(prior),
                        byrow=TRUE)

    validity <- (y < max_range) & (y > min_range) # features with abnormal values are not considered

    post_density <- 1 / postse * exp(-(y-postmean)^2 / (2*postse^2)) * prior_mat # not normalized yet

    normalized_density <- post_density / rowSums(post_density)

    return(list(post_density=normalized_density, mask=validity))

}

#' Given posterior densities fo all the valid features, generate a consensus posterior density and use it for final prediction
#' @param post_density posterior density matrix, nfeature * nlabels, each row sum to one
#' @param feature_significance a vector with length nfeature representing the usefulness of each feature, based on F statistic of cubic spline fit
#' @param validity a binary vector with length equal to number of features
#' @param labels label values
#' @return posterior probability, MAP, posterior mean and the importance of each feature
#' @export
eigen_posterior <- function(post_density, feature_significance=rep(1, nrow(post_density)), validity, labels){

    # weigh each feature by simpson index - 1/length(labels) (based on Jensen inequality)
    # take first eigen vector

    simpson_indices <- rowSums(post_density^2)
    adjusted_simpson_indices <- simpson_indices - 1/length(labels)
    weight_vec <- validity / simpson_indices * feature_significance

    names(weight_vec) <- rownames(post_density)

    density_square <- t(post_density) %*% (weight_vec * post_density)

    eigen_result <- eigen(density_square, symmetric=TRUE)

    evector <- eigen_result$vectors[, 1] # guarantee to be of the same sign given Perron Frobenius theorem
    stopifnot(all(evector >= 0) | all(evector <= 0))
    consensus_posterior <- evector/sum(evector)

    max_a_posterior <- labels[which.max(consensus_posterior)]
    mean_posterior <- sum(labels * consensus_posterior)

    output <- list(MAP=max_a_posterior, postmean=mean_posterior,
                   postdensity=consensus_posterior,
                   weights=weight_vec)

}



