

#' deconvolution based on l2 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param weights vector of weights for each feature, if null then is equally one
#' @param lambda1 penalty parameter for smoothness
#' @param lambda2 penalty parameter for LASSO shrinkage
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
l2_solve <- function(X, Y, labels, weights=NULL, lambda1=0, lambda2=0, verbose=FALSE){

    nsamples <- nrow(X)
    ngenes <- ncol(X)

    if (is.null(weights)){
        weights_mat <- diag(length(Y))
    } else{
        weights_mat <- diag(weights)
    }

    XXt <- X %*% weights_mat %*% t(X)
    XY <- X %*% weights_mat  %*% Y

    # penalty based on differences between ages
    omega <- penalty_smooth(labels=labels, lambda=lambda1)

    Pmat <- Matrix(XXt + omega, sparse=T)
    qvec <- -XY + lambda2/2


    # linear constraints
    ## larger than zero
    A_0 <- diag(nsamples)
    l_0 <- rep(0, nsamples)
    u_0 <- rep(Inf, nsamples)

    ## sum to one
    # A_1 <- rep(1, nsamples)
    # l_1 <- 1
    # u_1 <- 1

    A_mat <- Matrix(A_0, sparse=T)
    l_vec <- l_0
    u_vec <- u_0

    settings <- osqpSettings(alpha = 1.0, verbose=verbose)
    model <- osqp(P=Pmat, q=qvec, A=A_mat, l=l_vec, u=u_vec,
                  pars=settings)
    res <- model$Solve()

    solution <- res$x
    solution[solution < 0] <- 0

    fitted_Y <- as.vector(t(X) %*% solution)
    resids <- fitted_Y - Y

    normalized_solution <- solution/sum(solution)
    # trim very small weights to zero
    # trim_indices <- which(normalized_solution < 1/nsamples/10)
    # solution[trim_indices] <- 0
    # normalized_solution <- solution/sum(solution)

    label_estim <- sum(normalized_solution*labels)

    return(list(weights=solution, normalized_weights=normalized_solution,
                estim=label_estim, resids=resids))

}



#' deconvolution based on l1 or l2 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param weights vector of weights for each feature, if null then is equally one
#' @param lambda1 penalty parameter for smoothness
#' @param lambda2 penalty parameter for LASSO shrinkage
#' @param log whether to take log of the values
#' @param min_scale features whose IQR value smaller than min_scale will be removed
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
deconvolution <- function(X, Y, labels, weights=NULL, lambda1=0, lambda2=0, log=T,
                          min_scale=0.01,
                          verbose=FALSE){

    preprocessed_data <- preprocess(X=X, Y=Y,
                                    takelog=log,
                                    min_scale=min_scale)

    standardized_X <- preprocessed_data$normalized_X
    standardized_Y <- preprocessed_data$normalized_Y

    type <- match.arg(type)
    # if (type == "l2"){
    result <- l2_solve(X=standardized_X, Y=standardized_Y,
                       weights = weights,
                       labels=labels, lambda1=lambda1, lambda2=lambda2,
                       verbose=verbose)
    # } else{
    #     result <- l1_solve(X=standardized_X, Y=standardized_Y,
    #                        labels=labels, p0=p0, p1=p1,
    #                        verbose=verbose)
    # }
    return(result)
}


#' Iterative deconvolution based on l1 or l2 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param lambda1 penalty parameter for smoothness
#' @param lambda2 penalty parameter for LASSO shrinkage
#' @param nu stablizing parameter for reweighing features for deconvolution
#' @param epsilon tolerance value for stoppage of iterations
#' @param log whether to take log of the values
#' @param min_scale features whose IQR value smaller than min_scale will be removed
#' @param max_iter maximum number of iterations
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @export
iter_deconv <- function(X, Y, labels, lambda1=0, lambda2=0, log=T, min_scale=0.01,
                        nu=1e-3, epsilon=1e-2,
                        max_iter=20, verbose=FALSE){


    weights <- rep(0, length(labels))
    normalized_weights <- rep(0, length(labels))

    preprocessed_data <- preprocess(X=X, Y=Y,
                                    takelog=log,
                                    min_scale=min_scale)
    standardized_X <- preprocessed_data$normalized_X
    standardized_Y <- preprocessed_data$normalized_Y
    feature_weights <- rep(1, length(standardized_Y))

    j <- 0
    while(j < max_iter){
        print(j)
        result <- l2_solve(X=standardized_X, Y=standardized_Y, labels=labels, weights=feature_weights,
                                lambda1=lambda1, lambda2=lambda2)
        j <- j+1
        new_weights <- result$weights
        mask <- (new_weights > 0) | (weights > 0)
        new_normalized_weights <- result$normalized_weights
        new_top_indices <- order(new_weights, decreasing=T)[1:min(5, length(labels))]

        total_change <- sum(abs(new_normalized_weights - normalized_weights))
        print(total_change)
        print(new_top_indices)
        if (total_change < epsilon) break
        weights <- new_weights
        normalized_weights <- new_normalized_weights
        top_indices <- new_top_indices

        resids <- result$resids
        feature_weights <- 1/(nu + resids^2)
        feature_weights <- feature_weights / median(feature_weights)
    }

    label_estim <- sum(normalized_weights * labels)

    output <- list(weights=weights, normalized_weights=normalized_weights,
                   estim=label_estim, num_iter=j)

}




# deconvolution based on l1 loss
#
# @param X matrix of training data (nsample * nfeature)
# @param Y vector of test data
# @param labels labels of all the training samples
# @param p0 penalty parameter for similarity of weights for samples with the same label
# @param p1 penalty parameter for smoothness of average weights among three adjacent different labels
# @param verbose Whether to print out details of osqp function, default False
#
# @returns Estimated weights and the weighted sum of training labels
#
# @importFrom Matrix Matrix
# @importFrom osqp osqp osqpSettings
# @export

# l1_solve <- function(X, Y, labels, p0=0, p1=0, verbose=FALSE){
#
#     nsamples <- nrow(X)
#     ngenes <- ncol(X)
#
#     # penalty based on differences between ages
#     penalty_mat <- penalty_smooth(labels=labels,
#                                   p0=p0, p1=p1)
#
#     Pmat <- matrix(0, nrow=nsamples+ngenes, ncol=nsamples+ngenes)
#     Pmat[1:nsamples, 1:nsamples] <- penalty_mat
#     Pmat <- Matrix(Pmat, sparse=T)
#
#     qvec <- c(rep(0, nsamples), rep(1, ngenes))
#
#     # linear constraints
#     ## sum to one
#     A_0 <- c(rep(1, nsamples), rep(0, ngenes))
#     l_0 <- 1
#     u_0 <- 1
#
#     ## nonnegative weights
#     A_1 <- cbind(diag(nsamples), matrix(0, nrow=nsamples, ncol=ngenes))
#     l_1 <- rep(0, nsamples)
#     u_1 <- rep(Inf, nsamples)
#
#     ## absolute value 1
#     A_2 <- cbind(t(X), diag(ngenes))
#     l_2 <- Y
#     u_2 <- rep(Inf, ngenes)
#
#     ## absolute value 2
#     A_3 <- cbind(-t(X), diag(ngenes))
#     l_3 <- -Y
#     u_3 <- rep(Inf, ngenes)
#
#     A_mat <- Matrix(rbind(A_0, A_1, A_2, A_3), sparse=T)
#     l_vec <- c(l_0, l_1, l_2, l_3)
#     u_vec <- c(u_0, u_1, u_2, u_3)
#
#     settings <- osqpSettings(alpha = 1.0, verbose=verbose)
#     model <- osqp(P=Pmat, q=qvec, A=A_mat, l=l_vec, u=u_vec,
#                   pars=settings)
#
#     res <- model$Solve()
#
#     solution <- res$x[1:nsamples]
#     solution[solution < 1/nsamples/10] <- 0
#     solution <- solution/sum(solution)
#
#     label_estim <- sum(solution*labels)
#
#     return(list(weights=solution, estim=label_estim))
#
# }




