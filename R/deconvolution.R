

#' deconvolution based on l2 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param weights vector of weights for each feature, if null then is equally one
#' @param p0 penalty parameter for similarity of weights for samples with the same label
#' @param p1 penalty parameter for smoothness of average weights among three adjacent different labels
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
l2_solve <- function(X, Y, labels, weights=NULL, p0=0, p1=0, verbose=FALSE){

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
    penalty_mat <- penalty_smooth(labels=labels,
                                  p0=p0, p1=p1)

    Pmat <- Matrix(XXt + penalty_mat, sparse=T)
    qvec <- -XY


    # linear constraints
    ## larger than zero
    A_0 <- diag(nsamples)
    l_0 <- rep(0, nsamples)
    u_0 <- rep(Inf, nsamples)

    ## sum to one
    A_1 <- rep(1, nsamples)
    l_1 <- 1
    u_1 <- 1

    A_mat <- Matrix(rbind(A_0, A_1), sparse=T)
    l_vec <- c(l_0, l_1)
    u_vec <- c(u_0, u_1)

    settings <- osqpSettings(alpha = 1.0, verbose=verbose)
    model <- osqp(P=Pmat, q=qvec, A=A_mat, l=l_vec, u=u_vec,
                  pars=settings)

    res <- model$Solve()

    solution <- res$x
    solution[solution < 1/nsamples/10] <- 0
    solution <- solution/sum(solution)

    label_estim <- sum(solution*labels)

    return(list(weights=solution, estim=label_estim))

}



#' deconvolution based on l1 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param p0 penalty parameter for similarity of weights for samples with the same label
#' @param p1 penalty parameter for smoothness of average weights among three adjacent different labels
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
l1_solve <- function(X, Y, labels, p0=0, p1=0, verbose=FALSE){

    nsamples <- nrow(X)
    ngenes <- ncol(X)

    # penalty based on differences between ages
    penalty_mat <- penalty_smooth(labels=labels,
                                  p0=p0, p1=p1)

    Pmat <- matrix(0, nrow=nsamples+ngenes, ncol=nsamples+ngenes)
    Pmat[1:nsamples, 1:nsamples] <- penalty_mat
    Pmat <- Matrix(Pmat, sparse=T)

    qvec <- c(rep(0, nsamples), rep(1, ngenes))

    # linear constraints
    ## sum to one
    A_0 <- c(rep(1, nsamples), rep(0, ngenes))
    l_0 <- 1
    u_0 <- 1

    ## nonnegative weights
    A_1 <- cbind(diag(nsamples), matrix(0, nrow=nsamples, ncol=ngenes))
    l_1 <- rep(0, nsamples)
    u_1 <- rep(Inf, nsamples)

    ## absolute value 1
    A_2 <- cbind(t(X), diag(ngenes))
    l_2 <- Y
    u_2 <- rep(Inf, ngenes)

    ## absolute value 2
    A_3 <- cbind(-t(X), diag(ngenes))
    l_3 <- -Y
    u_3 <- rep(Inf, ngenes)

    A_mat <- Matrix(rbind(A_0, A_1, A_2, A_3), sparse=T)
    l_vec <- c(l_0, l_1, l_2, l_3)
    u_vec <- c(u_0, u_1, u_2, u_3)

    settings <- osqpSettings(alpha = 1.0, verbose=verbose)
    model <- osqp(P=Pmat, q=qvec, A=A_mat, l=l_vec, u=u_vec,
                  pars=settings)

    res <- model$Solve()

    solution <- res$x[1:nsamples]
    solution[solution < 1/nsamples/10] <- 0
    solution <- solution/sum(solution)

    label_estim <- sum(solution*labels)

    return(list(weights=solution, estim=label_estim))

}




#' deconvolution based on l1 or l2 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param weights vector of weights for each feature, if null then is equally one
#' @param p0 penalty parameter for similarity of weights for samples with the same label
#' @param p1 penalty parameter for smoothness of average weights among three adjacent different labels
#' @param log whether to take log of the values
#' @param type l1 loss or l2 loss
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
deconvolution <- function(X, Y, labels, weights=NULL, p0=0, p1=0, log=T, type=c("l1", "l2"),
                          verbose=FALSE){

    preprocessed_data <- preprocess(X=X, Y=Y,
                                    takelog=log)

    standardized_X <- preprocessed_data$normalized_X
    standardized_Y <- preprocessed_data$normalized_Y

    type <- match.arg(type)
    if (type == "l2"){
        result <- l2_solve(X=standardized_X, Y=standardized_Y,
                           weights = weights,
                           labels=labels, p0=p0, p1=p1,
                           verbose=verbose)
    } else{
        result <- l1_solve(X=standardized_X, Y=standardized_Y,
                           labels=labels, p0=p0, p1=p1,
                           verbose=verbose)
    }
    return(result)
}

