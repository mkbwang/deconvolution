

#' deconvolution based on l2 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param penalty penalty coefficient to force weights of similar samples to be similar
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
l2_solve <- function(X, Y, labels, penalty=0, verbose=FALSE){

    nsamples <- nrow(X)
    ngenes <- ncol(X)

    XXt <- X %*% t(X)
    XY <- X %*% Y

    # penalty based on differences between ages
    labels_mat <- matrix(labels, nrow=length(labels), ncol=length(labels))
    inverse_abs_distance <- penalty * ngenes / (abs(labels_mat - t(labels_mat)) + 1)
    penalty_mat <- -inverse_abs_distance
    diag(penalty_mat) <- diag(penalty_mat) + rowSums(inverse_abs_distance)

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
#' @param penalty penalty coefficient to force weights of similar samples to be similar
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
l1_solve <- function(X, Y, labels, penalty=0, verbose=FALSE){

    nsamples <- nrow(X)
    ngenes <- ncol(X)

    # penalty based on differences between ages
    labels_mat <- matrix(labels, nrow=length(labels), ncol=length(labels))
    inverse_abs_distance <- penalty * ngenes / (abs(labels_mat - t(labels_mat)) + 1)
    penalty_mat <- -inverse_abs_distance
    diag(penalty_mat) <- diag(penalty_mat) + rowSums(inverse_abs_distance)

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




#' deconvolution based on l1 loss
#'
#' @param X matrix of training data (nsample * nfeature)
#' @param Y vector of test data
#' @param labels labels of all the training samples
#' @param penalty penalty coefficient to force weights of similar samples to be similar
#' @param log whether to take log of the values
#' @param type l1 loss or l2 loss
#' @param verbose Whether to print out details of osqp function, default False
#'
#' @returns Estimated weights and the weighted sum of training labels
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
deconvolution <- function(X, Y, labels, penalty=0, log=T, type=c("l1", "l2"),
                          verbose=FALSE){

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

    type <- match.arg(type)
    if (type == "l2"){
        result <- l2_solve(X=standardized_X, Y=standardized_Y,
                           labels=labels, penalty=penalty,
                           verbose=verbose)
    } else{
        result <- l1_solve(X=standardized_X, Y=standardized_Y,
                           labels=labels, penalty=penalty,
                           verbose=verbose)
    }
    return(result)
}

