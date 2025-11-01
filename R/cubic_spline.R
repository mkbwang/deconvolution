


#' generate utility matrices for fitting penalized cubic spline regressions
#' @param knots unique knots for the cubic spline
#' @returns a list with knots (x), gap (h), Q matrix, R matrix and R^-1Q^t
#' @export
cspline_utils <- function(knots){

    knots <- sort(knots)
    K <- length(knots) #  number of knots
    h <- diff(knots)

    # set up Q matrix
    Q <- matrix(0, nrow=K, ncol=K-2)
    Q1 <- 1/(h[1:(K-2)])
    Q3 <- 1/(h[2:(K-1)])
    Q2 <- -Q1-Q3
    Q[cbind(1:(K-2), 1:(K-2))] <- Q1
    Q[cbind(2:(K-1), 1:(K-2))] <- Q2
    Q[cbind(3:K, 1:(K-2))] <- Q3

    # set up R matrix
    R <- diag(1/3*(h[1:(K-2)] + h[2:(K-1)]))
    R[row(R) == col(R)+1] <- R[row(R)+1 == col(R)] <- h[2:(K-2)]/6

    # R^-1Q^T
    RQt <- solve(R, t(Q))

    # Kmat
    K_mat <- Q %*% RQt

    # B mat
    B <- rbind(diag(K), RQt)

    return(list(x=knots, h=h, Q=Q, R=R, RQt=RQt, B=B, K_mat=K_mat))

}



# TODO: construct design matrix based on the input data
#' Design matrix for regression given the knots
#'
#' @param t input independent variable values vector
#' @param x vector of knots
#' @export
dmat_utils <- function(t, x){

    stopifnot(all(t >= min(x) & t <= max(x)))
    x <- sort(x)
    h <- diff(x)
    K <- length(x)

    # first construct the vector for values at each knot
    t_mat_1 <- matrix(t, nrow=length(t), ncol=K-1)
    mask_1 <- (t(t_mat_1) >= x[1:(K-1)] & t(t_mat_1) <= x[2:K])
    matg1 <- rbind( (1 - (t(t_mat_1) - x[1:(K-1)])/h) * mask_1, 0)
    matg2 <- rbind(0, (t(t_mat_1) - x[1:(K-1)]) / h * mask_1)
    mat_g <- pmax(matg1, matg2) |> t()

    # then construct the vector for second order gradients at each knot
    t_mat_2 <- t_mat_1[, 1:(K-2)]
    mat_gamma<- mask_1[1:(K-2), ] * ((t(t_mat_2) - x[1:(K-2)])^3 / (6 * h[1:(K-2)]) - (t(t_mat_2) - x[1:(K-2)]) * h[1:(K-2)] / 6) +
        mask_1[2:(K-1), ] * ( -(t(t_mat_2) - x[2:(K-1)])^3 / (6 * h[2:(K-1)]) + (t(t_mat_2) - x[2:(K-1)])^2 / 2 - (t(t_mat_2) - x[2:(K-1)]) * h[2:(K-1)] / 3)

    design_mat <- cbind(mat_g, t(mat_gamma))

    return(design_mat)

}



#' Fit cubic spline given input training dataframe and knots
#'
#' @param t independent variable values
#' @param y observed values
#' @param x knots
#' @param alphas penalty parameters to choose from
fit_cspline <- function(t, y, x, alphas=c(2^seq(-4, 4, 1))){

    N <- length(t)
    knots_utils <-cspline_utils(knots=x)
    B_mat <- knots_utils$B # transformation from g to vector of g and gamma
    K_mat <- knots_utils$K_mat # penalty component for smoothness
    RQt <- knots_utils$RQt # R^-1*Qt
    T_mat <- dmat_utils(t=t, x=x) # design matrix given input independent variable value
    U_mat <- T_mat %*% B_mat
    UtU <- t(U_mat) %*% U_mat
    eigen_UtU <- eigen(UtU)
    # stopifnot(all(eigen_UtU$values >= 0))
    eigen_K <- eigen(K_mat)
    # print(eigen_K$values)
    # stopifnot(all(eigen_K$values >= 0))

    eigen_ratio  <- sum(eigen_UtU$values) / sum(eigen_K$values)

    # choose the best alpha
    gcvs <- rep(0, length(alphas))
    for (j in 1:length(alphas)){

        alpha <- alphas[j]
        UtU_pen_inv <- chol2inv(chol(UtU + eigen_ratio * alpha*K_mat))

        g_hat <- UtU_pen_inv %*% t(U_mat) %*% y
        H_mat <- U_mat %*% UtU_pen_inv %*% t(U_mat)
        y_hat <- U_mat %*% g_hat
        resids <- y_hat - y

        EDF <- sum(diag(H_mat))
        sigma_hat2 <- sum(resids^2) / (N-EDF)
        gcvs[j] <- sum(resids^2) / (N-EDF)^2 * N

    }

    # refit with the best alpha
    best_alpha <- alphas[which.min(gcvs)]
    UtU_pen_inv <- chol2inv(chol(UtU + eigen_ratio * best_alpha*K_mat))

    proj_mat <- UtU_pen_inv %*% t(U_mat)
    g_hat <- proj_mat %*% y
    gamma_hat <- RQt %*% g_hat
    H_mat <- U_mat %*% UtU_pen_inv %*% t(U_mat)
    y_hat <- U_mat %*% g_hat
    resids <- y_hat - y

    EDF <- sum(diag(H_mat))
    sigma_hat2 <- sum(resids^2) / (N-EDF)

    pred_se <- diag(sigma_hat2 * (H_mat %*% H_mat))

    output <- list(x=x, g_hat=g_hat, gamma_hat=gamma_hat, sigma_hat2=sigma_hat2, EDF=EDF,
                   y_hat=y_hat, se=pred_se, proj_mat=proj_mat, alpha=best_alpha)

    # TODO: prediction error is E(Var(Y)) + Var(E(Y))
    # TODO:

}



