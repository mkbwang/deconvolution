


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

    return(list(x=knots, h=h, Q=Q, R=R, RQt=RQt))

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
    return(list(mat_g=mat_g, mat_gamma=t(mat_gamma)))

}



