

#' Ispline value for a four-knot equally spaced segment
#'
#' @param x input x vector value
#' @param theta value of the smallest knot
#' @param delta value of the gap
#' @param multiplicity maximium multiplicity at either end of the knots
#' @param loc left or right end knot to have more than zero multiplicity
#' @param normalize binary indicator, whether to normalize so that maximum value is one
#' @returns return a list that contains all the inputs and the corresponding y values
#'
#' @export
#'
ispline_basis_func <- function(x, theta, delta, multiplicity=c(0,1,2),
                               loc=c("left", "right"), normalize=FALSE){

    y <- rep(0, length(x))

    if (multiplicity == 0){
        knots <- theta + seq(0, 3) * delta
        mask1 <- (x >= knots[1]) & (x < knots[2])
        y[mask1] <- 1/(6*delta^2)*(x[mask1] - theta)^3
        mask2 <- (x >= knots[2]) & (x < knots[3])
        y[mask2] <- delta/6 + 1 / delta^2 * (x[mask2]-theta-delta) *
            (-1/3*x[mask2]^2 + (2/3*theta + 7/6*delta)*x[mask2] - 1/3*theta^2 - 7/6*theta*delta - 1/3*delta^2)
        mask3 <- (x >= knots[3]) & (x < knots[4])
        y[mask3] <- delta + 1/(6*delta^2)*(x[mask3]-theta-3*delta)^3
        mask4 <- (x >= knots[4])
        y[mask4] <- delta
        if (normalize){
            y <- y/delta
        }
    } else if (multiplicity == 1){
        if (loc == "left") {
            knots <- c(theta, theta, theta+delta, theta+2*delta)
            mask2 <- (x >= knots[2]) & (x < knots[3])
            y[mask2] <- 1 / (2 * delta^2) * (x[mask2]-theta) * (-x[mask2]^2 + (2*theta+2*delta) * x[mask2] - theta^2 - 2*theta*delta)
            mask3 <- (x >= knots[3]) & (x < knots[4])
            y[mask3] <- 2/3*delta + (x[mask3] - theta - 2 * delta)^3 / (6*delta^2)
            mask4 <- (x >= knots[4])
            y[mask4] <- 2/3*delta
        } else{
            knots <- c(theta, theta+delta, theta+2*delta, theta+2*delta)
            mask1 <- (x >= knots[1]) & (x < knots[2])
            y[mask1] <- 1/(6*delta^2)*(x[mask1]-theta)^3
            mask2 <- (x >= knots[2]) & (x < knots[3])
            y[mask2] <- delta/6 + 1/(2*delta^2)*(x[mask2]-theta-delta)*
                (-x[mask2]^2+(2*theta+3*delta)*x[mask2]-theta^2-3*theta*delta-delta^2)
            mask3 <- (x>=knots[3])
            y[mask3] <- 2/3*delta
        }
        if (normalize){
            y <- y/(2/3*delta)
        }
    } else{# multiplicity equal to 2
        if (loc == "left") {
            knots <- c(theta, theta, theta, theta+delta)
            mask3 <- (x >= knots[3]) & (x < knots[4])
            y[mask3] <- delta / 3 + (x[mask3] - theta - delta)^3 / (3 * delta^2)
            mask4 <- (x >= knots[4])
            y[mask4] <- delta / 3
        } else{
            knots <- c(theta, theta+delta, theta+delta, theta+delta)
            mask1 <- (x >= knots[1]) & (x < knots[2])
            y[mask1] <- (x[mask1] - theta)^3/(3*delta^2)
            mask2 <- (x >= knots[2])
            y[mask2] <- delta / 3
        }
        if (normalize){
            y <- y/(1/3*delta)
        }
    }
    return(list(x=x, y=y, theta=theta, delta=delta, multiplicity=multiplicity, loc=loc))

}



#' generate I-spline basis components
#'
#' @param minval smallest value of the knots
#' @param maxval largest value of the knots
#' @param xvals x values to be evalutated for spline values
#' @param nknots number of knots including the two end points
#'
#' @returns a matrix with nknots+1 basis values (each basis has length equal to xvals)
#'
#' @export
#'
ispline_basis_setup <- function(minval, maxval, xvals, nknots=11){

    gap <- (maxval - minval) / (nknots - 1)
    knots <- seq(minval, maxval, gap)
    knots <- c(rep(min(knots), 2), knots, rep(max(knots), 2))


    multiplicities <- rep(0, length(knots))
    multiplicities[c(1,2)] <- c(2,1)
    multiplicities[c(nknots, nknots+1)] <- c(1,2)

    mult_end <- rep("right", length(knots))
    mult_end[c(1,2)] <- rep("left", 2)


    ispline_y_vals <- list()
    basis_components <- matrix(0, nrow=length(xvals), ncol=nknots+1)

    for (j in 1:(nknots+1)){

        output <- ispline_basis_func(x=xvals,
                                     theta=knots[j],
                                     delta=gap,
                                     multiplicity=multiplicities[j],
                                     loc=mult_end[j],
                                     normalize = T)

        basis_components[, j] <- output$y

    }
    return(basis_components)

}




#' Deconvolution where the position with maximum weight is predetermined
#'
#' @param X training data matrix (nsample * nfeature)
#' @param Y test data vector (nfeature)
#' @param weights vector of weights for each feature, if null then is equally one
#' @param labels labels for all the training samples, assumed to be unique
#' @param log whether to take log before fitting deconvolution
#' @param min_scale features whose IQR value smaller than min_scale will be removed
#' @param peak the label value with maximum weight
#' @param verbose Whether to print out details of osqp function, default False
#' @returns weight bases
#'
#' @importFrom Matrix Matrix
#' @importFrom osqp osqp osqpSettings
#' @export
#'
deconvolution_fixed_peak <- function(X, Y, weights=NULL, labels, log=T, min_scale=0.01, peak,
                                     verbose=FALSE){

    # X is the expression matrix
    # Y is the expression vector
    # I assume that the labels are unique and already sorted
    # peak must be between the smallest and the largest label

    X <- X[order(labels), ]
    labels <- labels[order(labels)]

    # normalize X and Y
    preprocessed_data <- preprocess(X=X, Y=Y,
                                    takelog=log,
                                    min_scale=min_scale)
    standardized_X <- preprocessed_data$normalized_X
    standardized_Y <- preprocessed_data$normalized_Y

    if (is.null(weights)){
        weights_mat <- diag(length(standardized_Y))
    } else{
        weights_mat <- diag(weights)
    }

    standardized_X_withpeak <- rbind(standardized_X, rep(0, ncol(standardized_X)))
    labels_withpeak <- c(labels, peak)


    # first set up bases for labels smaller than or equal to peak
    n_labels_l <- sum(labels_withpeak <= peak)
    if (n_labels_l < 10){
        # use step functions
        unique_labels_l <- unique(labels_withpeak[labels_withpeak <= peak]) |> sort()
        l_bases <- matrix(0, ncol=length(unique_labels_l), nrow=length(labels_withpeak))
        # step function bases
        for (j in 1:length(unique_labels_l)){
            l_bases[, j] <- as.integer((labels_withpeak >= unique_labels_l[j]) & (labels_withpeak <= peak))
        }
    } else{
        # ispline bases and a constant basis
        nknots <- floor(n_labels_l/3)+1
        min_knot <- min(labels_withpeak)
        max_knot <- peak
        l_bases <- ispline_basis_setup(minval=min_knot, maxval=max_knot,
                                       xvals=labels_withpeak, nknots=nknots)
        l_bases <- cbind(rep(1, length(labels_withpeak)), l_bases)
        l_bases[labels_withpeak > peak, ] <- 0
    }
    l_bases <- l_bases - 0.5 * (labels_withpeak == peak)
    n_l_bases <- ncol(l_bases)

    # set up the bases for labels bigger than or equal to peak
    n_labels_g <- sum(labels_withpeak >= peak)
    if (n_labels_g < 10){
        # use step functions
        unique_labels_g <- unique(labels_withpeak[labels_withpeak >= peak]) |> sort()
        g_bases <- matrix(0, ncol=length(unique_labels_g), nrow=length(labels_withpeak))
        # step function bases
        for (j in 1:length(unique_labels_g)){
            g_bases[, j] <- as.integer((labels_withpeak >= peak) & (labels_withpeak <= unique_labels_g[j]))
        }
        # values at the labels equal to peak needs to be halved so that
    } else{
        # ispline bases and a constant basis
        nknots <- floor(n_labels_g / 3)+1
        min_knot <- peak
        max_knot <- max(labels_withpeak)
        g_bases <- 1 - ispline_basis_setup(minval=min_knot, maxval=max_knot,
                                       xvals=labels_withpeak, nknots=nknots)
        g_bases <- cbind(rep(1, length(labels_withpeak)), g_bases)
        g_bases[labels_withpeak < peak, ] <- 0
    }
    g_bases <- g_bases - 0.5 * (labels_withpeak == peak)
    n_g_bases <- ncol(g_bases)

    equality_constraints <- c(rep(1, n_l_bases), rep(-1, n_g_bases))


    Theta <-cbind(l_bases, g_bases) # all the bases
    Xt_Theta <- t(standardized_X_withpeak) %*% Theta
    XWXt <- t(Xt_Theta) %*% weights_mat %*% Xt_Theta
    XWY <- t(Xt_Theta) %*% standardized_Y

    Pmat <- XWXt
    qvec <- -XWY


    # linear constraints
    ## larger than zero
    A_0 <- diag(n_l_bases + n_g_bases)
    l_0 <- rep(0, n_l_bases + n_g_bases)
    u_0 <- rep(Inf, n_l_bases + n_g_bases)

    ## equal sum of coefficients for l_bases and g_bases
    A_1 <- equality_constraints
    l_1 <- 0
    u_1 <- 0

    A_combined <- Matrix(rbind(A_0, A_1), sparse=T)
    l_combined <- c(l_0, l_1)
    u_combined <- c(u_0, u_1)

    settings <- osqpSettings(alpha = 1.0, verbose=verbose)
    model <- osqp(P=Pmat, q=qvec, A=A_combined, l=l_combined, u=u_combined,
                  pars=settings)
    res <- model$Solve()

    solution <- res$x
    solution[solution < 0] <- 0

    estimated_betas <- Theta %*% solution

    return(list(labels=labels_withpeak, l_bases=l_bases, g_bases=g_bases,
                solution=solution, beta=estimated_betas))

}


