


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
#' @export
fit_cspline <- function(t, y, x, alphas=c(2^seq(-5, 4, 1))){

    N <- length(t)
    knots_utils <-cspline_utils(knots=x)
    B_mat <- knots_utils$B # transformation from g to vector of g and gamma
    K_mat <- knots_utils$K_mat # penalty component for smoothness
    RQt <- knots_utils$RQt # R^-1*Qt
    T_mat <- dmat_utils(t=t, x=x) # design matrix given input independent variable value
    U_mat <- T_mat %*% B_mat

    original_t <- t
    original_y <- y
    t_excluded <- c()
    y_excluded <- c()

    # for loop starts from here
    while(1){
        UtU <- t(U_mat) %*% U_mat
        eigen_ratio  <- sum(diag(UtU)) / sum(diag(K_mat))

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
        RSS <- sum(resids^2)
        TSS <- sum((y - mean(y))^2)
        ESS <- TSS - RSS
        sigma_hat2 <- RSS / (N-EDF)
        F_stat <- (ESS/(EDF-1)) / sigma_hat2
        P_value <- 1 - pf(q=F_stat, df1=EDF-1, df2=N-EDF)

        pred_se <- sqrt(diag(sigma_hat2 * (H_mat %*% H_mat)) + sigma_hat2)

        outlier_mask <- abs(resids) > 3*pred_se
        if (any(outlier_mask)){ # remove outliers and refit if necessary
            t_excluded <- c(t_excluded, t[outlier_mask])
            y_excluded <- c(y_excluded, y[outlier_mask])
            t <- t[!outlier_mask]
            N <- length(t)
            y <- y[!outlier_mask]
            U_mat <- U_mat[!outlier_mask, ]
        } else{
            break
        }
    }
    if (length(t_excluded) > 0){
        exclusion_mask <- (original_t %in% t_excluded) & (original_y %in% y_excluded)
    } else{
        exclusion_mask <- rep(FALSE, length(original_t))
    }



    output <- list(x=x, t=original_t, y=original_y, exclusion=exclusion_mask,
                   g_hat=g_hat, gamma_hat=gamma_hat, sigma_hat2=sigma_hat2, EDF=EDF,
                   Fstat=F_stat, pval=P_value, y_hat=y_hat, se=pred_se, proj_mat=proj_mat, alpha=best_alpha)

    return(output)

}


#' predict mean values and standard errors for given input independent variable values
#'
#' @param t vector of input values
#' @param fitted_result a list that contains x (knots of spline), g_hat (values at each knot), sigma_hat2 (estimated random error variance),
#' proj_mat (used for calculating prediction error)
#' @export
predict_cspline <- function(t, fitted_result){

    x <- fitted_result$x
    g_hat <- fitted_result$g_hat
    sigma_hat2 <- fitted_result$sigma_hat2
    proj_mat <- fitted_result$proj_mat

    T_mat <- dmat_utils(t=t, x=x)
    B_mat <- cspline_utils(knots=x)$B
    U_mat <- T_mat %*% B_mat

    y_hat <- U_mat %*% g_hat # predicted value

    H_mat <- U_mat %*% proj_mat
    sigma2_pred <- diag(sigma_hat2 * (H_mat %*% t(H_mat))) + sigma_hat2
    se_pred <- sqrt(sigma2_pred) # standard error
    return(list(y_hat=y_hat, se_pred=se_pred))

}






#' Fit cubic spline to all the features
#'
#' @param input_pheno continuous phenotype values of training data
#' @param input_marker_value marker values with dimension nfeature * nsample
#' @param output_pheno selected phenotype values for prediction of mean and standard error
#'
#' @returns list of spline fitting result and predicted mean/standard error
#' @export
fit_data <- function(input_pheno, input_marker_value, output_pheno){

    # confirm that the phenotype values for prediction fall in the range of the training data values
    stopifnot(all(output_pheno >= min(input_pheno)) & all(output_pheno <= max(input_pheno)))

    feature_names <- rownames(input_marker_value)

    # get EDF, p values and F stats
    feature_spline_summary <- data.frame(Feature=feature_names,
                                         EDF=0,
                                         Fstat=0,
                                         Pval=0,
                                         NumOutlier=0)

    input_outlier_mat <- matrix(0, nrow=nrow(input_marker_value),
                          ncol=ncol(input_marker_value))
    rownames(input_outlier_mat) <- rownames(input_marker_value)
    colnames(input_outlier_mat) <- colnames(input_marker_value)

    # estimated mean and standard error values for each feature for each output pheno
    output_mean_mat <- matrix(0, nrow=length(feature_names), ncol=length(output_pheno))
    output_se_mat <- matrix(0, nrow=length(feature_names), ncol=length(output_pheno))
    rownames(output_mean_mat) <- rownames(output_se_mat) <- rownames(input_marker_value)

    # knot selection
    unique_pheno <- sort(unique(input_pheno))
    indices <- c(seq(1, length(unique_pheno), 2), length(unique_pheno)) |> unique() # pick half of the unique values
    knots <- unique_pheno[indices] |> round() |> unique() # round to integers

    # fit csplines for each feature
    for (j in 1:length(feature_names)){

        feature_values <- input_marker_value[j, ]
        result_fit <- fit_cspline(t=input_pheno, y=feature_values, x=knots)
        input_outlier_mat[j, ] <- result_fit$exclusion
        feature_spline_summary$EDF[j] <- result_fit$EDF
        feature_spline_summary$Fstat[j] <- result_fit$Fstat
        feature_spline_summary$Pval[j] <- result_fit$pval
        feature_spline_summary$NumOutlier[j] <- sum(result_fit$exclusion)

        prediction_output <- predict_cspline(t=output_pheno,
                                       fitted_result=result_fit)
        output_mean_mat[j, ] <- prediction_output$y_hat
        output_se_mat[j, ] <- prediction_output$se_pred

    }

    return(list(spline_summary=feature_spline_summary,
                input_outlier=input_outlier_mat,
                output_mean=output_mean_mat,
                output_se=output_se_mat,
                output_pheno=output_pheno))

}



