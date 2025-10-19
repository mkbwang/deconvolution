

#' Fit penalized cubic regression splines
#'
#' @param df input dataframe with two columns "Age" and "Feature"
#' @param fname feature name to be shown in the plot title
#' @param k number of spline bases, default 10
#'
#'
#' @returns effective degrees of freedom, F statistics, p value, visualization
#'
#' @importFrom mgcv gam s
#' @importFrom stats pf predict
#' @importFrom ggplot2 ggplot geom_point aes geom_line labs
#' @export
gam_fit <- function(df, fname, k=10){

    # df has two columns, Age and Feature
    min_age <- min(df$Age)
    max_age <- max(df$Age)
    gam_model <- gam(Feature ~ s(Age, k=k, bs="cr", m=2), # cubic splines
                     data=df)

    RSS_gam <- sum(gam_model$residuals^2)
    DOF_gam <- sum(gam_model$edf1) # effective degrees of freedom
    DOF_rss <- nrow(df) - DOF_gam # degrees of freedom of residuals

    TSS <- sum((df$Feature - mean(df$Feature))^2) # total sum of squares
    ESS <- TSS - RSS_gam # explained sum of squares by gam
    F_stat <- (ESS/(DOF_gam - 1)) / (RSS_gam/DOF_rss)

    pval <- 1-pf(q=F_stat, df1=DOF_gam-1, df2=DOF_rss) # ANOVA F test

    pred_df <- data.frame(Age=seq(min_age, max_age, 0.1))
    pred_df$Feature <- predict(gam_model, newdata=pred_df)

    viz <- ggplot() + geom_point(data=df, aes(x=Age, y=Feature), size=1, alpha=0.6)+
        geom_line(data=pred_df, aes(x=Age, y=Feature), linewidth=1, alpha=0.8, color="blue") +
        labs(x="Age", y=fname, title=sprintf("EDF=%.2f, pval=%.4f", DOF_gam, pval))+
        theme(text = element_text(size = 10))
    return(list(EDF=DOF_gam, Fstat=F_stat, pval=pval, viz=viz))

}

