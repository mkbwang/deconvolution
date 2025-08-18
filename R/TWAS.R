


#' linear association test
#'
#' @param ages ages of all samples
#' @param expression_df gene expression data frame (sampe * gene)
#' @param takelog whether to take log, by default TRUE
#' @returns Regression coefficient and p values
#'
#' @importFrom stats lm
#' @export
twas <- function(ages, expression_df, takelog=T){
    if (takelog){
        expression_df <- log(expression_df)
    }
    gene_names <- colnames(expression_df)

    result <- data.frame(Gene=gene_names, effsize=0, pval=0)
    for (j in 1:length(gene_names)){
        expr <- expression_df[, gene_names[j]]
        lr_result <- lm(ages ~ expr) |> summary()
        result$effsize[j] <- lr_result$coefficients[2, 1]
        result$pval[j] <- lr_result$coefficients[2, 4]
    }
    return(result)
}



#' check monotonicity
#' @param ages ages of all samples
#' @param expression_df gene expression data frame (sampe * gene)
#' @param takelog whether to take log, by default TRUE
#' @returns Regression coefficient and p values
#' @importFrom stats median
#' @importFrom dplyr group_by summarise_all across everything %>%
#' @export
monotone_check <- function(ages, expression_df, takelog=T){

    if (takelog){
        expression_df <- log(expression_df)
    }
    gene_names <- colnames(expression_df)

    # take median by age
    expression_df$Age <- ages
    expression_df_meanbyage <- expression_df %>% group_by(Age) %>%
        summarise_all(median)

    unique_ages <- expression_df_meanbyage$Age
    expressions_mean_mat <- expression_df_meanbyage[, -1] |> as.matrix()
    changes <- expressions_mean_mat[1:(length(unique_ages)-1), ] < expressions_mean_mat[2:length(unique_ages), ]
    changes_df <- t(changes) |> as.data.frame()
    colnames(changes_df) <- sprintf("Age_%d_%d", unique_ages[1:(length(unique_ages)-1)],
                                    unique_ages[2:length(unique_ages)])

    rownames(changes_df) <- gene_names
    return(changes_df)

}
