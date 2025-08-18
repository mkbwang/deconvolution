

#' Calculate the mean values (maybe grouped by a column representing groups)
#'
#' @param expression_df gene expression data frame
#' @param group the variable name to group the samples
#' @returns a vector of not grouped by any variable and a matrix if grouped by a variable
#'
#' @importFrom dplyr group_by summarise across everything %>%
#' @importFrom rlang sym !!
#'
meanexp <- function(expression_df, group=NULL){

    if(is.null(group)){
        expression_mat <- as.matrix(expression_df)
        mean_expression <- colMeans(expression_mat)
    } else{
        mean_expression_df <- expression_df %>% group_by(!!sym(group)) %>%
            summarise(across(everything(), mean))
        mean_expression_df[[group]] <- NULL
        mean_expression <- as.matrix(mean_expression_df)
    }

    return(mean_expression)
}



#' Count the number of samples in each group
#' @param expression_df gene expression data frame
#' @param group the variable name to group the samples
#' @return a vector representing the counts in each group
#'
#' @importFrom dplyr group_by %>% n summarise
#' @importFrom rlang sym !!
#'
group_counts <- function(expression_df, group){

    count_df <- expression_df %>% group_by(!!sym(group)) %>%
        summarise(count=n())
    count_vec <- count_df$count
    return(count_vec)

}



#' Calculate sum of squares of error for each column
#' @param mat input matrix
#' @param vec mean vector
#' @param weights weights for different rows
#'
#' @returns a vector showing the sum of squares for each column
#'
sum_squares <- function(mat, vec, weights=NULL){

    nsamples <- nrow(mat)
    repeat_means <- matrix(vec, nrow=nsamples, ncol=length(vec),
                           byrow=TRUE)

    sq_difference <- (mat-repeat_means)^2
    if (is.null(weights)){
        SSE <- colSums(sq_difference)
    } else{
        SSE <- colSums(sq_difference * weights)
    }

    return(SSE)

}





#' ANOVA of individual gene expressions
#'
#' @param expression_df gene expression data frame
#' @param groups a vector of integers indicating group
#'
#' @returns a dataframe reflecting the total sum of squares, residual sum of squares, explained sum of squares and F statistic
#' @export
var_comp <- function(expression_df, groups){

    genenames <- colnames(expression_df)
    nsamples <- nrow(expression_df)
    expression_mat <- as.matrix(expression_df)

    geneexp_mean <- meanexp(expression_df) # overall mean
    TSS <- sum_squares(mat=expression_mat, vec=geneexp_mean)


    expression_df$Group <- groups
    geneexp_meanbygroup <- meanexp(expression_df, group="Group") # mean by group
    count_bygroup <- group_counts(expression_df, group="Group")
    ngroups <- length(count_bygroup)

    ESS <- sum_squares(geneexp_meanbygroup, geneexp_mean, weights=count_bygroup)
    RSS <- TSS - ESS
    Rsquare <- ESS / TSS
    F_stat <- (ESS/(ngroups - 1)) / (RSS/(nsamples - ngroups))

    SSdf <- data.frame(Gene=genenames,
                       TSS=TSS,
                       RSS=RSS,
                       ESS=ESS,
                       R2=Rsquare,
                       F_stat=F_stat)

    return(SSdf)

}


#' Filter genes that varies a lot between samples with different labels based on ANOVA
#'
#' @param geneexp_df gene expression dataframe with samples on rows and features on columns
#' @param labels labels of all the samples
#' @param min_label minimum label value to be considered for a sample to be included
#' @param max_label maximum label value to be considered for a sample to be included
#' @returns a list with the selected potential marker genes and the corresponding ANOVA results
#' @export
#' @importFrom stats pf p.adjust
f_test_filter <- function(geneexp_df, labels, min_label, max_label){

    mask <- (labels >= min_label) & (labels <= max_label)
    labels_subset <- labels[(labels >= min_label) & (labels <= max_label)]
    geneexp_subset <- geneexp_df[mask, ]
    geneexp_subset_ranks <- apply(geneexp_subset, 2, rank) |>
        as.data.frame() # calculate ranks

    # ANOVA result
    variance_df <-  var_comp(expression_df=geneexp_subset_ranks,
                             groups=labels_subset)

    df1 <- length(unique(labels_subset)) - 1
    df2 <- nrow(geneexp_subset_ranks) - length(unique(labels_subset))

    pvals <- 1 - pf(variance_df$F_stat, df1=df1, df2=df2)
    variance_df$pval <- pvals
    adjusted_pvals <- p.adjust(pvals, method="BH")
    variance_df$adjusted_pvals <- adjusted_pvals

    subset_genes <- variance_df$Gene[adjusted_pvals < 0.05]

    list(subset_genes=subset_genes, variance=variance_df)

}



#' Hierarchical clustering of genes based on standardized log gene expressions across different labels
#' @param geneexp_df gene expression dataframe with samples on rows and features on columns
#' @param labels labels of all samples
#' @param subset_genes a subset of genes to carry out hierarchical clustering
#' @param min_label only include samples whose labels are larger than or equal to min_label
#' @param max_label only include samples whose labels are smaller or equal to max_label
#'
#' @returns dendrogram of hierarchical clustering
#' @importFrom stats dist hclust sd median
#' @importFrom dplyr group_by summarise across everything %>%
#' @export
hcluster_genes <- function(geneexp_df, labels, subset_genes, min_label, max_label){

    mask <- (labels >= min_label) & (labels <= max_label)

    labels_subset <- labels[mask]
    geneexp_df_subset <- geneexp_df[mask, subset_genes]
    log_geneexp_df_subset <- log(geneexp_df_subset)
    log_geneexp_df_subset$Group <- labels_subset

    # calculate mean expression by group
    log_geneexp_bygroup <- log_geneexp_df_subset %>% group_by(.data$Group) %>%
        summarise(across(everything(), median))
    log_geneexp_bygroup$Group <- NULL

    # standardize each column
    mean_loggeneexp <- colMeans(log_geneexp_bygroup)
    repeat_means <- matrix(mean_loggeneexp, nrow=length(unique(labels_subset)),
                           ncol=length(mean_loggeneexp),
                           byrow=TRUE)
    sd_loggeneexp <- apply(log_geneexp_bygroup, 2, sd)
    repeat_sd <- matrix(sd_loggeneexp, nrow=length(unique(labels_subset)),
                           ncol=length(sd_loggeneexp),
                           byrow=TRUE)
    standardized_log_geneexp <- (log_geneexp_bygroup - repeat_means) / repeat_sd


    # calculate distance
    distance <- dist(t(standardized_log_geneexp), method="euclidean")

    # hierarchical clustering
    hc <- hclust(distance, method="ward.D2")

    return(hc)

}




