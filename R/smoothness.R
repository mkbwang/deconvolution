

#' Generate penalty matrix
#' @param labels phenotypes of all the samples
#' @param lambda penalty parameter for similarity of weights for samples with the same label
#' @returns penalty quadratic matrix
#' @export
penalty_smooth <- function(labels, lambda=1){

    # I already assume that the labels are unique already
    sorted_labels <- sort(labels)
    # penalty for samples with different but adjacent labels
    penalty_mat <- matrix(0, nrow=length(labels), ncol=length(labels))

    for (j in 2:(length(sorted_labels)-1)){

        weight_vec <- rep(0, length(labels))
        label_1 <- sorted_labels[j-1]
        index_1 <- which(labels == label_1)
        weight_vec[index_1] <- 1

        label_2 <- sorted_labels[j]
        index_2 <- which(labels == label_2)
        weight_vec[index_2] <- -2

        label_3 <- sorted_labels[j+1]
        index_3 <- which(labels == label_3)
        weight_vec[index_3] <- 1
        weight_mat <- as.matrix(weight_vec, ncol=1)
        penalty_mat <- penalty_mat + weight_mat %*% t(weight_mat)

    }

    penalty_mat <- penalty_mat * lambda

    return(penalty_mat)
}


