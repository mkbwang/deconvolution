

#' Generate penalty matrix
#' @param labels phenotypes of all the samples
#' @param p0 penalty parameter for similarity of weights for samples with the same label
#' @param p1 penalty parameter for smoothness of average weights among three adjacent different labels
#' @returns penalty quadratic matrix
#'@export
penalty_smooth <- function(labels, p0=1, p1=10){

    unique_labels <- unique(labels) |> sort()

    # penalty for samples with the same labels
    labels_mat <- do.call(cbind, replicate(length(labels), labels, simplify = FALSE))
    penalty_mat_0 <- -1*(labels_mat == t(labels_mat))
    diag_0 <- rep(0, length(labels))
    for (unique_label in unique_labels){
        indices <- which(labels == unique_label)
        diag_0[indices] <- sum(labels == unique_label) - 1
    }
    diag(penalty_mat_0) <- diag_0


    # penalty for samples with different but adjacent labels
    penalty_mat_1 <- matrix(0, nrow=length(labels), ncol=length(labels))

    for (j in 2:(length(unique_labels)-1)){

        weight_vec <- rep(0, length(labels))
        label_1 <- unique_labels[j-1]
        indices_1 <- which(labels == label_1)
        weight_vec[indices_1] <- -1/length(indices_1)

        label_2 <- unique_labels[j]
        indices_2 <- which(labels == label_2)
        weight_vec[indices_2] <- 2/length(indices_2)

        label_3 <- unique_labels[j+1]
        indices_3 <- which(labels == label_3)
        weight_vec[indices_3] <- -1/length(indices_3)
        weight_mat <- do.call(cbind, replicate(length(labels), weight_vec, simplify = FALSE))
        penalty_mat_1 <- penalty_mat_1 + weight_mat * t(weight_mat)

    }

    penalty_mat <- penalty_mat_0 *p0 + penalty_mat_1 *p1

    return(penalty_mat)
}


