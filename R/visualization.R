

#' Plot the predicted age against the true age
#'
#' @param truth truth vector
#' @param predicted vector of predicted values
#' @param diagonal whether to add an identity line, default TRUE
#' @param title plot title, optional
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_y_continuous geom_abline ggtitle theme element_text
#'
#' @export
#'
viz_predict <- function(truth, predicted, diagonal=TRUE, title=NULL){

    predictions_df <- data.frame(Truth=truth, Predicted=predicted)
    unique_ages <- unique(truth)
    min_age <- min(truth)
    max_age <- max(truth)

    myplot <- ggplot(predictions_df, aes(x=.data$Truth, y=.data$Predicted)) +
        geom_point() + scale_x_continuous(breaks=seq(min_age, max_age, 5), limits=c(min_age-1, max_age+1))+
        scale_y_continuous(breaks=seq(min_age, max_age, 5), limits=c(min_age-1, max_age+1))+
        theme(text=element_text(size=10))

    if(diagonal){
        myplot <- myplot + geom_abline(intercept = 0, slope = 1, color = "blue",                # Change line color
                                      linetype = "dashed",           # Change line style (options include "solid", "dashed", "dotted", etc.)
                                      size = 0.8)
    }

    if (!is.null(title)){
        myplot <- myplot + ggtitle(title)
    }

    return(myplot)

}




#' Visualize the deconvolution weights
#'
#' @param labels labels of training data
#' @param weights estimated weights
#' @param truth true value
#' @param exclude which training sample to remove, optional
#' @param title title of the ggplot
#'
#' @importFrom dplyr group_by %>% arrange summarise
#' @importFrom ggplot2 ggplot geom_line geom_point scale_x_continuous geom_vline xlab ylab ggtitle theme element_text
#' @export
viz_weights <- function(labels, weights, truth=NULL, exclude=NULL, title=NULL){

    weights_df <- data.frame(Labels = labels,
                             Weights= weights) %>%
        arrange(.data$Labels)

    if (!is.null(exclude)){
        weights_df <- weights_df[-exclude, ]
    }
    weights_df$ID <- seq(1, nrow(weights_df))

    # indices where I assign xticks
    starting_point <- weights_df %>% group_by(.data$Labels) %>%
        summarise(index=min(.data$ID))


    myplot <- ggplot(weights_df, aes(x=.data$ID, y=.data$Weights)) +
        geom_line(alpha=0.6) + geom_point(size=1.4, alpha=0.6) +
        scale_x_continuous(breaks=starting_point$index,
                           labels=as.character(starting_point$Labels)) +
        xlab("Age") + ylab("Weights")+
        theme(text=element_text(size=10))

    if (!is.null(truth)){

        # highlight ranges of labels with the same value as truth
        id_ranges <- weights_df$ID[weights_df$Labels == truth]
        myplot <- myplot + geom_vline(xintercept=min(id_ranges)-0.5, linetype="dashed", color="blue")+
            geom_vline(xintercept=max(id_ranges)+0.5, linetype="dashed", color="blue")

    }

    if(!is.null(title)){
        myplot <- myplot + ggtitle(title)
    }
    return(myplot)

}



#' Draw gene expression over age
#'
#' @param ages labels of training data
#' @param values estimated weights
#' @param title title of the ggplot
#'
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous xlab ylab ggtitle theme element_text
#' @export
viz_gene <- function(ages, values, title="Gene Plot"){

    unique_ages <- unique(ages)
    df <- data.frame(Age=ages,
                     Value=values)

    geneexp_plot <- ggplot(df, aes(x=.data$Age, y=.data$Value)) +
        geom_point() + scale_x_continuous(breaks=unique_ages) +
        xlab("Age") +ylab("Expression") + ggtitle(title)

    return(geneexp_plot)

}

