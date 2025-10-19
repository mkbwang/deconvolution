



#' simulate from trignometric functions
#'
#'@param knots x values
#'@param period period of the trignometric function
#'@param amplitude amplitude of thesinusoid function
#'@param phase phase of the trignometric function
#'@param shift intercept
#'@param noise standard deviation of noise which is added to the signal
#'@param type sine function or cosine function
#'
#'@export
#'@importFrom stats rnorm
sinusoid_generator <- function(knots, period, amplitude=1, phase=0, shift=0,
                               noise=0.1, type=c("sine", "cosine")){

    if (type == "sine"){
        myfunc <- sin
    } else{
        myfunc <- cos
    }

    signal <- myfunc(2 * pi / period * (knots - phase)) * amplitude + shift
    observed <- signal + rnorm(length(knots),mean=0, sd=noise)

    return(list(xvals=knots, yvals=observed))

}



#' Simulate from a linear function
#'
#' @param knots x values
#' @param slope slope
#' @param intercept intercept
#' @param noise standard deviation of noise which is added to the signal
#'
#' @export
#' @importFrom stats rnorm
linear_simulator <- function(knots, slope=0, intercept=0, noise=0.1){


    signal <- slope * knots + intercept
    observed <- signal + rnorm(length(knots),mean=0, sd=noise)

    return(list(xvals=knots, yvals=observed))

}


