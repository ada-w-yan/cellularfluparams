# complete doc

#' Analogue for dlnorm but in base 10
#'
#' \code{dl10norm} is an analogue for dlnorm but in base 10
#' 
#' @param x a numeric vector of quantiles.
#' @param mean a vector of means.
#' @param sd a vector of standard deviations.
#' @param log.p a logical vector of length 1.  
#' TRUE: return log of density, 
#' FALSE: return density
#' @return a numeric vector with the density of x in a log10normal distribution.
#' The log10normal distribution in log10 space has mu =  mean and sigma = sd.

dl10norm <- function(x, mean = 1, sd = 1, log.p = FALSE){
  dnorm(log10(x),mean,sd,log.p)
}

#' Analogue for plnorm but in base 10
#'
#' \code{pl10norm} is an analogue for plnorm but in base 10
#' 
#' @param q a numeric vector of quantiles.
#' @param mean a vector of means.
#' @param sd a vector of standard deviations.
#' @param lower.tail TRUE: calculate the probability of drawing a number less than q.
#' FALSE: = calculate the probability of drawing a number larger than q.
#' @param log a logical vector of length 1.  
#' TRUE: return log of probability, 
#' FALSE: return probability
#' @return a numeric vector with the probability of drawing a number less than 
#' (or greater than, depending on the value of \code{lower.tail})
#' in a log10normal distribution.
#' The log10normal distribution in log10 space has mu =  mean and sigma = sd.

pl10norm <- function(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE){
  pnorm(log10(q),mean,sd,lower.tail,log.p)
}

#' Analogue for rlnorm but in base 10
#'
#' \code{rl10norm} is an analogue for rlnorm but in base 10
#' 
#' @param n a numeric vector of length 1: number of observations
#' @param mean a vector of means.
#' @param sd a vector of standard deviations.
#' @return a numeric vector of length n, distributed according to the log10normal distribution.
#' The log10normal distribution in log10 space has mu =  mean and sigma = sd.

rl10norm <- function(n, mean = 1, sd = 1){
  10^rnorm(n,mean,sd)
}

#' Calculate the log likehood of data given true model outputs if the data is log10normal distributed
#' and there is an observation threshold
#'
#' \code{calc_LL_log10normal_threshold} calculates the log likehood of data given
#' true model outputs if the data is log10normal distributed
#' and there is an observation threshold
#' 
#' @param data numeric vector of observed data (where some number below \code{threshold}
#' is used to denote a measurement below the threshold)
#' @param true_values numeric vector of the same length as \code{data} which contains
#' the true model outputs
#' @param sd_noise numeric vector of length 1:
#' standard deviation of the noise distribution in log10 space
#' @param threshold numeric vector of length 1: observation threshold
#' @return numeric vector of length 1: the log likehood of the data given true 
#' model outputs
calc_LL_log10normal_threshold <- function(data,true_values,sd_noise,threshold){
  if(sd_noise == Inf) {
    return(0)
  }
  log10_true_values <- log10(true_values)
  above_threshold <- data >= threshold

  sum(dl10norm(data[above_threshold],log10_true_values[above_threshold],sd_noise,TRUE)) +
    # treat data below threshold as censored
    sum(pl10norm(threshold,log10_true_values[!above_threshold],sd_noise,
                 lower.tail = TRUE, log.p = TRUE))
}