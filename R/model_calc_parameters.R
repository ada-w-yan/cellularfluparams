# complete doc

#' calculate the mean generation time given model parameters
#' 
#' @param c1 rate of loss of free virus (due to loss of infectivity + entering target cell)
#' @param tau_E mean latent period
#' @param tau_I mean infectious period
#' @param n_L number of latent stages
#' @param n_I number of infectious stages
#' @param p_over_time "constant" for constant rate of virus production over
#' infected cell lifetime, "linear" for linearly increasing rate of virus production
#' over infected cell lifetime
#' @return numeric vector of length 1. mean generation time
calc_mean_gen_time <- function(c1, tau_E, tau_I, n_L, n_I,
                               p_over_time = "constant") {
  
  virus_mean_time <- 1 / c1
  mean_latent_period <- tau_E
  if(n_I == Inf) {
    mean_virion_production_time <- tau_I / 2
  } else if(p_over_time == "constant") {
    mean_virion_production_time <- (n_I + 1) / 2 / n_I * tau_I
  } else if(p_over_time == "linear") {
    mean_virion_production_time <- (2 * n_I + 1) / 3 / n_I * tau_I
  } else {
    stop("unknown value for p_over_time")
  }
  
  virus_mean_time + mean_latent_period + mean_virion_production_time
  
}

#' calculate the mean generation time given named vector of parameter values
#' 
#' @param pars named vector of parameter values
#' @return mean generation time
calc_mean_gen_time_wrapper <- function(pars) {
  calc_mean_gen_time(pars[["c_inf"]] + pars[["beta_inf"]] * pars[["T_0"]], 1 / pars[["k1"]], 1 / pars[["delta"]],
                     pars[["n_L"]], pars[["n_I"]])
}

#' calculates r for the TEIV model with Erlang-distributed latent and infectious periods
#' 
#' @inheritParams calc_mean_gen_time_wrapper
#' @inheritParams calc_mean_gen_time
#' @param warn_if_neg logical. if TRUE, throw warning if a negative value of r results
#' @return r: numeric vector of length 1
calc_r_teiv_stages <- function(pars, 
                               p_over_time = "constant", 
                               warn_if_neg = FALSE) {
  
  beta_T_0 <- pars[["beta"]] * pars[["T_0"]]
  p <- pars[["p_inf"]]
  c_inf <- pars[["c_inf"]] + pars[["beta_inf"]] * pars[["T_0"]]
  k <- pars[["k1"]]
  delta <- pars[["delta"]]
  n_L <- pars[["n_L"]]
  n_I <- pars[["n_I"]]

  make_LL_II_mat <- function(n_row, value) {
    if(!is.numeric(n_row)) {
      browser()
    }
    stopifnot(n_row > 0, isTRUE(all.equal(round(n_row), n_row)))
    
    if(n_row == 1) {
      matrix(-n_row * value)
    } else {
      vec <- double(n_row)
      vec[seq_len(2)] <- c(-n_row * value, n_row * value)
      mat <- toeplitz(vec)
      mat[upper.tri(mat)] <- 0
      mat
    }
  }
  
  LL_mat <- make_LL_II_mat(n_L, k)
  
  LI_mat <- matrix(0, nrow = n_L, ncol = n_I)
  
  IL_mat <- t(LI_mat)
  IL_mat[1,n_L] <- n_L * k
  
  II_mat <- make_LL_II_mat(n_I, delta)
  
  p_vec <- rep(p, n_I)
  if(p_over_time == "linear") {
    p_vec <- p_vec * seq_len(n_I) / n_I
  }
  V_mat_row <- matrix(c(double(n_L), p_vec, -c_inf), nrow = 1)
  V_mat_col <- matrix(c(beta_T_0, double(n_L + n_I - 1)), ncol = 1)
  mat1 <- rbind(LL_mat, IL_mat)
  mat2 <- rbind(LI_mat, II_mat)
  mat3 <- cbind(mat1, mat2, V_mat_col)
  
  eigenvalue_mat <- rbind(mat3, V_mat_row)
  ev <- eigen(eigenvalue_mat)
  ev <- ev$values
  ev <- ev[Im(ev) == 0]
  stopifnot(length(ev) > 0)
  r <- max(as.double(ev))
  if(warn_if_neg && r <= 0) {
    warning("negative value of r")
  }
  r
  
}

#' calculates R_0 for the TEIV model with Erlang-distributed latent and infectious periods
#' 
#' @inheritParams calc_mean_gen_time_wrapper
#' @inheritParams calc_mean_gen_time
#' @return R_0: numeric vector of length 1
calc_R_0 <- function(pars, p_over_time = "constant") {
  n_I <- pars$n_I
  p_vec <- rep(pars$p_inf, n_I)
  if(p_over_time == "linear") {
    p_vec <- p_vec * seq_len(n_I) / n_I
  }
  p <- mean(p_vec)
  pars$beta * p * pars$T_0 / pars$delta /
          (pars$c_inf + pars$beta_inf * pars$T_0)
}