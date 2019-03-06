# complete doc 

#' calculate how the cellular infection parameters change when a model parameter
#' is changed from the value for sH1N1-WT to that for H7N9-WT
#' 
#' @param rescale_par_name character:
#' "baseline", "p_inf", "delta", "c_inf", "k1", or "beta_inf".  If "baseline",
#' don't change any parameters and just calculate the value of the cellular 
#' infection parmaeter for sH1N1.  Otherwise, change the value of the named
#' parameter to that for H7N9, and calculate the value of the cellular
#' infectin parameter.
#' @param sample_idxs a vector of iteration numbers of the MCMC chains.
#' Each row of sample_idxs corresponds to an sH1N1 parameter set and a 
#' H7N9 parameter set.  So for each row of sample_idxs, we get one value for
#' each cellular infection parameter.
#' @return a named vector with elements "log10_R_0", "gen_time" and "r"
change_single_par_sH1N1_H7N9 <- function(rescale_par_name, sample_idxs) {
  transform_pars <- transform_pars_wrapper()
  sensitivity <- "main"
  dirs <- set_dirs(sensitivity)
  strains <- c("sH1N1-WT", "H7N9-WT")
  get_pars <- function(strain_idx) {
    strain <- strains[strain_idx]
    strain_dir <- get_strain_dir(dirs, strain)
    filename_in <- paste0(strain_dir, "1.RData")
    pars <- get_MCMC_chain(filename_in)
    pars <- as.data.frame(pars)
    pars <- pars[sample_idxs[strain_idx],]
    pars
  }
  pars <- lapply(seq_along(strains), get_pars)
  par_names <- colnames(pars[[1]])
  sH1N1_pars <- transform_pars(pars[[1]])
  H7N9_pars <- transform_pars(pars[[2]])
  
  if(rescale_par_name != "baseline") {
    sH1N1_pars[[rescale_par_name]] <- H7N9_pars[[rescale_par_name]]
    if(rescale_par_name == "beta_inf") {
      sH1N1_pars[["beta"]] <- H7N9_pars[["beta"]]
    }
  }
  
  log10_R_0 <- log10(calc_R_0(sH1N1_pars))
  gen_time <- calc_mean_gen_time_wrapper(sH1N1_pars)
  r <- calc_r_teiv_stages(sH1N1_pars)
  c("log10_R_0" = log10_R_0, "gen_time" = gen_time, "r" = r)
}

#' calculate how the cellular infection parameters change when model parameters
#' are changed from the value for sH1N1-WT to that for H7N9-WT
#' 
#' @return a data frame.
#' sum_stat is the name of the cellular infection parameter: "log10_R_0", "gen_time" or "r".
#' par_name is the name of the model parameter that is changed: 
#' "p_inf", "delta", "c_inf", "k1", or "beta_inf"
#' 2.5%, 50% and 97.5% are the quantiles of the percentage change in the cellular
#' infection parameter when the value of the model parameter is changed from
#' that of sH1N1-WT to that of H7N9-WT.
#' @importFrom magrittr %>%
calc_change_all_pars_sH1N1_H7N9_CI <- function() {
  calc <- TRUE
  sample_size <- 10
  sensitivity <- "main"
  rescale_par_names <- c("baseline", "p_inf", "delta", "c_inf", "k1", "beta_inf")
  
  strains <- c("sH1N1-WT", "H7N9-WT")
  dirs <- set_dirs(sensitivity)
  
  gen_strain_sample_idx <- function(strain) {
    strain_dir <- get_strain_dir(dirs, strain)
    filename_in <- paste0(strain_dir, "1.RData")
    pars <- get_MCMC_chain(filename_in)
    n_iter <- nrow(pars)
    sample.int(n_iter, size = sample_size)
  }
  
  sample_idxs <- vapply(strains, gen_strain_sample_idx, numeric(sample_size))
  
  
  change_single_par_sH1N1_H7N9_wrapper <- function(rescale_par_name) {
    apply(sample_idxs, 1, function(x) change_single_par_sH1N1_H7N9(rescale_par_name, x))
  }
  results_dir <- get_strain_dir(dirs, strains[1])
  sum_stats <- lapply(rescale_par_names,
                      change_single_par_sH1N1_H7N9_wrapper)
  
  par_names_plot <- c("baseline", "p_inf", "tau_I", "c_inf", "tau_L", "beta_inf")
  names(sum_stats) <- par_names_plot
  sum_stats <- lapply(sum_stats, function(x) x / sum_stats[["baseline"]]) 
  sum_stats <- sum_stats[-1]
  quantile_wrapper <- function(sum_stats, par_name) {
    sum_stats <- as.data.frame(t(apply(sum_stats, 1, quantile, probs = c(0.025, 0.5, 0.975))))
    sum_stats <- (sum_stats - 1) * 100
    sum_stats$sum_stat <- rownames(sum_stats)
    sum_stats$par_name <- par_name
    sum_stats
  }
  sum_stats <- Map(quantile_wrapper, sum_stats, par_names_plot[-1]) %>%
    do.call(rbind, .)
  filename_out <- paste0(results_dir, "sensitivity_sum_stats.rds")
  saveRDS(sum_stats, filename_out)
  sum_stats
}