# complete doc

#' calculate the mean generation number at the time of peak infectious viral load,
#' for a given set of parameters
#' 
#' @param pars a named vector of parmaeter values
#' @param wash logical.  If FALSE, exclude wash step for multi-cycle experiment.
#' If TRUE, do whatever the original multi-cycle experiment did.
#' @return a numeric vector of length 1: the mean generation number for that parameter set
calc_mean_gen_no <- function(pars, wash) {
  calc_generations <- calc_generations_closure(wash = wash)
  generations <- calc_generations(pars)
  generations <- generations/sum(generations)
  gen_nos <- seq_along(generations)
  mean_gen_no <- calc_moment_from_pmf(gen_nos, generations, 1)
  mean_gen_no
}

#' calculate the distribution of the mean generation number over the posterior
#' for a given strain
#' @inheritParams run_exp
#' @return a vector of samples from the distribution of the mean generation number
#' @export 
calc_mean_gen_no_by_strain <- function(strain, sensitivity) {
  dirs <- set_dirs(sensitivity)
  data_source <- determine_data_source_for_strain(strain)
  # take the wash step out for PR8 strains (wash step already absent for WT strains)
  wash <- (data_source == "WT")
  strain_dir <- get_strain_dir(dirs, strain)
  filename_in <- paste0(strain_dir, "1.RData")
  pars <- get_MCMC_chain(filename_in)
  parnames <- c("log10_c_inf", "obs_threshold")
  par_idx <- vnapply(parnames, function(x) which(x == colnames(pars)))
  pars <- as.data.frame(pars)
  pars <- pars[,seq(par_idx[1], par_idx[2])]
  mean_gen_no <- apply(pars, 1, function(x) calc_mean_gen_no(x, wash))
  saveRDS(mean_gen_no, paste0(strain_dir, "mean_gen_no.rds"))
  invisible(mean_gen_no)
}

#' calculate the generation number distribution at the time of peak viral load
#' for the maximum likelihood parameters for a given strain
#' @inheritParams run_exp
#' @return a data frame with three columns.  strain contains the name of the strain;
#' gen has generation numbers 1:20; and V has the proportion of infectious virions
#' belonging to that generation at the time of peak viral load.
#' @export 
calc_generation_strain <- function(strain, sensitivity){
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strain)
  filename_in <- paste0(strain_dir, "max_LL_params.rds")
  pars <- readRDS(filename_in)
  data_source <- determine_data_source_for_strain(strain)
  # take the wash step out for PR8 strains (wash step already absent for WT strains)
  wash <- (data_source == "PR8")
  p_over_time <- get_p_over_time_from_sensitivity(sensitivity)
  calc_generations <- calc_generations_closure(wash)
  generations <- calc_generations(pars)
  plot_df <- data.frame(gen = seq_along(generations), V = generations/sum(generations))
  plot_df$strain <- strain
  plot_df
}

#' make a function to calculate the generation number distribution at the time 
#' of peak infectious viral load, for a given set of parameters
#' @inheritParams calc_mean_gen_no
#' @useDynLib cellularfluparams
#' @return a function
calc_generations_closure <- function(wash) {
  G <- 20
  tol = 1e-6
  # path1 <- paste0("model_teiv_stages_gen_G_inoculum.R")
  
  # gen <- compile_model(path1)
  gen <- model_teiv_stages_gen_G_inoculum
  
  transform_pars <- transform_pars_wrapper(G = G)
  
  # function to calculate the viral load at the start of the incubation period
  calc_V_inf_0 <- function(pars_incubation) {
    pars_incubation[["MOI"]] * pars_incubation[["T_0"]] / pars_incubation[["supernatant_vol"]]
  }
  
  ## function to extract the number of latent and infectious cells at t = 0
  extract_iv <- function(sol, col_name) {
    iv <- unname(sol[2,grep(col_name,colnames(sol))])
    if(length(iv) > 1) {
      iv <- matrix(iv, ncol = G)
    }
    iv[iv <= tol] <- 0
    iv
  }
  
  calc_generations <- function(pars) {
    pars <- transform_pars(pars)
    pars$G <- G
    
    mod <- gen(user=pars)
    
    pars_incubation <- pars
    ## calculate starting viral load for incubation (t = -1)
    pars_incubation$V_inf_0 <- calc_V_inf_0(pars_incubation)
    do.call(mod$set_user, pars_incubation)
    t_end <- 168
    if(wash) {
      sol_incubation <- solve_ODE_error_handling(mod, c(-pars_incubation$incubation_period, 0), tol)
      
      ## extract number of latent and infectious cells at t = 0
      pars$L_0 <- extract_iv(sol_incubation, "L")
      pars$I_0 <- extract_iv(sol_incubation, "I")
      pars$T_0 <- extract_iv(sol_incubation, "T")
      
      do.call(mod$set_user,pars)
      # solve ODE at these times
      t1 <- seq(0, t_end)
    } else {
      t1 <- seq(-1, t_end)
    }
    
    # solve ODE which tracks generations
    sol <- as.data.frame(solve_ODE_error_handling(mod, t1, tol = 1e-10))
    
    # sum viral load across generations to determine the overall peak time
    viral_load <- sum_stages(sol, "V")
    time_idx <- which.max(viral_load)
    # get viral load in each generation at the peak time
    generations <- as.numeric(sol[time_idx, grep("V", colnames(sol))])
    return(generations)
  }
  calc_generations
}

#' sum time series of odin solution across compartments defined as vectors
#' 
#' @param df data frame containing odin solution
#' @param compartment_names character vector. Names of compartments to sum over
#' @return a data frame with one column for each element of compartment_names,
#' containing the solution summed over that vector compartment
sum_stages <- function(df, compartment_names){
  extract_and_sum_stages <- function(df, compartment_name) {
    extracted_df <- df[,grep(compartment_name, colnames(df)), drop = FALSE]
    if(ncol(extracted_df) == 0) {
      double(nrow(df))
    } else {
      apply(extracted_df, 1, sum)
    }
  }
  
  summed_df <- vapply(compartment_names, extract_and_sum_stages, numeric(nrow(df)), df = df)
  colnames(summed_df) <- compartment_names
  summed_df
}

#' calculate the generation number distribution as cellular infection parameters are varied
#' 
#' @param constant_param character vector: "log10_R_0", "gen_time" or "r"
#' @return a list with three objects.
#' plot_df is a data frame with three columns. gen is the generation number (1:20),
#' V is the proportion of infectious virions belonging to that generation at 
#' the time of peak viral load, and infection_par_value is the value of the
#' independent variable (which depends on constant_param).
#' all_pars is a list with three elements.
#' independent is a vector of values for the independent parameter,
#' dependent is a vector of values for the dependent parameter,
#' and constant is a vector of (repeating) values for the constant parameter.
#' Which parameter out of log10_R_0, gen_time and r is independent, dependent
#' or constant depends on constant_param.
#' sf is the number of significant figures to which the values in all_pars are
#' given (used for plotting later)
#' @export
calc_gen_dist_vary_params <- function(constant_param) {
  strain <- "pH1N1-PR8"
  sensitivity <- "main"
  data_source <- determine_data_source_for_strain(strain)
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strain)
  filename_in <- paste0(strain_dir, "max_LL_params.rds")
  pars <- readRDS(filename_in)
  
  transform_pars <- transform_pars_wrapper()
  
  transformed_pars <- transform_pars(pars)
  
  if(constant_param == "gen_time") {
    target_log10_R_0 <- seq(2, 5)
    current_log10_R_0 <- log10(calc_R_0(transformed_pars))
    
    change_p <- function(pars, current_log10_R_0, target_log10_R_0) {
      pars[["log10_p_inf_mc"]] <- pars[["log10_p_inf_mc"]] + target_log10_R_0 - 
        current_log10_R_0
      pars
    }
    
    pars <- lapply(target_log10_R_0, change_p, 
                   pars = pars, 
                   current_log10_R_0 = current_log10_R_0)
    names(pars) <- target_log10_R_0
    dependent_pars <- vnapply(pars, 
                              function(x) calc_r_teiv_stages(transform_pars(x)))
    constant_par_value <- calc_mean_gen_time_wrapper(transformed_pars)
    
  } else if(constant_param == "log10_R_0") {
    target_gen_time <- seq(20, 50, by = 10)
    
    current_gen_time <- calc_mean_gen_time_wrapper(transformed_pars)
    
    change_burst_size <- function(pars, current_gen_time, target_gen_time) {
      tau_I <- 10^pars[["log10_tau_I"]]
      target_tau_I <- tau_I + target_gen_time - current_gen_time
      pars[["log10_p_inf_mc"]] <- log10(10^pars[["log10_p_inf_mc"]] / target_tau_I * tau_I)
      pars[["log10_tau_I"]] <- log10(target_tau_I)
      pars
    }
    
    pars <- lapply(target_gen_time, change_burst_size, 
                   pars = pars, 
                   current_gen_time = current_gen_time)
    names(pars) <- target_gen_time
    
    dependent_pars <- vnapply(pars, 
                              function(x) calc_r_teiv_stages(transform_pars(x)))
    constant_par_value <- log10(calc_R_0(transformed_pars))
    
  } else if(constant_param == "r") {
    target_gen_time <- seq(20, 50, by = 10)
    current_gen_time <- calc_mean_gen_time_wrapper(transformed_pars)
    
    current_r <- calc_r_teiv_stages(transformed_pars)
    
    calc_new_p_inf <- function(pars) {
      calc_r_diff_given_pars_and_new_value <- function(new_value) {
        pars[["log10_p_inf_mc"]] <- new_value
        
        r_diff <- calc_r_teiv_stages(transform_pars(pars)) - current_r
        r_diff
      }
      
      log10_par_range <- c(-2, 3)
      log10_par_vec <- seq(log10_par_range[1], log10_par_range[2], length.out = 100)
      diff_range <- vnapply(log10_par_vec, calc_r_diff_given_pars_and_new_value)
      
      sign_change_ind <- find_sign_change_ind(diff_range)
      guess <- log10_par_vec[sign_change_ind + c(0, 1)]
      pars_new <- uniroot(calc_r_diff_given_pars_and_new_value, guess)$root
      pars[["log10_p_inf_mc"]] <- pars_new
      pars
    }
    change_p_delta <- function(pars, target_gen_time, current_r) {
      tau_I <- 10^pars[["log10_tau_I"]]
      target_tau_I <- tau_I + target_gen_time - current_gen_time
      pars[["log10_tau_I"]] <- log10(target_tau_I)
      pars <- calc_new_p_inf(pars)
      pars
    }
    
    pars <- lapply(target_gen_time, change_p_delta, 
                   pars = pars, 
                   current_r = current_r)
    names(pars) <- target_gen_time
    
    dependent_pars <- vnapply(pars, 
                              function(x) log10(calc_R_0(transform_pars(x))))
    constant_par_value <- calc_r_teiv_stages(transformed_pars)
  } else {
    stop("unknown option for constant_param")
  }
  
  calc_generations_wrapper <- function(infection_par_value, pars) {
    calc_generations <- calc_generations_closure(FALSE)
    generations <- calc_generations(pars)
    
    plot_df <- data.frame(gen = seq_along(generations), V = generations/sum(generations))
    
    plot_df$infection_par_value <- as.numeric(infection_par_value)
    plot_df
  }
  
  plot_df <- lapply_w_name(pars, calc_generations_wrapper)
  plot_df <- do.call(rbind, plot_df)
  plot_df$infection_par_value <- as.factor(plot_df$infection_par_value)
  
  sf <- 3
  all_pars <- list(independent = names(pars),
                   dependent = signif(dependent_pars, sf),
                   constant = rep(as.character(signif(constant_par_value, sf)), 4))
  list(plot_df = plot_df, all_pars = all_pars, sf = sf)
}

#' calculation for Fig. S5
#' for either the model where the viral load exponentially increases then stays constant,
#' or the model where the viral load exponentially increases then exponentially
#' decreases at the same rate, calculate the number of virions in each generation at a given time
#' @param time numeric vector of length 1.  time at which to calculate the number of virions in each generation
#' @param gen_time numeric vector of length 1.  fixed generation time
#' @param viral_load_function function.  
#' only calc_generations_fixed_growth_until_I_max is implemented currently
#' @inheritParams calc_generations_fixed_growth_until_I_max
#' @return numeric vector containing the number of virions in each generation at a given time.
#' @export
calc_generations_at_time_fixed <- function(time, gen_time = 1, V_0, R_0, I_max, viral_load_fn) {
  # work out generation number at the time
  n_gens <- floor(time / gen_time + 1)
  # calculate viral load at that generation number
  V <- viral_load_fn(n_gens, V_0, R_0, I_max)
  # make a vector n_gens long, and put the viral load at the end.  
  # all other entries are zero
  V_vec <- double(n_gens)
  V_vec[n_gens] <- V[n_gens]
  return(V_vec)
}

# calculation for Fig. S5

#' for a model where the viral load exponentially increases then stays constant, 
#' calculate the number of virions in each generation
#' @param n_gens numeric vector of length 1. number of generations to calculate up until
#' @param V_0 numeric vector of length 1.  Inoculum size
#' @param R_0 numeric vector of length 1.  Basic reproduction number
#' @param I_max numeric vector of length 1.  Maximum number of infected cells
#' @return numeric vector of length n_gens.  Number of virions in each generation
calc_generations_fixed_growth_until_I_max <- function(n_gens, V_0, R_0, I_max) {
  V <- double(n_gens) # vector to store viral load
  V[1] <- V_0 # inoculum
  if(n_gens > 1) {
    for(i in seq(2, n_gens)) { # in each generation
      I <- min(V[i - 1], I_max) # each virion from the previous generation 
      # goes into a cell.  Infected cells magically renew during each generation...
      V[i] <- I * R_0 # each infected cell produces R_0 virions
    }
  }
  V
}