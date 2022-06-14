# complete doc

#' postprocessing of MCMC output
#' 
#' @param filename character vector of length 1: name of .RData file outputted by MCMC
#' @param run_parallel logical. Whether to run calculations in parallel
#' @param gen_summary_statistics closure.  NULL or function to calculate summary statistics
#' @param calc_residuals closure. NULL or function to calculate residuals
#' @param get_prediction_ci_df closure. NULL or function to calculate credible
#' intervals for model predictions
postprocess_for_plotting <-function(filename, 
                                    run_parallel = TRUE, 
                                    gen_summary_statistics = NULL, 
                                    calc_residuals = NULL, 
                                    get_prediction_ci_df = NULL){
  
  # careful loading of MCMC output
  loaded_output <- load_into_list(filename)
  loaded_output <- loaded_output[c("parTab", "data", "CREATE_POSTERIOR_FUNC", "CREATE_PRIOR_FUNC", "startTab", "diagnostics", "output", "filenames", "mcmcPars")]
  list2here(loaded_output)
  remove(loaded_output)
  
  PRIOR_FUNC <- CREATE_PRIOR_FUNC(parTab)
  f <- CREATE_POSTERIOR_FUNC(parTab, data, PRIOR_FUNC)
  
  prior_chain <- sample_from_prior(PRIOR_FUNC, parTab, startTab, data, filenames,
                                   CREATE_PRIOR_FUNC, run_parallel)
  
  ## read in chain
  
  current_dir <- dirname(filename)
  filenames <- paste(current_dir, basename(filenames), sep = "/")
  chain_filenames <- lapply(output, function(x) paste(current_dir, basename(x$file), sep = "/"))
  
    chain <- lapply(chain_filenames, read_csv_or_rds)
    chain <- lapply(chain,data.table::as.data.table)
    
    # discard values from adaptive period and burn-in
    chain <- Map(function(x, y) x[sampno > (y$adaptive_period +
                                              diagnostics$burn_in * mcmcPars[["thin"]]),],
                 chain, output)

  
  # discard starting values from higher temperatures
  parallel_tempering_flag <- "temperature" %in% names(mcmcPars) && 
    length(mcmcPars[["temperature"]]) > 1
  
  if(parallel_tempering_flag && !is.data.frame(startTab[[1]])){
    startTab <- lapply(startTab, function(x) x[[1]])
  }
  
  ## trace plots need to happen before converged chains are combined
  
  # calc summary statistics
  if(!is.null(gen_summary_statistics)){
    
    summary_list <- gen_summary_statistics(parTab, startTab)
    calc_summary <- summary_list$calc_summary

    augment_Tab_with_summary <- function(parTab, summary_list) {
      augment_values <- calc_summary(parTab$values)
      augment_frame <- data.frame("values" = as.numeric(augment_values),
                                  "names" = names(augment_values),
                                  "fixed" = numeric(length(augment_values)),
                                  "lower_bound" = summary_list$lower_bound,
                                  "upper_bound" = summary_list$upper_bound,
                                  "steps" = numeric(length(augment_values)),
                                  "names_plot" = summary_list$par_names_plot,
                                  stringsAsFactors = FALSE)
      rbind(parTab, augment_frame)
    }
    augment_Tab_with_summary_partial <- function(x) augment_Tab_with_summary(x, summary_list = summary_list)
    
    # if list of data frame, lapply, otherwise just do on single data frame
    if(!is.data.frame(startTab)){
      startTab <- lapply(startTab, augment_Tab_with_summary_partial)
    } else{
      startTab <- augment_Tab_with_summary_partial(startTab)
    }
    
    original_par_names <- parTab$names
    parTab <- augment_Tab_with_summary_partial(parTab)
    
    bind_summary_to_chain <- function(chain, original_par_names){
      # determine number of summary statistics; if this is 1, don't transpose

      summary_stats <- apply(chain, 1,
                             function(x) calc_summary(as.numeric(x[original_par_names])))
      
      if(!is.matrix(summary_stats)) {
        summary_temp <- calc_summary(as.numeric(chain[1, original_par_names]))
        summary_stats <- t(as.matrix(summary_stats))
        rownames(summary_stats) <- names(summary_temp)
      }
      cbind(chain,t(summary_stats))
    }
    
    bind_summary_to_chain_partial <- function(x) bind_summary_to_chain(x,
                                                                       original_par_names = original_par_names)
    
    # calculate priors for summary statistics
    prior_chain <- bind_summary_to_chain_partial(prior_chain)
    
  } else {
    bind_summary_to_chain_partial <- NULL
  }
  
  ## if converged, combine chains
  if(diagnostics$converged && length(filenames) > 1){
    
      chain <- do.call(rbind,chain)

    invisible(postprocess_chain(chain, parTab, data, bind_summary_to_chain_partial,
                                calc_residuals, get_prediction_ci_df, filenames[1], combined_flag = TRUE))
  } else {
    chain_residual_list <- Map(function(x, y) postprocess_chain(x, parTab, data, bind_summary_to_chain_partial,
                                                                calc_residuals, get_prediction_ci_df, y, combined_flag = FALSE),
                               chain, filenames)
  }
  
  invisible(NULL)
}

#' more postprocessing of MCMC output
#' 
#' @param chain data frame containing MCMC chain
#' @param parTab data frame of model parameters created by specify_parameters
#' @param data data frame of data used for fitting in MCMC
#' @param bind_summary_to_chain_partial closure created by postprocess_for_plotting
#' used to calculate summary statistics and bind them to chain
#' @inheritParams postprocess_for_plotting
#' @param combined_flag logical.  INdicates whether this chain is combined from 
#' parallel chains that have converged
postprocess_chain <- function(chain, parTab, data, bind_summary_to_chain_partial,
                              calc_residuals, get_prediction_ci_df, filename, combined_flag){
  
  if(combined_flag){
    combined_str <- "_combined"
  } else {
    combined_str <- ""
  }  
  
  if(!is.null(bind_summary_to_chain_partial)) {
    chain <- bind_summary_to_chain_partial(chain)
  }
  
  if(!is.null(calc_residuals)) {
    residuals <- calc_residuals(data, chain, parTab)
  } else {
    residuals <- NULL
  }
  
  if(!is.null(get_prediction_ci_df)) {
    prediction_ci_df <- get_prediction_ci_df(chain, data)
    saveRDS(prediction_ci_df, paste0(filename, "ci.rds"))
  }

  unfixed_pars <- parTab$fixed == 0 & !is.na(parTab$values)
  par_names_unfixed <- parTab$names[unfixed_pars]
        # this is the effective sample size of the concatenated chain, which is
        # different to the sum of the effective sample sizes of each chain, 
        # which was calculated by lazymcmc
  effective_size_summary <- coda::effectiveSize(mcmc(chain[,par_names_unfixed,with = FALSE]))
  write(effective_size_summary, paste0(filename,".diagnostics"), append=TRUE, ncolumns=1000, sep = ",")

  par_names_plot_unfixed <- parTab$names_plot[unfixed_pars]
  prctile_table <- print_prctiles(chain[,colnames(chain) %in% par_names_unfixed, with = FALSE],
                                  par_names_plot = par_names_plot_unfixed)
  print(xtable::xtable(prctile_table),
  sanitize.text.function = identity,
  file = paste0(filename,"_prctile",combined_str,".tex"))

  invisible(NULL)
}

#' samples from joint prior distribution
#' 
#' @param PRIOR_FUNC NULL for uniform prior, otherwise closure specifying prior
#' @inheritParams postprocess_for_plotting
#' @param startTab data frame of smae format as parTab, storing starting
#' parameter values
#' @param filenames vector of filenames in which to temporarily store samples
#' @param CREATE_PRIOR_FUNC NULL or closure to create PRIOR_FUNC -- redundant
#' @return data frame of samples from the prior distribution
sample_from_prior <- function(PRIOR_FUNC, parTab, startTab, data, filenames,
                              CREATE_PRIOR_FUNC, run_parallel){
  if(is.null(PRIOR_FUNC)){ # if uniform prior, sample with LHS
    
    prior_chain_length <- 1e4
    prior_chain <- matrix(rep(parTab$values,prior_chain_length),
                          nrow = prior_chain_length,
                          byrow = TRUE)
    prior_chain <- data.frame(prior_chain)
    colnames(prior_chain) <- parTab$names
    prior_chain[,parTab$fixed == 0] <- 
      lhs::randomLHS(prior_chain_length,sum(parTab$fixed == 0))
    prior_chain[,parTab$fixed == 0] <- 
      lapply(which(parTab$fixed == 0),
             function(x) prior_chain[,x] * 
               (parTab[x,"upper_bound"] - parTab[x,"lower_bound"]) + parTab[x,"lower_bound"])
    
  } else { # else sample using MCMC
    
    mcmcPars_prior <- gen_default_mcmcPars(iterations=1e4,opt_freq=1e3,adaptive_period=0)
    filenames_prior <- paste0(filenames,"_prior")
    
    prior_list <- run_MCMC_loop(startTab, data, mcmcPars_prior, filenames_prior,
                                CREATE_PRIOR_FUNC, PRIOR_FUNC,
                                run_parallel = run_parallel)
    prior_output <- prior_list$output
    rm(prior_list)
    
    prior_chain <- lapply(prior_output, function(x) read_csv_or_rds(x$file))
    prior_chain <- do.call(rbind,prior_chain)
    lapply(prior_output, function(x) file.remove(x$file))
  }
  prior_chain <- data.table::as.data.table(prior_chain)
  prior_chain
}