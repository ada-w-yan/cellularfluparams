# compelete doc

#' Generate random starting values from a parameter table
#'
#' \code{generate_random_start} generates random starting values for all parameters to be fitted
#' 
#' @param parTab a parameter table data frame of the form needed by the \code{run_mcmc} function
#' @param random_seed random seed for generation
#' The info used here are the columns
#' parTab$fixed: 1 indicates the parameter is fixed, 0 indicates the parameter is to be fitted
#' parTab$lower_bound: lower bounds for parameters allowed for proposal steps
#' parTab$upper_bound: lower bounds for parameters allowed for proposal steps
#' @return a data frame of the same for as parTab, but the values of the fitted 
#' parameters are chosen from a uniform distribution with bounds [lower_bound, upper_bound]
generate_random_start <- function(parTab,random_seed){
    set.seed(random_seed)
    startTab <- parTab
    for(i in which(parTab$fixed == 0)){
        startTab[i,"values"] <- runif(1,parTab[i,"lower_bound"],parTab[i,"upper_bound"])
    }
    startTab
}

#' set default random number generator seed for each
set_default_seeds <- function(n_replicates, n_temperatures) {
    seeds <- lapply(seq_len(n_replicates), 
                    function(x) n_temperatures * (x -1) + seq_len(n_temperatures))
    seeds
    
}
#' Generate vector of times at which to solve ODE for fitting
#' 
#' \code{gen_time_list} generates a vector of times at which to solve ODE for fitting
#' 
#' @param sampling_time a numeric vector of length n: times at which we have data
#' @return a list with the elements
#' mcmc_solving_time: a numeric vector: if sampling time includes 0, mcmc_solving_time = sampling_time,
#' otherwise append 0
#' sampling_ind: a numeric vector containing the indices of mcmc_solving_time which are in sampling_time
gen_time_list <- function(sampling_time){
    if(sampling_time[1] != 0){
        mcmc_solving_time <- c(0,sampling_time) # solve for 0 <= t <= 11, 100 timesteps
    } else{
        mcmc_solving_time <- sampling_time
    }
    sampling_ind <- mcmc_solving_time %in% sampling_time
    list(mcmc_solving_time = mcmc_solving_time, sampling_ind = sampling_ind)
}



#' generate default MCMC parameters
#' 
#' \code{gen_default_mcmcPars} supplies some default MCMC parameters.
#' 
#' All arguments are optional (to overwrite defaults)
#' 
#' @param iterations numeric vector of length 1: number of iterations for which to run MCMC (excluding adaptive period)
#' @param popt: numeric vector of length 1: optimum acceptance ratio for adapting step size
#' @param opt_freq: numeric vector of length 1: adjust step size per this many number of iterations during adaptive period
#' @param thin: numeric vector of length 1: save state of chain and parameter values to file every this many iterations
#' @param adaptive_period: numeric vector of length 1: number of iterations for which to adapt MCMC
#' @param save_block: numeric vector of length 1: number of iterations to hold in memory
#' @return a named vector with the above elements
#' @export
gen_default_mcmcPars <- function(iterations=10000,
                                 popt=0.44,
                                 opt_freq=1000,
                                 thin=1,
                                 adaptive_period=5000,
                                 save_block=100){
    list("iterations"=iterations,"popt"=popt,"opt_freq"=opt_freq,"thin"=thin,"adaptive_period"=adaptive_period,"save_block"=save_block)
}

#' generate default MCMC parameter values
#' 
#' @param temperature numeric vector of temperatures for parallel tempering
#' @param parallel_tempering_iter numeric vector of length 1: number of
#' iterations before attempting to swap adjacent chains
#' @return list of default parameter values
#' @export
gen_mcmcPars <- function(temperature = 1, parallel_tempering_iter){
    mcmcPars <- gen_default_mcmcPars(iterations=2e5,opt_freq=1e4,adaptive_period=2e5,thin = 100)
    mcmcPars<- c(mcmcPars, list("max_adaptive_period" = 2e5, "adaptiveLeeway" = 0.2))
    mcmcPars <- c(mcmcPars, list("max_total_iterations" = 1e7, "temperature" = temperature))
    if(!missing(parallel_tempering_iter)){
        mcmcPars <- c(mcmcPars, list("parallel_tempering_iter" = parallel_tempering_iter))
    }
    mcmcPars
}

#' generate default MCMC parameter values for a short chain
#' 
#' @inheritParams gen_mcmcPars
#' @return list of default parameter values for a short chain
#' @export
gen_mcmcPars_short <- function(temperature = 1, parallel_tempering_iter){
    mcmcPars <- gen_default_mcmcPars(iterations=1e4,opt_freq=1e3,adaptive_period=2e3,thin = 10)
    mcmcPars<- c(mcmcPars, list("max_adaptive_period" = 2e3, "adaptiveLeeway" = 0.2))
    mcmcPars <- c(mcmcPars, list("max_total_iterations" = 1e4, "temperature" = temperature))
    if(!missing(parallel_tempering_iter)){
        mcmcPars <- c(mcmcPars, list("parallel_tempering_iter" = parallel_tempering_iter))
    }
    mcmcPars
}

#' produce confidence intervals for model predictions 
#' 
#' \code{gen_prediction_ci} produces confidence intervals for model predictions 
#' 
#' @param data_df data frame inputted into MCMC inference
#' @param chain: data frame resulting from reading in csv file outputted from \code{run_mcmc}
#' should discard burn-in before running this function.
#' @param prediction_names: character vector: column names in \code{chain} which
#' contain the model predictions
#' @param sampling_time: numeric vector: times at which data is sampled
#' @return a data frame which consists of columns containing the CIs
#' for the model predictions at the sampling times (padded with NAs if necessary)
gen_prediction_ci <- function(data_df, chain, prediction_names, 
             sampling_time){
        data_time  <- unique(data_df$t)
        # extract model predictions from chain data frame
        
        chain <- chain[, prediction_names, with = FALSE] 
        # calculate CIs
        dynamics <- apply(chain,2,function(y) quantile(y,probs = c(.025,.975)))
        # pad with NAs and bind to data_df
        dynamics_df <- matrix(NA, nrow = nrow(data_df), ncol = nrow(dynamics))
        
        # if n_measurements conducted at each time point, plot once only at each time point
        match_idx <- unique(match(sampling_time,data_df$t))
        dynamics_df[match_idx,] <- t(dynamics)
        colnames(dynamics_df) <- rownames(dynamics)
        dynamics_df
}


#' reads summary stats from a tex file generated by running MCMC
#' 
#' @param filename character vector of length 1
#' @param param_name name of parameter for which to retrieve summary stats
#' @return double vector containing summary stats
read_summary_stats <- function(filename, param_name) {
    text_in <- readLines(con = filename)
    param_idx <- grep(param_name,text_in, fixed = TRUE)
    param_text <- text_in[param_idx]
    param_values <- strsplit(param_text, "&")
    param_values <- tryCatch(param_values[[1]][-1], error = function(e) NA)
    param_values <- as.numeric(gsub("\\\\", "", param_values))
    param_values
}

#' given an RData file, load the associated MCMC chains and find the iterations with maximum likelihood
#' 
#' @param filename_in character vector of length 1: name of .RData file outputted by MCMC
#' @param n_par_sets numeric vector of length 1: find the n_par_sets parameter sets
#' with the highest likelihood
#' @return data frame with iterations as rows
extract_max_LL_iters <- function(filename_in, n_par_sets) {
    chain <- get_MCMC_chain(filename_in)
    if(n_par_sets == 1) {
        max_LL_iters <- chain[which.max(chain$lnlike),]
    } else {
        max_LL_iters <- chain[order(chain$lnlike, decreasing=TRUE)[seq_len(n_par_sets)],]
    }
    max_LL_iters
}

#' given an RData file, load the associated MCMC chains and find the parameters with maximum likelihood
#' 
#' @inheritParams extract_max_LL_iters
#' @return vector of parameter values
get_max_LL_params <- function(filename_in) {
    load(filename_in)
    par_names <- parTab$names
    loaded_vars <- ls()
    loaded_vars <- loaded_vars[!(loaded_vars %in% c("par_names", "filename_in"))]
    rm(list = loaded_vars)
    max_LL_iters <- extract_max_LL_iters(filename_in, n_par_sets = 1)
    colnames_max_LL_iters <- colnames(max_LL_iters)
    max_LL_iters <- as.numeric(max_LL_iters)
    par_names <- intersect(colnames_max_LL_iters, par_names)
    max_LL_params <- max_LL_iters[colnames_max_LL_iters %in% par_names]
    names(max_LL_params) <- par_names
    max_LL_params
}

#' extract chain from MCMC output
#' 
#' \code{get_MCMC_chain} extracts chain from MCMC output, discarding burn-in
#' @inheritParams extract_max_LL_iters
#' @return chain data table containing MCMC runs
#' @import data.table
#' @export
get_MCMC_chain <- function(filename_in) {
    load(filename_in)
    data_dir <- dirname(filename_in)
    
    chain_filenames <- lapply(output, function(x) paste0(data_dir, "/", basename(x$file)))
    chain <- lapply(chain_filenames, read_csv_or_rds)
    chain <- lapply(chain,data.table::as.data.table)
    chain <- lapply(seq_len(n_replicates),
                    function(x) chain[[x]]
                    [sampno > (output[[x]]$adaptive_period +
                                   diagnostics$burn_in * mcmcPars[["thin"]]),])
    chain <- do.call(rbind, chain)
}

#' given an RData file, load the associated MCMC chains, find the parameters with maximum likelihood, and save
#' 
#' @inheritParams extract_max_LL_iters
#' @return vector of parameter values
#' @export
get_and_save_max_LL_params <- function(filename_in) {
    params <- get_max_LL_params(filename_in)
    current_dir <- dirname(filename_in)
    saveRDS(params, file = paste0(current_dir, "/max_LL_params.rds"))
    params
}

#' select evenly spaced iterations from MCMC chain
#' 
#' @param chain data table containing MCMC iterations
#' @param n_samples numeric vector of length 1. number of samples to select
#' @return data table of selected iterations
thin_chain <- function(chain,n_samples){
    if(n_samples == 0 || n_samples >= nrow(chain)){
        chain
    } else {
        chain[round(seq(1,nrow(chain),length.out = n_samples)),]
    }
}