#' setup and run MCMC
#' 
#' @inheritParams run_all
setup_and_run_MCMC <-function(starting_pars = NULL,
                              starting_point_seed, 
                              filenames, 
                              run_parallel = TRUE, 
                              mvr = TRUE,
                              mvr_init = NULL,
                              specify_parameters,
                              generate_data_and_posterior,
                              mcmcPars = gen_mcmcPars(),
                              mcmc_seed = NULL,
                              generate_start_values = function(x, y, z) generate_random_start(x, y)){
  
    lapply(filenames, mkdir_from_filename)
  ###############################################################################
  ## setting up random number generation
  ###############################################################################

  ## provide some seeds to generate starting points, if not given
  n_temperatures <- length(mcmcPars[["temperature"]])
  n_replicates <- length(filenames)
  parallel_tempering_flag <- (n_temperatures > 1)
  
  if(missing(starting_point_seed)){
    starting_point_seed <- make_seed(n_replicates, n_temperatures)
  }
  
  ## if seeds manually given, check that the number of seeds is the same as the
  ## number of filenames
  if(parallel_tempering_flag){
    
    tryCatch(invisible(vapply(starting_point_seed, identity, double(n_temperatures))),
             error = function(c) "length of temperatures not equal to length of seeds"
    )
    
  }

  if(is.null(starting_pars) && n_replicates != length(starting_point_seed)){
    stop("length of filenames not equal to length of seeds")
  }  
  
  ###############################################################################
  ## specifying model parameters
  ###############################################################################
  
  parTab <- specify_parameters()
  
  ################################################################################
  ## generating data and defining likelihood function
  ################################################################################
  
  temp_list <- generate_data_and_posterior(parTab)
  list2here(temp_list, var_names = c("data", "CREATE_POSTERIOR_FUNC", "CREATE_PRIOR_FUNC"))
  rm(temp_list)

  PRIOR_FUNC <- CREATE_PRIOR_FUNC(parTab)
  f <- CREATE_POSTERIOR_FUNC(parTab, data, PRIOR_FUNC)

  # check that likelihood function works
  
  check_likelihood <- function(f, values) {
    f_output <- f(parTab$values)
    lik <- f_output$lik
    
    stopifnot(is.numeric(lik) && length(lik) == 1 && !is.na(lik) && lik != -Inf)
  }

  check_likelihood(f, parTab$values)

  ################################################################################
  ## generating random starting values
  ################################################################################
  
  if(is.null(starting_pars)) {
    if(parallel_tempering_flag){
      startTab <- lapply(starting_point_seed, 
                         function(y) lapply(y, function(x) generate_start_values(parTab, x, data)))
    } else {
      startTab <- lapply(starting_point_seed, function(x) generate_start_values(parTab, x, data))
    }
  } else {
    
    assign_starting_pars <- function(parTab, starting_pars) {
      startTab <- parTab
      startTab$values <- starting_pars
      startTab
    }
    
    if(parallel_tempering_flag){
      startTab <- lapply(starting_pars, 
                         function(y) lapply(y, function(x) assign_starting_pars(parTab, x)))
    } else {
      startTab <- lapply(starting_pars, function(x) assign_starting_pars(parTab, x))
    }
  }

  if(mvr){
    
    if(is.null(mvr_init)) {
      define_initial_mvrPars <- function(parTab) {
        n_row_covMat <- sum(parTab$fixed == 0)
        covMat <- diag(nrow(parTab))
        list(covMat,2.38/sqrt(n_row_covMat),w=0.8)
      }
      
      mvrPars <- define_initial_mvrPars(parTab)
      
      if(parallel_tempering_flag){
        mvrPars <- rep(list(mvrPars), n_temperatures)
      }
    } else {
      mvrPars <- mvr_init
    }
    
    mvrPars <- rep(list(mvrPars), n_replicates)

  } else {
    mvrPars <- NULL
  }

  save_all(filenames[1])

  ################################################################################
  ## running MCMC
  ################################################################################
  message("setup complete")
  posterior_list <- run_MCMC_loop(startTab, data, mcmcPars, filenames,
                                  CREATE_POSTERIOR_FUNC, mvrPars, PRIOR_FUNC,
                                  seed = mcmc_seed,
                                  run_parallel = run_parallel)
  
  list2here(posterior_list, var_names = c("diagnostics", "output"))
  rm(posterior_list)
  
  save_all(filenames[1])
  invisible(NULL)
}