#' estimate parameters for the different models in the study
#' 
#' @param strain character: sH1N1, pH1N1, H5N1, H7N9 for the WT strains;
#' ENG195, Turkey05 for pH1N1-PR8 and H5N1-PR8 respectively
#' @param sensitivity character: choice of "main" (model in main paper),
#' "n_L_n_I" (changing the number of latent and infectious stages), 
#' "loss" (ignoring loss of virus due to entry into target cells),
#' "mc" (use multi-cycle data only),
#' "linear_p" (production rate of virions increases linearly with time)
#' @param short logical.  if TRUE, run short MCMC chain for testing purposes;
#' otherwise, run until convergence
#' @return NULL (results saved to file)
#' @export
run_exp <- function(strain, sensitivity, short) {
  
  data_source <- determine_data_source_for_strain(strain)
  strain_short_name <- shorten_strain_name(strain)
  if(sensitivity == "mc") {
    # exclude single-cycle and mock-yield data
    inc_data_logical <- c(F,F,T,T,F)
  } else {
    inc_data_logical <- c(T,T,T,T,T)
  }
  if(data_source == "PR8") {
    # exclude total viral load data
    inc_data_logical[c(2,4)] <- FALSE
  }
  
  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)
  
  # set directories to save outputs
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strain)
  filenames = paste0(strain_dir, vapply(seeds, function(x) x[1], double(1)))
  
  # arguments to pass to setup functions
  if(sensitivity == "linear_p") {
    p_over_time <- "linear"
  } else {
    p_over_time <- "constant"
  }
  
  loss <- sensitivity != "loss" # TRUE unless ignoring loss of virus due to entry into target cells
  n_L <- n_I <- ifelse(sensitivity == "n_L_n_I", 60, 10)
  
  # define function to specify priors for model parameters
  specify_parameters <- function()
    specify_parameters_fn(data_source = data_source,
                          inc_data_logical = inc_data_logical,
                          loss = loss,
                          n_L = n_L,
                          n_I = n_I,
                          p_over_time = p_over_time)
  # define function to read in data and construct function to calculate likelihood
  generate_data_and_posterior <- function(x)
    get_data_and_generate_posterior_fn(x, strain,
                                       inc_data_logical)
  # define function to calculate summary parameters using model parameters
  gen_summary_statistics <- function(x, y) gen_summary_statistics_fn(x, y, p_over_time)
  # define function to find CIs for model predictions
  get_prediction_ci_df <- get_prediction_ci_df_all
  
  # define parameters for MCMC
  if(short) {
    mcmcPars <- gen_mcmcPars_short()
  } else {
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  }
  
  # define function to generate random starting values
  generate_start_values <- function(x, y, z) generate_random_start(x, y)
  
  # run MCMC
  run_all(starting_point_seed = seeds, 
          filenames = filenames, 
          run_parallel = TRUE, 
          mvr = TRUE,
          specify_parameters = specify_parameters,
          generate_data_and_posterior = generate_data_and_posterior,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = gen_summary_statistics, 
          calc_residuals = NULL, 
          get_prediction_ci_df = get_prediction_ci_df,
          run = c("run" = TRUE, "process" = TRUE))
}

#' wrapper for running and postprocessing MCMC chain
#' 
#' @param starting_pars NULL for randomly generated start values, or
#' n_replicates by n_temperatures list of vectors of length n_pars,
#' where n_pars is the number of parameters, n_replicates is the number of parallel
#' chains to run to assess convergence, and n_temperatures is the number of 
#' temperatures used for parallel tempering.  If no parallel tempering, use
#' an n_replicates  list of vectors of length n_pars
#' @param starting_point_seed list of n_replicates vectors, each of length 
#' n_temperatures, used for random number generator for generating starting
#' point for each chain.
#' @param filenames vector of n_replicates filenames in which to store results
#'(without extension)
#' @param run_flags logical vector with elements
#' run: run MCMC if TRUE
#' process: read in a process MCMC chain if TRUE
#' @param mvr logical.  if TRUE, use multivariate proposal, else use univariate proposal
#' @param mvr_init list of parameters for multivariate proposal.  if NULL, use defaults
#' @param specify_parameters closure to specify parameters by generating the parTab
#' data frame.
#' @param generate_data_and_posterior closure to load data and create the function
#' used to evaluate the likelihood.
#' @param mcmcPars list of tuning parameters for MCMC
#' @param mcmc_seed list of length n_replicates, each element being an integer.
#' Used to seed RNG for MCMC.  if NULL, use default
#' @param generate_start_values closure to generate starting parameter values
#' if starting_pars = NULL
#' @inheritParams postprocess_for_plotting
#' @export
run_all <- function(starting_pars = NULL,
                    starting_point_seed, 
                    filenames, 
                    run_flags = c(run = TRUE, process = TRUE),
                    run_parallel = TRUE, 
                    mvr = TRUE,
                    mvr_init = NULL,
                    specify_parameters,
                    generate_data_and_posterior,
                    mcmcPars = gen_mcmcPars(),
                    mcmc_seed = NULL,
                    generate_start_values = function(x, y, z) generate_random_start(x, y),
                    gen_summary_statistics = NULL, 
                    calc_residuals = NULL, 
                    get_prediction_ci_df = NULL) {
  
  if(run_flags[["run"]]) {
    setup_and_run_MCMC(starting_pars,
                       starting_point_seed,
                       filenames,
                       run_parallel,
                       mvr,
                       mvr_init,
                       specify_parameters,
                       generate_data_and_posterior,
                       mcmcPars,
                       mcmc_seed,
                       generate_start_values)
  }
  
  if(run_flags[["process"]]) {
    filename <- paste0(filenames[1], ".RData")
    postprocessed_list <- postprocess_for_plotting(filename = filename, 
                                                   run_parallel = run_parallel, 
                                                   gen_summary_statistics = gen_summary_statistics, 
                                                   calc_residuals = calc_residuals, 
                                                   get_prediction_ci_df = get_prediction_ci_df)
    
  }
}