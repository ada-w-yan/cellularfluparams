# complete doc

#'creates a closure to calculate the likelihood and model predictions for given
#' parameter values
#'
#' @inheritParams postprocess_chain
#' @inheritParams sample_from_prior
#' @inheritParams calc_mean_gen_time
#' @useDynLib cellularfluparams
#' @return closure to calculate likelihood and model predictions for given
#' parameter values
CREATE_POSTERIOR_FUNC <-
  function(parTab, data, PRIOR_FUNC, p_over_time = "constant") {
    predict_data <- names(data)
    
    ## function to extract the number of latent and infectious cells at t = 0
    extract_iv <- function(sol, col_name) {
      unname(sol[2, grep(col_name, colnames(sol))])
    }
    
    # function to calculate the viral load at the start of the incubation period
    calc_V_inf_0 <- function(pars_incubation) {
      pars_incubation[["MOI"]] * pars_incubation[["T_0"]] / pars_incubation[["supernatant_vol"]]
    }
    
    has_data <- function(x)
      has_data_outer(x, predict_data)
    in_data <- vapply(c("sc", "mc", "mock"), has_data, logical(1))
    
    # compile ODEs

    if (in_data[["sc"]] || in_data[["mc"]]) {
      if (parTab[parTab$names == "n_I", "values"] == 1) {
        # gen <- compile_model("model_teiv_stages_e_dual.R")
        gen <- model_teiv_stages_e_dual
      } else if(p_over_time == "constant") {
        # gen <- model_teiv_stages_dual
        gen <- compile_model("model_teiv_stages_dual.R")
      } else {
        gen <- model_teiv_stages_dual_linear_p
      }
    }
    
    if (in_data[["mock"]]) {
      gen_mock <- compile_model("model_mock_yield.R")
      # gen_mock <- model_mock_yield
    }
    
    # make function to transform a vector of parameter values into the
    # form needed for ODE solving
    transform_pars <- transform_pars_wrapper(parTab, predict_data)
    
    # get default parameters
    parameter_list <- transform_pars(parTab$values)
    
    gen_mod <- function(parameter_list) {
      parameter_list[["L_0"]] <- double(parameter_list[["n_L"]])
      parameter_list[["I_0"]] <- double(parameter_list[["n_I"]])
      gen(user = parameter_list)
    }
    if (in_data[["sc"]]) {
      mod <- gen_mod(parameter_list$sc)
    }
    if (in_data[["mc"]]) {
      mod <- gen_mod(parameter_list$mc)
    }
    
    if (in_data[["mock"]]) {
      mod_mock <- gen_mock(user = parameter_list$mock)
    }
    
    # find times at which to solve ODEs
    
      solving_time_list <- list(
        sc = c(data$sc_inf$t, data$sc_tot$t),
        mc = c(data$mc_inf$t, data$mc_tot$t),
        mock = data$mock$t
      )
      solving_time_list <-
        lapply(solving_time_list, make_solving_time)

    experiment <-
      c(
        sc_inf = "sc",
        sc_tot = "sc",
        mc_inf = "mc",
        mc_tot = "mc",
        mock = "mock"
      )
    experiment <- experiment[predict_data]
    
    # find indices of solving times at which samples are taken
      find_sampling_ind <- function(data, solving_time) {
        sampling_times <- data[!is.na(data$V_noise), "t"]
        match(sampling_times, solving_time)
      }
      
      # find indices fo solving times at which to make predictions
      find_prediction_ind <- function(data_name, solving_time) {
        current_exp <- experiment[[data_name]]
        # predict at all solving times
        seq_along(solving_time[[current_exp]])
      }
      
      sampling_ind <- lapply(predict_data,
                             function(x)
                               find_sampling_ind(data[[x]],
                                                 solving_time_list[[experiment[x]]]))
      names(sampling_ind) <- predict_data
      prediction_ind <- lapply(predict_data,
                               function(x)
                                 find_prediction_ind(x, solving_time_list))
      names(prediction_ind) <- predict_data
      prediction_sampling_indices <- list(prediction = prediction_ind,
                                          sampling = sampling_ind)

    
    # function to name a vector of predictions
    name_predictions <- function(prediction) {
      names(prediction) <- as.character(seq_along(prediction))
      prediction
    }
    
    compartment_names <- c(
      sc_inf = "V_inf",
      sc_tot = "V_tot",
      mc_inf = "V_inf",
      mc_tot = "V_tot",
      mock = "V"
    )
    compartment_names <- compartment_names[predict_data]
    
    extract_values_from_sol <- function(experiment, sol, ind) {
      lapply(predict_data[grepl(experiment, predict_data)],
             function(x)
               sol[ind[[x]], compartment_names[[x]]])
    }
    
    calc_sol_sc_mc <-
      function(transformed_pars,
               experiment,
               solving_time_list,
               prediction_sampling_indices,
               mod,
               tol) {
        pars_incubation <- transformed_pars[[experiment]]
        ## calculate starting viral load for incubation (t = -1)
        pars_incubation$V_inf_0 <- calc_V_inf_0(pars_incubation)
        pars_incubation$L_0 <- double(pars_incubation$n_L)
        pars_incubation$I_0 <- double(pars_incubation$n_I)
        do.call(mod$set_user, pars_incubation)
        sol_incubation <-
          solve_ODE_error_handling(mod, c(-pars_incubation$incubation_period, 0), tol)
        
        
        ## extract number of latent and infectious cells at t = 0
        pars <- transformed_pars[[experiment]]
        pars$L_0 <- extract_iv(sol_incubation, "L")
        pars$I_0 <- extract_iv(sol_incubation, "I")
        pars$T_0 <- extract_iv(sol_incubation, "T")
        
        do.call(mod$set_user, pars)
        solving_time_subset <- solving_time_list[[experiment]]
        solving_time_subset <-
          solving_time_subset[solving_time_subset >= 0]
        sol <- solve_ODE_error_handling(mod, solving_time_subset, tol)
        sol <- rbind(sol_incubation[1, ], sol)
        # find viral load at prediction and sampling time points
        
        lapply(prediction_sampling_indices,
               function(x)
                 extract_values_from_sol(experiment, sol, x))
      }
    
    # function to calculate log likelihood
    f <- function(pars) {
      transformed_pars <- transform_pars(pars)
      tol <- 1e-6
      
        if (in_data[["mc"]]) {
          ## solve for mc experiment
          V_mc <- calc_sol_sc_mc(
            transformed_pars,
            "mc",
            solving_time_list,
            prediction_sampling_indices,
            mod,
            tol
          )
        } else {
          V_mc <- list(NULL, NULL)
        }
        
        if (in_data[["mock"]]) {
          ## solve for mock experiment
          do.call(mod_mock$set_user, transformed_pars$mock)
          sol <- solve_ODE_error_handling(mod_mock,
                                          solving_time_list$mock,
                                          tol)
          
          # find viral load at prediction and sampling time points
          V_mock <- lapply(prediction_sampling_indices,
                           function(x)
                             extract_values_from_sol("mock", sol, x))
        } else {
          V_mock <- list(NULL, NULL)
        }
        
        if (in_data[["sc"]]) {
          V_sc <- calc_sol_sc_mc(
            transformed_pars,
            "sc",
            solving_time_list,
            prediction_sampling_indices,
            mod,
            tol
          )
        } else {
          V_sc <- list(NULL, NULL)
        }
        
        # concatenate predicted and sampled viral loads
        # from different experiments
        V <- lapply(names(prediction_sampling_indices),
                    function(x)
                      c(V_sc[[x]], V_mc[[x]], V_mock[[x]]))
        names(V) <- names(prediction_sampling_indices)
        
        # extract observation error parameters
        sigma_pars <- transformed_pars$sigma
        sigma_vec <-
          sigma_pars[c("sigma_inf",
                       "sigma_tot",
                       "sigma_inf",
                       "sigma_tot",
                       "sigma_mock")]
        sigma_vec <-
          sigma_vec[c(in_data[["sc"]], in_data[["sc"]], in_data[["mc"]], in_data[["mc"]], in_data[["mock"]])]
        
        # calculate log likelihood of data given model
        individual_lik <-
          Map_vapply(
            function(x, y, z)
              calc_LL_log10normal_threshold(x[!is.na(x$V_noise), "V_noise"], y,
                                            z, threshold = sigma_pars[["obs_threshold"]]),
            double(1),
            data,
            V[["sampling"]],
            sigma_vec
          )
        lik <- sum(individual_lik)
        
        # give predictions names
        V[["prediction"]] <-
          lapply(V[["prediction"]], name_predictions)
        names(V[["prediction"]]) <- predict_data
        misc <- unlist(V[["prediction"]])
        list(lik = lik, misc = misc)
    }
    f
  }


#' create closure to calculate summary statistics and provide bounds for them
#' 
#' @inheritParams postprocess_chain
#' @inheritParams sample_from_prior
#' @inheritParams calc_mean_gen_time
#' @return a list with elements
#' "calc_summary": closure to calculate
gen_summary_statistics_fn <-
  function(parTab, startTab, p_over_time) {
    sum_stat_names <- c("log10_R_0",
                        "r",
                        "gen_time")
    # if called without arguments, output names of what the outputs would have been
    if (missing(parTab) && missing(startTab)) {
      return(sum_stat_names)
    }
    
    # function to transform parameters to calculate summary statistics
    transform_pars <- transform_pars_wrapper(parTab,
                                             predict_data = "mc_inf")
    
    # functions to calculate summary statistics given parameter values
    
    calc_log10_R_0 <- function(pars) {
      log10(calc_R_0(pars, p_over_time))
    }
    
    calc_gen_time <- function(pars) {
      calc_mean_gen_time(
        c1 = pars[["c_inf"]] + pars[["beta_inf"]] * pars[["T_0"]],
        tau_E = 1 / pars[["k1"]],
        tau_I = 1 / pars[["delta"]],
        n_L = pars[["n_L"]],
        n_I = pars[["n_I"]],
        p_over_time = p_over_time
      )
    }
    calc_r <- function(pars) {
      calc_r_teiv_stages(pars, p_over_time)
    }
    
    # closure to generate function to calculate summary statistics given parameter values
    calc_summary <- function(values) {
      values <- transform_pars(values)
      
      calc_funcs <- paste0("calc_", sum_stat_names)
      calc_funcs <- lapply(calc_funcs, match.fun)
      

        if (any(grepl("mc_V_inf_0", parTab$names))) {
          values <- values[["mc"]]
        } else {
          values <- values[["sc"]]
        }
        sum_stats <-
          vapply(calc_funcs, function(f)
            f(values), double(1))
        names(sum_stats) <- sum_stat_names
        sum_stats
    }
    
    # nice parameter names for plotting
    
    par_names_plot <- c("log10_R_0" = "$\\log_{10}R_{0}$",
                        "r" = "$r$",
                        "gen_time" = "$T_G$")
    
    lower_bound <- c("log10_R_0" = 2,
                     "r" = 0,
                     "gen_time" = 0)
    upper_bound <- c("log10_R_0" = 6,
                     "r" = 1,
                     "gen_time" = 100)
    
    list(
      calc_summary = calc_summary,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      par_names_plot = par_names_plot
    )
  }

#' read data for experiment, and define functions to calculate likelihood
#'
#' @inheritParams postprocess_for_plotting
#' @inheritParams determine_data_source_for_strain
#' @param inc_data_logical logical vector of length 5.  indicates which of
#' \code{c("sc_inf", "sc_tot", "mc_inf", "mc_tot", "mock")} we have data for
#' @inheritParams calc_mean_gen_time
#' @return list with four elements:
#' data: list of data frames containing the data for each experiment
#' CREATE_POSTERIOR_FUNC: wrapper to create function to calculate
#' probability of data given parameters
#' CREATE_PRIOR_FUNC: wrapper to create function to calculate
#' probability of data given prior
get_data_and_generate_posterior_fn <- function(parTab,
                                                  strain = "sH1N1-WT",
                                                  inc_data_logical,
                                                  p_over_time = "constant") {
  # make character vector of experiments to predict
  # if we have any single-cycle data, includes all SC experiments;
  # if we have any multi-cycle data, includes all MC experiments;
  # otherwise equal to teh data sets we have
  predict_data <- make_inc_data(inc_data_logical, TRUE)$inc_data
  
  data_list <-
    read_and_process_all_data(strain = strain, predict_data)
  
  list(
    data = data_list,
    CREATE_POSTERIOR_FUNC = function(x, y, z)
      CREATE_POSTERIOR_FUNC(x, y, z, p_over_time),
    CREATE_PRIOR_FUNC = function(x)
      NULL
  )
}

#' checks if any string in the vector inc_data contains
#' any of the strings in the vector str
#' 
#' @param str character vector
#' @param inc_data character vector
#' @return logical vector of length 1
has_data_outer <- function(str, inc_data) {
  any(vapply(str, function(x)
    any(grepl(x, inc_data)), logical(1)))
}

#' make character vector of data to infer/predict given inc_data_logical vector
#' 
#' @inheritParams get_data_and_generate_posterior_fn
#' @param predict logical.  If TRUE, if say data for "sc_inf" is included,
#' "sc_tot" will be included in the output also.  If FALSE, directly
#' translate inc_data_logical into character vector
#' @return character vector
make_inc_data <-
  function(inc_data_logical = rep(T, 5),
           predict = FALSE) {
    data_names <- c("sc_inf", "sc_tot", "mc_inf", "mc_tot", "mock")
    stopifnot(length(inc_data_logical) == length(data_names))
    
    if (!is.null(names(inc_data_logical))) {
      inc_data_logical <- inc_data_logical[data_names]
    }
    
    if (predict) {
      update_inc_data_logical <- function(inc_data_logical, idx) {
        if (any(inc_data_logical[idx])) {
          inc_data_logical[idx] <- T
        }
        inc_data_logical
      }
      
      sc <- c(1, 2)
      mc <- c(3, 4)
      inc_data_logical <-
        update_inc_data_logical(inc_data_logical, sc)
      inc_data_logical <-
        update_inc_data_logical(inc_data_logical, mc)
    }
    
    inc_data <- data_names[inc_data_logical]
    stopifnot(!any(is.na(inc_data)))
    list(data_names = data_names, inc_data = inc_data)
  }

# make a vector of prediction times from a vector of sampling times
make_prediction_time <-
  function(sampling_time,
           min_prediction_period,
           t_end) {
    sort(unique(c(
      sampling_time, seq(0, t_end, by = min_prediction_period)
    )))
  }

#' make a vector of solving times from a vector of times
#'
#' @param times_in numeric vector of times
#' @return numeric vector: sorted unique times, with t = 0
make_solving_time <- function(times_in) {
  sort(unique(c(0, times_in)))
}


#' read PR8 data from prism file
#'
#' @return list of data frame with all data in file
#' @export
read_PR8_data <- function() {
  # extract all data from file
  data(PR8_data)
  data_lists <- PR8_data
  # for each experiment
  
  process_list <- function(data_list) {
    n_replicates <- vapply(data_list[-1], ncol, double(1))
    timepoints <- as.vector(data_list[[1]])
    viral_load <- unlist(lapply(data_list[-1], as.vector))
    strains <- names(data_list[-1])
    strain_indicator <-
      rep(strains, times = n_replicates * length(timepoints))
    replicate_indicator <-
      rep(unlist(lapply(n_replicates, seq_len)),
          each = length(timepoints))
    # observation threshold
    threshold <- 10
    # make columns for time, strain, replicate number, viral load and whether viral load is below threshold
    out_df <-
      data.frame(
        "t" = rep(timepoints, length.out = length(strain_indicator)),
        "strain" = strain_indicator,
        "replicate" = replicate_indicator,
        "V" = viral_load,
        "below_threshold" = !is.na(viral_load) &
          viral_load < threshold
      )
    # for plotting purposes, set below threshold data poitns to some arbitrary value
    out_df[out_df$below_threshold, "V"] <- threshold / 10
    rownames(out_df) <- NULL
    out_df
  }
  data_lists <- lapply(data_lists, process_list)
}

#' read and format data for single-cycle, multi-cycle and mock yield experiments
#'
#' @inheritParams determine_data_source_for_strain
#' @param inc_data character vector of up to length 5.  a subset of
#' \code{c("sc_inf", "sc_tot", "mc_inf", "mc_tot", "mock")} indicating which data
#' sets we're reading in
#' @return list of data frames.  Each data frame is a data set, and has columns
#' "t" and "V_noise".  Each value of V_noise is a measurement corresponding to the time
#' in the same row.  If the value of V_noise is NA, that row is a placeholder
#' indicating that later on, we want to predict the viral load at that time.
read_and_process_all_data <-
  function(strain = "sH1N1-WT", inc_data) {
    # for scraped or provided data, we make functions to read the data, then apply
    # them to different files depending on whether the data is scraped or provided
    
    # make filenames containing data
    
      tmax <- c("mc" = 150,
                "sc" = 24,
                "mock" = 72)
    data_source <- determine_data_source_for_strain(strain)
    if (data_source == "WT") {
      filename_short <-
        c("SC_TCID", "SC_RNA", "MC_TCID", "MC_RNA", "MY_TCID")
      experiment <- c("sc", "sc", "mc", "mc", "mock")
      ext <- ".dat"
      strain_short_name <- shorten_strain_name(strain)
      read_and_process_data <- function(filename, experiment) {
        data("WT_data")
        my_data <- WT_data[[filename]]
        colnames(my_data) <- c("t", "V_noise")
        # include zero time point for mock yield data only
        if (experiment == "mock") {
          my_data <- my_data[my_data$t >= 0, ]
        }  else {
          my_data <- my_data[my_data$t > 0, ]
          my_data <-
            rbind(data.frame(t = c(-1, 0), V_noise = c(NA, NA)), my_data)
        }
        my_data <-
          rbind(my_data, data.frame(t = tmax[[experiment]], V_noise = NA))
        
        my_data
      }
      data_dir <- paste0("../../../data/", data_source, "/")
      
      filename <- paste0(strain_short_name, "_", filename_short)
      data_names <- make_inc_data()$data_names
      names(filename) <- data_names
      
      data_list <- Map(read_and_process_data, filename, experiment)
      
    } else if (data_source == "PR8") {
      data_list_unformatted <- read_PR8_data()
      
      extract_strain <- function(data_df, strain, experiment) {
        # predict every 6 hPR8 for multi-cycle data; every hour for single-cycle data;
        # and every 24 hPR8 for mock yield data
        min_prediction_period <- c("sc" = 1,
                                   "mc" = 6,
                                   "mock" = 24)
        # if single-cycle data, need to choose between 1st and 2nd data set
        second_in_table <- grepl("2nd", data_df$strain)
        # if we're trying to get the second data set
        if (!grepl("old", strain)) {
          # and this is the single-cycle data (as indicated by there being a second run)
          if (any(second_in_table)) {
            # then use the second run
            data_df <- data_df[second_in_table, ]
          }
        } else {
          # else use the first run
          data_df <- data_df[!second_in_table, ]
          strain <- gsub("_old", "", strain)
        }
        data_df <- data_df[grep(strain, data_df$strain), ]
        # censor values below threshold
        if (length(data_df$below_threshold) > 0) {
          data_df[data_df$below_threshold, "V"] <- -1
        }
        data_df <- data_df[c("t", "V")]
        colnames(data_df) <- c("t", "V_noise")
        # deal with back titration data points
        supernatant_vol <- 3
        if (experiment == "mock") {
          inoculum_vol <- 3
        } else {
          inoculum_vol <- .5
        }
        
        data_df[data_df$t < 0, "V_noise"] <- NA
        
        data_df <- data_df[!is.na(data_df$t), ]
        if (experiment == "mock") {
          data_df[data_df$t < 0, "t"] <- 0
        }
        
        t_vec <- data_df$t
        # pad data frame with prediction time points if necessary
        t_vec <-
          unique(c(seq(0, tmax[[experiment]], by = min_prediction_period[[experiment]]), tmax[[experiment]]))
        t_vec <- t_vec[!(t_vec %in% data_df$t)]
        if (length(t_vec) > 0) {
          data_df <- rbind(data_df, data.frame(t = t_vec, V_noise = NA))
        }
        data_df <- data_df[order(data_df$t), ]
        data_df
      }
      
      experiment <- c("mc", "sc", "mock")
      
      old_strain_names <- c("pH1N1-PR8" = "ENG195", "H5N1-PR8" = "Turkey05")
      strain <- old_strain_names[[strain]]
      
      # read data for the given strain
      data_list <- Map(function(x, y, z)
        extract_strain(x, strain, y),
        data_list_unformatted,
        experiment)
      
      names(data_list) <- c("mc_inf", "sc_inf", "mock")
      # re-order data for consistency, and construct empty total virus data frames for prediction
      data_list <-
        data_list[c("sc_inf", "sc_inf", "mc_inf", "mc_inf", "mock")]
      names(data_list) <- make_inc_data()$data_names
      data_list$sc_tot$V_noise <- NA
      data_list$mc_tot$V_noise <- NA
    } else {
      stop("unknown data source")
    }
    
    # output only the required data sets
    data_list <- data_list[inc_data]
    data_list
    
  }

#' make parameter table 
#'
#' @param data_source character string: "WT", "PR8"
#' @param n_L numeric vector of length 1. number of latent stages
#' @param n_I numer vector of length 1.  number of infectious stages
#' @inheritParams get_data_and_generate_posterior_fn
#' @param loss logical vector of length 1.  if TRUE, fit model with loss of virus
#' due to entry into target cells. if FALSE, fit model without loss of virus
#' due to entry into target cells.
#' @inheritParams calc_mean_gen_time
#' @return data frame with columns
#' values: default parameter value
#' names: internal name for parameter
#' fixed:  if 0, fit parameter, otherwise fix parameter
#' lower_bound: lower bound of parameter for fitting
#' upper_bound: upper bound of parameter for fitting
#' steps: width of proposal distribution in parameter space
#' names_plot: character vector which which to label parameter distirbution plot
specify_parameters_fn <-
  function(data_source,
           n_L = 10,
           n_I = 10,
           inc_data_logical = rep(T, 5),
           loss = TRUE,
           p_over_time = "constant") {
    # make vectors for the names of the full set of data sets; the data sets we have; and
    # the experiments we want to predict the viral load for, respectively
    temp_list <-
      lapply(c(FALSE, TRUE), function(x)
        make_inc_data(inc_data_logical, x))
    data_names <- temp_list[[1]]$data_names
    inc_data <- temp_list[[1]]$inc_data
    predict_data <- temp_list[[2]]$inc_data
    
    default_step <- 0.1
    
    # inefficient but easier to read
    # check lower and upper bounds for everything
    
    # initial number of target cells
    T_0 <- switch(
      data_source,
      "PR8" = 2.5e6,
      "WT" = 1e6
    )
    
    # supernatant volume differs between experiments
    supernatant_vol <- switch(
      data_source,
      "PR8" = 3,
      "WT" = 10
    )
    
    # make a function that greps predict_data
    need_to_predict_data <-
      function(x)
        has_data_outer(x, predict_data)
    # make a function that greps inc_data
    has_data <- function(x)
      has_data_outer(x, inc_data)
    
    # the above setup is because sometimes we want to fit to total + mock yield data, but not the infectious
    # SC and MC data
    
    if (p_over_time != "constant") {
      p_tot_upper_bound <- p_inf_upper_bound <- 10
    } else {
      p_tot_upper_bound <- 6
      p_inf_upper_bound <- 3
    }
    
    parTab <- data.frame(
      values = -1.5,
      names = "log10_c_inf",
      fixed = 0,
      lower_bound = -3,
      upper_bound = 0,
      steps = default_step,
      names_plot = "$\\log_{10}c_{inf}$",
      stringsAsFactors = FALSE
    )
    
    parTab <-
      rbind(
        parTab,
        data.frame(
          values = supernatant_vol,
          names = "supernatant_vol",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "supernatant_vol"
        )
      )
    # if there is any single-cycle or multi-cycle data, we can fit beta
    if (need_to_predict_data(c("sc", "mc"))) {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = -5,
            names = "log10_beta_inf",
            fixed = 0,
            lower_bound = -10,
            upper_bound = -1,
            steps = default_step,
            names_plot = "$\\log_{10}\\beta_{inf}$"
          )
        )
      
      # whether loss = 1 or 0 controls inclusion or exclusion of teh loss term in the ODEs
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = as.numeric(loss),
            names = "loss",
            fixed = 1,
            lower_bound = -Inf,
            upper_bound = Inf,
            steps = default_step,
            names_plot = "loss"
          )
        )
      
      if (loss) {
        parTab <-
          rbind(
            parTab,
            data.frame(
              values = -1,
              names = "log10_pfu_per_inf",
              fixed = 0,
              lower_bound = -7,
              upper_bound = 0,
              steps = default_step,
              names_plot = "$\\log_{10} pfu per inf$"
            )
          )
      } else {
        parTab <-
          rbind(
            parTab,
            data.frame(
              values = 0,
              names = "log10_pfu_per_inf",
              fixed = 1,
              lower_bound = -Inf,
              upper_bound = Inf,
              steps = default_step,
              names_plot = "$\\log_{10} pfu per inf$"
            )
          )
      }
      
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0,
          names = "log10_p_tot",
          fixed = as.numeric(!has_data("tot")),
          lower_bound = 0,
          upper_bound = p_tot_upper_bound,
          steps = default_step,
          names_plot = "$\\log_{10}p_{tot}$"
        )
      )
      
      
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = 4,
            names = "tau_E",
            fixed = 0,
            lower_bound = 0.1,
            upper_bound = 12,
            steps = default_step,
            names_plot = "$\\tau_E$"
          )
        )
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = 1.5,
            names = "log10_tau_I",
            fixed = 0,
            lower_bound = 0,
            upper_bound = 3,
            steps = default_step,
            names_plot = "$\\log_{10}\\tau_I$"
          )
        )
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = 1e-3,
            names = "c_tot",
            fixed = 1,
            lower_bound = -Inf,
            upper_bound = Inf,
            steps = default_step,
            names_plot = "$c_{tot}$"
          )
        )
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = n_L,
            names = "n_L",
            fixed = 1,
            lower_bound = -Inf,
            upper_bound = Inf,
            steps = default_step,
            names_plot = "$n_L$"
          )
        )
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = n_I,
            names = "n_I",
            fixed = 1,
            lower_bound = -Inf,
            upper_bound = Inf,
            steps = default_step,
            names_plot = "$n_I$"
          )
        )
    }
    
    # if we have any single-cycle data, we need to fit p_inf_sc
    # (regardless of whether we have the single-cycle infectious viral load data)
    if (need_to_predict_data("sc")) {
      parTab <- rbind(
        parTab,
        data.frame(
          values = -1,
          names = "log10_p_inf_sc_on_p_inf_mc",
          fixed = 0,
          lower_bound = -5,
          upper_bound = 0,
          steps = default_step,
          names_plot = "$\\log_{10}(p_{inf,SC}/p_{inf,MC})$"
        )
      )
    }
    
    # if we have any multi-cycle data, we need to fit p_inf_mc
    # (regardless of whether we have the multi-cycle infectious viral load data)
    # actually now that we have express p_tot and p_inf_sc in terms of p_inf_mc,
    # need some value of p_inf_mc even if no mc data.
    # this value needs to be fitted even if we have no mc data, because it is the
    # lower bound of p_tot and the upper bound of p_inf_sc -- artificially restrict
    # range if we don't fit it
    if (need_to_predict_data(c("sc", "mc"))) {
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0,
          names = "log10_p_inf_mc",
          fixed = 0,
          lower_bound = -2,
          upper_bound = p_inf_upper_bound,
          steps = default_step,
          names_plot = "$\\log_{10}p_{inf,MC}$"
        )
      )
    }
    
    ## initial conditions
    
    if (need_to_predict_data("sc")) {
      # different initial multiplicites of infection for our data and for Simon et al. data
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = log10(5),
            names = "log10_sc_MOI",
            fixed = 0,
            lower_bound = log10(supernatant_vol /
                                  T_0),
            upper_bound = log10(10),
            steps = default_step,
            names_plot = "$\\log_{10} SC MOI$"
          )
        )
      
      # if we are using the Simon et al. data, the measurement frequency for
      # the single-cycle data is very high, so we can fix the initial SC
      # viral load to the geomentric mean for the initial measurements.
      # otherwise infer the initial SC viral load
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0,
          names = "log10_sc_V_inf_0",
          fixed = 0,
          lower_bound = -2,
          upper_bound = 7,
          steps = default_step,
          names_plot = "$\\log_{10}V_{inf,0,SC}$"
        )
      )
      
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0,
          names = "log10_sc_V_tot_0_on_sc_V_inf_0",
          fixed = as.numeric(!has_data("sc_tot")),
          lower_bound = 0,
          upper_bound = 6,
          steps = default_step,
          names_plot = "$\\log_{10}(V_{tot,0,SC}/V_{inf,0,SC})$"
        )
      )
    }
    
    # if we are using the Simon et al. data, the measurement frequency for
    # the multi-cycle data is very high, so we can fix the initial MC
    # viral load to the geomentric mean for the initial measurements.
    # otherwise infer the initial MC viral load
    
    if (need_to_predict_data("mc")) {
      parTab <- rbind(
        parTab,
        data.frame(
          values = -10,
          names = "log10_mc_MOI",
          fixed = as.numeric(data_source == "WT"),
          lower_bound = log10(supernatant_vol /
                                T_0),
          upper_bound = -1,
          steps = default_step,
          names_plot = "log10 MC MOI"
        )
      )
      
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0,
          names = "log10_mc_V_inf_0",
          fixed = 0,
          lower_bound = -2,
          upper_bound = 4,
          steps = default_step,
          names_plot = "$\\log_{10}V_{inf,0,MC}$"
        )
      )
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0,
          names = "log10_mc_V_tot_0_on_mc_V_inf_0",
          fixed = as.numeric(!has_data("mc_tot")),
          lower_bound = 0,
          upper_bound = 6,
          steps = default_step,
          names_plot = "$\\log_{10}(V_{tot,0,MC}/V_{inf,0,MC})$"
        )
      )
    }
    
    # need initial target cell number for SC and MC experients
    if (need_to_predict_data(c("sc", "mc"))) {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = 1,
            names = "sc_incubation_period",
            fixed = 1,
            lower_bound = -Inf,
            upper_bound = Inf,
            steps = default_step,
            names_plot = "SC incubation period"
          )
        )
      
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = T_0,
            names = "T_0",
            fixed = 1,
            lower_bound = -Inf,
            upper_bound = Inf,
            steps = default_step,
            names_plot = "$T_0$"
          )
        )
    }
    
    # different priors for initial mock yield viral load concentrion between our experiments
    # and Simon et al.
    if (need_to_predict_data("mock")) {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = 7.5,
            names = "mock_log10_V_0",
            fixed = 0,
            lower_bound = 5,
            upper_bound = 10,
            steps = default_step,
            names_plot = "$\\log_{10}V_{0,MY}$"
          )
        )
    }
    
    ## observation parameters
    parTab <-
      rbind(
        parTab,
        data.frame(
          values = 0.4,
          names = "sigma_mock",
          fixed = 0,
          lower_bound = 0,
          upper_bound = 2,
          steps = default_step,
          names_plot = "$\\sigma_{mock}$"
        )
      )
    
    parTab <- rbind(
      parTab,
      data.frame(
        values = 1,
        names = "sigma_inf_on_sigma_mock",
        fixed = 1,
        lower_bound = 0,
        upper_bound = 2,
        steps = default_step,
        names_plot = "$\\sigma_{mock}$"
      )
    )
    
    parTab <-
      rbind(
        parTab,
        data.frame(
          values = 0.6,
          names = "sigma_tot",
          fixed = 0,
          lower_bound = 0,
          upper_bound = 2,
          steps = default_step,
          names_plot = "$\\sigma_{tot}$"
        )
      )
    
    if (data_source == "WT") {
      obs_threshold <- 0
    } else {
      obs_threshold <- 10
    }
    # osebrvation threshold of 10 pfu.  treat values below threshold as censored
    parTab <-
      rbind(
        parTab,
        data.frame(
          values = obs_threshold,
          names = "obs_threshold",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "obs threshold"
        )
      )
    
    # if we don't have any total viral load data
    if (!has_data("tot")) {
      # if sigma_tot is set to infinity, the likelihood term associated with the
      # total viral load is set to zero
      parTab[parTab$names == "sigma_tot",  c("fixed", "values")] <-
        c(1, Inf)
      
    }
    if (!has_data("inf")) {
      # if sigma_inf is set to infinity, the likelihood term associated with the
      # infectious viral load for the single-cycle or multi-cycle experiments is set to zero
      parTab[parTab$names == "sigma_inf_on_sigma_mock",  "values"] <-
        Inf
    }
    
    # convert logicals to numerics
    parTab$fixed <- as.numeric(parTab$fixed)
    # set bounds for fixed parameters to c(-Inf Inf)
    parTab[parTab$fixed == 1, "lower_bound"] <- -Inf
    parTab[parTab$fixed == 1, "upper_bound"] <- Inf
    
    parTab
  }

#' create a closure which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
#'
#' @inheritParams postprocess_for_plotting
#' @param predict_data character vector up to length 5.  subset of
#' \code{c("sc_inf", "sc_tot", "mc_inf", "mc_tot", "mock")}. indicates the
#' experiments about which we want to make predictions
#' @param G numeric vector: number of generations to track
#' @return a function which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
transform_pars_wrapper <- function(parTab, predict_data, G) {
  if(missing(parTab)) {
    par_names <- NULL
    need_to_predict_data <- function(x) "mc" %in% x
  } else {
    par_names <- parTab$names
    log10_ind <- grep("log10", par_names)
    need_to_predict_data <-
      function(x)
        has_data_outer(x, predict_data)
  }
  if(missing(G)) {
    G <- NULL
  }
  
  transform_pars <- function(pars) {
    mc_only <- FALSE
    if(is.null(par_names)) {
      par_names <- names(pars)
      log10_ind <- grep("log10", par_names)
      mc_only <- TRUE
    }
    pars[log10_ind] <- sapply(pars[log10_ind], function(x)
      10 ^ (x))
    names(pars) <- par_names
    
      if (need_to_predict_data(c("sc", "mc"))) {
        if (pars[["loss"]] == 0) {
          beta1 <- pars[["log10_beta_inf"]]
          inf_per_virion <-
            beta_tot <- pars[["log10_beta_inf"]] <- 0
        } else {
          inf_per_virion <- pars[["log10_pfu_per_inf"]]
          beta1 <-
            pars[["log10_beta_inf"]] * pars[["supernatant_vol"]] / inf_per_virion
          beta_tot <- pars[["log10_beta_inf"]] / inf_per_virion
        }
      }
      
      # list parameters needed to predict outcome of single-cycle experiment
      if (need_to_predict_data("sc")) {
        sc_pars <- list(
          beta = beta1,
          beta_inf = pars[["log10_beta_inf"]],
          beta_tot = beta_tot,
          inf_per_virion = inf_per_virion,
          p_inf = pars[["log10_p_inf_sc_on_p_inf_mc"]] * pars[["log10_p_inf_mc"]],
          p_tot = pars[["log10_p_tot"]],
          k1 = 1 / pars[["tau_E"]],
          delta = 1 / pars[["log10_tau_I"]],
          c_inf = pars[["log10_c_inf"]],
          c_tot = pars[["c_tot"]],
          T_0 = pars[["T_0"]],
          MOI = pars[["log10_sc_MOI"]],
          supernatant_vol = pars[["supernatant_vol"]],
          incubation_period = pars[["sc_incubation_period"]],
          V_inf_0 = pars[["log10_sc_V_inf_0"]],
          V_tot_0 = pars[["log10_sc_V_tot_0_on_sc_V_inf_0"]] * pars[["log10_sc_V_inf_0"]],
          n_L = pars[["n_L"]],
          n_I = pars[["n_I"]]
        )
      } else {
        sc_pars <- NULL
      }
      
      # list parameters needed to predict outcome of multi-cycle experiment
      if (need_to_predict_data("mc")) {
        mc_pars <- list(
          beta = beta1,
          beta_inf = pars[["log10_beta_inf"]],
          beta_tot = beta_tot,
          inf_per_virion = inf_per_virion,
          p_inf = pars[["log10_p_inf_mc"]],
          p_tot = pars[["log10_p_tot"]],
          k1 = 1 / pars[["tau_E"]],
          delta = 1 / pars[["log10_tau_I"]],
          c_inf = pars[["log10_c_inf"]],
          c_tot = pars[["c_tot"]],
          T_0 = pars[["T_0"]],
          MOI = pars[["log10_mc_MOI"]],
          supernatant_vol = pars[["supernatant_vol"]],
          incubation_period = pars[["sc_incubation_period"]],
          V_inf_0 = pars[["log10_mc_V_inf_0"]],
          V_tot_0 = pars[["log10_mc_V_tot_0_on_mc_V_inf_0"]] * pars[["log10_mc_V_inf_0"]],
          n_L = pars[["n_L"]],
          n_I = pars[["n_I"]],
          L_0 = double(pars[["n_L"]]),
          I_0 = double(pars[["n_I"]])
        )
        if(!is.null(G)) {
          mc_pars$L_0 <- matrix(0, nrow = mc_pars$n_L, ncol = G)
          mc_pars$I_0 <- matrix(0, nrow = mc_pars$n_I, ncol = G)
        }
      } else {
        mc_pars <- NULL
      }
      
      # list parameters needed to predict outcome of mock yield experiment
      if (need_to_predict_data("mock")) {
        mock_pars <- list(c = pars[["log10_c_inf"]],
                          V_0 = pars[["mock_log10_V_0"]])
      } else {
        mock_pars <- NULL
      }
      
      # list measurement error parameters
      sigma_pars <-
        list(
          sigma_inf = pars[["sigma_inf_on_sigma_mock"]] *
            pars[["sigma_mock"]],
          sigma_tot = pars[["sigma_tot"]],
          sigma_mock = pars[["sigma_mock"]],
          obs_threshold = pars[["obs_threshold"]]
        )
      
      if(mc_only) {
        mc_pars
      } else {
        list(
          "sc" = sc_pars,
          "mc" = mc_pars,
          "mock" = mock_pars,
          "sigma" = sigma_pars
        )
      }
  }
  transform_pars
}

#' extract model predictions form MCMC chain and calculate CIs
#' 
#' @param chain MCMC chain
#' @param data_list data used to fit model
#' @return list of same lenght as data_list.  Each element contains:
#' data_df: the corresponding element of data_list
#' V_ci: the predicted viral load corresponding to that data
get_prediction_ci_df_all <- function(chain, data_list) {
  get_prediction_ci_df_single <- function(data_name) {
    # extract experiment from data list
    data_df <- data_list[[data_name]]
    col_names <- colnames(data_df)
    noise_var_str <- col_names[grep("noise", col_names)]
    # retrieve prediction times from data frame
    prediction_time <- make_solving_time(data_df$t)
    # make a t = 0 data point if there isn't one, just for axis limit purposes
    
    if (!any(data_df$t == 0)) {
      zero_df <- data.frame(t = 0, noise_var_str = NA)
      # if the following line is absent, the colname will be "noise_var_str"
      # rather than the value of noise_var_str
      colnames(zero_df) <- c("t", noise_var_str)
      data_df <- rbind(zero_df, data_df)
    }
    
    var_str <- paste0(data_name, ".")
    
    # extract the credible intervals for the experiment being plotted
    chain_names <- names(chain)
    
    prediction_names <-
      chain_names[grepl(var_str, chain_names, fixed = TRUE)]
    
    ci <-
      gen_prediction_ci(data_df, chain, prediction_names, unique(data_df$t))
    list(data_df = data_df, V_ci = ci)
  }
  
  prediction_ci_df <-
    lapply(names(data_list), get_prediction_ci_df_single)
  names(prediction_ci_df) <- names(data_list)
  prediction_ci_df
}