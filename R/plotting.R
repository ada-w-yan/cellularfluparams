# complete doc

#' make viral load plots in paper: Fig 1...
#' 
#' @param data_source character: WT or PR8
#' @param experiment character: "sc", "mc" or "mock" for single-cycle, 
#' multi-cycle and mock-yield experiments respectively
#' @inheritParams run_exp
#' @return ggplot object
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
plot_viral_load <- function(data_source, experiment, sensitivity) {
  
  strains <- determine_strains_for_exp(data_source) %>%
    paste0(., "-", data_source)
  dirs <- set_dirs(sensitivity)
  
  get_data_and_ci_dfs <- function(strain) {
    strain_dir <- get_strain_dir(dirs, strain)
    df_list <- readRDS(paste0(strain_dir, "1ci.rds"))
    data_df <- df_list[[experiment]]$data_df
    ci_df <- cbind(data_df["t"],as.data.frame(df_list[[experiment]]$V_ci))
    colnames(ci_df) <- c("t","lower","upper")
    data_df$strain <- strain
    ci_df$strain <- strain
    list(data_df = data_df, ci_df = ci_df)
  }
  
  dfs <- lapply(strains, get_data_and_ci_dfs)
  
  bind_dfs <- function(dfs, df_type) {
    all_df <- lapply(dfs, function(x) x[[df_type]])
    all_df <- do.call(rbind, all_df)
    all_df$strain <- factor(all_df$strain, levels = strains)
    all_df
  }
  
  df_names <- names(dfs[[1]])[1:2]
  dfs <- lapply(df_names, function(x) bind_dfs(dfs, x))
  names(dfs) <- df_names
  
  y_lims <- switch(experiment,
                   sc_inf = c(1e2,1e8),
                   mc_inf = c(1, 1e9),
                   mock = c(1e4,1e9),
                   sc_tot = c(1e4,1e11),
                   mc_tot = c(1,1e12))
  axis_break_interval <- switch(experiment,
                                sc_inf = 2,
                                mc_inf = 3,
                                mock = 2,
                                sc_tot = 3,
                                mc_tot = 2)
  
  axis_breaks <- seq(log10(y_lims[1]), log10(y_lims[2]), by = axis_break_interval)
  axis_break_labels <- paste0("$10^{", axis_breaks, "}$")
  axis_breaks <- 10^axis_breaks
  
  point_shape <- 24
  point_size <- 3
  x_label <- "Time (hpi)"
  error_string <- c("lower", "upper")
  noise_var <- "V_noise"
  if(data_source == "PR8") {
    viral_load_str <- "Viral load (PFU/mL)"
  } else if(grepl("tot", experiment)) {
    viral_load_str <- "Viral load (copy number/mL)"
  } else {
    viral_load_str <- "Viral load (TCID$_{50}$/mL)"
  }
  strain_colours <- define_palette(data_source)
  
  measurement_times <- unique(dfs$data_df[!is.na(dfs$data_df[,"V_noise"]),"t"])
  
  time_df <- expand.grid(t = measurement_times, y = y_lims)
  time_df$group <- as.character(time_df$t)
  
  if(grepl("sc", experiment)) {
    interval <- 6
  } else {
    interval <- 24
  }
  
  x_max <- max(dfs$data_df[,"t"])
  
  x_breaks <- seq(0, x_max, by = interval)
  jitter <- position_jitter(width = 2, height = 0)
  
  g <- ggplot(dfs$data_df) + 
    geom_point(aes(x = t, y = V_noise,
                   fill = strain, group = strain),
               colour = "transparent",
               shape = point_shape, 
               size = point_size) +
    geom_ribbon(data = dfs$ci_df[!is.na(dfs$ci_df[,error_string[1]]),],
                mapping = aes_string(x = "t", ymin = error_string[1],
                                     ymax = error_string[2], 
                                     fill = "strain", 
                                     group = "strain"),
                alpha = .3) +
    scale_y_log10(latex2exp::TeX(viral_load_str), expand = c(0,0), breaks = axis_breaks, labels = latex2exp::TeX(axis_break_labels)) +
    theme_bw() +
    xlab(x_label) + 
    ggtitle("") +
    scale_x_continuous(expand = c(0,0), breaks = x_breaks) +
    coord_cartesian(xlim = c(0, x_max), ylim = y_lims)
  g <- g + scale_fill_manual(values = strain_colours)
  if(experiment == "mock") {
    g <- g + theme(plot.title = element_text(hjust = 0.5),
                   legend.justification=c(1,0),
                   legend.position=c(1,0),
                   text = element_text(size = 24),
                   plot.margin = unit(c(0,1,0,0), "cm"))
  } else {
    g <- g + theme(plot.title = element_text(hjust = 0.5),
                   legend.position="none",
                   text = element_text(size = 24),
                   plot.margin = unit(c(0,1,0,0), "cm"))
  }
  
  ggsave(paste0(dirs[["figs"]], data_source, "_", experiment, ".pdf"), g,
         width = 15, height = 15, units = "cm")
  g
  
}

# Fig. 2A-C, Fig. 3C
#' plot summary statistics or parameters across strains, for a particular set of
#' experiments
#' 
#' @param sum_stat character vector of length 1.  if "cell_infect_params", plot generation
#' time, R_0 and initial growth rate.
#' if "mean_gen_no", plot the mean generation number.
#' @inheritParams run_exp
#' @return ggplot object
#' @import ggplot2
#' @export
plot_stats_across_strains <- function(sum_stat, sensitivity) {
  
  data_sources <- c("WT", "PR8")
  strains <-determine_strains_for_exp("all", TRUE)
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strains)
  if(sum_stat == "mean_gen_no") {
    filename_str <- paste0(sum_stat, "_prctile.tex")
  } else {
    filename_str <- "1_prctile_combined.tex"
  }
  filenames <- paste0(strain_dir, filename_str)
  
  param_names <- switch(sum_stat,
                        "cell_infect_params" = c("$T_G$", "$\\log_{10}R_{0}$", "$r$"),
                        "mean_gen_no" = "mean_gen_no")
  
  file_param_grid <- expand.grid("filenames" = filenames, "param_names" = param_names)
  file_param_grid <- cbind(file_param_grid, data.frame("strain" = factor(rep(strains, length(param_names)), levels = strains)))
  
  sum_stats <- Map(function(x,y) read_summary_stats(as.character(x), as.character(y)), file_param_grid$filenames, file_param_grid$param_names)
  file_param_grid <- cbind(file_param_grid,do.call(rbind, sum_stats))
  file_param_grid <- file_param_grid[,-1]
  colnames(file_param_grid)[3:5] <- c("lower", "median", "upper")
  
  if(sum_stat == "cell_infect_params") {
    bound_vec <- rep(c(100, 6, 1), each = length(strains))
    bound_vec[seq(1, length(file_param_grid$param_names), length(strains))] <- 0
  } else {
    bound_vec <- c(0, 20)
  }
  
  # sensible parameter names for plotting
  levels(file_param_grid$param_names) <- switch(sum_stat,
                                                "cell_infect_params" = c("Mean generation time (h)",
                                                                         "log10(Reproduction number)",
                                                                         "Initial growth rate (/h)"),
                                                "mean_gen_no" = "Mean generation no.")
  
  levels(file_param_grid$strain) <- unlist(lapply(data_sources, determine_strains_for_exp, w_data_source = TRUE))
  
  dummy <- data.frame(strain = file_param_grid$strain, median = bound_vec,
                      param_names = file_param_grid$param_names)  
  # determine colour palette for plot
  
  strain_colours <- unlist(lapply(data_sources, define_palette))
  names(strain_colours) <- levels(file_param_grid$strain)
  
  # make plots
  g <- ggplot(file_param_grid, aes(x = strain, y = median)) + geom_point(aes(color = strain)) +
    geom_errorbar(aes(x = strain, ymin = lower, ymax = upper, color = strain)) +
    facet_wrap(~param_names, nrow = 1, 
               scales = "free", labeller = label_wrap_gen(width=30), drop = TRUE) +
    scale_y_continuous("", expand = c(0,0)) + 
    scale_x_discrete("") + 
    scale_colour_manual(values = strain_colours) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none",
          strip.text.x = element_text(size = 8)) + 
    geom_blank(data = dummy)
  
  width_each <- 5
  width <- width_each * length(levels(file_param_grid$param_names))
  ggsave(paste0(dirs[["figs"]], sum_stat, ".pdf"), g, width = width, height = 7.5, units = "cm")
  g
}

# Fig. 2D-E
#' plot the bivariate density plot for R_0 and gen_time for multiple strains on the same plot
#' 
#' @inheritParams plot_viral_load
#' @inheritParams run_exp
#' @return a ggplot object. the bivariate density plot for R_0 and gen_time for multiple strains on the same plot
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
plot_bivariate_density_multiple_strains <- function(data_source, sensitivity) {
  dirs <- set_dirs(sensitivity)
  strains <- determine_strains_for_exp(data_source, TRUE)
  strain_dir <- get_strain_dir(dirs, strains)
  filenames <- paste0(strain_dir, "1cell_infect_params.rds")
  strain_colours <- define_palette(data_source)
  # get the samples and probability density function for eadch strain
  read_and_process_data <- function(strain, filename) {
    # read in samples
    original_sum_stats_names <- rownames(readRDS(filename))
    my_summary_stats_names <- c("log10_R_0", "gen_time", "r")
    data_df <- t(subset_summary_stats_from_rds(filename, my_summary_stats_names, n_samples = 0))
    data_df <- as.data.frame(data_df)
    data_df$R_0 <- log10(data_df$R_0)
    names(data_df) <- my_summary_stats_names
    data_df$strain <- strain
    data_df
  }
  
  data_df <- Map(read_and_process_data, strains, filenames)
  
  # put the samples for all the strains in one data frame
  data_df <- do.call(rbind, data_df)
  data_df$strain <- factor(data_df$strain, levels = strains)
  
  # find the median of each summary stat for each strain
  summary_df <- data_df %>% 
    dplyr::group_by(strain) %>% 
    dplyr::summarise(log10_R_0 = median(log10_R_0),
                     gen_time  = median(gen_time),
                     r = median(r))
  
  # make density plot
  make_plot <- function(sum_stats) {
    axis_labels <- c("log10_R_0" = "Basic reproduction number",
                     "gen_time" = "Mean generation time (h)",
                     "r" = "Initial growth rate (h$^{-1})")
    axis_lims <- list("log10_R_0" = c(0, 6),
                      "gen_time" = c(0, 100),
                      "r" = c(0, 1))
    axis_break_interval <- c("log10_R_0" = 1,
                             "gen_time" = 25,
                             "r" = .25)
    
    make_axis_breaks <- function(axis_lims, axis_break_interval) {
      seq(axis_lims[1], axis_lims[2], by = axis_break_interval)
    }
    
    axis_breaks <- Map(make_axis_breaks, axis_lims, axis_break_interval)
    axis_break_labels <- axis_breaks
    axis_break_labels[["log10_R_0"]] <- paste0("10^", 
                                               axis_break_labels[["log10_R_0"]])
    g <- ggplot(data_df, aes_string(x = sum_stats[1],
                                    y = sum_stats[2],
                                    group = "strain")) +
      stat_density2d(aes(alpha=..level.., fill = strain),geom="polygon",lwd=0,bins=5)+
      geom_point(data = summary_df, aes(colour = strain)) +
      scale_x_continuous(latex2exp::TeX(axis_labels[[sum_stats[1]]]), 
                         limits = axis_lims[[sum_stats[1]]],
                         breaks = axis_breaks[[sum_stats[1]]],
                         labels = latex2exp::TeX(axis_break_labels[[sum_stats[1]]])) + 
      scale_y_continuous(latex2exp::TeX(axis_labels[[sum_stats[2]]]), 
                         limits = axis_lims[[sum_stats[2]]],
                         breaks = axis_breaks[[sum_stats[2]]],
                         labels = latex2exp::TeX(axis_break_labels[[sum_stats[2]]])) + 
      guides(alpha=FALSE, colour = FALSE) +
      theme_bw() +
      theme(text = element_text(size = 20),
            legend.position = c(0.8, 0.2))
    strain_names <- determine_strains_for_exp(data_source, TRUE)
    g <- g +
      scale_fill_manual(values = strain_colours,
                        labels = strain_names) +
      scale_color_manual(values = darken(strain_colours),
                         labels = strain_names)
    ggsave(paste0(dirs[["figs"]], sum_stats[1] ,"_", sum_stats[2], "_", data_source,".pdf"),
           g, width = 15, height = 15, units = "cm")
    g
  }
  
  sum_stats_list <- matrix(c("gen_time", "log10_R_0", 
                             "log10_R_0", "r",
                             "gen_time", "r"), 
                           ncol = 2, 
                           byrow = TRUE)
  apply(sum_stats_list, 1, make_plot)
  
}

#' make Fig 3A-B
#' 
#' @inheritParams run_exp
#' @return a list of ggplot objects
#' @import ggplot2
#' @export
make_Fig3 <- function(sensitivity) {
  
  dirs <- set_dirs(sensitivity)
  
  data_source <- c("PR8", "WT")
  strains <- lapply(data_source, determine_strains_for_exp, w_data_source = TRUE)
  plot_df <- lapply(strains, function(x) lapply(x, calc_generation_strain,
                                                sensitivity = sensitivity))
  plot_peak_dist_all_strains <- function(strains, plot_df) {
    data_source <- determine_data_source_for_strain(strains[1])
    strain_colours <- define_palette(data_source)
    gen_df <- do.call(rbind, plot_df)
    gen_df$strain <- factor(gen_df$strain, levels = strains)
    gen_df$gen <- as.numeric(gen_df$gen)
    
    if(data_source == "PR8") {
      lty_values <- c(1,2)
    } else {
      lty_values <- c(3,1,2,4)
    }
    y_label <- "Viral load (proportion)"
    g <- ggplot(gen_df, aes(x = gen, y = V, group = strain)) +
      geom_line(aes(lty = strain, color = strain)) +
      scale_y_continuous(latex2exp::TeX(y_label), expand = c(0,0)) +
      theme_bw() +
      theme(legend.justification=c(1,1),
            legend.position=c(.9,.9),
            text = element_text(size = 24)) +
      
      geom_point(aes(color = strain), size = 3) +
      scale_x_continuous("Generation number", 
                         limits = c(1, 16),
                         breaks = seq(2, 20, by = 2)) +
      coord_cartesian(ylim = c(0, 1))
    
    
    strain_names <- determine_strains_for_exp(data_source, TRUE)
    g <- g + scale_color_manual(values = strain_colours,
                                labels = strain_names) +
      scale_linetype_manual(values = lty_values,
                            labels = strain_names)
    g
  }
  
  peak_gen_all <- Map(plot_peak_dist_all_strains, strains, plot_df)
  Map(function(x, y) ggsave_wch(x, y, width = 15, height = 15, units = "cm"),
      paste0(dirs[["figs"]], data_source, "_gen_dist_peak.pdf"), peak_gen_all)
  peak_gen_all
}

#' Make Fig. 4
#' 
#' @param constant_param character: "gen_time" (A), "log10_R_0" (B), or "r" (C).
#' The parameter to hold constant.
#' @return a ggplot object
#' @import ggplot2
#' @export
plot_gen_dist_vary_params <- function(constant_param) {
  list2here(calc_gen_dist_vary_params(constant_param))
  
  par_lookup <- switch(constant_param,
                       "gen_time" = c("log10_R_0" = "independent",
                                      "r" = "dependent"),
                       "log10_R_0" = c("gen_time" = "independent",
                                       "r" = "dependent"),
                       "r" = c("gen_time" = "independent",
                               "log10_R_0" = "dependent"))
  par_lookup <- c(par_lookup, c("constant"))
  names(par_lookup)[3] <- constant_param
  
  legend_labels <- paste0("[", all_pars[[par_lookup[["gen_time"]]]], ", ", 
                          all_pars[[par_lookup[["log10_R_0"]]]], ", ",
                          all_pars[[par_lookup[["r"]]]],"]")
  
  legend_title <- "$T_G$ (h), $\\log_{10}R_0$, $r$ (h$^{-1}$)"
  r_range <- range(as.numeric(all_pars[[par_lookup[["r"]]]]))
  if(constant_param == "r") {
    plot_title <- paste0("r = ", signif(r_range[1], sf))
  } else {
    plot_title <- paste0(signif(r_range[1], sf), "$\\leq r \\leq$", signif(r_range[2], sf))
  }
  
  g <- ggplot(plot_df, aes(x = gen, 
                           y = V, 
                           color = infection_par_value, 
                           group = infection_par_value)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    coord_cartesian(xlim = range(plot_df$gen), ylim = c(0, 1), expand = FALSE) +
    xlab("Generation no.") +
    ylab("Viral load (proportion)") +
    ggtitle(latex2exp::TeX(plot_title)) +
    scale_color_manual(values = gg_color_hue(length(legend_labels)), 
                       labels = legend_labels,
                       name = latex2exp::TeX(legend_title)) +
    theme(text = element_text(size = 16),
          legend.position = c(0.5, 0.8),
          legend.text = element_text(size = 16))
  dirs <- set_dirs("main")
  ggsave(paste0(dirs[["figs"]], constant_param, 
                "_peak.pdf"),
         g, width = 10, height = 10, units = "cm")
  g
}

#' make Fig. S3
#' 
#' @inheritParams run_exp
#' @return a list of ggplot objects, each with one bivariate plot
#' @import ggplot2
#' @export
plot_bivariate_density_single_strain <- function(strain) {
  sensitivity <- "main"
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strain)
  filename <- paste0(strain_dir, "/1cell_infect_params.rds")
  # get the samples and probability density function for eadch strain
  
  # read in samples
  cell_infect_param_names <- c("log10_R_0", "gen_time", "r")
  data_df <- t(subset_summary_stats_from_rds(filename, cell_infect_param_names, n_samples = 0))
  data_df <- as.data.frame(data_df)
  data_df$R_0 <- log10(data_df$R_0)
  names(data_df) <- cell_infect_param_names
  
  make_plot <- function(sum_stats) {
    axis_labels <- c("gen_time" = "Generation time (h)",
                     "log10_R_0" = "log10(Reproduction number)",
                     "r" = "Growth rate (/h)")
    g <- ggplot(data_df, aes_string(x = sum_stats[1],
                                    y = sum_stats[2])) +
      stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
      guides(alpha=FALSE, colour = FALSE) +
      theme_bw() +
      xlab(axis_labels[sum_stats[1]]) +
      ylab(axis_labels[sum_stats[2]]) +
      theme(text = element_text(size = 28),
            legend.position = c(0.8, 0.2))
    ggsave(paste0(dirs[["figs"]], sum_stats[1] ,"_", sum_stats[2], "_", strain,".pdf"),
           g, width = 15, height = 15, units = "cm")
    g
  }
  
  
  sum_stats_list <- matrix(c("gen_time", "log10_R_0", 
                             "log10_R_0", "r",
                             "gen_time", "r"), 
                           ncol = 2, 
                           byrow = TRUE)
  apply(sum_stats_list, 1, make_plot)
  
  
}

#' Fig S3: plot histogram of cellular infection parameter
#' @param cell_infect_param character: "gen_time", "log10_R_0" or "r".
#' @inheritParams run_exp
#' @return a ggplot object
#' @import ggplot2
#' @export
plot_cell_infect_param_histogram_by_strain <- function(strain, sum_stat) {
  sensitivity <- "main"
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strain)
  filename <- paste0(strain_dir, "1cell_infect_params.rds")
  chain <- readRDS(filename)
  chain <- chain[sum_stat,]
  x_lab <- switch(sum_stat,
                  "gen_time" = "Generation time (h)",
                  "log10_R_0" = "log10(Reproduction number)",
                  "r" = "Growth rate (/h)")
  
  g <- plot_histogram_general(vec = chain, plot_xlab = x_lab, limits = range(chain))
  ggsave(paste0(dirs[["figs"]], strain, "_", sum_stat, ".pdf"), g,
         width = 15, height = 15, units = "cm")
  g
}

#' plot Fig S4
#' 
#' @return a ggplot object
#' @import ggplot2
#' @export
plot_change_all_pars_sH1N1_H7N9_CI <- function() {
  sensitivity <- "main"
  strains <- c("sH1N1-WT", "H7N9-WT")
  dirs <- set_dirs(sensitivity)
  results_dir <- get_strain_dir(dirs, strains[1])
  filename_out <- paste0(results_dir, "sensitivity_sum_stats.rds")
  
  sum_stats <- readRDS(filename_out)
  
  
  sum_stats[,"sum_stat"] <- as.factor(sum_stats[,"sum_stat"])
  levels(sum_stats[,"sum_stat"]) <- c("Generation time (h)", "log10(Reproduction number)", "Growth rate (/h)")
  par_labels <- c("beta_inf" = parse(text = latex2exp::TeX('$\\beta_{inf}$')),
                  "c_inf" = parse(text = latex2exp::TeX('$c_{inf}$')),
                  "p_inf" = parse(text = latex2exp::TeX('$p_{inf,MC}$')),
                  "tau_I" = parse(text = latex2exp::TeX('$\\tau_I$')),
                  "tau_L" = parse(text = latex2exp::TeX('$\\tau_L$')))
  g <- ggplot(sum_stats, aes(x = par_name, color = par_name)) +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    facet_wrap(~sum_stat) +
    scale_y_continuous("Change from baseline (%)", limits = c(-100,305)) +
    scale_x_discrete("Parameter", labels = par_labels) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size=20),
          axis.text.x = element_text(angle=90, hjust=1),
          strip.text.x = element_text(size = 12))
  ggsave(paste0(dirs[["figs"]], "change_single_par.pdf"), g, width = 20, height = 10, units = "cm")
  g
}

#' plot Fig. S5
#' 
#' @inheritParams calc_gen_dist_vary_params
#' @return a ggplot object
#' @import ggplot2
#' @export
plot_gen_dist_fixed_gen_time <- function(constant_param) {
  viral_load_fn <- calc_generations_fixed_growth_until_I_max
  I_max <- Inf
  r <- 1
  gen_time <- 12
  n_strains <- 4
  V_0 <- 1
  R_0 <- 1e4
  V_max <- 1e14
  
  if(constant_param == "r") {
    constant_par_value <- r
    gen_time <- seq_len(n_strains) * 6
    R_0 <- exp(r * gen_time)
    r <- rep(r, n_strains)
    infection_pars <- log10(R_0)
    dependent_pars <- gen_time
  } else if (constant_param == "gen_time") {
    constant_par_value <- gen_time
    R_0 <- 10^(seq_len(n_strains) + 2)
    r <- log(R_0) / gen_time
    gen_time <- rep(gen_time, n_strains)
    infection_pars <- log10(R_0)
    dependent_pars <- r
  } else {   
    constant_par_value <- log10(R_0)
    gen_time <- seq_len(n_strains) * 6
    r <- log(R_0) / gen_time
    R_0 <- rep(R_0, n_strains)
    infection_pars <- gen_time
    dependent_pars <- r
  }
  
  n_gens <- ceiling(log(V_max) / log(R_0))
  time <- n_gens * gen_time
  
  calc_generations_wrapper <- function(gen_time, R_0, time) {
    gen_dist <- calc_generations_at_time_fixed(time, gen_time, V_0, R_0, I_max, viral_load_fn)
    length(gen_dist)
  }
  
  n_gens <- unlist(Map(calc_generations_wrapper, gen_time, R_0, time))
  sf <- 3
  plot_df <- data.frame(gen = n_gens, 
                        group = signif(infection_pars, sf))
  
  
  all_pars <- list(independent = signif(infection_pars, sf),
                   dependent = signif(dependent_pars, sf),
                   constant = rep(as.character(signif(constant_par_value, sf)), 4))
  par_lookup <- switch(constant_param,
                       "gen_time" = c("log10_R_0" = "independent",
                                      "r" = "dependent"),
                       "log10_R_0" = c("gen_time" = "independent",
                                       "r" = "dependent"),
                       "r" = c("log10_R_0" = "independent",
                               "gen_time" = "dependent"))
  par_lookup <- c(par_lookup, c("constant"))
  names(par_lookup)[3] <- constant_param
  
  legend_labels <- paste0("[", all_pars[[par_lookup[["gen_time"]]]], ", ", 
                          all_pars[[par_lookup[["log10_R_0"]]]], ", ",
                          all_pars[[par_lookup[["r"]]]],"]")
  
  legend_title <- "$T_G$ (h), $\\log_{10}R_0$, $r$ (h$^{-1}$)"
  plot_df$group <- as.factor(plot_df$group)
  independent_par_name <- names(par_lookup)[par_lookup == "independent"]
  
  x_lab <- switch(independent_par_name,
                  "gen_time" = "Generation time (h)",
                  "log10_R_0" = "log10(Reproduction number)")
  
  R_0_range <- range(as.numeric(all_pars[[par_lookup[["log10_R_0"]]]]))
  if(constant_param == "log10_R_0") {
    plot_title <- paste0("$\\log_{10}R_0$ = ", signif(R_0_range[1], sf))
  } else {
    plot_title <- paste0(signif(R_0_range[1], sf), "$\\leq \\log_{10}R_0 \\leq$", signif(R_0_range[2], sf))
  }
  
  g <- ggplot(plot_df, aes(x = group, 
                           y = gen,
                           fill = group)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    coord_cartesian(ylim = c(0, 10), expand = FALSE) +
    xlab(x_lab) +
    ylab("Generation no.") +
    scale_x_discrete(breaks = as.numeric(infection_pars), labels = infection_pars) +
    scale_y_continuous(breaks = seq(2,10,2))+
    ggtitle(latex2exp::TeX(plot_title)) +
    scale_fill_manual(values = gg_color_hue(length(legend_labels)),
                      labels = legend_labels,
                      name = latex2exp::TeX(legend_title)) +
    theme(text = element_text(size = 16),
          legend.position = c(0.5, 0.8),
          legend.text = element_text(size = 16))
  
  dirs <- set_dirs("main")
  ggsave(paste0(dirs[["figs"]], constant_param, "_fixed_peak.pdf"),
         g, width = 10, height = 10, units = "cm")
  g
}

#' make correlation table for supplementary csv
#' 
#' @inheritParams run_exp
#' @return a matrix
#' @export
make_corr_table <- function(strain) {
  sensitivity <- "main"
  dirs <- set_dirs(sensitivity)
  strain_dir <- get_strain_dir(dirs, strain)
  filename <- paste0(strain_dir, "1.RData")
  
  chain <- get_MCMC_chain(filename)
  par_names_old <- c("log10_c_inf",
                     "log10_beta_inf",
                     "tau_E",
                     "log10_tau_I",
                     "log10_p_inf_sc_on_p_inf_mc",
                     "log10_p_inf_mc",
                     "log10_p_tot",
                     "log10_pfu_per_inf",
                     "log10_sc_MOI",
                     "log10_sc_V_inf_0",
                     "log10_sc_V_tot_0_on_sc_V_inf_0",
                     "log10_mc_MOI",
                     "log10_mc_V_inf_0",
                     "log10_mc_V_tot_0_on_mc_V_inf_0",
                     "mock_log10_V_0",
                     "sigma_mock",
                     "sigma_tot")
  par_names_new <- c("log_10(c)",
                     "log_10(beta_inf)",
                     "tau_L",
                     "log_10(tau_I)",
                     "log_10(p_{SC,pfu/mL}/p_{MC,pfu/mL})",
                     "log_10(p_{MC,pfu/mL})",
                     "log_10(p_{RNA/mL})",
                     "log_10(n)",
                     "log_10(V_{SC,pfu/mL,-1})",
                     "log_10(V_{SC,pfu/mL,0})",
                     "log_10(V_{SC,RNA/mL,0}/V_{SC,pfu/mL,0})",
                     "log_10(V_{MC,pfu/mL,-1})",
                     "log_10(V_{MC,pfu/mL,0})",
                     "log_10(V_{MC,RNA/mL,0}/V_{MC,pfu/mL,0})",
                     "log_10(V_{mock,pfu/mL,0})",
                     "sigma_pfu",
                     "sigma_RNA")
  experiment <- determine_data_source_for_strain(strain)
  mc_incub_param_old <- "log10_mc_MOI"
  mc_incub_param_new <- "log_10(V_{MC,pfu/mL,-1})"
  RNA_params_old <- c("log10_p_tot", 
                      "log10_sc_V_tot_0_on_sc_V_inf_0", 
                      "log10_mc_V_tot_0_on_mc_V_inf_0",
                      "sigma_tot")
  RNA_params_new <- c("log_10(p_{RNA/mL})", 
                      "log_10(V_{SC,RNA/mL,0}/V_{SC,pfu/mL,0})", 
                      "log_10(V_{MC,RNA/mL,0}/V_{MC,pfu/mL,0})",
                      "sigma_RNA")
  
  if(experiment == "PR8") {
    par_names_old <- par_names_old[!(par_names_old %in% RNA_params_old)]
    par_names_new <- par_names_new[!(par_names_new %in% RNA_params_new)]
    MOI_params <- c("log_10(V_{SC,pfu/mL,-1})", "log_10(V_{MC,pfu/mL,-1})")
  } else {
    par_names_old <- par_names_old[!(par_names_old %in% mc_incub_param_old)]
    par_names_new <- par_names_new[!(par_names_new %in% mc_incub_param_new)]
    MOI_params <- "log_10(V_{SC,pfu/mL,-1})"
  }
  
  chain <- as.data.frame(chain)
  T_0 <- chain[1, "T_0"]
  supernatant_vol <- chain[1, "supernatant_vol"]
  MOI_conversion_factor <- log10(T_0 / supernatant_vol)
  chain <- chain[,par_names_old]
  colnames(chain) <- par_names_new
  chain[,"log_10(n)"] <- -chain[,"log_10(n)"]
  chain[,MOI_params] <- chain[,MOI_params] + MOI_conversion_factor
  corrs <- cor(chain)
  quantiles <- t(apply(chain, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  max_LL <- readRDS(paste0(strain_dir, "max_LL_params.rds"))
  max_LL <- max_LL[par_names_old]
  names(max_LL) <- par_names_new
  
  all_mat <- cbind(data.frame(strain = strain, parameter = par_names_new),
                   quantiles, data.frame("max likelihood" = max_LL), corrs)
  write.table(all_mat, paste0(dirs[["figs"]], "pars.csv"), sep = ",",
              row.names = FALSE, append = TRUE)
  invisible(all_mat)
  
}

#' Emulates the ggplot colour palette.
#'
#' \code{gg_color_hue} emulates the ggplot colour palette.
#' 
#' @param n numeric vector of length 1. The number of colours required.
#' @return a character vector of length \code{n}. Each element is a hex colour code.
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Specifies the \code{y} scale of a \code{ggplot} object.
#'
#' \code{add_y_scale} specifies the \code{y} scale of a plot.
#' 
#' Called by \code{\link{plot.results.pyro}} to unify plots.
#' @param in_g_obj \code{ggplot} object to specify scale for.
#' @param log_scale logical vector of length 1. If true, use log scale on 
#' \code{y}-axis; use linear scale otherwise.
#' @param y_label character vector of length 1. Y-axis label.
#' @param y_max numeric vector of length 1. Upper limit of y-axis.
#' @return a \code{ggplot} object: \code{in_g_obj} with specified \code{y}-scale.
#' @import ggplot2
add_y_scale <- function(in_g_obj, log_scale = TRUE, y_label, y_max){
    if(log_scale) {
        scale_fn <- "scale_y_log10"
        y_min <- 1
    } else {
        scale_fn <- "scale_y_continuous"
        y_min <- 0
    }
    in_g_obj + do.call(scale_fn, list(y_label, expand = c(0,0))) +
        coord_cartesian(ylim = c(y_min,y_max))
}

#' construct plot filename
#' 
#' @param filename character vector of length 1. experiment/chain identifier
#' @param par_name character vector of length 1. parameter name
#' @param description character vector of length 1. description of plot
#' @return character vector of length 1.
make_filename <- function(filename, par_name, description) {
  
    paste0(filename,"_", description, par_name,".pdf")
}

#' darken a vector of hex codes for colours
#' 
#' @param colour_in vector of hex codes for colours
#' @return vector of hex codes for darkened colours
darken <- function(colour_in) {
    offset_const <- -.3
    grDevices::adjustcolor(colour_in, offset = c(rep(offset_const, 3), 0))
}

#' lighten a vector of hex codes for colours
#' 
#' @param colour_in vector of hex codes for colours
#' @return vector of hex codes for lighened colours
lighten <- function(colour_in) {
    offset_const <- .3
    grDevices::adjustcolor(colour_in, offset = c(rep(offset_const, 3), 0))
}

#' define colour palette for plotting
#' @inheritParams plot_viral_load
#' @return character vector of colour hexes
define_palette <- function(data_source) {
    if(data_source == "PR8") {
        palette <- rev(gg_color_hue(2))
    } else if(data_source == "WT") {
        colours <- c("green" = "#497f1f", "blue" = "#2139f7", "red" = "#db2a05", "yellow" = "#e59026")
        palette  <- c("sH1N1-WT" = colours[["green"]], 
                      "pH1N1-WT" = colours[["blue"]], 
                      "H5N1-WT" = colours[["red"]], 
                      "H7N9-WT" = colours[["yellow"]])
        
    } else {
        stop("unknown data source")
    }
}

#' format a numeric vector of breaks nicely for printing
#' 
#' @param breaks numeric vector
#' @param limits numeric vector of length 2. If present, breaks must be within limits
#' @return list with two elements:
#' breaks: numeric vector of same length
#' break_names: character vector of same length
#' order_magnitude: highest order of magnitude of breaks

gen_nice_breaks <- function(breaks, limits){
    
    # find highest order of magnitude
    order_magnitude <- floor(log10(max(abs(breaks))))
    # divide breaks by that order of magnitude
    breaks <- breaks / (10 ^ order_magnitude)
    # find minimum difference between breaks
    min_diff <- min(abs(diff(breaks)))
    # find order of minimum difference between breaks
    order_min_diff <- floor(log10(min_diff))
    # round breaks accordingly
    breaks <- round(breaks / 10 ^ (order_min_diff - 1)) * 10 ^ (order_min_diff - 1)
    
    # keep breaks within limits
    if(!missing(limits)){
        # divide limits by that order of magnitude
        limits <- limits / (10 ^ order_magnitude)
        limits_temp <- double(2)
        
        while(limits_temp[1] == limits_temp[2]){
            limits_temp[1] <- ceiling(limits[1] / 10 ^ order_min_diff) * 10 ^ order_min_diff
            limits_temp[2] <- floor(limits[2] / 10 ^ order_min_diff) * 10 ^ order_min_diff
            order_min_diff <- order_min_diff - 1
        }
        
        limits <- limits_temp
        
        replace_breaks_under_limit <- function(breaks, limit){
            if(any(breaks < limit)){
                breaks <- breaks[breaks >= limit]
                breaks <- c(breaks,limit)
            }
            breaks
        }
        
        breaks <- replace_breaks_under_limit(breaks,limits[1])
        breaks <- -replace_breaks_under_limit(-breaks,-limits[2])
        
    }
    
    break_names <- as.character(breaks)
    breaks <- breaks * (10 ^ order_magnitude)
    list(breaks = breaks, break_names = break_names, 
         order_magnitude = order_magnitude)
}

#' version of ggsave which suppresses warnings that rows are removed from the data frame
#' @import ggplot2
ggsave_wch <- function(...) {
    h <- function(w) {
        if(any(grepl("Removed", w))) {
            invokeRestart("muffleWarning")
        }
    }
    withCallingHandlers(ggsave(...), warning = h)
}

#' wrapper function to feed arguments to ggsave
#' 
#' @param g either a ggplot object, or a list of ggplot objects of length m,
#' or a nested list of ggplot objects (length n, each of length m)
#' @param filename_fn function which makes a filename by taking n arguments and
#' outputting a character vector of length 1. If g is a ggplot object or a list
#' of ggplot objects of length m, then n = 1
#' @param width numeric vector of length 1. width of figure in cm
#' @param height numeric vector of length 1. height of figure in cm
#' @param filename_args if g is a ggplot object, a charater vector of length 1.
#' if g is a list of ggplot objects of length m, a character vector of length m.
#' if g is a nested list of ggplot objects, a list of length 2, containing a 
#' character vector of length n and a character vector of length m respectively.
#' @return the original ggplot object
#' @import ggplot2
ggsave_wrapper <- function(g, filename_fn, width, height, filename_args) {
    
    ggsave_dims <- function(x, y) ggsave_wch(filename_fn(x), y, width = width, height = height, units = "cm")
    ggsave_dims2 <- function(x1, x2, y) ggsave_wch(filename_fn(x1, x2), y, width = width, height = height, units = "cm")
    
    if(is.ggplot(g[[1]][[1]])) {
        idx_matrix <- expand.grid(seq_along(g), seq_along(g[[1]]))
        apply(idx_matrix, 1, function(x) ggsave_dims2(filename_args[[1]][x[1]], filename_args[[2]][x[2]], g[[x[1]]][[x[2]]]))
    } else if(is.ggplot(g)){
        ggsave_dims(filename_args[[1]][1], g)
    } else {
        Map(ggsave_dims, filename_args[[1]], g)
    }
    invisible(g)
}

#' plot a nicer-looking histogram
#' 
#' @param vec numeric vector of valuse from which to make histogram
#' @param special_value optional parameter.  If specified, should be numeric 
#' vector of length 1. Plot the bar for this value in a different colour. Used, 
#' for example, to draw attention to the bar at x = 0.
#' @param breaks optional parameter.  If specified, should be numeric vector.
#' Specifies breaks on horizontal axis.  If not given, use horizontal axis 
#' limits as breaks.
#' @param vline_values optional parameter.  If specified, should be numeric
#' vector. Plot vertical lines at these values.
#' @param limits optional parameter. If specified, should be numeric vector of
#' length 2.  Use as horizontal axis limits. If not given use the range of vec
#' (if there are values below 0) or [0, max(vec)] (if there are no values below
#' zero). 
#' @param plot_title optional parameter. If specified, should be character
#' vector of length 1. Title of plot.
#' @param plot_xlab character vector of length 1. Horizontal axis label.
#' @return ggplot object
#' @import ggplot2
plot_histogram_general <- function(vec, special_value, breaks, vline_values, limits, plot_title, plot_xlab){
  
  if(missing(breaks)){
    breaks <- limits
  }
  
  if(missing(limits)){
    if(any(vec < 0)){
      limits <- range(vec)
    } else {
      limits <- c(0,max(vec))
    }
  }
  
  breaks_list <- gen_nice_breaks(breaks, limits)
  
  if(breaks_list$order_magnitude == 1){
    cat_xlab_string <- " ($\\times 10$)"
  } else if (breaks_list$order_magnitude == 0) {
    cat_xlab_string <- ""
  } else {
    cat_xlab_string <- paste0(" ($\\times 10^{",breaks_list$order_magnitude,"}$)")
  }
  
  if(missing(special_value)){
    plot_df <- data.frame(vec = vec)
  } else {
    plot_df <- data.frame(vec = vec, is_special = vec == special_value)
  }
  
  g <- ggplot(plot_df, aes(x = vec)) +
    theme_bw() +
    scale_x_continuous(latex2exp::TeX(paste0(plot_xlab, cat_xlab_string)),
                       limits = limits, 
                       breaks = breaks_list$breaks, 
                       labels = breaks_list$break_names,
                       expand = c(0,0)) +
    scale_y_continuous("",breaks = NULL)
  
  if(!missing(plot_title)){
    g <- g + ggtitle(plot_title)
  }
  
  if(missing(special_value)){
    g <- g + geom_histogram(bins = 30) +
      theme(text = element_text(size = 32), aspect.ratio = 1)
  } else {
    g <- g + geom_histogram(aes(fill = is_special), bins = 30) +
      theme(text = element_text(size = 32), aspect.ratio = 1, legend.position = "none")
  }
  
  if(!missing(vline_values)) {
    g <- g + geom_vline(xintercept = vline_values, linetype = "dotted")
  }
  
  g
}