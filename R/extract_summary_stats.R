# complete doc

#' calculate and save cellular infection parameters for MCMC chain
#' 
#' @param data_dir directory where 1.RData lives
#' @return 3xn matrix where each row corresponds to a cellular infection parameter and n is
#' the number of MCMC iterations
#' @export
extract_cell_infect_params <- function(data_dir) {
  filename <- paste0(data_dir, "1.RData")
  print(getwd())
  load(filename)
  chain <- get_MCMC_chain(filename)
    temp_list <- gen_summary_statistics_fn(parTab, parTab, "constant")
    calc_summary <- temp_list$calc_summary
    par_names <- parTab$names
    cell_infect_params <- apply(chain, 1,
                           function(x) calc_summary(as.numeric(x[par_names])))

  saveRDS(cell_infect_params, paste0(data_dir, "1cell_infect_params.rds"))
  cell_infect_params
}

#' test whether the ratio of a summary statistic between two strains
#' differs from 1
#' 
#' @param data_dir the n directories where the 1.RData files live
#' @inheritParams plot_stats_across_strains
#' @inheritParams plot_viral_load
#' @return a list of length n, where n is the number of pair combinations of data_dir.
#' each element contains a list with the elements prctile_ratios and p_value.
#' prctile_ratio is a matrix where each column corresponds to a summary statistic,
#' and the rows give the 2.5, 50 and 97.5th percentiles of the ratio of each
#' statistic between the two strains.
#' p_value is a named vector where each element is twice the proportion of the
#' sampled ratios greater than 1, or the proportion of the sampled ratios smaller
#' than 1 (whichever is the smaller)
#' @export
pairwise_stats_tests <- function(data_dir, sum_stat, data_source) {
  strain <- determine_strains_for_exp(data_source)

  # calculate the summary statistics if haven't already
  if(sum_stat == "cell_infect_params") {
    filenames <- paste0(data_dir, "1", sum_stat, ".rds")
  } else {
    filenames <- paste0(data_dir, sum_stat, ".rds")
  }

  # read in the summary stats
  stats <- lapply(filenames, readRDS)
  if(sum_stat == "mean_gen_no") {
    convert_to_matrix <- function(x) {
      x <- matrix(x, nrow = 1)
      rownames(x) <- "mean_gen_no"
      x
    }
    stats <- lapply(stats, convert_to_matrix)
  }

  names(stats) <- strain
  
  # undo log transforms of summary stats
  
  undo_log10 <- function(x) {
    log10_var <- grepl("log10", rownames(x))
    x[log10_var,] <- 10^x[log10_var,]
    rownames(x) <- sub("log10_", "", rownames(x))
    rownames(x) <- sub("log10", "", rownames(x))
    x
  }
  
  stats <- lapply(stats, undo_log10)

  # generate all pairwise combinations of strains
  pairs <- combn(strain, 2)
  pairs <- lapply(seq_len(ncol(pairs)), function(x) pairs[,x])
  
  n_samples <- 1e6
  
  sample_ratio <- function(matrix_x, matrix_y, samples_x, samples_y) {
    matrix_y[,samples_y] / matrix_x[,samples_x]
  }
  
  # take n_samples samples of each summary statistic for each strain
  samples <- lapply(stats, function(x) sample.int(ncol(x), size = n_samples, replace = T))

  # calculate the ratios of the sampled summary statistics between the strains
  ratio_samples <- lapply(pairs, function(x) sample_ratio(stats[[x[1]]], stats[[x[2]]],
                                                            samples[[x[1]]], samples[[x[2]]]))
  
  # calculate the percentiles and p values of the ratios
  prctile <- c(2.5, 50, 97.5)
  
  calc_prctile_and_p_value <- function(y) {
    if(!is.matrix(y)) {
      y <- matrix(y, nrow = 1)
      rownames(y) <- rownames(stats[[1]])
    }
    prctile_ratios <- apply(y, 1, function(x) quantile(x, probs = prctile / 100,
                                                       na.rm = TRUE))
    calc_p_value <- function(x) {
      x <- x[!is.na(x)]
      min(sum(x > 1), sum(x < 1)) / length(x) * 2
    }
    p_value <- apply(y, 1, calc_p_value)
    list(prctile_ratios = prctile_ratios, p_value = p_value)
  }
  lapply(ratio_samples, calc_prctile_and_p_value)
}

#' from a summary stats .rds file, get a subset of the summary stats, undo log10 and thin
#' 
#' needed
#' @param filename summary stats .rds file
#' @param my_summary_stats_names vector of summary stats to subset
#' @param n_samples number of samples to take.  If 0, take all samples
#' @param undo_log10 if TRUE, for summary stats beggining with log10, take 10^
#' change summary stat name to reflect. if FALSE, do nothing
#' @return subsetted and thinned matrix of summary stats
subset_summary_stats_from_rds <- function(filename, my_summary_stats_names,
                                          n_samples = 1e3, undo_log10 = TRUE) {
  # read summary stats file
  summary_stats <- readRDS(filename)
  
  # thin uniformly
  if(n_samples > 0 && ncol(summary_stats) > n_samples) {
    summary_stats <- summary_stats[,round(seq(1, ncol(summary_stats), length.out = n_samples))]
  }

  # undo log10
  summary_stats <- summary_stats[my_summary_stats_names,]
  if(undo_log10) {
    summary_stats[grep("log10", my_summary_stats_names),] <- 10^summary_stats[grep("log10", my_summary_stats_names),]
    my_summary_stats_names <- gsub("log10_", "", my_summary_stats_names)
  }
  rownames(summary_stats) <- my_summary_stats_names

  summary_stats
}