# complete doc

#' extract name of data source from strain name
#' 
#' @param strain name of strain
#' @return character vector of length 1
#' @export
determine_data_source_for_strain <- function(strain) {
  strsplit(strain, "-")[[1]][2]
}

#' remove data source from strain name
#' @inheritParams determine_data_source_for_strain
#' @return character vector of length 1
#' @export
shorten_strain_name <- function(strain) {
  strsplit(strain, "-")[[1]][1]
}

#' get directory names containing MCMC output for strains
#' 
#' @param dirs output of set_dirs
#' @param strains character vector of strain names
#' @return character vector of same length as strains
#' @export
get_strain_dir <- function(dirs, strains) {
  paste0(dirs[["base"]], strains, "/")
}

#' list all strains for a data source
#'
#' @param data_source character string: "WT", "PR8" or "all"
#' @param w_data_source logical. if TRUE, append data source to strain name
#' (compulsory if data_source == "all")
#' @return character vector containing strain names
#' @export
determine_strains_for_exp <- function(data_source, w_data_source = FALSE) {
  strains <- switch(
    data_source,
    "WT" = c("sH1N1", "pH1N1", "H5N1", "H7N9"),
    "PR8" = c("pH1N1", "H5N1"),
    "all" = unlist(lapply(c("WT", "PR8"), 
                          function(x) paste0(determine_strains_for_exp(x), "-", x)))
  )
  if(data_source %in% c("WT", "PR8") && w_data_source) {
    strains <- paste0(strains, "-", data_source)
  }
  strains
}

#' set directories to store computation results and plots
#' 
#' @inheritParams run_exp
#' @return list with element "base": character: directory name for computation results
#' "figs": directory name for figures
#' @export
set_dirs <- function(sensitivity) {
  paste_sensitivity_str <- function(base_dir, sensitivity) {
    paste0(base_dir, switch(
      sensitivity,
      "main" = "",
      "n_L_n_I" = "60/",
      paste0(sensitivity, "/")
    ))
  }
  # dirs <- c("temp/gen_time_paper/", "temp/gen_time_paper/figs/")
  dirs <- c("results/", "figs/")
  dirs <- paste_sensitivity_str(dirs, sensitivity)
  
  names(dirs) <- c("base", "figs")
  create_dir <- function(dir_name) {
    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
    }
  }
  invisible(lapply(dirs, create_dir))
  dirs
}

#' extract the value of p_over_time from the sensitivity input parameter
#' @inheritParams run_exp
#' @return "linear" or "constant"
get_p_over_time_from_sensitivity <- function(sensitivity) {
  switch(sensitivity,
         "linear_p" = "linear",
         "constant")
}