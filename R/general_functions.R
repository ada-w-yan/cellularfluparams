# complete doc

#' checks if file exists in rds or csv format, and reads rds file preferentially
#' 
#' checks if file exists in rds or csv format, and reads from original extension if given; 
#' else reads rds file preferentially
#' 
#' @param filename character vector of length 1: with or without extension
#' @return the read data frame
read_csv_or_rds <- function(filename,...){
  # check if extension given; if so, remove and keep track of original extension
  original_ext <- NULL
  if(as.logical(match(substr(filename,nchar(filename)-3,nchar(filename)),
                      c(".rds", ".csv"), nomatch = 0))){
    original_ext <- substr(filename,nchar(filename)-3,nchar(filename))
    filename <- substr(filename,1,nchar(filename)-4)
  }
  
  if(file.exists(paste0(filename,".rds"))){
    if(file.exists(paste0(filename,".csv"))){
      # if file exists in both extensions and original extension was not given, read rds
      if(is.null(original_ext)){
        warning(paste0(filename, " exists in both rds and csv -- reading rds"))
        try_read(filename,rds = TRUE,...)
      } else if(original_ext == ".rds"){
        # if file exists in both extensions and original extension was  given, 
        # read original extension
        try_read(filename,rds = TRUE,...)
      } else {
        try_read(filename,rds = FALSE,...)
      }
    } else {
      # if file only exists in one extension, read that extension
      try_read(filename,rds = TRUE,...)
    }
  } else if (file.exists(paste0(filename,".csv"))){
    try_read(filename,rds = FALSE,...)
    # throw error if file doesn't exist in either extension
  } else {
    stop(paste0(filename, " not found"))
  }
}

#' reads rds or csv with tryCatch
#' 
#' reads rds or csv with tryCatch
#' 
#' rationale for existence: if it's not known beforehand whether we're reading
#'  an rds or csv, may pass nonsensical arguments, so ignore if these are causing an error
#' 
#' @param filename character vector of length 1: with or without extension
#' @return the read data frame
try_read <- function(filename, rds, ...){
  tryCatch({
    if(rds){
      readRDS(paste0(filename,".rds"),...)
    } else {
      try_fread(paste0(filename,".csv"),...)
    }
  }, error = function(c){
    if(rds){
      readRDS(paste0(filename,".rds"))
    } else {
      try_fread(paste0(filename,".csv"))
    }
    warning(paste0("Additional arguments ignored when reading ", filename))
  })
}

#' data.table::fread with tryCatch
#' 
#' data.table::fread sometimes crashes.  If crashing, try read.csv then converting.
#' 
#' @param filename character vector of length 1 ending in .csv
#' @return the read data frame
try_fread <- function(filename, ...) {
  tryCatch ({
    data.table::fread(filename, ...)
  }, error = function(c) {
    data_table <- read.csv(filename, ...)
    data.table::as.data.table(data_table)
  })
}

#' thins a csv or rds file, saves it with new name, optionally deletes old file
#' 
#' thins a csv or rds file, saves it with new name, optionally deletes old file
#' 
#' @param filename character vector of length 1: with .csv or .rds extension
#' @param thin numeric vector of length 1 (must be in integer if not integer format):
#' thin every this many iterations
#' @param save_as_rds logical vector of length 1: if TRUE, save thinned chain as 
#' rds regardless of extension of filename; if FALSE, save as original extension
#' @param remove_after logical vector of length 1: remove old file if TRUE
#' @return NULL
thin_csv_or_rds <- function(filename, thin, save_as_rds, remove_after){
  
  # check if thin is an integer
  if(!all.equal(thin, as.integer(thin))){
    stop(paste0("thin = ", as.character(thin), ": should be integer"))
  }
  
  # read file
  filename_wo_ext <- substr(filename, 1, nchar(filename)-4)
  ext <- substr(filename, nchar(filename)-3, nchar(filename))
  if(ext == ".rds"){
    table <- readRDS(filename)
  } else if(ext == ".csv"){
    table <- data.table:: fread(filename)
  } else {
    stop(paste0(filename, " is neither csv nor rds"))
  }
  
  #thin
  table <- table[seq(1, nrow(table), by = thin),]
  
  #save thinned file
  if(save_as_rds || ext == ".rds"){
    saveRDS(table, paste0(filename_wo_ext, "_thin_", as.character(thin),".rds"))
  } else{
    write.table(table, paste0(filename_wo_ext, "_thin_", as.character(thin),ext),
                row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)
  }
  
  # remove old file
  if(remove_after){
    file.remove(filename)
  }
  invisible(NULL)
}

#' saves all variables in parent environment
#' 
#' same as save.image except for parent environment rather than global
#' 
#' @param filename character vector of length 1: no extension
#' @return NULL
save_all <- function(filename){
  save(list = ls(all.names = TRUE, envir = parent.frame()), file = paste0(filename,".RData"), envir = parent.frame())
}

#' makes a directory from a filename
#' 
#' makes a directory from a filename
#' 
#' @param filename character vector of length 1
#' @return NULL
mkdir_from_filename <- function(filename){
  dir_name <- dirname(filename)
  if(!dir.exists(dir_name)){
    dir.create(dir_name,recursive = TRUE)
  }
  invisible(NULL)
}

#' calculate geometric mean of a vector
#' 
#' @param x numeric vector
#' @param na.rm logical vector of length 1
#' @return numeric vector of length 1: geometric mean
gm_mean <- function(x, na.rm=TRUE){
  if(na.rm) {
    denom <- sum(!is.na(x))
  } else {
    denom <- length(x)
  }
  exp(sum(log(x[x > 0]), na.rm=na.rm) / denom)
}

#' a version of Map() where we specify the format of the output
#' 
#' @param FUN function
#' @param FUN.VALUE prespecified type of return value
#' @return vector with prespecified type of return value
Map_vapply <- function(FUN, FUN.VALUE, ...){
  Map_output <- Map(FUN, ...)
  vapply(Map_output, identity, FUN.VALUE)
}

#' solve an ODE with error handling
#' 
#' if an error occurs, solve with a smaller time interval
#' for solution elements except the time vector, set to tol if the solution is below tol
#' 
#' @param mod an ode_system object generated by odin
#' @param solving_time times at which to solve the ODE
#' @param tol numeric vector of length 1: tolerance
#' @return a deSolve object: solution of the ODE
solve_ODE_error_handling <- function(mod, solving_time, tol) {
  
  ## maximum number of times to try to solve with smaller time interval
  max_fail_iter <- 5
  
  sol <- tryCatch(mod$run(solving_time), error = function(e) "error") # solve ODEs
  
  ## if error occurred while solving ODEs, try smaller solving interval
  if(!is.matrix(sol)) {
    fail_iter <- 0
    n_div <- 10
    solving_time_temp <- solving_time
    while(fail_iter < max_fail_iter && !is.matrix(sol)) {
      # subdivide each time step into n_div steps
      solving_time_temp <- interpolate_vector(solving_time_temp, n_div)
      # resolve
      sol <- tryCatch(mod$run(solving_time_temp), error = function(e) "error")
      fail_iter <- fail_iter + 1
    }
    # retrive solution for original time steps
    sol <- sol[solving_time_temp %in% solving_time,]
  }
  
  # if ODE solver exits early, pad with last value
  sol <- sol[sol[,1] %in% solving_time,]
  if(nrow(sol) < length(solving_time)) {
    last_row <- sol[nrow(sol),]
    rep_last_rows <- rep(list(last_row), length(solving_time) - nrow(sol))
    sol <- do.call(rbind, c(list(sol), rep_last_rows))
    sol[,1] <- solving_time
  }
  
  # for solution elements except the time vector, set to tol if the solution is below tol
  non_time_elements <- sol[,-1]
  non_time_elements[non_time_elements < tol] <- tol
  sol[,-1] <- non_time_elements
  
  sol
}

#' interpolate ordered numeric vector into n_div intervals betwen successive values
#' 
#' @param vec numeric vector to interpolate
#' @param n_div numeric vector of length 1: number of intervals into which to subdivide
#' successive values of the vector
#' @return numeric vector of length (length(vec) - 1) * n_div + 1: interpolated vector
interpolate_vector <- function(vec, n_div) {
  x <- seq(1, length(vec), by = 1/n_div)
  temp <- approx(seq_along(vec), vec, xout = x)
  vec_out <- temp$y
  vec_out
}

#' compile an odin model
#' 
#' @param mod_path character vector of length 1: path to ODEs
#' @return compiled model
compile_model <- function(mod_path){
  dir_name <- "R/"
  if(substr(mod_path, 1, nchar(dir_name)) != dir_name) {
    mod_path <- paste0(dir_name, mod_path)
  }
  odin::odin(mod_path, verbose=FALSE) # compile model
}


#' a version of lapply() that supplies FUN with both the name and the value of each component
#' 
#' @param X see documentation for lapply
#' @param FUN function which takes a character vector of length 1 as a first argument and 
#' something else as a second argument
#' @return see documentation for lapply
lapply_w_name <- function(X, FUN){
  Map(FUN, names(X), unname(X))
}

#' a version of a function that suppresses warnings
#' 
#' @param f a function
#' @return the same function without warnings
no_warn <- function(f) {
  function(...) {
    suppressWarnings(f(...))
  }
}

#' make string from number in format acceptable for filename
#' 
#' \code{make_filename_from_number} makes a string from a number, replacing 
#' periods and plus signs
#' @param x numeric vector of length 1
#' @param decimal_points numeric vector of length 1. number of decimal points to
#'   round to
#' @param scientific logical vector of length 1.  whether to use scientific
#'   notation
#' @return character vector of lenth 1: string suitable for use as filename
make_filename_from_number <- function(x, decimal_points = 3, scientific = FALSE)
{
  if(scientific) {
    formatted_x <- format(signif(x, decimal_points + 1), scientific = TRUE)
  } else {
    formatted_x <- round(x, digits = decimal_points)
  }
  formatted_x <- sub(".", "point", formatted_x, fixed = TRUE)
  formatted_x <- sub("+", "", formatted_x, fixed = TRUE)
  formatted_x
}


#' load contents of RData file into a list
#' 
#' @param filename character vector of length 1: RData filename
#' @param var_names optional argument specifying variables to load.  If missing,
#' load all variables
#' @return list containing variables with names var_names
load_into_list <- function(filename, var_names) {
  load_env <- new.env()
  load(filename, envir = load_env)
  output <- as.list(load_env)
  if(!missing(var_names)) {
    output <- get_vars_from_list_with_check(output, var_names)
  }
  output
}

#' gather variables with given names from an environment into a list
#' 
#' @param var_names character vector of variable names to gather
#' @param envir environment in which to find parameters
#' @return list of variables with names var_names
list_vars_from_environment <- function(var_names, envir = parent.frame()) {
  env_vars <- as.list(envir)
  env_vars <- get_vars_from_list_with_check(env_vars, var_names)
  env_vars
}

#' dump variables from a list to the parent frame
#' 
#' @param x list containing variables
#' @param var_names character vector containing names of variables to dump into 
#' parent frame
#' @param overwrite optional parameter controlling the behaviour if any variables
#' with the same name(s) alraedy exist in the parent frame. If set to "warn",
#' throws a warning; if set to "error", throws an error
#' @return NULL
list2here <- function(x, var_names, overwrite) {
  if(!missing(var_names)) {
    x <- get_vars_from_list_with_check(x, var_names)
  }
  
  
  if(!missing(overwrite)) {
    parent_frame_vars <- ls(parent.frame)
    overwriting_vars <- intersect(names(x), parent_frame_vars)
    if(length(overwriting_vars) > 0) {
      if(overwrite == "warn") {
        lapply(overwriting_vars, function(x) warning(cat("overwriting", x)))
      } else if(overwrite == "error") {
        stop(cat("Attempting to overwrite", x[1]))
      }
    }
  }
  
  list2env(x, envir = parent.frame())
  invisible(NULL)
}

#' extract variables from list, throwing an error if they are not found
#' 
#' @param x list of variables
#' @param var_names character vector containing names of variables to extract
#' @return list of selected variables
get_vars_from_list_with_check <- function(x, var_names) {
  missing_vars <- var_names[!(var_names %in% names (x))]
  if(length(missing_vars) > 0) {
    stop(cat("variables missing from list: ", paste0(missing_vars, collapse = " ")))
  }
  x <- x[var_names]
  x
}

#' make factor out of a numeric vector
#' 
#' @param vec numeric vector
#' @return a character version of the vector, with levels in numeric order
make_factor <- function(vec) {
  factor(as.character(vec), levels = unique(as.character(sort(vec))))
}

#' return names of all functions in a .R file
#' 
#' @param filename location of .R file
#' @return character vector where each entry is the name of a function in the file
get_func_names <- function(filename) {
  suppress_final_line_warning <- function(w) {
    if(any(grepl("incomplete final line found on", w, fixed = TRUE))) {
      invokeRestart("muffleWarning")
    }
  }
  # read the .R file, suppressing warnings about incomplete final line
  code <- withCallingHandlers(readLines(filename), warning = suppress_final_line_warning)
  # remove spaces
  code <- gsub(" ", "", code)
  # identify lines which define functions
  func_str <- "<-function("
  potential_func_lines <- grep(func_str, code, fixed = TRUE)
  
  # determine whether each function is an outermost function
  is_outside_func <- function(code, potential_func_line) {
    # assume functions on first line are outside functions
    if(potential_func_line == 1) {
      return(TRUE)
    }
    # paste all code up to name of potential function
    pasted_code <- paste0(code[seq_len(potential_func_line - 1)], collapse = "")
    potential_func_subline <- strsplit(code, split = func_str, fixed = TRUE)
    potential_func_subline <- potential_func_subline[1]
    pasted_code <- paste0(pasted_code, potential_func_subline, collapse = "")
    # count number of open and close curly brackets up to potential function name
    count_characters <- function(pasted_code, char_in) {
      n <- gregexpr(char_in, pasted_code, fixed = TRUE)
      length(n[[1]])
    }
    n_brackets <- vapply(c("{", "}"), function(x) count_characters(pasted_code, x), numeric(1))
    # the functino is an outermost function if the number of open and close brackets is the same
    n_brackets[1] == n_brackets[2]
  }
  
  func_lines <- potential_func_lines[vapply(potential_func_lines, 
                                            function(x) is_outside_func(code, x),
                                            logical(1))]
  # split off the function names
  func_lines <- strsplit(code[func_lines], split = func_str, fixed = TRUE)
  func_lines <- vapply(func_lines, function(x) x[1], character(1))
  func_lines
}

#' a version of apply. if MARGIN = 1, the values in each column of X are passed to FUN
#' as named arguments according to the column name. if MARGIN = 2, the values in
#' each row of X are passed to FUN
#' as named arguments according to the row name
#' @param X see arguments for apply.  If MARGIN = 1, colnames(X) must match the named
#' arguments of FUN.
#' @param MARGIN see arguments for apply
#' @param FUN see arguments for apply.  Must have named arguments
#' @return see apply
apply_named_args <- function(X, MARGIN, FUN) {
  apply(X, MARGIN, function(X) do.call(FUN, as.list(X)))
}

#' make a default list of random number generator seeds for MCMC
#' 
#' @param n_replicates the number of parallel chains to run to assess convergence
#' @param n_temperatures the number of temperatures
#' @return an n_replicates by n_temperatures of integers
make_seed <- function(n_replicates, n_temperatures) {
  starting_point_seed <- lapply(seq_len(n_replicates),
                                function(x) (x - 1) * n_temperatures +
                                  seq_len(n_temperatures))
}

#' find the indices at which a vector changes sign
#' 
#' @param x a numeric vector
#' @return a logical vector of length length(x) - 1.  TRUE at element i indicates
#' that x[i] has a different sign from x[i+1].
find_sign_change_ind <- function(x) {
  which(abs(diff(sign(x))) == 2)
}

#' short for \code{vapply(X, FUN, logical(1), ...)}
#' 
#' short for \code{vapply(X, FUN, logical(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, logical(1), ...)}: logical vector of 
#' same length as X
vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}

#' short for \code{vapply(X, FUN, integer(1), ...)}
#' 
#' short for \code{vapply(X, FUN, integer(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, integer(1), ...)}: integer vector of 
#' same length as X
viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}

#' short for \code{vapply(X, FUN, numeric(1), ...)}
#' 
#' short for \code{vapply(X, FUN, numeric(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, numeric(1), ...)}: numeric vector of 
#' same length as X
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

#' short for \code{vapply(X, FUN, character(1), ...)}
#' 
#' short for \code{vapply(X, FUN, character(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, character(1), ...)}: character vector
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

#' check if x is within \code{tol} of \code{x}
#' 
#' @param x numeric vector of length 1
#' @param tol tolerance
#' @return logical vector of length 1
is_integer_like <- function(x, tol = sqrt(.Machine$double.eps)) {
  is.integer(x) || (is.numeric(x) && abs(x - round(x)) < tol)
}

#' find the file containing a function
#' 
#' @param func_name name of function
#' @param dir_name directory in which to search (only searches in that directory
#' + subdirectories one level down)
#' @return the name of the file which has the function
find_func_file <- function(func_name, dir_name = ".") {
  # list current directory + directories one level down
  dirs <- list.dirs(path = dir_name, recursive = FALSE)
  dirs <- c(dirs, dir_name)
  # exclude Rproj.user -- not a relevant directory
  dirs <- dirs[!grepl("Rproj.user", dirs, fixed = TRUE)]
  # list files in those directories
  filenames <- lapply(dirs, function(x) list.files(x, pattern="\\.R$", full.names=TRUE))
  # list functions in those files
  func_names <- lapply(filenames, function(x) lapply(x, get_func_names))
  # find index of file containing target function
  in_file <- lapply(func_names, function(x) vlapply(x, function(y) any(grepl(func_name, y))))
  # find directory of that file
  in_dir <- vlapply(in_file, any)
  # return filename
  func_file <- filenames[in_dir][[1]][in_file[in_dir][[1]]]
  func_file
}

#' loads a character vector of packages
#' 
#' loads a character vector of packages
#' 
#' @param package_vector: character vector to load
#' @return NULL
load_packages <- function(package_vector){
  lapply(package_vector,function(x) do.call(library,list(x)))
  invisible(NULL)
}

#' sources a character vector of files
#' 
#' sources a character vector of files
#' 
#' @param source_vector: character vector of files to source
#' @return NULL
source_files <- function(source_vector){
  lapply(source_vector,function(x) source(x))
  invisible(NULL)
}

#' calculate the nth moment from a probability mass function
#' 
#' @param x numeric vector. x-values of probability mass function.  Must include
#' all values where the probability mass function is non-zero.
#' @param y numeric vector of the same length as x.  Values of probability mass
#' function at the points in x.
#' @param n integer >= 0.  Calculate nth moment.
#' @param central logical with default value FALSE.  If TRUE, calculate central
#' moment, otherwise calculate raw moment.
#' @return numeric vector of length 1: the nth moment
calc_moment_from_pmf <- function(x, y, n, central = FALSE) {
  stopifnot(n >= 0 && is_integer_like(n))
  stopifnot(length(x) == length(y))
  stopifnot(all.equal(sum(y), 1))
  if(n == 0) {
    return(0)
  }
  
  if(central) {
    first_moment <- calc_moment_from_pmf(x, y, n = 1, central = FALSE)
    adjustment <- first_moment
  } else {
    adjustment <- 0
  }
  
  moment <- sum(y * (x - adjustment)^n)
  moment
}

#' calculate percentiles for fitted parameters
#' 
#' calculate percentiles for fitted parameters
#' 
#' @param chain the MCMC chain: a data frame with n columns, where n is the number
#' of fitted parameters
#' @param prctiles a numeric vector of length m containing the percentiles to be calculated
#' (between 0 and 1)
#' @param par_names_plot parameter names for the output data frame
#' @return a data frame with n rows and m columns containing the percentiles
#' @export 
print_prctiles <- function(chain, prctiles = c(.025,.5,.975), par_names_plot){
  prctile_table <- lapply(chain,function(x) quantile(x,prctiles, na.rm = TRUE))
  prctile_table <- t(as.data.frame(prctile_table))
  rownames(prctile_table) <- par_names_plot
  col_names <- as.character(prctiles*100)
  col_names <- trimws(format(col_names, digits = 3))
  col_names <- paste0(col_names,"\\%")
  colnames(prctile_table) <- col_names
  prctile_table
}