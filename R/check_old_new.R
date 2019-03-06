check_old_new <- function(filenames) {
  output_list <- lapply(filenames, load_into_list)

  for (i in seq(2, ncol(output_list[[1]]$parTab))) {
    print(all.equal(output_list[[1]]$parTab[,i], output_list[[2]]$parTab[,i]))
  }
  print(all.equal(output_list[[1]]$data, output_list[[2]]$data))
  chains <- lapply(filenames, get_MCMC_chain)
  print(all.equal(colnames(chains[[1]]), colnames(chains[[2]])))
  invisible(NULL)
}