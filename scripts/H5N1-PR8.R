devtools::load_all("~/git_repos/cellularfluparams/")
library(lazymcmc)
library(parallel)
run_exp("H5N1-PR8", "main", FALSE)